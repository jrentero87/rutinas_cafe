# -*- coding: utf-8 -*-
"""
@author: Jesús Rentero Bonilla
Rutina 2: Posición e intensidad del flat
Objetivo: Determinar la posición de los órdenes en ciertas columnas del CCD
          para controlar si ha habido algún desplazamiento en X-Y

"""

from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import ascii
import matplotlib.pyplot as plt
import os.path
from dateutil import parser
import astropy.time
import matplotlib.gridspec as gridspec # GRIDSPEC !
import datetime
from jdcal import gcal2jd

"""
Fichero que almacena las posiciones de cada uno de los ordenes medidas con el DS9 para la columna central
"""
INPUT_ORDEN="./ordenes_input.txt"

"""
Fichero que almacena el ajuste de las posiciones de cada orden para una imagen flat que tomamos como referencia.
"""
AJUSTE_INICIAL="./ordenes_inicial.txt"

"""
Constante donde se almacena el nombre del fichero Master para la rutina 01. En él se almacenarán las desviaciones medias de cada noche.
"""
FICH_MASTER="./Rut02_dat/ordenes_master.txt"

"""
Función que devuelve la gausiana con los parámetros:
- x: altura de la campana.
- a: amplitud de la curva.
- x0: centro de la gaussiana.
- sigma: anchura de la gausiana.
- zero: valor umbral.
"""
def gaus(x,a,x0,sigma,zero):
    return zero + a*np.exp(-(x-x0)**2/(2*sigma**2))


"""
Funcion que obtiene la matriz de datos a partir de una imagen de flat.
"""
def getMatrizDatos(arcoFits):
    # Abrimos el fichero de arco
    hdulist=fits.open(arcoFits);
    #Obtenemos la matriz con los datos
    tbdata = hdulist[0].data
    #cerramos el fichero
    hdulist.close();
    #Hallamos la matriz traspuesta, puesto que hdulist contiene la matriz traspuesta de la imagen
    tbdata_traspuesta=tbdata.transpose()
    return tbdata_traspuesta

"""
Funcion que devuelve un vector con las posiciones iniciales de cada orden.
Estas posiciones deben haberse calculado para el centro de la imagen y 
deberán estar en un fichero cuyo formato en cada linea será: Id_orden#PosY
donde:
- Id_orden: es el identificador de cada uno de los ordenes
- PosY: es la posición Y del orden en el centro de la imagen
"""
def getConfiguracion(fichero):
    # Abrimos el fichero que contiene la configuración inicial
    infile = open(fichero,'r')
     #Creamos un vector vacio para almacenar cada posicion
    ordenesPosY=[]
    for line in infile:
        #Troceamos la linea, almacenando en orden[0] el Id, y en orden[1] la posición
        orden=line.split(",")
        idOrden=orden[0]
        #Comprobamos que la linea no sea un comentario, es decir, que no comience por @
        if idOrden[0] != "@":
            # Almacenamos la posicion del orden en el vector
            ordenesPosY.append(int(orden[1]))
    # Cerramos el fichero
    infile.close()
    return ordenesPosY

"""
Funcion que a partir del fichero fits con cada orden y el fichero de configuración,
genera una matriz por cada coeficiente que se ajuste. En dicha matriz contendrá el
valor del coeficiente para cada columna en la imagen.
Para el cálculo, se ha cogido la columna central de la imagen y 17 columnas a la izquierda
y a la derecha de la columna central con la separación de 60 píxeles.
"""
def generarAjuste(fich_ordenes, fich_conf):
    # Obtenemos la matriz con los datos del fichero fits
    mat=getMatrizDatos(fich_ordenes)
    # Obtenemos las posiciones del fichero de configuración de cada uno de los órdenes
    posiciones=getConfiguracion(fich_conf)
    # Definimos el rango de los pixeles de la imagen
    XX = np.arange(0,2048)
    # Fijo la posición de la primera columna en el centro de la imagen
    posX=1024
    # Vector y matriz para almacenar las nuevas posiciones de los ordenes
    newPos=[]
    matPosY=[]
    matPosX=[]
    # Vector y matriz para almacenar los nuevos valores de sigma
    newSigma=[]
    matSigma=[]
    # Vector y matriz para almacenar los nuevos valores umbrales
    newUmbral=[]
    matUmbral=[]
    # Definimos el salto entre cada columna
    salto=60
    # Repetimos este proceso de ajuste para la columna central y para 17 columnas más a la izquierda de esta
    for i in range(17):
        # Sumo para la posición determinada 5 columnas y obtengo sus valores en un vector YY
        YY = np.sum(mat[posX:posX+5,:], axis=0)
        for y0 in posiciones:
            x=XX[y0-9:y0+9]
            y=YY[y0-9:y0+9]
            try:
                p0=[np.max(y)-y[0],y0,2.,y[0]]
                # Realizamos el ajuste
                coeff,pcov = curve_fit(gaus,x,y,p0)
                # Almacenamos los valores obtenidos en los vectores para cada orden
                newPos.append(coeff[1])
                newSigma.append(coeff[2])
                newUmbral.append(coeff[3])
            except:
                #Almacenamos los valores obtenidos en los vectores para cada orden
                newPos.append(y0)
                newSigma.append(0.0)
                newUmbral.append(0.0)
                print "Rutina 2 WARNING: el flat "+fich_ordenes+" no se ajustó correctamente"
            
        matPosX.insert(0,(posX+posX+5.)/2.)
        #Desplazamos la posicion de la columna "salto" pixeles a la izquierda
        posX=posX-salto
        #Actualizamos el nuevo vector de posiciones
        posiciones = np.array(newPos[:])- (np.array(newPos[:])*4./2048.)
        #Almacenamos los datos en las matrices correspondientes, al inicio de las mismas
        matPosY.insert(0,posiciones)
        matSigma.insert(0,newSigma)
        matUmbral.insert(0,newUmbral)
        
        #print matPosX
        # Vaciamos los vectores
        newPos=[]
        newSigma=[]
        newUmbral=[]
    
    #Repetimos el proceso de ajute para los valores que están a la derecha de la columna central
    #Inicializamos el vector de posiciones al valor obtenido para la columna central
    #posiciones=np.array(matPosY[len(matPosY)-1:])[0]
    posiciones=getConfiguracion(fich_conf)
    posX=1024+salto
    for i in range(17):
        YY = np.sum(mat[posX:posX+5,:], axis=0)
        for y0 in posiciones:
            x=XX[y0-9:y0+9]
            y=YY[y0-9:y0+9]
            try:
                p0=[np.max(y)-y[0],y0,2.,y[0]]
                # Realizamos el ajuste
                coeff,pcov = curve_fit(gaus,x,y,p0)
                # Almacenamos los valores obtenidos en los vectores para cada orden
                newPos.append(coeff[1])
                newSigma.append(coeff[2])
                newUmbral.append(coeff[3])
            except:
                #Almacenamos los valores obtenidos en los vectores para cada orden
                newPos.append(y0)
                newSigma.append(0.0)
                newUmbral.append(0.0)
                print "Rutina 2 WARNING: el flat "+fich_ordenes+" no se ajustó correctamente"
        
        matPosX.append((posX+posX+5.)/2.)
        #Desplazamos la posicion de la columna "salto" pixeles a la derecha
        posX=posX+salto
        #Actualizamos el nuevo vector de posiciones
        posiciones = np.array(newPos[:])+ (np.array(newPos[:])*4./2048.)
        #Almacenamos los datos en las matrices correspondientes, al final de las mismas
        matPosY.append(posiciones)
        matSigma.append(newSigma)
        matUmbral.append(newUmbral)
        
        # Vaciamos los vectores
        newPos=[]
        newSigma=[]
        newUmbral=[]
    
    return matPosY,matSigma,matUmbral,matPosX

"""
Función que se encarga de escribir el contenido de una matriz en un fichero
"""
def escribirMatriz(mat,nomFichero):
    matriz=np.array(mat).transpose()
    outfile = open(nomFichero,"w")
    for i in range(0,len(matriz)):
        for j in range(0,len(matriz[i])):
            dato=round(matriz[i][j],4)
            outfile.write(str(dato))
            if j<len(matriz[i])-1:
                outfile.write(",")
        outfile.write("\n")
    outfile.close()
            

"""
Esta funcion se encarga de generar el ajute para una lista de ficheros de flat.
Este listado de flats vendrá dado en un fichero que se le pasará a la función por parámetro.
Se generará un fichero con las posiciones de las órdenes de cada una de las imagenes flat
"""    
def rutina02Run(listaFlat):
    # Creamos una lista para almacenar todos los ajustes de una noche
    listaAjustes=[]
    # Abrimos el fichero con el listado de ficheros flat
    infile = open(listaFlat,'r')
    # Procesamos cada una de las lineas del fichero, y generamos las estadísticas para cada fichero de flat
    for line in infile:
        #Eliminamos de la linea el retorno de carro (\n)
        line=line.strip()
        #Comprobamos que la linea tenga información y no sea una linea en blanco
        if len(line)>0:
            #Obtenemos el ajuste de cada orden
            matPos,matSigma,matUmbral,matPosX = generarAjuste(line,INPUT_ORDEN)
            #Escribimos en un fichero el resultado
            nomFichero = line[0:len(line)-5]+"_"+line[0:6]+"_dat.txt"
            escribirMatriz(matPos,"./Rut02_dat/"+nomFichero[nomFichero.index('/')+1:])
            listaAjustes.append(np.array(matPos))
    # Realizamos el chequeo
    checkRutina02(listaAjustes, listaFlat)


"""
Función que se encarga de genera el fichero Master de la rutina y de chequear los datos
Se le proporciona una lista con el ajuste de todos los ficheros flat de una noche
"""
def checkRutina02(listaAjustes, listaFlat):
    # Obtenemos la matriz con el ajuste inicial con el que calcularemos la desviación de las posiciones de cada orden
    table = ascii.read(AJUSTE_INICIAL, format='csv')
    ajusteInicial=np.array(table)
    
    # Creamos un vector donde almacenaremos las desviaciones de los ordenes 10, 40 y 70
    desviacionO10=[]
    desviacionO40=[]
    desviacionO70=[]
    
    # Calculamos la desviación de la posición de cada orden con respecto al flat inicial para la columna central (17)
    for ajuste in listaAjustes:
        desv10=ajusteInicial[9][17]-ajuste[17][10]
        desviacionO10.append(desv10)
        desv40=ajusteInicial[39][17]-ajuste[17][40]
        desviacionO40.append(desv40)
        desv70=ajusteInicial[69][17]-ajuste[17][70]
        desviacionO70.append(desv70)
    # Calculamos las desviaciones medias para cada orden
    desvMedia10=np.mean(desviacionO10)
    desvMedia40=np.mean(desviacionO40)
    desvMedia70=np.mean(desviacionO70)
    
    # Añadimos el valor de la media de todas las desviaciones de los ordenes 10 40 70 de la noche al fichero master
    # Si existe el fichero lo abrimos en modo "a", sino lo creamos
    if os.path.exists(FICH_MASTER):
        file=open(FICH_MASTER,"a")
    else:
        file=open(FICH_MASTER,"w")
        file.write("@juldate,desv_Orden10,desv_Orden40,desv_Orden70\n")
    # Abrimos el fichero con el listado de ficheros flat
    infile = open(listaFlat,'r')
    # Obtenemos el dia juliano para uno de los ficheros flat de la noche
    imagen = infile.readline().strip()
    infile.close()
    juldate=getDiaJuliano(imagen)
    #Comprobamos que la entrada en el fichero no exista. En caso de no existir escribimos nueva entrada
    if not existeNoche(np.int(juldate)):
        file.write(str(np.int(juldate))+","+str(round(desvMedia10,4))+","+str(round(desvMedia40,4))+","+str(round(desvMedia70,4))+"\n")
    file.close()
    
    #Realizamos el checkeo para la rutina02
    # Si las desviciones medias de los ordenes son menores a 20 milipíxeles
    if desvMedia10 < 0.02 and desvMedia10 > -0.02:
        print "... Desviación media del orden 10: %.2f ... OK"%(desvMedia10)
    else:
        print "... Desviación media del orden 10: %.2f ... NO OK! - CHECK"%(desvMedia10)
    if desvMedia40 < 0.02 and desvMedia40 > -0.02:
        print "... Desviación media del orden 40: %.2f ... OK"%(desvMedia40)
    else:
        print "... Desviación media del orden 40: %.2f ... NO OK! - CHECK"%(desvMedia40)
    if desvMedia70 < 0.02 and desvMedia70 > -0.02:
        print "... Desviación media del orden 70: %.2f ... OK"%(desvMedia70)
    else:
        print "... Desviación media del orden 70: %.2f ... NO OK! - CHECK"%(desvMedia70)
    
    
    


"""
Esta función devolverá true en caso de que exista en el fichero bias_master.txt 
una entrada para la noche que se introduce por parámetro, y false en caso contrario.
El dia juliano introducido por parámetro debe ser un valor entero.
"""
def existeNoche(diaJuliano):
    # Abrimos el fichero en modo lectura
    infile=open(FICH_MASTER,'r')
    # Creamos una variable que inicialmente inicializamos a falso
    existe=False
    # Recorremos el fichero y comparamos cada una de las lineas
    for line in infile:
        #Ignoramos las lineas que comiencen por @, puesto que se trata de un comentario en el fichero
        if line[0]!='@':
            # Obtenemos el dia juliano almacenado en la primera posicion de la linea
            linea=line.split(',')
            juldate=np.int(linea[0])
            if juldate==diaJuliano:
                existe=True
    infile.close()
    return existe  

"""
Funcion que obtiene el día juliano a partir de una imagen fit que se le pasa por parámetro.
"""
def getDiaJuliano(imagenFit):
    # Abrimos el fichero
    f=fits.open(imagenFit);
    # Obtenemos la fecha y finalmente con astropy obtenemos el dia juliano
    date=f[0].header["DATE"]
    dt = parser.parse(date)
    time = astropy.time.Time(dt)
    juldate = time.jd
    return juldate
        
    
"""
Funcion que se encarga de generar el fichero inicial con el que se van a comparar los demas ficheros flats
Para ello, se toma una imagen FLAT de referencia (flat_ref)
"""
def cargarAjustes(flat_ref):
# Generamos el fichero inicial con el que se van a comparar los demás.
    matPos, matSigma, matUmbral, matPosX=generarAjuste(flat_ref, INPUT_ORDEN)
    escribirMatriz(matPos, AJUSTE_INICIAL)

"""
Funcion encargada de añadir pintar y añadir al historial los resultados obtenidos en la noche que se esta ejecutando
"""
def plotHistory():
    colnames = ('jd','desv10','desv40','desv70')
    table = ascii.read('Rut02_dat/ordenes_master.txt', format='csv', names=colnames, comment='@')	
    jd  = np.array(table["jd"])
    desv10   = np.array(table["desv10"])
    desv40   = np.array(table["desv40"])
    desv70   = np.array(table["desv70"])
    
    today = datetime.datetime.now()
    today = astropy.time.Time(today)
    jd_today = np.int(today.jd)
    jd_ini=jd_today-180
    
    plt.figure(figsize=(12,7))
    gs = gridspec.GridSpec(3,1)
    gs.update(left=0.08, right=0.95, bottom=0.08, top=0.93, wspace=0.2, hspace=0.1)
    
    ax = plt.subplot(gs[0,0])
    ax.set_ylabel(r'$\Delta y$ (pix) - Orden 10')
    ax.get_xaxis().set_ticks([])
    ax.set_ylim([-1,1])
    ax.set_xlim([0,180])
    arr = desv10
    plt.errorbar(jd-jd_ini,desv10,yerr=0,fmt='o',c='red')
    for year in range(10):
    	jdyear = gcal2jd(2011+year,1,1)
    	plt.axvline(jdyear[0]+jdyear[1]-jd_ini, ls='--', c='black')
    	begin = jdyear[0]+jdyear[1]-jd_ini
    	ax.annotate(np.str(2011+year), xy=(begin+150, 890), xycoords='data', fontsize=14)
     
    
    ax = plt.subplot(gs[1,0])
    ax.set_ylabel(r'$\Delta y$ (pix) - Orden 40')
    label=r'JD-'+str(jd_ini)+' (days)'
    ax.set_xlabel(label)
    ax.set_xlim([0,180])
    ax.set_ylim([-1,1])
    arr = desv40
    plt.errorbar(jd-jd_ini,desv40,yerr=0,fmt='o',c='red')
    for year in range(10):
    	jdyear = gcal2jd(2011+year,1,1)
    	plt.axvline(jdyear[0]+jdyear[1]-jd_ini, ls='--', c='black')
    	begin = jdyear[0]+jdyear[1]-jd_ini
    	ax.annotate(np.str(2011+year), xy=(begin+150, 890), xycoords='data', fontsize=14)
    
    ax = plt.subplot(gs[2,0])
    ax.set_ylabel(r'$\Delta y$ (pix) - Orden 70')
    label=r'JD-'+str(jd_ini)+' (days)'
    ax.set_xlabel(label)
    ax.set_xlim([0,180])
    ax.set_ylim([-1,1])
    plt.errorbar(jd-jd_ini,desv70,yerr=0,fmt='o',c='red')
    for year in range(10):
    	jdyear = gcal2jd(2011+year,1,1)
    	plt.axvline(jdyear[0]+jdyear[1]-jd_ini, ls='--', c='black')
    	begin = jdyear[0]+jdyear[1]-jd_ini
    	ax.annotate(np.str(2011+year), xy=(begin+150, 890), xycoords='data', fontsize=14)


"""  
i=getMatrizDatos("./flat_160106_evening.fits")
posiciones=getConfiguracion("./ordenes_input.txt")
matPos,matSigma,matUmbral,matPosX = generarAjuste("./flat_160106_evening.fits","./ordenes_input.txt")
escribirMatriz(matPos,"posicionesY.txt")
y0=150
XX = np.arange(0,2048)
x = XX[y0-9:y0+9]
#sumo las columnas.
YY = np.sum(i[1024:1029,:],axis=0)
y = YY[y0-9:y0+9]
p0 = [np.max(y)-y[0],y0,2.,y[0]]
coeff,pcov = curve_fit(gaus,x,y,p0)
print coeff
plt.plot(x,y)
ymodel = gaus(x,coeff[0],coeff[1],coeff[2],coeff[3])
plt.plot(x,ymodel)


print "%.2f"%(coeff[1])

m=[]
m.insert(0,[2,3,4,5])
m.insert(0,[6,6,7,8])
print "%.2f"%(m[0][0])
"""
