# -*- coding: utf-8 -*-
"""
@author: Jesús Rentero Bonilla
Rutina 1: ARC SPOTS
Objetivo: Medir las posiciones X-Y (en píxeles) de 200 spots en las imágenes
          arco crudas obtenidas durante la noche.
"""

from astropy.io import fits
from scipy import ndimage
import numpy as np
import os.path
import astropy.time
from astropy.io import ascii
import matplotlib.pyplot as plt
from dateutil import parser
from lmfit import  Model
from astroML.stats import sigmaG
from os import listdir
import matplotlib.gridspec as gridspec # GRIDSPEC !
import datetime
from jdcal import gcal2jd

# Para instalar ephem: pip install lmfit

"""
Constante para almacenar la distancia tomada desde el dentro del spot
hasta los extremos de la ventana. Por lo tanto, la ventana que tomaremos
para la muestra será de tamaño 2*TAM_VENTANA x 2*TAM_VENTANA
"""
TAM_VENTANA=10

"""
Constante para almacenar el nombre del fichero de spots que tomaremos como referencia
donde se almacena el centro del spot y la intensidad del mismo. A partir de este fichero
se generarán las estadísticas para los demás ficheros.
"""
INPUT_SPOT="input_spot.txt"

"""
Constante donde se almacena el nombre del fichero Master para la rutina 01. En él se almacenarán las desviaciones medias de cada noche.
"""
FICH_MASTER="./Rut01_dat/desviaciones_master.txt"

"""
Funcion que obtiene la matriz de datos a partir de una imagen de arco.
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
Función que obtiene el promedio de todos los elementos que contiene una matriz
que se le pasa por parámetro.
"""
def getPromedio(matriz):
    suma=0
    elementos=0
    # Recorremos la matriz y sumamos todos los elementos
    for i in range(0,len(matriz)):
        for j in range(0,len(matriz[i])):
            suma=suma+int(matriz[i][j])
            elementos=elementos+1
    # Calculamos el promedio con decimales
    promedio=float(suma)/float(elementos)
    return promedio

"""
Funcion que obtiene la submatriz con centro en las coordenadas (X,Y) y con
una anchura y altura 2*TAM_VENTANA
"""
def getSubMatriz(posX, posY, matriz):
    coordX=int(posX)
    coordY=int(posY)
    subMatriz=matriz[coordX-TAM_VENTANA:coordX+TAM_VENTANA,coordY-TAM_VENTANA:coordY+TAM_VENTANA] 
    return subMatriz
    
"""
Funcion que obtiene las coordenadas del centro del spot.
Para ello se recibe por parametros las coordenadas del spot y la matriz con los datos de la imagen
A continuacion se genera una submatriz con un tamaño de ventana 2*TAM_VENTANAx2*TAM_VENTANA
Se calcula el centro para dicha matriz con la funcion center_of_mass
Finalmente, se realiza un cambio de coordenadas de la submatriz a la matriz principal
La función devuelve el centro del spot.
"""
def getCentroSpot(coordX,coordY,matriz):
    #Obtenemos la submatriz
    subM=getSubMatriz(coordX,coordY,matriz)
    #Obtenemos el centro del spot
    centroSub=ndimage.measurements.center_of_mass(subM)
    #Realizamos un cambio de coordenadas para obtener el centro en la matriz de datos
    desvX=centroSub[0]-TAM_VENTANA
    centroX=coordX+desvX
    desvY=centroSub[1]-TAM_VENTANA
    centroY=coordY+desvY
    centro=[centroX,centroY]
    return centro

"""
Funcion que obtiene el centro del spot a partir de la esquina superior izquierda de la ventana
VERSION 0. A partir del centro de masas del spot

def getCentroVentana(venX, venY, matriz):
    #Obtenemos la submatriz
    subM=matriz[venX:venX+TAM_VENTANA*2,venY:venY+TAM_VENTANA*2]
    #Obtenemos el centro del spot contenido en la ventana
    centroSub=ndimage.measurements.center_of_mass(subM)
    #Realizamos un cambio de coordenadas para obtener el centro en la matriz de datos
    desvX=centroSub[0]-TAM_VENTANA
    centroX=venX+TAM_VENTANA+desvX
    desvY=centroSub[1]-TAM_VENTANA
    centroY=venY+TAM_VENTANA+desvY
    centro=[centroX,centroY]
    return centro
"""

"""
Funcion Gaussiana para ajustar los spots en una dimension
"""
def gaussian(x, amp, cen, wid, level):
    return (amp/(np.sqrt(2*np.pi)*wid)) * np.exp(-(x-cen)**2 /(2*wid**2))  + level

"""
Funcion que obtiene el centro del spot a partir de la esquina superior izquierda de la ventana
Se realizan dos ajustes de 1-dimensión. Para la coordenada X del centro del spot colapsaremos
la matriz en vertical y ajustaremos a la Gaussiana, cuyo centro será la coordenada X. 
Para la coordenada Y haremos lo mismo solo que en sentido horizontal.
"""
def getCentroVentana(venX, venY, matriz):
    #Obtenemos la submatriz
    subM=matriz[venX:venX+TAM_VENTANA*2,venY:venY+TAM_VENTANA*2]
    #Obtenemos el centro del spot contenido en la ventana en la dirección X
    x = np.arange(TAM_VENTANA*2)
    y = np.sum(subM,axis=1)
    gmod = Model(gaussian)
    resultX = gmod.fit(y, x=x, amp=np.max(y)-np.median(y), cen=TAM_VENTANA, wid=2, level=np.median(y))
    #Obtenemos el centro del spot contenido en la ventana en la dirección Y
    x2 = np.arange(TAM_VENTANA*2)
    y2 = np.sum(subM,axis=0)
    gmod2 = Model(gaussian)
    resultY = gmod2.fit(y2, x=x2, amp=np.max(y2)-np.median(y2), cen=TAM_VENTANA, wid=2, level=np.median(y2))
    #Centro del spot {X,Y}
    centroSub = np.array([resultX.best_values["cen"], resultY.best_values["cen"]])        
    #Realizamos un cambio de coordenadas para obtener el centro en la matriz de datos
    desvX=centroSub[0]-TAM_VENTANA
    centroX=venX+TAM_VENTANA+desvX
    desvY=centroSub[1]-TAM_VENTANA
    centroY=venY+TAM_VENTANA+desvY
    centro=[centroX,centroY]
    return centro

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
Funcion que genera el fichero input_spot.txt. Este fichero contendrá la siguiente 
información: Id_spot##coordX#coordY#valor_medio_spot, donde:
- IdSpot = identificador del spot
- posVenX = coordenada X de la esquina superior izquierda de la ventana que contiene al spot
- posVenY = coordenada Y de la esquina superior izquierda de la ventana que contiene al spot
- posX = coordenada X del centro del spot estudiado
- posY = coordenada Y del centro del spot estudiado
- Intensidad= Valor promedio de la submatriz que contiene los datos del spot
Para ello se recibe como parametro:
- ficheroSpot: contiene la información obtenida con el ds9 de una imagen inicial
- matriz: continee la información de los datos de la imagen
"""
def generarInputSpot(ficheroSpot,matriz):
    #Abrimos el fichero que contiene las coordenadas iniciales de los spots
    infile = open(ficheroSpot,'r')
    #Abrimos un fichero de escritura (input_spot.txt) donde escribimos los resultados
    outfile = open(INPUT_SPOT,"w")
    #Escribimos la primera linea del fichero con los comentarios
    outfile.write("@IdSpot,posVenX,posVenY,posX,posY,Intensidad\n")
    #Recorremos este fichero obteniendo cada una de las lineas
    for line in infile:
        #Troceamos la linea, almacenando en spot[0] el id, spot[1] PosX, spot[2] PosY
        spot=line.split(",")
        idSpot=spot[0]
        #Comprobamos que la linea no sea un comentario, es decir, que no comience por @
        if idSpot[0] != "@":
            posXSpot=int(spot[1])
            posYSpot=int(spot[2])
            #Obtenemos el centro del spot para calcular la ventana a partir de el centro
            centro=getCentroSpot(posXSpot,posYSpot,matriz)
            #Calculamos las coordenadas X e Y de la esquina superior izquierda de la ventana
            venX=int(centro[0]-TAM_VENTANA)
            venY=int(centro[1]-TAM_VENTANA)
            #Recalculamos el centro del spot contenida en la ventana a partir de sus coordenadas superior izquierda
            centroVen=getCentroVentana(venX,venY,matriz)
            #Obtenemos la submatriz para calcular la intensidad del spot realizando la suma
            subM=matriz[venX:venX+TAM_VENTANA*2,venY:venY+TAM_VENTANA*2]
            intensidad=np.sum(subM)
            #Tomamos una precisión de 4 decimales para el calculo del centro
            cenX=round(centroVen[0],4)
            cenY=round(centroVen[1],4)
            #Escribimos los datos en el fichero
            outfile.write(idSpot+","+str(venX)+","+str(venY)+","+str(cenX)+","+str(cenY)+","+str(intensidad)+"\n")
    #Cerramos ambos ficheros    
    infile.close()
    outfile.close() 
    
"""
Función que se encarga de generar un fichero de salida con nombre "arcoFits_dat.spot"
con la siguiente información en cada linea del fichero:
- IdSpot = identificador del spot de la imagen arcoFits
- posX = coordenada X del centro del spot
- posY = coordenada Y del centro del spot
- distX = distancia en el eje X del spot de la imagen arcoFits con respecto al fichero inputSpot que tenemos como referencia
- distY = distancia en el eje Y del spot de la imagen arcoFits con respecto al fichero inputSpot que tenemos como referencia
- Intensidad = Valor promedio de la submatriz que contiene los datos del spot
Para ello se recibe como parametro:
- inputSpot = fichero de muestra con el que se van a realizar las comparaciones
- arcoFits = imagen de arco a analizar.
"""
def generarEstadisticas(inputSpots, arcoFits):
    # Obtenemos la matriz de datos de la imagen arcoFits
    matrizDat=getMatrizDatos(arcoFits)
    # Abrimos el fichero con la información calculada previamente de los spots (posicion e intensidad)
    infile = open(inputSpots,'r')
    # Creamos el fichero de estadisticas
    nomFichero = arcoFits[0:len(arcoFits)-5]+"_"+arcoFits[0:6]+".spot"
    outfile = open("./Rut01_dat/"+nomFichero[nomFichero.index('/')+1:],"w")
    # Escribimos la primera linea del fichero con los comentarios
    outfile.write("@IdSpot,posX,posY,distX,distY,Intensidad,diaJuliano\n")
    # Obtenemos el dia juliano en el que se ha realizado la imagen arcoFits
    diaJuliano=getDiaJuliano(arcoFits)
    # Recorremos el fichero donde tenemos las coordenadas de los spots para generar las estadisticas
    for line in infile:
         #Troceamos la linea, almacenando en spot[0] el id, spot[1] venX, spot[2] venY, spot[3] posX, spot[4] posY
        spot=line.split(",")
        idSpot=spot[0]
        #Comprobamos que la linea no sea un comentario, es decir, que no comience por @
        if idSpot[0] != "@":
            #Obtenemos la posicion de la ventana y los centros de referencia
            venX=int(spot[1])
            venY=int(spot[2])
            posXSpot=float(spot[3])
            posYSpot=float(spot[4])
            #Obtenemos el centro del spot de la imagen que a analizar
            centro=getCentroVentana(venX,venY,matrizDat)
            #Calculamos las distancias de los respectivos centros
            distX=(posXSpot-centro[0])
            distY=(posYSpot-centro[1])
            #Obtenemos la submatriz para calcular la intensidad del spot
            subM=matrizDat[venX:venX+TAM_VENTANA*2,venY:venY+TAM_VENTANA*2]
            intensidad=np.sum(subM)
            #Tomamos una precisión de 4 decimales para el calculo del centro y las distancias
            cenX=round(centro[0],4)
            cenY=round(centro[1],4)
            distanciaX=round(distX,4)
            distanciaY=round(distY,4)
            diaJul=round(diaJuliano,6)
            #Escribimos los datos en el fichero
            outfile.write(idSpot+","+str(cenX)+","+str(cenY)+","+str(distanciaX)+","+str(distanciaY)+","+str(intensidad)+","+str(diaJul)+"\n")
    #Cerramos ambos ficheros    
    infile.close()
    outfile.close() 


"""
Funcion que se encarga de generar todas las estadísticas para el listado de imágenes de arco
que se encuentren en el fichero que se le introduce por parámetro.
IMPORTANTE: Estos ficheros de ARCO deberán estar en el mismo directorio que la rutina y 
que el fichero con dicho listado.
"""
def rutina01Run(listaArcos):
    # Abrimos el fichero con el listado de ficheros arco
    infile = open(listaArcos,'r')
    # Procesamos cada una de las lineas del fichero, y generamos las estadísticas para cada fichero de arco
    for line in infile:
        #Eliminamos de la linea el retorno de carro (\n)
        line=line.strip()
        #Comprobamos que la linea tenga información y no sea una linea en blanco
        if len(line)>0:
            generarEstadisticas(INPUT_SPOT,line)

"""
Función que realiza el promedio de las desviaciones de todos los spots para un fichero en concreto, y el promedio de las intensidades
"""
def getPromedioDesv(fichero):
    #Abrimos el fichero
    infile = open("./Rut01_dat/"+fichero[fichero.index('/')+1:],'r')
    # Inicializamos contadores
    sumaX=0.0
    sumaY=0.0
    sumaInt=0.0
    # Recorremos el fichero donde tenemos la informacion de los spots
    for line in infile:
         #Troceamos la linea, almacenando en spot[0] el id, spot[1] PosX, spot[2] PosY, spot[3] distX, spot[4] distY, spot[5] intensidad
        spot=line.split(",")
        idSpot=spot[0]
        #Comprobamos que la linea no sea un comentario, es decir, que no comience por @
        if idSpot[0] != "@":
            sumaX=sumaX+float(spot[3])
            sumaY=sumaY+float(spot[4])
            sumaInt=sumaInt+float(spot[5])
    infile.close()
    promedioX=sumaX/200
    promedioY=sumaY/200
    promedioInt=sumaInt/200
    return [promedioX,promedioY,promedioInt]
    #print "Fichero: %s. Anchura ventana: %d px. Desviacion en X: %.4f px. Desviacion en Y: %.4f px" %(fichero,TAM_VENTANA*2, promedioX, promedioY)
    
"""
Funcion que obtiene un array con las intensidades del fichero ARCO de referencia
"""
def getIntensidadReferencia():
    #Abrimos el fichero
    infile = open(INPUT_SPOT,'r')
    # Inicializamos contadores
    intensidadesRef=[]
    # Recorremos el fichero donde tenemos la informacion de los spots
    for line in infile:
         #Troceamos la linea, almacenando en spot[0] el id, spot[1] PosVenX, spot[2] PosVenY, spot[3] posX, spot[4] posY, spot[5] intensidad
        spot=line.split(",")
        idSpot=spot[0]
        #Comprobamos que la linea no sea un comentario, es decir, que no comience por @
        if idSpot[0] != "@":
            intensidadesRef.append(float(spot[5]))
    infile.close()
    return intensidadesRef

"""
Función que muestra el promedio de las desviaciones de todos los spots estudiados y de las intensidades
Para ello, hace uso del fichero que contiene el litado de imágenes de arco y accede
a los ficheros "arco_dat.txt" para mostrar el promedio de las desviaciones de los spots.
Es necesario haber ejecutado la rutina para poder ejecutar esta función, es decir, deben existir los ficheros
"""
def promedioDistancias(listaArcos):
    #Inicializamos vectores donde se almacenara las desviaciones medias de cada fichero asi como las intensidades
    desvX=[]
    desvY=[]
    intensidad=[]
    # Abrimos el fichero con el listado de ficheros arco
    infile = open(listaArcos,'r')
    # Procesamos cada una de las lineas del fichero, y generamos las estadísticas para cada fichero de arco
    for line in infile:
        #Comprobamos que la linea tenga información y no sea una linea en blanco
        if len(line)>0:
            #Eliminamos de la linea el retorno de carro (\n), eliminamos la extensión .fits y añadimos "_fecha.spot"
            line=line.strip()
            line=line[0:len(line)-5]+"_"+line[0:6]+".spot"
            fichDat=getPromedioDesv(line)
            desvX.append(fichDat[0])
            desvY.append(fichDat[1])
            intensidad.append(fichDat[2])
    #Calculamos las intensidades normalizadas de la noche
    intRef=getIntensidadReferencia()
    intNorm=np.array(intensidad[:])/np.mean(intRef)
    infile.close()
    return [np.mean(desvX), np.mean(desvY), np.mean(intNorm)]      

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
Función que se encarga de genera el fichero Master de la rutina y de chequear los datos
Se le pasa por parametro la lista de ficheros arco.
"""
def checkRutina01(listaArcos):
    #Obtenemos el promedio de las desviaciones en X y en Y de los spots y el promedio de las intensidades normalizadas
    [desvX, desvY, intNorm]=promedioDistancias(listaArcos)
    # Añadimos el valor de la media de todas las desviaciones de los spots de la noche al fichero master
    # Si existe el fichero lo abrimos en modo "a", sino lo creamos
    if os.path.exists(FICH_MASTER):
        file=open(FICH_MASTER,"a")
    else:
        file=open(FICH_MASTER,"w")
        file.write("@juldate,desvX_media,desvY_media,intensidad_Norm\n")
    # Abrimos el fichero con el listado de ficheros arco
    infile = open(listaArcos,'r')
    # Obtenemos el dia juliano para uno de los ficheros arco de la noche
    imagen = infile.readline().strip()
    infile.close()
    juldate=getDiaJuliano(imagen)
    #Comprobamos que la entrada en el fichero no exista. En caso de no existir escribimos nueva entrada
    if not existeNoche(np.int(juldate)):
        file.write(str(np.int(juldate))+","+str(round(desvX,4))+","+str(round(desvY,4))+","+str(intNorm)+"\n")
    file.close()
    
    #Realizamos el checkeo para la rutina01
    # Si las desviciones medias de los spots son menores a 20 milipíxeles y la intensidad normalizada esta entre el 0.99% y el 1.01%
    if desvX < 0.02 and desvX>-0.02:
        print "... Desviación media en eje X: %.2f ... OK"%(desvX)
    else:
        print "... Desviación media en eje X: %.2f ... NO OK! - CHECK"%(desvX)
    if desvY < 0.02 and desvY>-0.02:
        print "... Desviación media en eje Y: %.2f ... OK"%(desvY)
    else:
        print "... Desviación media en eje Y: %.2f ... NO OK! - CHECK"%(desvY)
    if intNorm > 0.99 and intNorm < 1.01:
        print "... Ruido de lectura medio: %.2f ... OK"%(intNorm)
    else:
        print "... Ruido de lectura medio: %.2f ... NO OK! - CHECK"%(intNorm)
        
"""
Funcion encargada de añadir pintar y añadir al historial los resultados obtenidos en la noche que se esta ejecutando
"""
def plotHistory():
    colnames = ('jd','desvX','desvY','intNorm')
    table = ascii.read('Rut01_dat/desviaciones_master.txt', format='csv', names=colnames, comment='@')	
    jd  = np.array(table["jd"])
    desvX   = np.array(table["desvX"])
    desvY   = np.array(table["desvY"])
    intNorm   = np.array(table["intNorm"])
    
    today = datetime.datetime.now()
    today = astropy.time.Time(today)
    jd_today = np.int(today.jd)
    jd_ini=jd_today-180
    
    plt.figure(figsize=(12,7))
    gs = gridspec.GridSpec(3,1)
    gs.update(left=0.08, right=0.95, bottom=0.08, top=0.93, wspace=0.2, hspace=0.1)
    
    ax = plt.subplot(gs[0,0])
    ax.set_ylabel(r'$\Delta x$ (pix)')
    ax.get_xaxis().set_ticks([])
    ax.set_ylim([-4,4])
    ax.set_xlim([0,180])
    arr = desvX
    plt.errorbar(jd-jd_ini,desvX,yerr=0,fmt='o',c='red')
    for year in range(10):
    	jdyear = gcal2jd(2011+year,1,1)
    	plt.axvline(jdyear[0]+jdyear[1]-jd_ini, ls='--', c='black')
    	begin = jdyear[0]+jdyear[1]-jd_ini
    	ax.annotate(np.str(2011+year), xy=(begin+150, 890), xycoords='data', fontsize=14)
     
    
    ax = plt.subplot(gs[1,0])
    ax.set_ylabel(r'$\Delta y$ (pix)')
    label=r'JD-'+str(jd_ini)+' (days)'
    ax.set_xlabel(label)
    ax.set_xlim([0,180])
    ax.set_ylim([-4,4])
    arr = desvY
    plt.errorbar(jd-jd_ini,desvY,yerr=0,fmt='o',c='red')
    for year in range(10):
    	jdyear = gcal2jd(2011+year,1,1)
    	plt.axvline(jdyear[0]+jdyear[1]-jd_ini, ls='--', c='black')
    	begin = jdyear[0]+jdyear[1]-jd_ini
    	ax.annotate(np.str(2011+year), xy=(begin+150, 890), xycoords='data', fontsize=14)
    
    ax = plt.subplot(gs[2,0])
    ax.set_ylabel('Norm. Intensity')
    label=r'JD-'+str(jd_ini)+' (days)'
    ax.set_xlabel(label)
    ax.set_xlim([0,180])
    ax.set_ylim([0.5,1.5])
    
    for year in range(10):
    	jdyear = gcal2jd(2011+year,1,1)
    	plt.axvline(jdyear[0]+jdyear[1]-jd_ini, ls='--', c='black')
    	begin = jdyear[0]+jdyear[1]-jd_ini
    	ax.annotate(np.str(2011+year), xy=(begin+150, 890), xycoords='data', fontsize=14)
    
    arr = intNorm
    
    plt.scatter(jd-jd_ini,intNorm,c=arr, cmap='winter',vmin=3.5, vmax=6)
    plt.savefig('spots_history_CAFE.pdf')


        
"""
Plot de los resultados de la noche que se esta ejecutando
"""
def Plot1night(night):
	#Directorio de trabajo y numero de ficheros
    DIR = './Rut01_dat'
    nfiles= len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])
	#Inicializamos algunas variables
    XX = np.zeros((nfiles, 199))
    YY = np.zeros((nfiles, 199))
    dX = np.zeros((nfiles, 199))
    dY = np.zeros((nfiles, 199))
    IN = np.zeros((nfiles, 199))
    JD = np.zeros((nfiles, 199))
	#Leemos todos los ficheros creados para cada arco y los almacenamos en las matrices correspondientes
    ii = 0
    for file in listdir('./Rut01_dat'):
        colnames = ('IdSpot','posX','posY','distX','distY','Intensidad','jd')
        table = ascii.read('Rut01_dat/'+file, format='csv', names=colnames, comment='@')
        jda  = np.array(table["jd"])
        x   = np.array(table["posX"])
        y   = np.array(table["posY"])
        dx  = np.array(table["distX"])
        dy  = np.array(table["distY"])
        I	= np.array(table["Intensidad"])
        XX[ii,:] = x
        YY[ii,:] = y
        dX[ii,:] = dx
        dY[ii,:] = dy
        IN[ii,:] = I
        JD[ii,:] = jda
        ii = ii+1
    
	# Calculamos los offsets de cada spot respecto a la mediana de ese spot en todos los arcos de la noche
    nXX = np.zeros((nfiles, 199))
    nYY = np.zeros((nfiles, 199))
    nIN = np.zeros((nfiles, 199))
    for i in range(199):
        nXX[:,i] = (XX[:,i]-np.median(XX[:,i]))#*0.037517/5500. * 299792458
        nYY[:,i] = (YY[:,i]-np.median(YY[:,i]))#*0.037517/5500. * 299792458
        nIN[:,i] = IN[:,i]/np.mean(IN[:,i])
    
    # Inicializamos el plot
    plt.figure(figsize=(12,7))
    gs = gridspec.GridSpec(3,1)
    gs.update(left=0.08, right=0.95, bottom=0.08, top=0.93, wspace=0.2, hspace=0.1)
    
    # Plot para los offsets relativos en la dirección X
    ax = plt.subplot(gs[0,0])
    ax.set_ylabel(r'$\Delta x$ (mpix)')
    ax.get_xaxis().set_ticks([])
    ax.set_ylim([-50,55])
    ax.set_xlim([5,20])
    for i in range(199):
        Xplot = (XX[:,i]-np.median(XX[:,i]))*1.e3
        plt.plot((JD[:,i]-2457594.)*24.,Xplot,'+',c='Silver',zorder=-1,alpha=0.6)
    for i in range(nfiles):
        Xplot = nXX[i,:]*1.e3
        plt.errorbar((JD[i,0]-2457594.)*24.,np.median(Xplot),yerr=sigmaG(Xplot),fmt='o',c='b',zorder=1)
    
    # Plot para los offsets relativos en la dirección Y
    ax = plt.subplot(gs[1,0])
    ax.set_ylabel(r'$\Delta y$ (mpix)')
    ax.get_xaxis().set_ticks([])
    ax.set_ylim([-50,50])
    ax.set_xlim([5,20])
    for i in range(199):
        Xplot = (YY[:,i]-np.median(YY[:,i]))*1.e3
        plt.plot((JD[:,i]-2457594.)*24.,Xplot,'+',c='Silver',zorder=-1,alpha=0.6)
    for i in range(nfiles):
        Xplot = nYY[i,:]*1.e3
        plt.errorbar((JD[i,0]-2457594.)*24.,np.median(Xplot),yerr=sigmaG(Xplot),fmt='o',c='r',zorder=1)
    # Plot para la intensidad
    ax = plt.subplot(gs[2,0])
    ax.set_ylabel('Norm. Intensity')
    ax.set_xlabel('JD-2457594 (h)')
    ax.set_ylim([0.95,1.02])
    ax.set_xlim([5,20])
    for i in range(199):
        plt.plot((JD[:,i]-2457594.)*24.,nIN[:,i],'+',c='Silver',zorder=-1,alpha=0.6)
    for i in range(nfiles):
        plt.errorbar((JD[i,0]-2457594.)*24.,np.median(nIN[i,:]),yerr=sigmaG(nIN[i,:]),fmt='o',c='forestgreen',zorder=1)
    
    plt.savefig("./Rut01_dat/Rutina01_plot_1night_"+night[0:6]+".pdf") 

"""
tbdata=getMatrizDatos("./cali_0075.fits")
promedio=np.sum(tbdata)
print promedio
subM=getSubMatriz(811,1840,tbdata)
print subM

print tbdata[811][1840]
centro=getCentroSpot(811,1840,tbdata)
print centro

generarInputSpot("./spots.txt",tbdata)
rutina01Run("./listaArcos.txt")
promedioDistancias("./listaArcos.txt")
"""