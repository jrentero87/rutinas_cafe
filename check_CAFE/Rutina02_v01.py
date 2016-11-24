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
from scipy import asarray as exp
import matplotlib.pyplot as plt
import sys

"""
Fichero que almacena las posiciones de cada uno de los ordenes
"""
INPUT_ORDEN="./ordenes_input.txt"

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
