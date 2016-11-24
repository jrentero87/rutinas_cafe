# -*- coding: utf-8 -*-
"""
@author: Jesús Rentero Bonilla
Rutina 4: Control del nivel de BIAS
Objetivo: Medir el nivel de BIAS de cada fichero. Se generará un único fichero
          llamado "nivel_bias_directorio.txt" donde se almacene la media, mediana y desv. típica
          del nivel de bias de cada una de las imágenes BIAS.
          Esta rutina añadirá una entrada en el fichero "bias_master.txt". Este fichero contendrá
          la mediana del nivel de bias de una noche específica, junto con el día juliano.
"""

from astropy.io import fits
import numpy as np
import astropy.time
from dateutil import parser
import os.path
from astroML.stats import sigmaG

"""
Definición de constantes:
- FICH_BIAS: fichero que almacena el listado de ficheros bias de la noche
- FICH_MASTER: fichero donde se almacenan las estadísticas de todas las noches.
"""
FICH_BIAS="biasFits.txt"
FICH_MASTER="./Rut04_dat/bias_master.txt"

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
        
        
    

def runRutina04(directorio):
     # Abrimos el fichero con el listado de ficheros bias
    infile = open(FICH_BIAS,'r')
    # Abrimos el fichero donde escribiremos los resultados
    outfile = open("./Rut04_dat/nivel_bias_"+directorio+".txt","w")
    outfile.write("@fichero, bias_medio, bias_mediana, bias_desvTipica, dia_juliano\n")
    # Variable para almacenar todos los valores de todos los bias de una noche
    biasNoche=[]
    # Procesamos cada una de las lineas del fichero
    for line in infile:
        #Eliminamos de la linea el retorno de carro (\n)
        line=line.strip()
        #Comprobamos que la linea tenga información y no sea una linea en blanco
        if len(line)>0:
            # Abrimos el fichero de bias
            hdulist=fits.open(line);
            #Obtenemos la matriz con los datos
            tbdata = hdulist[0].data
            #Obtenemos el dia juliano del bias
            date=hdulist[0].header["DATE"]
            dt = parser.parse(date)
            time = astropy.time.Time(dt)
            juldate = time.jd
            #cerramos el fichero
            hdulist.close();
            nombre=line[line.index("/")+1:]
            media=np.mean(tbdata)
            mediana=np.median(tbdata)
            desviacion=sigmaG(tbdata)
            biasNoche.append(tbdata)
            outfile.write(nombre+","+str(round(media,4))+","+str(mediana)+","+str(round(desviacion,4))+","+str(round(juldate,6))+"\n")     
    outfile.close()
    infile.close()
    
    # Añadimos el valor de la mediana de todos los bias de la noche al fichero master
    # Si existe el fichero lo abrimos en modo "a", sino lo creamos
    if os.path.exists(FICH_MASTER):
        file=open(FICH_MASTER,"a")
    else:
        file=open(FICH_MASTER,"w")
        file.write("@juldate,bias_mediana,bias_medio,bias_desvTipica\n")
    mediana_total=np.median(biasNoche)
    media_total=np.mean(biasNoche)
    desvTipica_total=sigmaG(biasNoche)
    #Comprobamos que la entrada en el fichero no exista. En caso de no existir escribimos nueva entrada
    if not existeNoche(np.int(juldate)):
        file.write(str(np.int(juldate))+","+str(round(mediana_total,4))+","+str(round(media_total,4))+str(round(desvTipica_total,4))+"\n")
    file.close()