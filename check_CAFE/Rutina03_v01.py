# -*- coding: utf-8 -*-
"""
@author: Jesús Rentero Bonilla
Rutina 3: Degradación del CCD
Objetivo: Determinar si el CCD está sufriendo algún tipo de degradación
         o pérdida de eficiencia.
         La perdida de eficiencia en el detector hace que para obtener una señal
         ruido determinada debamos integrar cada vez más tiempo. Por lo que mediremos
         la señal ruido para un mismo tiempo de integración y ver cómo evoluciona
         dicha señal ruido a lo largo del tiempo
         Esta rutina solo podrá ser ejecutada para espectros ya reducidos.
"""

#NOTA: Para instalar astroML: conda install -c astropy astroml=0.3

from astropy.io import fits
import numpy as np
import os.path
from astroML.stats import sigmaG
import astropy.time
from dateutil import parser
import os.path
from os import listdir

LOG_SNR="./Rut03_dat/log_snr.txt"

"""
Esta función devolverá true en caso de que exista en el fichero log_snr.txt 
una entrada para la noche que se introduce por parámetro, y false en caso contrario.
"""
def existeNoche(diaJuliano):
    # Abrimos el fichero en modo lectura
    infile=open(LOG_SNR,'r')
    # Creamos una variable que inicialmente inicializamos a falso
    existe=False
    # Recorremos el fichero y comparamos cada una de las lineas
    for line in infile:
        #Ignoramos las lineas que comiencen por @, puesto que se trata de un comentario en el fichero
        if line[0]!='@':
            # Obtenemos el dia juliano almacenado en la primera posicion de la linea
            linea=line.split(',')
            if linea[0]==diaJuliano:
                existe=True
    infile.close()
    return existe

"""
Función que se encarga de calcular la señal ruido para un fichero determinado.
Este fichero debe contener el espectro ya reducido. Almacena los resultados en un fichero
log_snr.txt
"""
def procesar(fichero):
    # Abrimos el fichero
    f=fits.open(fichero)
    # Obtenemos la matriz con los datos
    im = f[0].data
    
    # Calculamos la señal ruido del orden 42, correspondiente a 5500 Angstrom
    # Utilizaremos el rango de píxeles entre 1500 y 1700, puesto que en este rango solo hay continuo
    # Haremos el cociente de la mediana de los datos entre la desviación típica
    snr = np.median(im[42,1500:1700]) / sigmaG(im[42,1500:1700])
    
    # Obtenemos el tiempo de exposicion
    exptime = np.float(f[0].header["EXPTIME"])
    
    # Obtenemos el día juliano y lo pasamos a entero
    date=f[0].header["DATE"]
    dt = parser.parse(date)
    time = astropy.time.Time(dt)
    juldate = time.jd
    
    # Dividimos el tiempo de exposicion entre 10 para hallar la relación Señal-Ruido/Tiempo-exposicion
    time_exp10=exptime/10
    
    # Calculamos snr/time_exp10
    snr_time=snr/time_exp10
    
    # Obtenemos el nombre del objeto observado
    name = f[0].header["OBJECT"]
    
    # Escribimos en el fichero de registro
    if os.path.exists(LOG_SNR):
        file=open(LOG_SNR,"a")
    else:
        file=open(LOG_SNR,"w")
        file.write("@juldate,snr/exptime,object\n")
    # Comprobamos que ese mismo fichero no se haya procesado anteriormente, para ello comparamos con el dia juliano
    if not existeNoche(str(round(juldate,6))):
        file.write(str(round(juldate,6))+","+str(round(snr_time,4))+","+name+"\n")
    file.close()


"""
Esta función procesa todos los ficheros de un directorio, y ejecutar la función procesar
para aquellos ficheros que estén reducidos, es decir, cuya extensión es .disp_cor.fits
"""
def runRutina03(directorio):
    # Recorremos el directorio
    for fichero in listdir(directorio):
        # Comprobamos que exista el fichero y que se trata de un fichero reducido
        if os.path.isfile(directorio+"/"+fichero) and fichero.endswith(".disp_cor.fits"):
            rutaFich=directorio+"/"+fichero
            procesar(rutaFich)
    
    
    
runRutina03("./prueba")

"""
fichero="HD109358_120706_0002.disp_cor.fits"
rutina03run(fichero)
"""