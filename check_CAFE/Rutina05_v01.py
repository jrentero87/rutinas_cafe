# -*- coding: utf-8 -*-
"""
@author: Jesús Rentero Bonilla
Rutina 5: Tiempos de observación. Aprovechamiento de la noche.
Objetivo: Obtener la eficiencia de la noche de observación. 
          Para ello se mediran todos los tiempos de exposición de cada imagen 
          en toda la noche y se dividirá entre la duración de la noche, para calcular así la eficiencia.
          También obtendrá el tiempo dedicado a ficheros arco y el tiempo usado para ficheros de ciencia.
"""


import ephem
import os.path
from os import listdir
from astropy.io import fits
import numpy as np
import astropy.time
from dateutil import parser

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
Función que se encarga de lanzar la rutina y generar las estadísticas a partir del directorio
que contiene todos los ficheros de observación de una noche
"""
def runRutina05(directorio):
    observatorio=ephem.Observer()
    #Obtenemos la fecha de observacion a partir del directorio    
    anio="20"+directorio[0:2]
    mes=directorio[2:4]
    dia=directorio[4:6]
    # Fijamos la fecha a las 12:00 para que se calcule correctamente el próximo ocaso y crepúsculo
    fecha=anio+"/"+mes+"/"+dia+" 12:00"
    # Definimos posicion del telescopio y fecha
    observatorio.lat='37.2300'
    observatorio.lon='357.4537'
    observatorio.date=fecha
    # Definimos astronomical twilight
    observatorio.horizon='-18'
    sol=ephem.Sun()
    # Hallamos el twilight para asegurarnos que todos los ficheros se han realizado dentro de ese tiempo
    inicio=observatorio.next_setting(sol, use_center=True)
    fin=observatorio.next_rising(sol, use_center=True)

    #Hallamos el dia juliano para el twilight
    inicioTw = ephem.julian_date(inicio)
    finTw = ephem.julian_date(fin)
    
    #Inicializamos los tiempos
    tiempoArco=0.0
    tiempoTotal=0.0   
    
    #Calculamos la duracion de la noche en horas, minutos y segundos
    inicioDT=inicio.datetime()
    finDT=fin.datetime()
    duracionNoche=finDT-inicioDT
    segundosNoche=np.float(duracionNoche.seconds)
    minutosNoche=np.float(segundosNoche)/60.0
    horasNoche=np.float(segundosNoche)/3600.0
    
    #Calculamos el numero de ficheros arco, flats y bias
    numArcos=0
    numFlats=0
    numBias=0
    
    # Recorremos el directorio y obtenemos los tiempos
    for fichero in listdir(directorio):
        if os.path.isfile(directorio+"/"+fichero) and fichero.endswith(".fits"):
            rutaFich=directorio+"/"+fichero
            # Abrimos el fichero
            f=fits.open(rutaFich)
            # Obtenemos el tiempo de exposicion del fichero y la fecha
            tiempo=f[0].header["EXPTIME"]
            tExposicion=np.float(tiempo)
            fechaJul=getDiaJuliano(rutaFich)
            # Obtenemos el tipo de fichero que estamos tratando
            objeto=f[0].header["OBJECT"]
            if objeto.startswith("[arc]"):
                numArcos=numArcos+1
            if objeto.startswith("[flat]"):
                numFlats=numFlats+1
            if objeto.startswith("[Bias]"):
                numBias=numBias+1
            # Comprobamos que la fecha del fichero este dentro de los limites del twilight
            if inicioTw<fechaJul and fechaJul<finTw:
                # Comprobamos si es de tipo arco y si es así sumamos su tiempo de exposicion
                if objeto.startswith("[arc]"):
                    tiempoArco=tiempoArco+tExposicion
                # Sumamos el tiempo total de exposicion de todos los ficheros de la noche
                tiempoTotal=tiempoTotal+tExposicion
    
    eficiencia=(tiempoTotal/segundosNoche)*100.0
    tiempoCiencia=tiempoTotal-tiempoArco
    print "Numero de ficheros arco: %d"%(numArcos)
    print "Numero de ficheros flat: %d"%(numFlats)
    print "Numero de ficheros BIAS: %d"%(numBias)
    print "Tiempo total de exposicion: %.2f horas"%(tiempoTotal/3600.0)
    print "Tiempo total para ficheros ARCO: %.2f horas"%(tiempoArco/3600.0)
    print "Tiempo total para ciencia: %.2f horas"%(tiempoCiencia/3600.0)
    print "EL APROVECHAMIENTO HA SIDO DEL: %.2f por ciento"%(eficiencia)
    
    # Almacenamos los resultados en un fichero: eficiencia_fecha.txt
    # Abrimos el fichero donde escribiremos los resultados
    outfile = open("./Rut05_dat/eficiencia_"+directorio+".txt","w")
    outfile.write("Numero de ficheros arco: "+str(numArcos)+"\n")
    outfile.write("Numero de ficheros flat: "+str(numFlats)+"\n")
    outfile.write("Numero de ficheros BIAS: "+str(numBias)+"\n")
    outfile.write("Tiempo total de exposicion: "+str(tiempoTotal/3600.0)+" horas\n")
    outfile.write("Tiempo total para ficheros ARCO: "+str(tiempoArco/3600.0)+" horas\n")
    outfile.write("Tiempo total para ciencia: "+str(tiempoCiencia/3600.0)+" horas\n")
    outfile.write("EL APROVECHAMIENTO HA SIDO DEL: "+str(eficiencia)+" por ciento\n")
    outfile.close()
    

