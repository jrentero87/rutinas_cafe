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
from astropy.io import ascii
import datetime
from jdcal import gcal2jd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec # GRIDSPEC !

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
        
"""
Funcion encargada de llevar a cabo la ejecucion de la rutina 4
"""
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
        file.write(str(np.int(juldate))+","+str(round(mediana_total,4))+","+str(round(media_total,4))+","+str(round(desvTipica_total,4))+"\n")
    file.close()
    # Realizamos el checkeo de valores umbrales. 
    # Si el bias medio está entre 810 y 830 es correcto, y si el ruido de lectura es menor que 6 será también correcto.
    if media_total >= 810 and media_total <=830:
        print "... Nivel BIAS medio: %.2f ADUs ... OK"%(media_total)
    else:
        print "... Nivel BIAS medio: %.2f ADUs ... NO OK! - CHECK"%(media_total)
    if desvTipica_total < 6:
        print "... Ruido de lectura medio: %.2f ADUs ... OK"%(desvTipica_total)
    else:
        print "... Ruido de lectura medio: %.2f ADUs ... NO OK! - CHECK"%(desvTipica_total)
        
"""
Funcion encargada de añadir pintar y añadir al historial los resultados obtenidos en la noche que se esta ejecutando
"""
def plotHistory():
    colnames = ('jd','bias','noise','std')
    table = ascii.read('Rut04_dat/bias_master.txt', format='csv', names=colnames, comment='@')	
    jd  = np.array(table["jd"])
    bias   = np.array(table["bias"])
    std   = np.array(table["std"])
    
    today = datetime.datetime.now()
    today = astropy.time.Time(today)
    jd_today = np.int(today.jd)
    jd_ini=jd_today-180
    
    plt.figure(figsize=(12,7))
    gs = gridspec.GridSpec(2,1)
    gs.update(left=0.08, right=0.95, bottom=0.08, top=0.93, wspace=0.2, hspace=0.1)
    
    ax = plt.subplot(gs[0,0])
    ax.set_ylabel(r'Bias (ADUs)')
    ax.get_xaxis().set_ticks([])
    ax.set_ylim([800,900])
    ax.set_xlim([0,180])
    arr = bias
    plt.errorbar(jd-jd_ini,bias,yerr=0,fmt='o',c='red')
    for year in range(10):
    	jdyear = gcal2jd(2011+year,1,1)
    	plt.axvline(jdyear[0]+jdyear[1]-jd_ini, ls=':', c='gray')
    	begin = jdyear[0]+jdyear[1]-jd_ini
    	ax.annotate(np.str(2011+year), xy=(begin+150, 890), xycoords='data', fontsize=14)
    plt.grid(ls=':',c='gray')
    plt.axhline(810,ls='--',c='red')
    plt.axhline(830,ls='--',c='red')
    
    ax = plt.subplot(gs[1,0])
    ax.set_ylabel(r'Ruido de lectura (ADUs)')
    label=r'JD-'+str(jd_ini)+' (days)'
    ax.set_xlabel(label)
    ax.set_xlim([0,180])
    ax.set_ylim([2,7])
    
    for year in range(10):
    	jdyear = gcal2jd(2011+year,1,1)
    	plt.axvline(jdyear[0]+jdyear[1]-jd_ini, ls=':', c='gray')
    	begin = jdyear[0]+jdyear[1]-jd_ini
    	ax.annotate(np.str(2011+year), xy=(begin+150, 890), xycoords='data', fontsize=14)
    plt.grid(ls=':',c='gray')
    plt.axhline(6,ls='--',c='red')
    
    arr = std
    
    plt.scatter(jd-jd_ini,std,c=arr, cmap='winter',vmin=3.5, vmax=6)
    plt.savefig('bias_history_CAFE.pdf')