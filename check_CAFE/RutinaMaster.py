# -*- coding: utf-8 -*-
"""
@author: Jesús Rentero Bonilla
Rutina Master.
Objetivo: Ejecutar cada una de las rutinas de manera automática.
          Esto se realizará sobre el directorio que la rutina Master reciba
          por parámetros. La rutina master se encargará de ejecutar las
          siguientes rutinas:
          - Rutina 01: ARC-Spots. Medir las posiciones de los spots de una imagen arco.
          - Rutina 02: Posición e intensidad del flat.
                       Determinará la posición de las órdenes por columnas en el CCD.
          - Rutina 03: Degradación del CCD (No se ejecutará para espectros no reducidos)
          - Rutina 04: Nivel de BIAS.
          - Rutina 05: Eficiencia de la noche. (tiempo exposicion/tiempo empleado)
"""
# Para instalar ephem: pip install pyephem
import sys
from astropy.io import fits
import os.path
from os import system
from os import listdir
import Rutina01_v01
import Rutina02_v01
import Rutina04_v01
import Rutina05_v01

"""
Constantes para almacenar la ruta de los ficheros arco y flats que tomamos como referencia
"""
ARCO_REF="./cali_0061.fits"
FLAT_REF="./cali_0032.fits"

"""
Constantes donde almacenamos los nombres de los ficheros que contienen el listado
de ficheros arco,flat y bias 
"""
FICH_ARCO="arcoFits.txt"
FICH_FLAT="flatFits.txt"
FICH_BIAS="biasFits.txt"

"""
Funcion que se encarga de generar las listas de ficheros para arco, flats y bias
del directorio que se recibe por parámetro 
"""
def generarListaFicheros():
    # Obtenemos el directorio donde se encuentran todos los ficheros
    direct=sys.argv[1]
    # Definimos un directorio auxiliar de trabajo
    directAux=direct+'_aux'
    
    # Realizamos una copia del directorio con el que vamos a trabajar
    # En linux
    #system('cp -r '+direct+' '+directAux)
    directAux=direct
    
    #Creamos ficheros para arco, flat y bias
    arcoFits=open(FICH_ARCO,"w")
    flatFits=open(FICH_FLAT,"w")
    biasFits=open(FICH_BIAS,"w")
    
    # Recorremos el directorio
    for fichero in listdir(directAux):
        if os.path.isfile(directAux+"/"+fichero) and fichero.endswith(".fits"):
            rutaFich=directAux+"/"+fichero
            # Abrimos el fichero
            f=fits.open(rutaFich)
            # Obtenemos el tipo de fichero que estamos tratando
            objeto=f[0].header["OBJECT"]
            if objeto[0]=='[':
                tipo=objeto[:objeto.index(']')+1]
            else:
                tipo='[science]'
            #print "%s - %s"%(rutaFich,tipo)
            #Clasificamos los ficheros segun su tipo y creamos una lista de ficheros para cada tipo
            if tipo=='[arc]':
                arcoFits.write(rutaFich+"\n")
            elif tipo=='[flat]':
                flatFits.write(rutaFich+"\n")
            elif tipo=='[Bias]':
                biasFits.write(rutaFich+"\n")
                    
            f.close()
    
    arcoFits.close()
    flatFits.close()
    biasFits.close()


"""
Función que se encarga de lanzar las rutinas 1 y 2
"""
def run_Rutina01_Rutina02(directorio):
    #Arrancamos la rutina 01. 
    print "EJECUTANDO RUTINA 01: ARC-SPOTS ..."
    print "==================================="
    # Obtenemos la matriz de datos del fichero que cogemos como referencia
    tbdata=Rutina01_v01.getMatrizDatos(ARCO_REF)
    # Generamos el fichero input_spot.txt que utilizaremos para el estudio
    Rutina01_v01.generarInputSpot("./spots.txt",tbdata)
    # Lanzamos la rutina generando para cada fichero arco un fichero de datos con los resultados
    Rutina01_v01.rutina01Run(FICH_ARCO)
    #Rutina01_v01.promedioDistancias(FICH_ARCO)
    # Generamos el fichero Master de la primera rutina:
    Rutina01_v01.checkRutina01(FICH_ARCO)
     # Si no existe el fichero pdf, generamos el plot para la rutina 01
#    if not os.path.exists("./Rut01_dat/Rutina01_plot_1night_"+directorio[0:6]+".pdf"):
    Rutina01_v01.Plot1night(directorio)
    
    # Cargamos ajustes de la rutina02
    Rutina02_v01.cargarAjustes(FLAT_REF)
    # Lanzamos la rutina 02.
    print "EJECUTANDO RUTINA 02: Posición e intensidad del flat ..."
    print "========================================================"
    Rutina02_v01.rutina02Run(FICH_FLAT)
    


#Comprobamos que se ha introducido un parámetro al programa y que sea un directorio
if len(sys.argv)==2:
    if os.path.exists(sys.argv[1]) and not os.path.isfile(sys.argv[1]):
        generarListaFicheros()
        run_Rutina01_Rutina02(sys.argv[1])
        print "EJECUTANDO RUTINA 04: Control del nivel de BIAS ..."
        print "==================================================="
        Rutina04_v01.runRutina04(sys.argv[1])
        print "EJECUTANDO RUTINA 05: Calculando tiempos de observación ..."
        print "==================================================="
        Rutina05_v01.runRutina05(sys.argv[1])
        # Hacemos los plots
        Rutina01_v01.plotHistory()
        Rutina02_v01.plotHistory()
        Rutina04_v01.plotHistory()
    else:
        print "El directorio introducido no existe"
else:
    print "El numero de parámetros es incorrecto."
    print "Debes introducir el directorio de trabajo:"
    print "SINTAXIS: python RutinaMaster [directorio]"

