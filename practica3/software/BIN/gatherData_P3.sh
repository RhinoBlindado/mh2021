#!/bin/bash

# [CASTELLANO]
# Función: Realiza las ejecuciones de los algoritmos y almacenar los datos de salida.
# Practica 3.b - Técnicas de Búsqueda basadas en Poblaciones para el Problema del Agrupamiento con Restricciones
# Asignatura: Metaheuristicas
# Autor: Valentino Lugli (Github: @RhinoBlindado)
# Mayo, Junio 2021
# 

# 
# [ENGLISH]
# Function: Execute all the algorithms and save the outputs.
# Practice 3.b - Population-based search techniques for the Clustering Problem with Restrictions
# Course: Metaheuristics
# Author: Valentino Lugli (Github: @RhinoBlindado)
# May, June 2021


# Datasets
dataName=("ZOO" "GLASS" "BUPA")

# Algos
algorithm=("ES")

# Variables
instances=(101 214 345)
dimensions=(16 9 5)
clusters=(7 7 16)
inputSet=("./input/zoo_set.dat" "./input/glass_set.dat" "./input/bupa_set.dat")
inputR10=("./input/zoo_set_const_10.const" "./input/glass_set_const_10.const" "./input/bupa_set_const_10.const")
inputR20=("./input/zoo_set_const_20.const" "./input/glass_set_const_20.const" "./input/bupa_set_const_20.const")
outputDir="./output"

# Seeds obtained from random.org randon integer generator
seed=(860681 980472 206894 954919 426969)

len=${#algorithm[@]}

# For each algorithm...
for (( i=0; i<len; i+=1 )); do

    echo ">>${algorithm[$i]}"
    
    # For each dataset...
    for (( j=0; j<3; j+=1)); do
    
        echo ">>${dataName[$j]}"
        echo "${dataName[$j]}">> $outputDir/${algorithm[$i]}_10.txt
        echo "${dataName[$j]}">> $outputDir/${algorithm[$i]}_20.txt
        # Loop 5 times for each data set, 10% & 20%.
        for (( k=0; k<5; k+=1 )); do
        
            echo -n "[$(($k + 1))/5]..."
            ./${algorithm[$i]}_exe ${instances[$j]} ${dimensions[$j]} ${clusters[$j]} ${inputSet[$j]} ${inputR10[$j]} ${seed[$k]} >> $outputDir/${algorithm[$i]}_10.txt
            ./${algorithm[$i]}_exe ${instances[$j]} ${dimensions[$j]} ${clusters[$j]} ${inputSet[$j]} ${inputR20[$j]} ${seed[$k]} >> $outputDir/${algorithm[$i]}_20.txt
            echo "Done"
            
        done
            
    done

done
