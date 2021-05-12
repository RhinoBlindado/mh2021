#!/bin/bash

# [CASTELLANO]
# Función: Realiza las ejecuciones de los algoritmos y almacenar los datos de salida.
# Practica 2.b - Algoritmos Greedy y Busqueda Local para PAR
# Asignatura: Metaheuristicas
# Autor: Valentino Lugli (Github: @RhinoBlindado)
# Abril, Mayo 2021
# 

# 
# [ENGLISH]
# Function: Execute all the algorithms and save the outputs.
# Practice 2.b - Greedy and Local Search Algorithms for Clustering with Restrictions
# Course: Metaheuristics
# Author: Valentino Lugli (Github: @RhinoBlindado)
# April, May 2021


# Datasets
dataName=("ZOO" "GLASS" "BUPA")

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

for (( i=0; i<3; i+=1 )); do
    echo "---${dataName[$i]}---"
    echo "  GENE"
    echo "  - 10% Restrictions"
    for (( j=0; j<5; j+=1 )); do
        echo -n "[$(($j + 1))/5]..."
        echo "[$(($j + 1))/5]" >> $outputDir/${dataName[$i]}_GENE_10.txt
        ./GENE_exe ${instances[$i]} ${dimensions[$i]} ${clusters[$i]} ${inputSet[$i]} ${inputR10[$i]} ${seed[$j]} >> $outputDir/${dataName[$i]}_GENE_10.txt
        echo "Done"
    done
    echo "  - 20% Restrictions"
    for (( j=0; j<5; j+=1 )); do
        echo -n "[$(($j + 1))/5]..."
        echo "[$(($j + 1))/5]" >> $outputDir/${dataName[$i]}_GENE_20.txt
        ./GENE_exe ${instances[$i]} ${dimensions[$i]} ${clusters[$i]} ${inputSet[$i]} ${inputR20[$i]} ${seed[$j]} >> $outputDir/${dataName[$i]}_GENE_20.txt
        echo "Done"
    done
    echo "  MEME"
    echo "  - 10% Restrictions"
        for (( j=0; j<5; j+=1 )); do
        echo -n "[$(($j + 1))/5]..."
        echo "[$(($j + 1))/5]" >> $outputDir/${dataName[$i]}_MEME_10.txt
        ./MEME_exe ${instances[$i]} ${dimensions[$i]} ${clusters[$i]} ${inputSet[$i]} ${inputR10[$i]} ${seed[$j]} >> $outputDir/${dataName[$i]}_MEME_10.txt
        echo "Done"
    done
    echo "  - 20% Restrictions"
    for (( j=0; j<5; j+=1 )); do
        echo -n "[$(($j + 1))/5]..."
        echo "[$(($j + 1))/5]" >> $outputDir/${dataName[$i]}_MEME_20.txt
        ./MEME_exe ${instances[$i]} ${dimensions[$i]} ${clusters[$i]} ${inputSet[$i]} ${inputR20[$i]} ${seed[$j]} >> $outputDir/${dataName[$i]}_MEME_20.txt
        echo "Done"
    done
done
