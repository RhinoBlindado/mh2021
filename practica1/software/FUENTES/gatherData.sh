#!/bin/bash

if [ ! -f "./COPKM_exe" ]; then
    echo " COPKM_exe not found, compiling."
    make COPKM_exe
fi

if [ ! -f "./BL_exe" ]; then
    echo " BL_exe not found, compiling."
    make BL_exe
fi

# Seeds obtained from random.org random integer generator
seed=(860681 980472 206894 161195 420696)

echo "---ZOO---"
# Setting up the variables
instances=101
dimensions=16
clusters=7
inputSet="./input/zoo_set.dat"
inputR="./input/zoo_set_const_10.const"
outputDir="./output/"

echo "--10% Restrictions"
echo "--COPKM"
for (( i=0; i<5; i+=1 )); do
	echo -n "[$(($i + 1))/5]..."
	./COPKM_exe $instances $dimensions $clusters $inputSet $inputR ${seed[$i]} >> $outputDir/zoo_COPKM_10.txt
	echo "Done"
done

echo "--BL"
for (( i=0; i<5; i+=1 )); do
	echo -n "[$(($i + 1))/5]..."
	./BL_exe $instances $dimensions $clusters $inputSet $inputR ${seed[$i]} >> $outputDir/zoo_BL_10.txt
	echo "Done"
done

echo "--20% Restrictions"
echo "--COPKM"
inputR="./input/zoo_set_const_20.const"
for (( i=0; i<5; i+=1 )); do
    echo -n "[$(($i + 1))/5]..."
	./COPKM_exe $instances $dimensions $clusters $inputSet $inputR ${seed[$i]} >> $outputDir/zoo_COPKM_20.txt
	echo "Done"
done

echo "--BL"
for (( i=0; i<5; i+=1 )); do
	echo -n "[$(($i + 1))/5]..."
	./BL_exe $instances $dimensions $clusters $inputSet $inputR ${seed[$i]} >> $outputDir/zoo_BL_20.txt
	echo "Done"
done

echo "---GLASS---"
# Setting up the variables
instances=214
dimensions=9
clusters=7
inputSet="./input/glass_set.dat"
inputR="./input/glass_set_const_10.const"
outputDir="./output/"

echo "--10% Restrictions"
echo "--COPKM"
for (( i=0; i<5; i+=1 )); do
	echo -n "[$(($i + 1))/5]..."
	./COPKM_exe $instances $dimensions $clusters $inputSet $inputR ${seed[$i]} >> $outputDir/glass_COPKM_10.txt
	echo "Done"
done

echo "--BL"
for (( i=0; i<5; i+=1 )); do
	echo -n "[$(($i + 1))/5]..."
	./BL_exe $instances $dimensions $clusters $inputSet $inputR ${seed[$i]} >> $outputDir/glass_BL_10.txt
	echo "Done"
done

echo "--20% Restrictions"
echo "--COPKM"
inputR="./input/glass_set_const_20.const"
for (( i=0; i<5; i+=1 )); do
    echo -n "[$(($i + 1))/5]..."
	./COPKM_exe $instances $dimensions $clusters $inputSet $inputR ${seed[$i]} >> $outputDir/glass_COPKM_20.txt
	echo "Done"
done

echo "--BL"
for (( i=0; i<5; i+=1 )); do
	echo -n "[$(($i + 1))/5]..."
	./BL_exe $instances $dimensions $clusters $inputSet $inputR ${seed[$i]} >> $outputDir/glass_BL_20.txt
	echo "Done"
done

echo "---BUPA---"
# Setting up the variables
instances=345
dimensions=5
clusters=16
inputSet="./input/bupa_set.dat"
inputR="./input/bupa_set_const_10.const"
outputDir="./output/"

echo "--10% Restrictions"
echo "--COPKM"
for (( i=0; i<5; i+=1 )); do
	echo -n "[$(($i + 1))/5]..."
	./COPKM_exe $instances $dimensions $clusters $inputSet $inputR ${seed[$i]} >> $outputDir/bupa_COPKM_10.txt
	echo "Done"
done

echo "--BL"
for (( i=0; i<5; i+=1 )); do
	echo -n "[$(($i + 1))/5]..."
	./BL_exe $instances $dimensions $clusters $inputSet $inputR ${seed[$i]} >> $outputDir/bupa_BL_10.txt
	echo "Done"
done

echo "--20% Restrictions"
echo "--COPKM"
inputR="./input/bupa_set_const_20.const"
for (( i=0; i<5; i+=1 )); do
    echo -n "[$(($i + 1))/5]..."
	./COPKM_exe $instances $dimensions $clusters $inputSet $inputR ${seed[$i]} >> $outputDir/bupa_COPKM_20.txt
	echo "Done"
done

echo "--BL"
for (( i=0; i<5; i+=1 )); do
	echo -n "[$(($i + 1))/5]..."
	./BL_exe $instances $dimensions $clusters $inputSet $inputR ${seed[$i]} >> $outputDir/bupa_BL_20.txt
	echo "Done"
done

