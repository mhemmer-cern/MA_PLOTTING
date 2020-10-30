#!/bin/bash

path="/mnt/wwn-0x50000395b2b85149-part1/PreparedData/2020_10_29_pp13TeV_background_JJMC/" # mnt/wwn-0x50000395b2b85149-part1/PreparedData/2020_09_03_pp_13TeV_Pi0Rotation/MC/
pathdata="/mnt/wwn-0x50000395b2b85149-part1/PreparedData/2020_10_29_pp13TeV_background_JJMC/"
# Declare a string array with type
declare -a OutputFiles=("OmegaToPiZeroGamma_2074" "OmegaToPiZeroGamma_2084")
declare -a EventCut=("0008e113" "0008d113")
PhotonConvCut="00200009327000008250400000"
ClusterCut="411791206f032230000"
PionCut="01631036000000d0"
declare -a OmegaCut=("0r631031000000d0" "0v631031000000d0" "0x631031000000d0")
declare -a Mode=("R")

for (( i = 0; i < 2; i++ )); do
  for val in "${OmegaCut[@]}"; do
    rm -rf "JJ/"${OutputFiles[$i]}"/"$val
    mkdir -p "JJ/"${OutputFiles[$i]}"/"$val
    rm -rf "Data/"${OutputFiles[$i]}"/"$val
    mkdir -p "Data/"${OutputFiles[$i]}"/"$val
    rm -rf "Comp/"${OutputFiles[$i]}"/"$val
    mkdir -p "Comp/"${OutputFiles[$i]}"/"$val
    for mode in "${Mode[@]}"; do
      ./install/plot $pathdata $path ${OutputFiles[$i]} ${EventCut[$i]} $PhotonConvCut $ClusterCut $PionCut $val 'E' 'J' $mode
    done
  done
done
