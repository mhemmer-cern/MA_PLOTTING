#!/bin/bash

## with PhotonSelection UPDATED##
path="/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/"
pathdata="/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/"
declare -a OutputFiles=("OmegaToPiZeroGamma_2074" "OmegaToPiZeroGamma_2084")
declare -a EventCut=("0008e113" "0008d113")
PhotonConvCut="00200009327000008250400000"
ClusterCut="411791206f032230000"
PionCut="01631031000000d0"
declare -a OmegaCut=("0r631031000000d0" "0v631031000000d0" "0x631031000000d0")

## with PhotonSelection Pi0 Rotation##
# path="/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/"
# pathdata="/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/"
# declare -a OutputFiles=("OmegaToPiZeroGamma_2078" "OmegaToPiZeroGamma_2088")
# declare -a EventCut=("0008e113" "0008d113")
# PhotonConvCut="00200009327000008250400000"
# ClusterCut="411791206f032230000"
# PionCut="01631031000000d0"
# declare -a OmegaCut=("0x631031000000d0" "0y631031000000d0" "0z631031000000d0")

## with PhotonSelection##
# path="/mnt/wwn-0x50000395b2b85149-part1/PreparedData/2020_10_29_pp13TeV_background_JJMC/" # mnt/wwn-0x50000395b2b85149-part1/PreparedData/2020_09_03_pp_13TeV_Pi0Rotation/MC/
# pathdata="/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_10_29_pp13TeV/"
# declare -a OutputFiles=("OmegaToPiZeroGamma_2074" "OmegaToPiZeroGamma_2084")
# declare -a EventCut=("0008e113" "0008d113")
# PhotonConvCut="00200009327000008250400000"
# ClusterCut="411791206f032230000"
# PionCut="01631036000000d0"
# declare -a OmegaCut=("0r631031000000d0" "0v631031000000d0" "0x631031000000d0")

## without PhotonSelection##
# path="/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/"
# pathdata="/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/"
# declare -a OutputFiles=("OmegaToPiZeroGamma_6074" "OmegaToPiZeroGamma_6084")
# declare -a EventCut=("0008e113" "0008d113")
# PhotonConvCut="00200009327000008250400000"
# ClusterCut="411791206f032230000"
# PionCut="01631031000000d0"
# declare -a OmegaCut=("0r631031000000d0" "0v631031000000d0" "0x631031000000d0")

## with AP like cut##
# path="/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_04_pp13TeV/"
# pathdata="/mnt/wwn-0x50000395b2b85149-part1/PreparedData2/2020_12_03_pp13TeV/"
# declare -a OutputFiles=("OmegaToPiZeroGamma_2077" "OmegaToPiZeroGamma_2087")
# declare -a EventCut=("0008e113" "0008d113")
# PhotonConvCut="00200009327000008250400000"
# ClusterCut="411791206f032230000"
# PionCut="01631031000000d0"
# declare -a OmegaCut=("0x631031000000d0" "0x631031010000d0" "0x631031020000d0" "0x631031030000d0")


# Declare a string array with type
declare -a modenames=("Standard" "LowerRebin" "HigherRebin" "SmallerFitRange" "HigherFitRange" "OneSigma" "ThreeSigma" "GausAndPol")
rm -rf "Test/"
mkdir -p "Test/"
for (( i = 0; i < 2; i++ )); do
  for val in "${OmegaCut[@]}"; do
    rm -rf "Comp/"${OutputFiles[$i]}"/"$val
    mkdir -p "Comp/"${OutputFiles[$i]}"/"$val
    for (( j = 0; j < 8; j++ )) do
      rm -rf "JJ/"${OutputFiles[$i]}"/"${modenames[$j]}"/"$val
      mkdir -p "JJ/"${OutputFiles[$i]}"/"${modenames[$j]}"/"$val
      rm -rf "Data/"${OutputFiles[$i]}"/"${modenames[$j]}"/"$val
      mkdir -p "Data/"${OutputFiles[$i]}"/"${modenames[$j]}"/"$val
      rm -rf "Comp/"${OutputFiles[$i]}"/"${modenames[$j]}"/"$val
      mkdir -p "Comp/"${OutputFiles[$i]}"/"${modenames[$j]}"/"$val
      ./install/plot $pathdata $path ${OutputFiles[$i]} ${EventCut[$i]} $PhotonConvCut $ClusterCut $PionCut $val $j
    done
  done
done
