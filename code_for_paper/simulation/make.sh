#!/bin/bash

## example: sh make.sh script_dt

echo ----Creating result folder----
echo results/res_$(date +"%Y-%m-%d")  # results/res_2020-12-16
mkdir results/res_$(date +"%Y-%m-%d")

args=("$@")
echo ----Using script: ${args[0]}.R----  # script_dt

echo ----Creating code folder----
echo results/res_$(date +"%Y-%m-%d")/codes_$(date +"%H-%M")_${args[0]}
mkdir results/res_$(date +"%Y-%m-%d")/codes_$(date +"%H-%M")_${args[0]}

echo ----Copying codes to result folder----
cp ${args[0]}.R results/res_$(date +"%Y-%m-%d")/codes_$(date +"%H-%M")_${args[0]}/
cp -r R results/res_$(date +"%Y-%m-%d")/codes_$(date +"%H-%M")_${args[0]}/

echo ----Running script...----
Rscript ${args[0]}.R
echo ----Finished all!----
