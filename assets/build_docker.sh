#!/bin/bash

dockerfiles=$(ls *.dockerfile)
for dockerfile in $dockerfiles;
do
    name=$(echo $dockerfile | cut -d. -f1)
    docker build -t stavisvols/psp_$name -f $dockerfile ./
    docker push stavisvols/psp_$name
done

