#!/usr/bin/env bash

images=$(ls ../assets/*.dockerfile | sed 's|../assets/\([^\.]*\).dockerfile|\1|g')
for image in $images
do
	apptainer pull stavisvols-psp_$image-latest.img docker://stavisvols/$image:latest
done


