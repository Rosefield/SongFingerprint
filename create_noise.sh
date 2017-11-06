#!/bin/bash

IS=(1 2 5 8 10)
PREFIX=_noisy_

for i in "${IS[@]}"; do
    ffmpeg -i $1.wav -i whitenoise_$i.wav -filter_complex amix=duration=shortest $1$PREFIX$i.wav    

done
