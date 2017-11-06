#!/bin/bash

IS=(1 2 5 8 10)
PREFIX=_echo_

for i in "${IS[@]}"; do
    ffmpeg -i $1.wav -filter_complex aecho=0.8:0.88:$((i*10)):0.4 $1$PREFIX$i.wav    

done
