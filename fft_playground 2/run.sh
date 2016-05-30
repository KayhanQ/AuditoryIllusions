#!/bin/sh

ffmpeg -y -i "$1" -ac 1 -c:a pcm_f32le -f data -map 0:0 ./in.pcm
./foo > ./out.pcm
play -r 44100 -e float -b 32 -L -c 1 -t raw ./out.pcm