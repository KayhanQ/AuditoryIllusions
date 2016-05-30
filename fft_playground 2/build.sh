#!/bin/sh

g++ -std=c++11 -Wall -O3 -o foo ./*.cpp fft/*.cpp -I. -Ifft
