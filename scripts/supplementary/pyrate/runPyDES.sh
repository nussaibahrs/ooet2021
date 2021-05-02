#!/bin/sh

DIRECTORY="pyrate_output"

for filename in $DIRECTORY/*rep1.txt; do
    echo Running PyRateDES2 for $filename
    python3 PyRateDES2.py -d $filename -A 2 -q 66 61.6000 56.0000 47.8000 41.3000 38.0000 33.9000 28.1000 23.0300 15.9700 11.6200  5.3330  2.5880  0.0117  0.0000 -DivdE -DivdD -n 1000000
done
