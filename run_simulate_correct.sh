#!/bin/bash

PRO="protanopia"
DEU="deuteranopia"
TRI="tritanopia"

blindtype=$TRI

input="./images/ishihara_23.png"

output_nocorrect="./nocorrect"
output_preprocess="./corrected"
output_final="./final"

sens=".1"

# Copy original image for comparison
cp $input "./original"

# Run simulating on original as a control test.
python SimulateColorBlind.py --type $blindtype --sensitivity $sens --yes --out $output_nocorrect $input

# Run the correct then simulate.
python CorrectColorBlind.py --type $blindtype --sensitivity $sens --yes --out $output_preprocess $input
python SimulateColorBlind.py --type $blindtype --sensitivity $sens --yes --out $output_final $output_preprocess".png"
