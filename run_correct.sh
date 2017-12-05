#!/bin/bash

PRO="protanopia"
DEU="deuteranopia"
TRI="tritanopia"

blindtype=$PRO

input="./images/Plate13.gif"

output="./run_result"

sens=".3"

# Run the CorrectColorBlind.py script
python CorrectColorBlind.py --type $blindtype --sensitivity $sens --show --yes --out $output $input
