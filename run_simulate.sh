#!/bin/bash

PRO="protanopia"
DEU="deuteranopia"
TRI="tritanopia"

blindtype=$TRI

input="./images/flowers.jpg"

output="./run_result"

sens=".5"

# Run the SimulateColorBlind.py script
python SimulateColorBlind.py --type $blindtype --sensitivity $sens --show --yes --out $output $input
