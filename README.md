# vaspdos
 Atom-spesific Projected LDOS plot for VASP output
 
 **WARNING:** This code is still in development. Many bugs and wrong results may happen. Please use with care.

## Requirements
 1. python3
 2. numpy
 3. matplotlib

## How to use
 In VASP job directory, run

 `python3 vaspdos.py -ion <atom number>`

 Make sure you have already run VASP for lm-decomposed DOS calculation. This code reads the PROCAR file to plot the projected LDOS for a spesific atom in the system.
