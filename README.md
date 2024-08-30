# FlyGo

## How to Use
The **gene_main.py** file is the main file to run. It is almost entirely set up! You'll just need to 
edit the lines indicated with your directories. In the example code, there are 2 directories- one
for the Hyperoxia .csvs and another for the Mutant .csvs (that Timothy made). 

## Dependencies
This uses pandas (any version should be fine), and standard python libraries (requests, os, warnings, and glob). This
should work on any system where Pandas can be run.

## Structure
The `gene_module.py` script is structured as a module, so all the functions can be found there. The relevant functions
are imported 

## How it works
Moreorless, this loads in the flybase data, then the obo data, matches them up, and does Timothy's preprocessing to link
them all together

## Contact
rschanta@udel.edu