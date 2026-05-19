# PIKOEGUI
GUI of PIKOE code

provide
(1) User freindly interface and tooltips.
(2) Simple plot tool for quick analysis. 
(3) Global optical potentials for the input to pikoe.

## How to use 
In pikoe1.1/qt_files/ folder,  run "python3 qt_pikoe.py". 

## Requirements 
0. Python3 environment. 
   
1. PIKOE executable is included in pikoe1.1/bin/ folder.
In case of linux or other OS, one have to compile and configure for the executable location. 

2. omget executable is included in pikoe1.1/qt_files/omget_RIPL3
In case of linux or other OS, one have to compile and configure for the executable location.

Note: When loading/saving 'cnt' file directly, 
      it is assumed that pikoe is executed at the location of cnt file. 
      (default: input 'cnt' file is located in ../cnt/ )
      In other words, the file path in the input must be relative to the cnt file location. 
       
