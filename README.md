# BLSthermalFit
This repo is Matlab code for fitting detection function of the BLS. It is used for fitting maximal detecttable k

Use with Matlab CFtool (curve fitting toolbox)

before opening CFtool import the simulation data and make variables globals 

disable "Auto fit" (fitting takes too much time - start it only manually) .On x-axis put the detected TFP frequency and on y-axis. 
Change to custom equation. - y = f(f)
and model:
fitBLSthermal(f,FWTM, FWTM, A,dc) - symmetric detection function
fitBLSthermal(f,FWTMx, FWTMy, A,dc) - assymetric detection function
fitBLSthermal(f+df,FWTM, FWTM, A,dc) - compansate for inaccuracy in the external field
