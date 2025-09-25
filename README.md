# RNO_G_balloon_project

This part of the code read all the files into python form and plot it respectively, that include:  
1. Balloon trace
2. theory results
3. experimental results
4. error measured by pulser data  (not included in this repositories)

## Balloon trace
**Balloon_class.ipynb**  
  
use the gps data to extract the balloon trace, and have the useful function like:     
search the time period that balloon is within certain radius of certain radius    
generate the position table based on the start and end time period (used in theory results)  

## Theory results
**theory/**  
  
generated theorical raytracing result for dt with a given source using radiopropa  
it have 2 mode:  
1. given the file array and start and end time to generate theory result by time (take result from balloon trace)  
2. given the fixed postion of balloon and station to generate it  
in this experiment, greenland_simple, greenland_poly5 ice model is used, but it is easily swap to any ice model supported by radiopropa  
note: sphere_sizes_i need to be manually tuned if switched between far away object and near object  

## Experimental results
**balloon_data_analyzer.ipynb**  
supported by:   
**event_v2.py** (signal event process class)  
**sampler2_0.py** (data reading module)  
  
the analyzer read in the whole run of data, make it the format easy to filter, process and write out  
and provide the optimal method to plot things like, waveform, fft, spectrogram etc  
note, in this experiment, the method balloon_phase_extracter() is used to extract the phase  
which is include in Bob's disscertation at: https://biblio.ugent.be/publication/01JQNZR41J5WFNT2GNNBC3K174  

## Summary of results
**Balloon_project_plotter_ver2.0.ipynb**  
  
to summarize all the results from all three modules above, a sperate module, plotter is introduced,   
This read in the specific info from balloon, theory and experiment  
and allow user to adjust the specific parameter of the plot   
This module is less general, it is frequently changed based on the new plotting needs  
