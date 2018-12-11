# Cepheid Experiment

This repository collects all python programs that were used for the Cepheid experiment performed in winter of 2018 at ETH.

The programs have to be called in a certain order and assume a certain directory structure.

Here is the basic directory structure, for the moment without the data:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

-o Data Redution
 |-> day1.sh
 |-> day2.sh
 |-> day3.sh
 ...
 |-> day_n.sh if more measuring days are available
 |-> data_reduction.py
 
-o Day1
 | Here comes the data for Day 1
 
-o Day2
 | Here comes the data for Day 2
 
-o Day3
 | Here comes the data for Day 3
 
 ...
 
 -o Day_n
  | If more measuring days are available
  
  cleaning.py
  find_flux.py <- This one is obsolete, but was kept to be sure
  find_relative_flux.py
  general_fit.py
  periodogram.py
  periodogram_reader.py
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  The data for each measuring day has to be structured as follows:
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Day_n -o
         |-o Dark_For_Flats
         | |-> Put here the darks for the flats
         |
         |-o Flats
         | |-> Put here the flats
         |
         |-o Star_Name_Dayn  <- So for example FF_Aql_Day1 is valid. A different directory has to be created for each star.
         | |-o Dark_For_Science
         | | |-> Put here the darks for the science images
         | |
         | |-o Science
             |-> Put here the science images
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 The files themselves do not need to be named in a certain way, but do need to be fits files with a correct header.
 
 To process the data (python2.7 and python3 have to be installed, ideally through anaconda, together with the astropy, numpy, scipy, matplotlib, astroscrappy and photutils packages), do as follows:
 1. Go into Data Reduction and in a terminal there use
   bash day1.sh
 and similarly for each other day. This is the most time-consuming program. This will create for each star for the given day three folders:
 
Star_Name_Dayn_Processed contains all processed fits files, where the background has been set to zero. Those are typically the cleanest images, but do behave strangely when scaling is applied to them.
Star_Name_Dayn_With_Background contains all processed fits files, where the background has not been set to zero. Using certain scales, they can be made to look 'natural' compared to the ones with no background, but cannot be used efficiently for the next computations.
Star_Name_Dayn_Errors contains fits files of the same dimensions as the original ones, except each pixel contains the error that was made during the whole correction process for that one pixel.
 
 2. Go back to the main folder and use
   python2 periodogram.py && cleaning.py
This will create the folder Individual_Fits which contains for each day and each star the Gaussian fits near the peak for both the Cepheid and the secondary standard star (both in 2D along the x- and y-axis starting at the peak, and in 3D near the peak), as well as additional fits files on which the secondary standard star that was used is marked.
Also, a pickle file that will be used file will be created, and cleaned by cleaning.py. cleaning.py was very specific to the experiment we made, so might not result in corrected data for other data that was set as input.

3. Finally, use
  python3 periodogram_reader.py
This will output the mean apparent magnitudes, the periods and false alarm probabilites in the console for each star, and a folder called Final_Results will be created with the light curves, the periodograms and the relative phase diagrams for the fitting frequency (NOT the actual frequency) inside for each star.
  
  
