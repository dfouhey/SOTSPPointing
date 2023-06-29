# SOTSPPointing
Code for using/loading updated Hinode SOT-SP Pointing.

Written by David Fouhey to accompany results from:

``Large-Scale Spatial Cross-Calibration of Hinode/SOT-SP and SDO/HMI'' 
by David F. Fouhey, Richard E.L. Higgins, Spiro K. Antiochos, Graham
Barnes, Marc DeRosa, J. Todd Hoeksema, K.D. Leka, Yang Liu, Peter W. Schuck,
Tamas I.  Gombosi. ApJS 264 49, 2023. 

Some of the pointing data, in particular for polar scans, is from an upcoming
paper by Ruoyu Wang, Richard Higgins, David Fouhey, et al. (most of the above)

The updated pointing data can be found at:
https://fouheylab.eecs.umich.edu/~fouhey/HinodeSOTSPLevel2Update.tgz

In order of likely utility, there are three scripts:


## Scripts

### minimalCorrect.py

This provides an easy, minimal piece of code that will make an estimate of what
the updated pointing should be at a particular date. The code has close to no
dependencies and should work with nearly any version of python. If the pointing
has a known correction (i.e., it was fit), that is used. Otherwise, it makes a 
prediction that should be substantially better than the current pointing.

For instance, here is the prediction at the time of writing:
```
> python minimalCorrect.py 20230608_232500
20230608_232500 23.438624 61.849458
```

### visualizePointingUpdate.py

This loads the updated pointing as well as a HMI field observation and makes
pictures that show that the new pointing better aligns the data. You'll need
to download hmi data from JSOC (http://jsoc.stanford.edu/ajax/exportdata.html)
in order to use this script.

### plotPointingUpdate.py

This script can take a folder full of SOT-SP Level 2 observations and the
updated data and: (a) dump a table that shows pointing and updates for the
shared observations; (b) fit models that predict the pointing update; (c) plot
plots quickly verifying the pointing update models; (d) plot the update as a
function of various properties. Since this script requires a complete set of
SOT-SP observations, it's less likely of interest.


## Pointing information


### Fits files for update

The pointing data tgz provides four groups of updated fits files. These are put
in four groups that are of varying interest. They are as follows:
- Main: primarily equatorial scans and based on co-alignment with SDO/HMI following Fouhey ApJS 2023
- Pole: primarily polar scans and based on co-alignment with SDO/HMI following a forthcoming paper
- MDIHighConfidence: primarily equatorial scans, based on co-alignment with SOHO/MDI. These are highly confident alignments and most ought to be correct.
- MDILowConfidence: primarily equatorial scans, based on co-alignment with SOHO/MDI. These include less confident alignments, and these are a bit noisier.


## Structure of each fits file




### pointingTableSOTSP.txt

This is a csv file containing all pointing information

```
datestr,XCENU,YCENU,XCENO,YCENO,DXCEN,DYCEN,totalTime,timeOfYear,timeOfDay,T_SPCCD,T_SPCEB
20110101_065127,-598.687317,-385.336823,-619.500305,-391.698090,20.812988,6.361267,4.274841,0.000782,0.285729,-35.485500,-2.910680
```

The columns are:
- datestr: the nominal date of the Level2 Scan
- XCENU: the *u*pdated XCEN
- YCENU: the *u*pdated YCEN
- XCENO: the *o*riginal XCEN. Notably, this is re-computed from the X_COORDINATE because the XCEN in the Level 2 header isn't always right
- YCENO: the *o*riginal YCEN, computed like same as XCENO.
- DXCEN: the change in XCEN (XCENU - XCENO)
- DYCEN: the change in YCEN (YCENU - YCENO)
- totalTime: the time since Hinode Launch (in years)
- timeOfYear: the time of year (in years), i.e., Jan 1 is 0, July 1 is ~0.5, Dec 31 is 1
- timeOfDay: the time of day (in days), i.e., Midnight is 0, 06:00 is 0.25, 18:00 is 0.75, and 11:59 is 1
- T_SPCCD: the temperature at the CCD (i.e., camera chip), as reported in SOT-SP Level 1 data. This is pulled from a somewhat arbitrary level 1 scan that went into the Level 2 data.
- T_SPCEB: the temperature at the CEB, camera electronics box, as reported in SOT-SP Level 1 data.


