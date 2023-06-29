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

## minimalCorrect.py

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

## visualizePointingUpdate.py

This loads the updated pointing as well as a HMI field observation and makes
pictures that show that the new pointing better aligns the data. You'll need
to download hmi data from JSOC (http://jsoc.stanford.edu/ajax/exportdata.html)
in order to use this script.

## plotPointingUpdate.py

This script can take a folder full of SOT-SP Level 2 observations and the
updated data and: (a) dump a table that shows pointing and updates for the
shared observations; (b) fit models that predict the pointing update; (c) plot
plots quickly verifying the pointing update models; (d) plot the update as a
function of various properties. Since this script requires a complete set of
SOT-SP observations, it's less likely of interest.

