# SOTSPPointing
Code for using/loading updated Hinode SOT-SP Pointing.

Written by David Fouhey to accompany results from:

``Large-Scale Spatial Cross-Calibration of Hinode/SOT-SP and SDO/HMI'' 
by David F. Fouhey, Richard E.L. Higgins, Spiro K. Antiochos, Graham
Barnes, Marc DeRosa, J. Todd Hoeksema, K.D. Leka, Yang Liu, Peter W. Schuck,
Tamas I.  Gombosi. ApJS 264 49, 2023. 

Some of the pointing data, in particular for polar scans, is from an upcoming
paper by Ruoyu Wang, Richard Higgins, David Fouhey, et al. (most of the above)

This repository contains scripts that can work with and without an updated pointing
information cache. This updated information can be found at:
https://fouheylab.eecs.umich.edu/~fouhey/HinodeSOTSPLevel2Update.tgz

There are three scripts in the repository. In order of likely usefulness, they
are:

# Scripts

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


# Pointing information

We provide an updated pointing data cache at :
https://fouheylab.eecs.umich.edu/~fouhey/HinodeSOTSPLevel2Update.tgz
as well a table of pointing updates in ``pointingTableSOTSP.txt``.



## pointingTableSOTSP.txt

This is a csv file containing all pointing information. The first two lines are: 

```
datestr,XCENU,YCENU,XCENO,YCENO,DXCEN,DYCEN,totalTime,timeOfYear,timeOfDay,T_SPCCD,T_SPCEB
20110101_065127,-598.687317,-385.336823,-619.500305,-391.698090,20.812988,6.361267,4.274841,0.000782,0.285729,-35.485500,-2.910680
```

The columns are:
- (1) datestr: the nominal date of the Level2 Scan
- (2) XCENU: the *u*pdated XCEN
- (3) YCENU: the *u*pdated YCEN
- (4) XCENO: the *o*riginal XCEN. Notably, this is re-computed from the X_COORDINATE because the XCEN in the Level 2 header isn't always right
- (5) YCENO: the *o*riginal YCEN, computed like same as XCENO.
- (6) DXCEN: the change in XCEN (XCENU - XCENO)
- (7) DYCEN: the change in YCEN (YCENU - YCENO)
- (8) totalTime: the time since Hinode Launch (in years)
- (9) timeOfYear: the time of year (in years), i.e., Jan 1 is 0, July 1 is ~0.5, Dec 31 is 1
- (10) timeOfDay: the time of day (in days), i.e., Midnight is 0, 06:00 is 0.25, 18:00 is 0.75, and 11:59 is 1
- (11) T_SPCCD: the temperature at the CCD (i.e., camera chip), as reported in SOT-SP Level 1 data. This is pulled from a somewhat arbitrary level 1 scan that went into the Level 2 data.
- (12) T_SPCEB: the temperature at the CEB, camera electronics box, as reported in SOT-SP Level 1 data.

The information for the pointing update comes from: co-alignment with SDO/HMI
as well as confident co-alignments to SOHO/MDI.


## Fits files for update

The pointing data tgz provides four groups of updated fits files. These are put
in four groups that are of varying interest. They are as follows:
- Main: primarily equatorial scans and based on co-alignment with SDO/HMI following Fouhey ApJS 2023
- Pole: primarily polar scans and based on co-alignment with SDO/HMI following a forthcoming paper
- MDIHighConfidence: primarily equatorial scans, based on co-alignment with SOHO/MDI. These are highly confident alignments and most ought to be correct.
- MDILowConfidence: primarily equatorial scans, based on co-alignment with SOHO/MDI. These include less confident alignments, and these are a bit noisier.


## Structure of each fits file

Each updated fits file is renamed identically to the old fits file. It has the
same number of extensions. The only extensions that are filled with data are 
the 38 and 39th, which were ``X_COORDINATE`` and ``Y_COORDINATE`` respectively. The remaining
fits data is dummy data, which is put in so that the same extension numbers are the same.

Here's an example loading a fits file (with some lines removed for brevity). The script 
visualizePointingUpdate.py shows how to use the fits file.

```
>>> import astropy.io.fits as fits
>>> X = fits.open("updateLevel2/Main/20141118_160506.fits")
>>> X[0].header
SIMPLE  =                    T / conforms to FITS standard
BITPIX  =                    8 / array data type
NAXIS   =                    0 / number of array dimensions
EXTEND  =                    T
XCEN    =   -334.5798034667969 / XCEN, updated
YCEN    =   -287.9953308105469 / YCEN, updated
XSCALE  =      0.2952414721443 / XSCALE with scale correction
YSCALE  =      0.3162061153349 / YSCALE with scale correction
...
T_SPCCD =             -32.5464 / L1 T_SPCCD: Temp at the CCD
T_SPCEB =             -1.82974 / L1 T_SPCEB: Temp at the Camera Electronics Box
L1DATE  = '2014-11-18T16:05:06.231' / Timestamp of data for Level 1 readings
PNTDATE = '20141118_163600'    / Timestamp of SDO/HMI data providing pointing
WARP00  =   0.5854536324953824 / Affine warp 0,0 entry; see comments
WARP01  = 0.004877281726691565 / Affine warp 0,1 entry; see comments
WARP02  =     1138.56124039986 / Affine warp 0,2 entry; see comments
WARP10  = -0.00455391520662339 / Affine warp 1,0 entry; see comments
WARP11  =   0.6270257951755208 / Affine warp 1,1 entry; see comments
WARP12  =    1351.306456012081 / Affine warp 1,2 entry; see comments
BNDMINX =                 1139 / Minimum X of warped region
BNDMINY =                 1348 / Minimum Y of warped region
BNDMAXX =                 1639 / Maximum X of warped region
BNDMAXY =                 1606 / Maximum Y of warped region
PACKDATE= '20230627'           / Date data was stored
COMMENT
COMMENT Pointing has been updated from SDO/MDI following
COMMENT Fouhey et al 2023 ApJS 264 49
COMMENT Important information on affine warp. The warp assumes:
COMMENT 1) Hinode/SOT-SP has been expanded using the slit position
COMMENT 2) X (left/right) is coordinate 1, Y (top/down) is coordinate 2
COMMENT 3) SDO/HMI has been flipped/rotated exactly 180 degrees
...
>>> X[38].data
array([[-464.1023 , -463.807  , -463.51172, ..., -206.61899, -206.3237 , -206.02843],
       [-464.09988, -463.80463, -463.50934, ..., -206.6166 , -206.32132, -206.02605],
        ...,
       [-463.13364, -462.83835, -462.54306, ..., -205.65031, -205.35504, -205.05977],
       [-463.13123, -462.83597, -462.54068, ..., -205.64793, -205.35266, -205.05737]], dtype=float32)

>>> X[39].data
array([[-351.38245, -351.38467, -351.3869 , ..., -353.3257 , -353.32794, -353.33017],
       [-351.0662 , -351.06842, -351.07065, ..., -353.00946, -353.0117 , -353.01392],
        ...,
       [-222.98666, -222.98889, -222.99112, ..., -224.92966, -224.93188, -224.93411],
       [-222.67043, -222.67265, -222.67488, ..., -224.61342, -224.61565, -224.61787]], dtype=float32)

```




