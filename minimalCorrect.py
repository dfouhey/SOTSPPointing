"""
Minimal piece of python code to get updated SOT-SP pointing information for a
given datetime Should work on nearly any version of python, including even
python2. This is meant to be easy to replicate in another language e.g., in
IDL, Julia, MATLAB, x86_64 assembly, etc.
"""
import os
import sys
import pdb
import datetime
from math import exp

#launch date; used to compute a feature for scans
launchHinode = datetime.datetime(year=2006,month=9,day=22,hour=21,minute=36,second=0)

def getTimeOfYear(date):
    """Given a datetime, return the seconds since the start of the year"""
    startOfYear = datetime.datetime(year=date.year, month=1, day=1, hour=0, minute=0, second=0)
    return (date-startOfYear).total_seconds() 

def loadLUTPlusLinear(filename):
    """Given a filename containing the linear+lookup table model, return the 
    linear coefficients and the lookup table. 
    
    Model storage format:

    (linear term), (constant term)
    (X for LUT Entry 0), (Y for LUT Entry 0)
    (X for LUT Entry 1), (Y for LUT Entry 1)
    ...
    (X for LUT Entry N), (Y for LUT Entry N)
    """
    data = open(filename).read().strip().split("\n")
    linearLine, lutLines = data[0].split(","), data[1:]
    linearCoeff = (float(linearLine[0]), float(linearLine[1]))
    LUT = []
    for line in lutLines:
        line = line.split(",")
        LUT.append((float(line[0]), float(line[1])))
    
    return linearCoeff, LUT

def applyLUT(LUT,xe,bandwidth=3.0/365):
    """Apply the soft-lookup table (basic Nadaraya Watson kernel regression).
    LUT: list of tuples (x,y) providing the LUT
    xe: evaluation location
    bandwidth: bandwidth 

    Given N input/outputs {x_i}/{y_} and an evaluation location xe, calculate 
        (\sum_{i=1}^N w(x_i,xe) * y_i) / (\sum_{i=1}^N w(w_i,xe))
    where w(x,q) is the Gaussian PDF centered at w, evaluated at xe with a
    standard deviation of bandwidth. This is a local, weighted average.
    """
    totalNumerator, totalDenominator = 0, 0
    bandwidthSquared = bandwidth**2
    for x, y in LUT:
        #don't need the scaling factor; this is constant across the samples
        scaledPdf = exp( -0.5 * ((x-xe)**2) / bandwidthSquared )

        #add to numerator and denominator
        totalNumerator += scaledPdf * y
        totalDenominator += scaledPdf
    return totalNumerator / totalDenominator


if __name__ == "__main__":
    if len(sys.argv) < 2:
        #if we have no arguments, just print an error message
        print("%s datetime" % (sys.argv[0]))
        print("Prints datetime, dx, dy for each datetime in the arguments")
        print("Datetime format: %Y%m%d_%H%M%S, or 2020 March 15, 13:14:04 is 20200315_131404")
        sys.exit(1)


    #load hard table (if they've been fit, don't predict)
    LUT2DXDY = {}
    data = open("pointingTableSOTSP.txt").read().strip().split("\n")[1:]
    for line in data:
        lineSplit = line.split(",")
        timestamp = lineSplit[0]
        dx = float(lineSplit[5])
        dy = float(lineSplit[6])
        LUT2DXDY[timestamp] = (dx,dy)

    #load the models (useful for unalignable deep scans, etc.)
    linearDx, LUTDx = loadLUTPlusLinear("dxModel.txt")
    linearDy, LUTDy = loadLUTPlusLinear("dyModel.txt")

    for evalDate in sys.argv[1:]:

        #try to help the user out; date formats are hard and annoying
        evalDate = evalDate.replace(" ","_")
        evalDate = evalDate.replace("/","").replace(".","").replace(":","")

        if evalDate in LUT2DXDY:
            #if the evaluation date is known, just return it
            dx, dy = LUT2DXDY[evalDate]
        else:
            #else, use the predictive model
            date = datetime.datetime.strptime(evalDate,"%Y%m%d_%H%M%S")

            totalTime = (date-launchHinode).total_seconds()/(3600*24*365.25)
            timeOfYear = getTimeOfYear(date)/(3600*24*365.25)

            dx = applyLUT(LUTDx, timeOfYear) + totalTime*linearDx[0] + linearDx[1]
            dy = applyLUT(LUTDy, timeOfYear) + totalTime*linearDy[0] + linearDy[1]

        print("%s %f %f" % (evalDate, dx, dy))
        

