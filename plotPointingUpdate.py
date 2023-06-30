import os
import random
import sys
import multiprocessing
import datetime
import numpy as np
import scipy.stats
import sunpy
import matplotlib.pyplot as plt
import astropy.io.fits as fits

#launch date; used to compute a feature for scans
launchHinode = datetime.datetime(year=2006,month=9,day=22,hour=21,minute=36,second=0)

def timeOfDay(date):
    """Given a datetime, return the seconds since the start of the day"""
    startOfDay = datetime.datetime(year=date.year, month=date.month, day=date.day, hour=0, minute=0, second=0)
    return (date-startOfDay).total_seconds()

def timeOfYear(date):
    """Given a datetime, return the seconds since the start of the year"""
    startOfYear = datetime.datetime(year=date.year, month=1, day=1, hour=0, minute=0, second=0)
    return (date-startOfYear).total_seconds() 

def kernelReg(x,y,sigma,xEval):
    """Nadaraya-Watson kernel regression

    Given N input/outputs {x_i}/{y_} and an evaluation location xe, calculate 
        (\sum_{i=1}^N w(x_i,xe) * y_i) / (\sum_{i=1}^N w(w_i,xe))
    where w(x,q) is the Gaussian PDF centered at w, evaluated at xe with a
    standard deviation of sigma. This is a local, weighted average.
    """
    rv = scipy.stats.norm(scale=sigma)

    if np.isscalar(xEval):
        xmmu = rv.pdf(xEval-x)
        return np.sum(xmmu*y)/np.sum(xmmu)
    else:
        yEval = np.zeros(xEval.shape)
        for xi in range(xEval.size):
            xmmu = rv.pdf(xEval[xi]-x)
            yEval[xi] = np.sum(xmmu*y)/np.sum(xmmu)
        return yEval


class LinearPlusLUT:
    """Class for a model that's a lookup table (LUT) with a linear residual"""

    def __init__(self, bandwidth=3.0/365, yOutlierQuantile=0.01):
        """
        bandwidth -- the bandwidth of the LUT kernel regression
        yOutlierQuantile -- ignore the top/bottom yOutlierQuantile in the fit
            this is important since some of the fits are really quite large and
            are outliers (often due to radically incorrect pointing in SOT-SP)
        """
        self.X = None
        self.y = None
        self.yOutlierQuantile = yOutlierQuantile
        self.bandwidth = bandwidth

    def fit(self, XLinear, XLUT, y):
        """Fit the model to predict y using a linear model in XLinear and LUT
            model in XLUT.
            
        This is done by doing a cross-validation to find      
        """
        NInit = XLUT.shape[0]
     
        #if we're filtering, filter 
        if self.yOutlierQuantile is not None:
            qLow = np.nanquantile(y, self.yOutlierQuantile)
            qHigh = np.nanquantile(y, 1-self.yOutlierQuantile)
            #keep if it's between the two quantiles
            k = (qLow < y) &  (y < qHigh)

            #subsample once as if the data were never passed in
            y = y[k]
            XLUT = XLUT[k]
            XLinear = XLinear[k]

        self.y = y
        self.XLUT = XLUT
        XLinear = XLinear
        N = self.XLUT.shape[0]

        numFolds = 5

        #assign data to folds
        foldId = np.floor(np.arange(N).astype(np.float)/N*numFolds)

        #cross-validated predictions
        cvLUT = np.zeros_like(self.y)
        for fi in range(numFolds):
            tr = foldId!=fi
            te = foldId==fi
            cvLUT[te] = kernelReg(self.XLUT[tr], self.y[tr], self.bandwidth, self.XLUT[te])

        #figure out the error after the LUT has been applied
        LUTResidual = self.y-cvLUT

        #predict it with a linear model; the predictions will get added
        XLinear1 = np.hstack([XLinear.reshape(-1,1), np.ones((N,1))])
        self.w, _, _, _ = np.linalg.lstsq(XLinear1, LUTResidual,rcond=None)

    def predict(self, XLinear, XLUT):
        """Predict a model given the linear feature in XLinear and the LUT
        feature in XLUT"""
        XLUTPred = kernelReg(self.XLUT, self.y, self.bandwidth, XLUT) 
        LinearPred = self.w[0]*XLinear+self.w[1]
        return XLUTPred + LinearPred

    def savePlaintext(self, target):
        """Save the model as a plaintext

        Model storage format:

        (linear term), (constant term)
        (X for LUT Entry 0), (Y for LUT Entry 0)
        (X for LUT Entry 1), (Y for LUT Entry 1)
        ...
        (X for LUT Entry N), (Y for LUT Entry N)
        """

        with open(target,"w") as fh:
            fh.write("%f,%f\n" % (self.w[0],self.w[1]))
            for i in range(self.XLUT.shape[0]):
                fh.write("%f,%f\n" %  (self.XLUT[i], self.y[i]))
        
    def loadPlaintext(self, filen):
        """Load the model from plaintext in the format above"""
        data = open(filen).read().strip().split("\n")
        linearModel = data[0].split(",")
        self.w = np.array([float(linearModel[0]), float(linearModel[1])])

        #there are N lines, but the first one is the linear term
        N = len(data)-1

        self.XLUT, self.y = np.zeros((N,)), np.zeros((N,))
        for i in range(1,len(data)):
            line = data[i].split(",")
            self.XLUT[i-1], self.y[i-1] = float(line[0]), float(line[1])


def handle(t):
    """Compute the pointing update and compute variable to correlate against
    Put in a global function to enable multiprocessing"""

    #unpack the arguments
    scanI, scan, origSrc, updateSrc = t

    print(scanI, scan)

    #open the original and upate
    scanOrig = fits.open(os.path.join(origSrc, scan))
    scanUpdate = fits.open(os.path.join(updateSrc, scan))

    #the original XCEN and YCEN; we recalculate this since the headers aren't
    #always quite right since XCEN gets loaded from Level1 in a suboptimal way
    #in some cases
    XCENO = np.mean(scanOrig[38].data[:,[0,-1]])
    YCENO = np.mean(scanOrig[39].data[[0,-1],:])

    #updated XCEN, YCEN
    XCENU = scanUpdate[0].header['XCEN']
    YCENU = scanUpdate[0].header['YCEN'] 

    #date of the scan
    dateStr = scan.replace(".fits","")
    date = datetime.datetime.strptime(scan,"%Y%m%d_%H%M%S.fits")

   
    #all the covariates as a tuple:
    #time since launch (y), time of year (y) , time of day (d), temp at CCD,
    #and temp at CEB
    covariates = (
         (date-launchHinode).total_seconds()/(3600*24*365.25),
         (timeOfYear(date)/(3600*24*365.25)),
         (timeOfDay(date)/(3600*24)),
         scanUpdate[0].header['T_SPCCD'],
         scanUpdate[0].header['T_SPCEB'],
        )

    #the original and upated XCEN, YCEN
    pointing = (XCENU, YCENU, XCENO, YCENO)

    return (dateStr, covariates, pointing)



if __name__ == "__main__":
    #
    #Code to do a few things; likely only of interest if you have all of the 
    #SOT-SP level 2 data. This will: (a) save the update table; (b) fit models
    #to predict dx/dy updates; (c) plot predictions from these models; and 
    #(e) plot the fit updates.

    #where the original SOTSP Level 2 data is
    origSrc = "SOTSPLevel2/"

    #where the data's been packed
    updateBase = "updateLevel2/"

    #where to dump the graphs
    visTarget = "plotGraphs/"

    #options to do 
    validOptions = ["savetable", "fitmodels", "plotpredictions", "plotfits"]
    todoList = [c.lower() for c in sys.argv[1:]]

    for t in todoList:
        if t not in validOptions:
            print("Unrecognized task %s" % t)
            print("Valid tasks: %s" % (" ".join(validOptions)))
            sys.exit(1)

    if len(todoList) == 0:
        print("No tasks specified")
        print("Valid tasks: %s" % (" ".join(validOptions)))
        sys.exit(1)

    toHandle = []

    #include alignments based on:
    #- Main -- all the SDO/HMI equatorial scans that could be aligned
    #- Pole -- SDO/HMI pole scans via a method from Wang et al. forthcoming
    #- MDIHighConfidence -- SOHO/MDI scans with confident alignments 
    #
    #The release also includes MIDLowConfidence, which includes ones that there
    #is less confidence in. These alignments show similar trends with more noise
    for sub in ["Main","Pole","MDIHighConfidence"]:
        updateSrc = os.path.join(updateBase, sub)
        toHandle += [(scanI+len(toHandle),scan,origSrc,updateSrc) for scanI, scan in enumerate(sorted(os.listdir(updateSrc)))]

    if not os.path.exists(visTarget):
        os.mkdir(visTarget)

    covariates = []
    pointing = []

    #Do this in multiprocessing; turn this up or down depending on your system
    P = multiprocessing.Pool(12)
    results = P.map(handle, toHandle)
    P.close()

    #stack the results into one numpy array
    dateStrs = [t[0] for t in results]
    covariates, pointing = np.vstack([t[1] for t in results]), np.vstack([t[2] for t in results])

    #dx := update in xcen (new XCEN - old XCEN), dy := update in ycen
    dx, dy = pointing[:,0]-pointing[:,2], pointing[:,1]-pointing[:,3]

    if "savetable" in todoList:
        with open("pointingTableSOTSP.txt","w") as fh:
            fh.write("datestr,XCENU,YCENU,XCENO,YCENO,DXCEN,DYCEN,totalTime,timeOfYear,timeOfDay,T_SPCCD,T_SPCEB\n")

            order = list(range(len(covariates)))
            order.sort(key=lambda i: dateStrs[i])
            for i in order:
                dxi, dyi = dx[i], dy[i] 
                fh.write("%s,%f,%f,%f,%f,%f,%f," % ((dateStrs[i],)+tuple(pointing[i,:])+(dxi,dyi)))
                fh.write("%f,%f,%f,%f,%f\n" % tuple(covariates[i,:]))


    if "fitmodels" in todoList:
        #fit the model with the linear part being total time and the
        #LUT being the time of year
        dxModelAll = LinearPlusLUT()
        dxModelAll.fit(covariates[:,0], covariates[:,1], dx)
        dxModelAll.savePlaintext("dxModel.txt")
   
        dyModelAll = LinearPlusLUT()
        dyModelAll.fit(covariates[:,0], covariates[:,1], dy)
        dyModelAll.savePlaintext("dyModel.txt")


    if "plotpredictions" in todoList:

        #Do a quick test to show predictions. This is a sanity check and
        #a demonstration of the code
        total = covariates[:,0]
        timeOfYear = covariates[:,1]

        N = dx.size

        #Test on every 5 just for this quick check
        te = np.arange(N) % 5 == 0
        tr = np.arange(N) % 5 != 0

        dxModel = LinearPlusLUT() 
        dxModel.fit(total[tr], timeOfYear[tr], dx[tr])
        dxPredict = dxModel.predict(total[te], timeOfYear[te])

        dyModel = LinearPlusLUT()
        dyModel.fit(total[tr], timeOfYear[tr], dy[tr])
        dyPredict = dyModel.predict(total[te], timeOfYear[te])

        #Plot a scatter
        plt.figure(figsize=(4,4))
        plt.scatter(dy[te],dyPredict,5)
        minV, maxV = np.nanquantile(np.hstack([dy[te], dyPredict]), [0.02, 0.98]) 
        plt.xlim(minV, maxV); plt.ylim(minV, maxV)
        plt.xlabel("Fit Change in YCEN (arcsec)")
        plt.ylabel("Predicted Change in YCEN (arcsec)")
        plt.savefig(visTarget+"/dyPredictScatter.pdf")
        plt.close()

        plt.figure(figsize=(4,4))
        plt.scatter(dx[te],dxPredict,5)
        minV, maxV = np.nanquantile(np.hstack([dx[te], dxPredict]), [0.02, 0.98]) 
        plt.xlim(minV, maxV); plt.ylim(minV, maxV)
        plt.xlabel("Fit Change in XCEN (arcsec)")
        plt.ylabel("Predicted Change in XCEN (arcsec)")
        plt.savefig(visTarget+"/dxPredictScatter.pdf")
        plt.close()


    if "plotfits" in todoList:
        #from now on for variable names, x = covariate, y = pointing update
        xLabel = {'tot': 'Total Time (y)','toy': 'Time of Year (y)', 'tod': 'Time of Day (d)',
                  't_spccd':'Temp at CCD (deg C)', 't_spceb':'Temp at CEB (deg C)'}
        yLabel = {'dx': 'Update in XCEN (arcsec)', 'dy': 'Update in YCEN (arcsec)'}

        yYLim = {'dx': [-5,50], 'dy': [-5,70]}
        xFigSize = {'tot': (8,4)}

        #generate plots following the paper
        for yName, yval in [("dx",dx),("dy",dy)]:
            for xName, xval in [("tot", covariates[:,0]), ("toy", covariates[:,1]), 
                                ("tod", covariates[:,2]), ("t_spccd", covariates[:,3]),
                                ("t_spceb", covariates[:,4])]:

                #compute a gaussian kde for showing density
                vals = np.hstack([xval.reshape(-1,1), yval.reshape(-1,1)]).T
                kernel = scipy.stats.gaussian_kde(vals)(vals)

                #default to (6,4)
                plt.figure(figsize=(xFigSize[xName] if xName in xFigSize else (4,4)))
                plt.scatter(xval, yval, 5, kernel, label='data')
                plt.ylim(yYLim[yName][0], yYLim[yName][1])
                plt.xlim(np.min(xval), np.max(xval))
                plt.xlabel(xLabel[xName])
                plt.ylabel(yLabel[yName])
                plt.tight_layout()
                plt.savefig(visTarget+"/"+xName+"_"+yName+".pdf")
                plt.close()


