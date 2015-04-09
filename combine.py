from numpy import load,save,histogramdd,arange,random,diff,array

datadir = '../data/'

unfacATLAS = load(datadir+'unfoldedac.npy')
ATLASsyst = {}
ATLASsyst['jer']  = load(datadir+'jer.npy')
ATLASsyst['LUMI'] = load(datadir+'LUMI.npy')

unfacCMS = load(datadir+'unfoldedac.npy')
CMSsyst = {}
CMSsyst['jer'] = load(datadir+'jer.npy')
CMSsyst['LUMI'] = load(datadir+'LUMI.npy')

binsac = arange(-0.020,0.024,0.004)
binssyst = arange(-1.5,1.8,0.3)

## bin likelihoods in multi-dimensional histograms
lhATLAS, binsATLAS = histogramdd([unfacATLAS]+[ATLASsyst[syst] for syst in ATLASsyst.keys()],
                                         bins=[binsac,binssyst,binssyst])
lhCMS, binsCMS = histogramdd([unfacCMS]+[CMSsyst[syst] for syst in CMSsyst.keys()],
                                     bins=[binsac,binssyst,binssyst])

## sampling 1000 points
sampledbins = [random.randint(0,len(binsac)-1,1000)] + [random.randint(0,len(binssyst)-1,1000)
                                                     for syst in ATLASsyst.keys()]
points = zip(*sampledbins)

## combine likelihoods
weights = [lhATLAS[point]*lhCMS[point]/(lhATLAS.sum()*lhCMS.sum()) for point in points]

## bin centers
midpoints_ac   = [binsac[:-1] + diff(binsac)/2]
midpoints_syst = [binssyst[:-1] + diff(binssyst)/2
                  for syst in ATLASsyst.keys()]
midpoints = array(midpoints_ac + midpoints_syst)

## combined traces
values = array([[mp[ibin] for mp,ibin in zip(midpoints,point)]
                for point in points])

import ROOT as r
combinedac = r.TH1D('combinedac','combinedac',len(binsac)-1,binsac[0],binsac[-1])
for val,w in zip(values,weights):
    combinedac.Fill(val[0],w)

combinedac.SaveAs("test.root")
