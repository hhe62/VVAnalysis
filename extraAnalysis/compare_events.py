import ROOT as r
import pdb,subprocess,math,array
import sys,json,os

if 'qqZZ' in sys.argv[2]:
    sample = 'zz4l-amcatnlo'
if 'ggZZ4e' in sys.argv[2]:
    sample = "ggZZ4e"

fa=r.TFile(sys.argv[1])
r.SetOwnership(fa,False)
sumweights_hist = fa.Get(str("/".join([sample, "sumweights"])))
r.SetOwnership(sumweights_hist, False)
totWgt = sumweights_hist.Integral(0,sumweights_hist.GetNbinsX()+1)
if sample == 'zz4l-amcatnlo':
    xsec = 1.218
    kfac = 1.0835
    lumi = 59.7*1000
if sample == 'ggZZ4e':
    xsec = 0.001586
    kfac = 1.7
    lumi = 59.7*1000
factor = xsec*kfac*lumi/totWgt

name='%s_eeee'%"LepPt1"
name2='%s_eeee'%"LepPt1Full"
varHist1 = fa.Get(sample+"/"+name).Clone()
varHist2 = fa.Get(sample+"/"+name2).Clone()
evt1 = varHist1.Integral(1,varHist1.GetNbinsX()+1)*factor #this stage hasn't included overflow
evt2 = varHist2.Integral(1,varHist2.GetNbinsX()+1)*factor
print("Weighted events, mZ 60,120: %s"%evt1)
print("Weighted events, m4l 80,110: %s"%evt2)