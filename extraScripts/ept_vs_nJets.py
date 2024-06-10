import ROOT 
import pdb 
import json
import array
import math

with open('varsFile.json') as var_json_file:
    myvar_dict = json.load(var_json_file)

_binning = {}
for key in myvar_dict.keys(): #key is the variable
    if "Mass" in key and not "Full" in key:
        _binning[key] = myvar_dict["MassAllj"]["_binning"]
    else:
        _binning[key] = myvar_dict[key]["_binning"]
_binning["Mass"] = _binning["MassAllj"]
#_binning["ept"] = [7.0,10.0,15.0,20.0,30.0,50.0,100.0,200.0]
_binning["ept"] = [7.0,20.0,50.0,100.0,150.0,200.0]

def rebin(hist,varName):
    ROOT.SetOwnership(hist, False)
    #No need to rebin certain variables but still might need overflow check
    if varName not in ['eta']:
        bins=array.array('d',_binning[varName])
        Nbins=len(bins)-1 
        hist=hist.Rebin(Nbins,hist.GetName()+"NoConfusion",bins)
    else:
        Nbins = hist.GetSize() - 2
    add_overflow = hist.GetBinContent(Nbins) + hist.GetBinContent(Nbins + 1)
    lastbin_error = math.sqrt(math.pow(hist.GetBinError(Nbins),2)+math.pow(hist.GetBinError(Nbins+1),2))
    hist.SetBinContent(Nbins, add_overflow)
    hist.SetBinError(Nbins, lastbin_error)
    hist.SetBinContent(Nbins+1,0)
    hist.SetBinError(Nbins+1,0)
    if not hist.GetSumw2(): hist.Sumw2()
    return hist

def listh(h): #return list
    
    list = [h.GetBinContent(i) for i in range(1,h.GetNbinsX()+1)]
    return list

def percentDiff(h1,h2):
    list = [abs(h1.GetBinContent(i)-h2.GetBinContent(i))/(h1.GetBinContent(i)+h2.GetBinContent(i))*2.*100. for i in range(1,h1.GetNbinsX()+1)]
    return list

def printr(l,ro):
    print([round(x,ro) for x in l])


pre_datasets = ["zz4l-amcatnlo","ggZZ4e","ggZZ2e2mu","ggZZ4m"]
xsecspre = [1.218*1.0835,0.001586*1.7,0.003194*1.7,0.001586*1.7]
checkgen = False
#if checkgen:
#    xsecsUL = xsecspre


#on-shell ZZ files, not UL but use previous naming

fUL = ROOT.TFile("TreeFile_ZZSelector_Hists07Jun2024-ZZ4l2018_MVAsel_Inclusive.root")
fhUL = ROOT.TFile("Hists07Jun2024-ZZ4l2018_MVA.root")

treetag = "_fTreeNtuple_"
channels = ["eeee","eemm","mmee","mmmm"]

lumitot = 59.74 #fb-1

hdict = {}
weightExpr = "weight"
if checkgen:
    weightExpr = "genWeight"

#=========================================
initBins = "400,0,200"
initBins2 = "50,0,5"
njmax = 3
variations = ["nominal", "RECO_up","RECO_dn","ID_up","ID_dn"]
histslist = []
#=========================================

for chan in channels:

    hdict[chan] = []
    #Only look at 4e channel for ept for now
    if chan != "eeee":
        continue

    for nj in range(0,njmax+1):
        nj_print = nj
        if nj==njmax:
            nj_print = "%sorMore"%njmax
        histnamesPt = ["ePt_with_%s_jets_"%nj_print+chan+"_"+vari for vari in variations]
        histnamesEta = ["eEta_with_%s_jets_"%nj_print+chan+"_"+vari for vari in variations]

        #Have to do this manually for now
        #0,1,2,3,4 for nominal and 4 up/dn shifts
        for var_ind,histname in enumerate(histnamesPt):
            exec("hPt%s = ROOT.TH1D(histname,histname,%s)"%(var_ind,initBins) )
            exec("histslist.append(hPt%s)"%var_ind) 
        for var_ind,histname in enumerate(histnamesEta):
            exec("hEta%s = ROOT.TH1D(histname,histname,%s)"%(var_ind,initBins2) )
            exec("histslist.append(hEta%s)"%var_ind) 

        #Only 1 histogram but keep the loop from previous codes
        for i,h in enumerate(histslist):
           
            datasets = pre_datasets
            fin = fUL
            fh = fhUL
            xsecs = xsecspre
            
            additional = ""
            lumi = lumitot
            tag = ""
            


                

            for j,ds in enumerate(datasets):
                treename = ds + tag + treetag + chan
                foldername = ds + tag
                tree = fin.Get(treename)
                h_sumw = fh.Get(foldername+"/"+"sumweights")
                sumweights = h_sumw.Integral(0,h_sumw.GetNbinsX()+1)
                if not tree:
                    print("something wrong getting tree %s"%treename)
                #pdb.set_trace()
                for evt in tree:
                    if evt.nJets == nj or (nj==njmax and evt.nJets > nj):
                        exec("h.Fill(evt.l1pt,evt.%s*xsecs[j]*lumi*1000/sumweights%s)"%(weightExpr,additional))
                        exec("h.Fill(evt.l2pt,evt.%s*xsecs[j]*lumi*1000/sumweights%s)"%(weightExpr,additional))
                        exec("h.Fill(evt.l3pt,evt.%s*xsecs[j]*lumi*1000/sumweights%s)"%(weightExpr,additional))
                        exec("h.Fill(evt.l4pt,evt.%s*xsecs[j]*lumi*1000/sumweights%s)"%(weightExpr,additional))
      
            h = rebin(h,"ept")
            hdict[chan].append(h)
            #print(h.GetName())
            
            #print(listh(h))

# Analyze results

#for k,h in enumerate(hdict[var]["eemm"]):
#    h.Add(hdict[var]["mmee"][k])

#*Printout*#
fout = ROOT.TFile("ePtEta_vs_nJets_HistOutput.root","RECREATE")
fout.cd()

for chan in channels:
    if chan == "mmee":
        continue
    
    #Only look at 4e channel for now
    if not chan == "eeee":
        continue
    for h in hdict[chan]:
        h.Write()
fout.Close()
print("Histograms written")
