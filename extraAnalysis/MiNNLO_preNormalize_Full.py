import ROOT 
import array
import pdb

def printTH1(h):
    print([h.GetBinContent(i) for i in range(1,h.GetNbinsX()+1)])

def printTH1N(h): #normalize
    hint = h.Integral(1,h.GetNbinsX())
    list = [h.GetBinContent(i)/hint for i in range(1,h.GetNbinsX()+1)]
    print(list)
    print(sum(list))

def listNh(h): #return normalized list
    hint = h.Integral(1,h.GetNbinsX())
    list = [h.GetBinContent(i)/hint for i in range(1,h.GetNbinsX()+1)]
    return list

def Ratiol(list1,list2):
    return [x/y for x,y in zip (list1,list2)]

def printr(l,ro):
    print([round(x,ro) for x in l])

#fqq = ROOT.TFile("qqOutputHistMiNNLO.root")
#fgg = ROOT.TFile("ggOutputHistMiNNLO.root")
fqq = ROOT.TFile("qqOutputHistPlotsMiNNLOFull.root")
fgg = ROOT.TFile("ggOutputHistPlotsMiNNLOFull.root")

vars = ["MassAllj","nJets","mjj","dEtajj","jetPt[0]","jetPt[1]","absjetEta[0]","absjetEta[1]"]
vars2 = ["m4l","CleanJet","DiJetMass","DeltaEtajj","LeadingJetPt","SubLeadingJetPt","LeadingJetEta","SubLeadingJetEta"]
vars = ["dEtajj"]
vars2 = ["DeltaEtajj"]
kfacs = [1.02244,0.98414,0.97058,0.95705,0.95456,0.92758,0.91712,0.87614,0.81093]

#xsec in fiducial region https://arxiv.org/pdf/2108.05337.pdf
qqxsec = 17.45
nNNLOxsec = 20.04
ggxsec = nNNLOxsec-qqxsec

#m4l distribution as reference for normalization
hqqm = fqq.Get("h_m4l")
hggm = fgg.Get("h_m4l")
hqqmNG = fqq.Get("UW_m4l")
hggmNG = fgg.Get("UW_m4l")
qqfac = qqxsec/hqqm.Integral(1,hqqm.GetNbinsX())
ggfac = ggxsec/hggm.Integral(1,hggm.GetNbinsX())
qqfacNG = qqxsec/hqqmNG.Integral(1,hqqmNG.GetNbinsX())
ggfacNG = ggxsec/hggmNG.Integral(1,hggmNG.GetNbinsX())


for i,var in enumerate(vars2):
    var1 = vars[i]
    if var != "m4l":
        hqq = fqq.Get("h_"+var)
        hgg = fgg.Get("h_"+var)
        hqqNoGen = fqq.Get("UW_"+var)
        hggNoGen= fgg.Get("UW_"+var)
    else:
        hqq = hqqm
        hgg = hggm
        hqqNoGen = hqqmNG
        hggNoGen= hggmNG
    hqqEW = fqq.Get("Wh_"+var)
    hggEW= fgg.Get("Wh_"+var)
    
    hqq.Scale(qqfac)
    hgg.Scale(ggfac)
    hqqEW.Scale(qqfac)
    hggEW.Scale(ggfac)
    hqqNoGen.Scale(qqfacNG)
    hggNoGen.Scale(ggfacNG)

    #hsum = hqq.Clone("tot_MassAllj_trueNNLO")
    hsum = hqq.Clone("tot_%s_trueNNLO"%var1)
    hsum.Add(hgg)

    hsumEW = hqqEW.Clone("tot_%s_trueEWC"%var1)
    hsumEW.Add(hggEW)

    hsumNoGen = hqqNoGen.Clone("tot_%s_trueNNLONoGenW"%var1)
    hsumNoGen.Add(hggNoGen)

    if i == 0:
        fout = ROOT.TFile("Sum_qqZZggZZ_new.root","RECREATE")
    else:
        fout = ROOT.TFile("Sum_qqZZggZZ_new.root","UPDATE")
    fout.cd()
    hsum.Write()
    hsumEW.Write()

    fout.Close()

    print(var1)
    #printTH1N(hqq)
    #printTH1N(hqqEW)
    #printTH1N(hgg)
    #printTH1N(hggEW)
    pdb.set_trace()
    printTH1N(hsum)
    #printTH1N(hsumNoGen)
    
    #Print group for checkin with paper
    '''
    print("Normalized hsum, hsumNoGen and paperNorm, followed by two ratios")
    paperNorm = [0.17590936266589538, 0.3781595634643841, 0.1963998830023097, 0.10069941959681805, 0.05568004727354278, 0.052217673625462646, 0.02160027834458957, 0.01459272741863249, 0.006666502098751068]    
    lhsum = listNh(hsum)
    lhsumNoGen = listNh(hsumNoGen)
    printr(lhsum,7)
    printr(lhsumNoGen,7)
    printr(paperNorm,7)
    printr(Ratiol(lhsum,paperNorm),7)
    printr(Ratiol(lhsumNoGen,paperNorm),7)
    '''
    printTH1N(hsumEW)
    print("")


#Temporary comparison
compare = False
if compare:
    fc1 = ROOT.TFile("Sum_qqZZggZZOld.root")
    hc1EW = fc1.Get("tot_MassAllj_trueEWC")
    hc1 = fc1.Get("tot_MassAllj_trueNNLO")
    #print([hc1EW.GetBinContent(i) for i in range(1,hsum.GetNbinsX()+1)])
    #print([hc1.GetBinContent(i) for i in range(1,hc1.GetNbinsX()+1)])

    fc2 = ROOT.TFile("Sum_qqZZggZZ.root")
    hc2EW = fc2.Get("tot_MassAllj_trueEWC")
    hc2 = fc2.Get("tot_MassAllj_trueNNLO")
    #pdb.set_trace()

    '''
    print("Previous:")
    print("With EWK")
    print([round(hc1EW.GetBinContent(i)/hc1EW.Integral(1,9),6) for i in range(1,hc1EW.GetNbinsX()+1)])
    print("Without EWK")
    print([round(hc1.GetBinContent(i)/hc1.Integral(1,9),6) for i in range(1,hc1.GetNbinsX()+1)])
    print("Ratio with/without")
    print([round(hc1EW.GetBinContent(i)/hc1EW.Integral(1,9)/hc1.GetBinContent(i)*hc1.Integral(1,9),6) for i in range(1,hc1.GetNbinsX()+1)])
    print("")
    print("New:")
    print("With EWK")
    print([round(hc2EW.GetBinContent(i)/hc2EW.Integral(1,9),6) for i in range(1,hc2EW.GetNbinsX()+1)])
    print(sum([hc2EW.GetBinContent(i)/hc2EW.Integral(1,9) for i in range(1,hc2EW.GetNbinsX()+1)]))
    #print([hsum.GetBinContent(i) for i in range(1,hc2EW.GetNbinsX()+1)])
    print("Without EWK")
    print([round(hc2.GetBinContent(i)/hc2.Integral(1,9),6) for i in range(1,hc2.GetNbinsX()+1)])
    print(sum([hc2.GetBinContent(i)/hc2.Integral(1,9) for i in range(1,hc2.GetNbinsX()+1)]))
    print("Ratio with/without")
    print([round(hc2EW.GetBinContent(i)/hc2EW.Integral(1,9)/hc2.GetBinContent(i)*hc2.Integral(1,9)/kfacs[i-1],9) for i in range(1,hc2.GetNbinsX()+1)])
    '''

