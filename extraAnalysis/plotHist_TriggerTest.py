import ROOT as r
import pdb,subprocess,math,array
import sys,json,os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-c", "--chan", dest="channel",
                  help="channel")

(options, args) = parser.parse_args()

#search DataMC and "#name by format" in this script to see input dependent part
def getTextBox(x,y,axisLabel,size=0.2,color=1,rotated=False):
    texS = r.TLatex(x,y,'#color[%s]{%s}'%(color,axisLabel))
    texS.SetNDC()
    #rotate for y-axis                                                                                                                                                                                             
    if rotated:
        texS.SetTextAngle(90)
    texS.SetTextFont(42)
    #texS.SetTextColor(ROOT.kBlack)                                                                                                                                                                                
    texS.SetTextSize(size)
    texS.Draw()
    return texS

def checkZeroBin(hist,label,histn):
    contents =[hist.GetBinContent(i) for i in range(1,hist.GetNbinsX()+1)]
    contentsn =[histn.GetBinContent(i) for i in range(1,histn.GetNbinsX()+1)]
    #print(contents)
    #print(contentsn)
    for i in range(1,hist.GetNbinsX()+1):
        if hist.GetBinContent(i)==0.:
            #print("WARNING: %s contains 0 in bin %s"%(label,i))
            if not histn.GetBinContent(i) == 0.:
                raise Exception('Zero denominator with nonzero numerator!')
            hist.SetBinContent(i,0.000001)
        if hist.GetBinContent(i)<0.:
            print("WARNING: %s contains negative value in bin %s"%(label,i))

def rebin(hist, binning):        
    bins = array.array('d', binning)
    hist = hist.Rebin(len(bins)-1, "", bins)
    #add overflow
    num_bins = hist.GetNbinsX()
    add_overflow = hist.GetBinContent(num_bins) + hist.GetBinContent(num_bins + 1)
    add_error = math.sqrt(math.pow(hist.GetBinError(num_bins),2)+math.pow(hist.GetBinError(num_bins+1),2))
    hist.SetBinContent(num_bins, add_overflow)
    hist.SetBinError(num_bins, add_error)
    hist.SetBinContent(num_bins+1, 0.)
    hist.SetBinError(num_bins+1, 0.)
    return hist


pdfcommand=['convert']
#with open('varsFile.json') as var_json_file:
#    myvar_dict = json.load(var_json_file)

nostat = False #no statistical error or not
samples = sys.argv[3].split(',') # results in ['zz4l-amcatnlo',"zz4l-powheg"] or ['ggZZ4e'] for comparison or single plotting
fullsamples = ['zz4l-amcatnlo','ggZZ4e','ggZZ4m','ggZZ4t','ggZZ2e2mu','ggZZ2e2tau']
fullkfac = [1.0835,1.7,1.7,1.7,1.7,1.7]
fullxsec = [1.218, 0.001586, 0.001586, 0.001586,0.003194,0.003194]
lumi = 59.7*1000 #name by format
if '17' in sys.argv[1]:
    lumi = 41.5*1000    
fullfac = []
channel = ''
if options.channel:
    channel = options.channel
if channel == 'mmee':
    ind0 = 1
    ind1 = 0
else:
    ind0 = 0
    ind1 = 1


#print("Plotting samples:",samples)
#print("=============================================================")


suffix = "_"+sys.argv[3].replace(',','_') #name by format
if 'Extra4eCut' in sys.argv[1]: 
    suffix = suffix+'_4eCut'
if '2e' in sys.argv[1]:
    suffix = suffix+'_2eCut'

colors = [4,2] if len(samples) ==2 else [3]
markers = [1,1] if len(samples) ==2 else [1]
binnings = [7.,10.,12.,20.,23.,30.,50.,100.,200.]
varstr="LepPt1Full LepPt2Full LepPt3Full LepPt4Full LepPt1 LepPt2 LepPt3 LepPt4"
#varstr2="e1PtSortedFull e2PtSortedFull e1PtSorted e2PtSorted"
varstr2="e1PtSortedFull e2PtSortedFull"


if 'DataMC' in samples or 'AllData' in samples: #name by format
    varlist = varstr2.split(' ')
else:
    varlist = varstr.split(' ')

#varlist=["nJets"]
for var in varlist:
    print("=======%s=========="%var)
    prettyVar = "Lepton%s p_{T} [GeV]"%(var.replace("LepPt","").replace("Full",""))
    if 'e1' in var or 'e2' in var: #name by format
        prettyVar = prettyVar.replace('Leptone','Electron').replace('PtSorted','')
    fnames=[sys.argv[1],sys.argv[2]] #first numerator, then denominator
    labels=[sname for sname in samples]
    for l,la in enumerate(labels):
        if la == 'DataMC':
            labels[l] = "qqZZ+ggZZ"
        if la == 'AllData':
            labels[l] = "Data"
    hists=[]

    #name by format
    if not 'AllData' in samples and not "DataMC" in samples:
        unfname=['%s_eeee'%var]
    else:
        unfname=['%s_eemm'%var,'%s_mmee'%var]
   

    
    num = 0.
    den = 0.

    fa=r.TFile(fnames[0])
    fb=r.TFile(fnames[1])
    r.SetOwnership(fa,False)
    r.SetOwnership(fb,False)

    if not 'Data' in samples[0]:
        sumweights_hist = fa.Get(str("/".join([samples[0], "sumweights"]))) #provided first hist is not data
    #sumweights_hist2 = fb.Get(str("/".join([samples[0], "sumweights"])))
    
        r.SetOwnership(sumweights_hist, False)
        totWgt = sumweights_hist.Integral(0,sumweights_hist.GetNbinsX()+1)
    #totWgt2 = sumweights_hist2.Integral(0,sumweights_hist.GetNbinsX()+1)
    #assert totWgt == totWgt2

    #print("==========Total Weight ===============")
    #print(totWgt,totWgt2)
        if 'zz4l-amcatnlo' in samples[0]:
            xsec = 1.218
            kfac = 1.0835
            #lumi = 59.7*1000
        if 'ggZZ4e' in samples[0]:
            xsec = 0.001586
            kfac = 1.7
            #lumi = 59.7*1000
        factor = xsec*kfac*lumi/totWgt

    if 'DataMC' in samples:
        
        for s,sample in enumerate(fullsamples):
            swgt_hist = fa.Get(str("/".join([fullsamples[s], "sumweights"])))
            r.SetOwnership(swgt_hist, False)
            swgt = swgt_hist.Integral(0,swgt_hist.GetNbinsX()+1)
            sfac = fullxsec[s]*fullkfac[s]* lumi/swgt
            fullfac.append(sfac)

    
    for i in range(len(samples)):
        if not 'DataMC' in samples[i]:
            if len(unfname) ==1:
                hunfa = fa.Get(samples[i]+"/"+unfname[0]).Clone()
                hunfb = fb.Get(samples[i]+"/"+unfname[0]).Clone()
            else:
                #pdb.set_trace()
                hunfa = fa.Get(samples[i]+"/"+unfname[ind0]).Clone()
                htmpa = fa.Get(samples[i]+"/"+unfname[ind1]).Clone()
                if not channel == 'eemm' and not channel == 'mmee':
                    hunfa.Add(htmpa)

                hunfb = fb.Get(samples[i]+"/"+unfname[ind0]).Clone()
                htmpb = fb.Get(samples[i]+"/"+unfname[ind1]).Clone()
                if not channel == 'eemm' and not channel == 'mmee':
                    hunfb.Add(htmpb)
        else: #this should correspond to eemm+mmee case only
            fullhistsa = []
            fullhistsb = []
            for s,sample in enumerate(fullsamples):
                hunfat = fa.Get(fullsamples[s]+"/"+unfname[ind0]).Clone()
                htmpat = fa.Get(fullsamples[s]+"/"+unfname[ind1]).Clone()
                if not channel == 'eemm' and not channel == 'mmee':
                    hunfat.Add(htmpat)
                fullhistsa.append(hunfat)

                hunfbt = fb.Get(fullsamples[s]+"/"+unfname[ind0]).Clone()
                htmpbt = fb.Get(fullsamples[s]+"/"+unfname[ind1]).Clone()
                if not channel == 'eemm' and not channel == 'mmee':
                    hunfbt.Add(htmpbt)
                fullhistsb.append(hunfbt)
            
            hunfa = fullhistsa[0]*fullfac[0]
            hunfb = fullhistsb[0]*fullfac[0]
            for s in range(1,len(fullsamples)):
                h_tmpa = fullhistsa[s]*fullfac[s]
                h_tmpb = fullhistsb[s]*fullfac[s]
                hunfa.Add(h_tmpa)
                hunfb.Add(h_tmpb)

        hunfa = rebin(hunfa,binnings)
        hunfb = rebin(hunfb,binnings)
        checkZeroBin(hunfb, 'denominator',hunfa)
        
        hunf_a_b = hunfa.Clone()
                                    #hunf_amb = hunfa.Clone()
        hunf_a_b.Divide(hunfb)
                                    #hunf_amb.Add(hunfb,-1) #this hist is not used
        hists.append(hunf_a_b) #append in the orders of labels
        
        if i==0:
            if not 'Data' in samples[0]:
                num = hunfa.Integral(1,hunfa.GetNbinsX()+1)*factor #only take amcnlo numerator and denominator for two MC comparison case
                print([round(hunfa.GetBinContent(j)*factor,3) for j in range(1,hunfa.GetNbinsX()+1)])
                den = hunfb.Integral(1,hunfb.GetNbinsX()+1)*factor
                print([round(hunfb.GetBinContent(j)*factor,3) for j in range(1,hunfb.GetNbinsX()+1)])
            else:
                num = hunfa.Integral(1,hunfa.GetNbinsX()+1)
                den = hunfb.Integral(1,hunfb.GetNbinsX()+1)
                print("==========MC bins=================")
                print([round(hunfa.GetBinContent(j),3) for j in range(1,hunfa.GetNbinsX()+1)])
                #print([round(hunfb.GetBinContent(j),3) for j in range(1,hunfb.GetNbinsX()+1)])
        
        if i == 1 and 'Data' in samples[i]:
            numd = hunfa.Integral(1,hunfa.GetNbinsX()+1)
            dend = hunfb.Integral(1,hunfb.GetNbinsX()+1)
            print("==========Data bins=================")
            print([round(hunfa.GetBinContent(j),3) for j in range(1,hunfa.GetNbinsX()+1)])
            #print([round(hunfb.GetBinContent(j),3) for j in range(1,hunfb.GetNbinsX()+1)])
        
        

        
    

    maxs = []

    c1=r.TCanvas("canvas")
    c1.cd()
    #c1.SetLogy()
    #if "Full" in var:
    #    c1.SetLogx()

    r.gStyle.SetOptDate(0)
    #pdb.set_trace()
    for i in range(len(hists)):
        hists[i].SetLineColor(colors[i])
        maxs.append(hists[i].GetMaximum())
        hists[i].SetTitle("")
        hists[i].GetYaxis().SetTitle("Trig Eff") #yt
        hists[i].GetYaxis().SetTitleOffset(1)
        hists[i].GetXaxis().SetTitle(prettyVar)
        hists[i].GetXaxis().SetTitleSize(0.04)
        hists[i].GetXaxis().SetTitleOffset(1.3)
        hists[i].SetStats(0)

    for i in range(len(hists)):
        hists[i].SetMaximum(1.2)
        hists[i].SetMinimum(0.)
        hists[i].SetMarkerStyle(1)
        if i == 0:   
            #hists[i].Draw("HIST P")
            if nostat:
                hists[i].Draw("HIST")
            else:
                hists[i].Draw()
            r.gStyle.SetLegendFont(42)
            r.gStyle.SetLegendTextSize(0.03)
            
            legend = r.TLegend (0.7 ,0.3 ,0.85 ,0.4)
            legend.SetFillStyle(0)
        else:
            #hists[i].SetMarkerStyle(markers[i])
            #hists[i].Draw("HIST P SAME")
            if nostat:
                hists[i].Draw("HIST SAME")
            else:
                hists[i].Draw("SAME")
        legend.AddEntry(hists[i],labels[i])

    legend.SetLineWidth (0)
    legend.Draw("same")
    latex = r.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    #if 'nonreg' in sys.argv[1]:
    textbox_num= getTextBox(0.74,0.28,"num. %s"%round(num,2),0.03,colors[0])
    textbox_den= getTextBox(0.74,0.25,"den. %s"%round(den,2),0.03,colors[0])
    if 'Data' in samples[1]:
        textbox_numd= getTextBox(0.74,0.22,"num. %s"%round(numd,2),0.03,colors[1])
        textbox_dend= getTextBox(0.74,0.19,"den. %s"%round(dend,2),0.03,colors[1])
    if "Full" in var:
        textbox= getTextBox(0.35,0.97,"With 80 GeV < m_{4l}< 110 GeV",0.03)
    else:
        textbox= getTextBox(0.35,0.97,"With 60 GeV < m_{Z1},m_{Z2} < 120 GeV",0.03)
    #latex.DrawLatex(0.74,0.83 ,"59.7fb^{-1}")
    dirName = 'TrigEffPlots'+ suffix
    if not os.path.isdir(dirName):
        os.mkdir(dirName)

    try:
        c1.SaveAs("%s/%s_TrigEff.png"%(dirName,var))
    except:
        print("Problem saving plot.")
    c1.Clear()
    pdfcommand.append("%s/%s_TrigEff.png"%(dirName,var))

pdfcommand.append('%s/TrigEff_plots%s.pdf'%(dirName,suffix))
subprocess.call(pdfcommand)
