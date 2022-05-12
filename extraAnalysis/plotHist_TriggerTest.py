import ROOT as r
import pdb,subprocess,math,array
import sys,json,os

def getTextBox(x,y,axisLabel,size=0.2,rotated=False):
    texS = r.TLatex(x,y,axisLabel)
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
    print(contents)
    print(contentsn)
    for i in range(1,hist.GetNbinsX()+1):
        if hist.GetBinContent(i)==0.:
            print("WARNING: %s contains 0 in bin %s"%(label,i))
            histn.SetBinContent(i,0.)
            hist.SetBinContent(i,0.1)
        if hist.GetBinContent(i)<0.:
            print("WARNING: %s contains negative value in bin %s"%(label,i))

def rebin(hist, binning):        
    bins = array.array('d', binning)
    hist = hist.Rebin(len(bins)-1, "", bins)
    #add overflow
    num_bins = hist.GetSize() - 2
    add_overflow = hist.GetBinContent(num_bins) + hist.GetBinContent(num_bins + 1)
    add_error = math.sqrt(math.pow(hist.GetBinError(num_bins),2)+math.pow(hist.GetBinError(num_bins+1),2))
    hist.SetBinContent(num_bins, add_overflow)
    hist.SetBinError(num_bins, add_error)
    return hist


pdfcommand=['convert']
#with open('varsFile.json') as var_json_file:
#    myvar_dict = json.load(var_json_file)

nostat = False #no statistical error or not
samples = sys.argv[3].split(',') # results in ['zz4l-amcatnlo',"zz4l-powheg"] or ['ggZZ4e'] for comparison or single plotting
print("Plotting samples:",samples)
print("=============================================================")


suffix = "_"+sys.argv[3].replace(',','_')
if 'Extra4eCut' in sys.argv[1]:
    suffix = suffix+'_4eCut'
colors = [3,2] if len(samples) ==2 else [3]
markers = [1,1] if len(samples) ==2 else [1]
binnings = [7.,10.,15.,20.,30.,50.,100.,200.]
varstr="LepPt1Full LepPt2Full LepPt3Full LepPt4Full LepPt1 LepPt2 LepPt3 LepPt4"
varlist = varstr.split(' ')
#varlist=["nJets"]
for var in varlist:
    prettyVar = "Lepton%s p_{T} [GeV]"%(var.replace("LepPt","").replace("Full",""))

    fnames=[sys.argv[1],sys.argv[2]] #first numerator, then denominator
    labels=[sname for sname in samples]
    hists=[]
    unfname='%s_eeee'%var
    num = 0.
    den = 0.

    fa=r.TFile(fnames[0])
    fb=r.TFile(fnames[1])
    r.SetOwnership(fa,False)
    r.SetOwnership(fb,False)
    sumweights_hist = fa.Get(str("/".join([samples[0], "sumweights"])))
    sumweights_hist2 = fb.Get(str("/".join([samples[0], "sumweights"])))
    
    r.SetOwnership(sumweights_hist, False)
    totWgt = sumweights_hist.Integral(0,sumweights_hist.GetNbinsX()+1)
    totWgt2 = sumweights_hist2.Integral(0,sumweights_hist.GetNbinsX()+1)
    assert totWgt == totWgt2
    print("==========Total Weight ===============")
    print(totWgt,totWgt2)
    if 'zz4l-amcatnlo' in samples:
        xsec = 1.218
        kfac = 1.0835
        lumi = 59.7*1000
    if 'ggZZ4e' in samples:
        xsec = 0.001586
        kfac = 1.7
        lumi = 59.7*1000
    factor = xsec*kfac*lumi/totWgt
    
    for i in range(len(samples)):
        hunfa = fa.Get(samples[i]+"/"+unfname).Clone()
        hunfb = fb.Get(samples[i]+"/"+unfname).Clone()
      
        hunfa = rebin(hunfa,binnings)
        hunfb = rebin(hunfb,binnings)
        checkZeroBin(hunfb, 'denominator',hunfa)
        
        hunf_a_b = hunfa.Clone()
        hunf_amb = hunfa.Clone()
        hunf_a_b.Divide(hunfb)
        hunf_amb.Add(hunfb,-1)
        if i==0:
            num = hunfa.Integral(1,hunf_amb.GetNbinsX())*factor #only take amcnlo numerator and denominator
            den = hunfb.Integral(1,hunf_amb.GetNbinsX())*factor

        hists.append(hunf_a_b) #append in the orders of labels
    

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
        hists[i].GetYaxis().SetTitleOffset(1.5)
        hists[i].GetXaxis().SetTitle(prettyVar)
        hists[i].GetXaxis().SetTitleSize(0.04)
        hists[i].GetXaxis().SetTitleOffset(1.3)
        hists[i].SetStats(0)

    for i in range(len(hists)):
        hists[i].SetMaximum(max(maxs)*1.2)
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
            
            legend = r.TLegend (0.7 ,0.35 ,0.85 ,0.5)
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
    textbox_num= getTextBox(0.73,0.35,"num. %s"%round(num,2),0.03)
    textbox_den= getTextBox(0.73,0.3,"den. %s"%round(den,2),0.03)
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

pdfcommand.append('%s/TrigEff_plots_%s.pdf'%(dirName,suffix))
subprocess.call(pdfcommand)
