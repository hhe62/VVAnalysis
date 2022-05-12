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

def checkZeroBin(hist,label):
    for i in range(1,hist.GetNbinsX()+1):
        if hist.GetBinContent(i)==0.:
            print("WARNING: %s contains 0 in bin %s"%(label,i))
        if hist.GetBinContent(i)<0.:
            print("WARNING: %s contains negative value in bin %s"%(label,i))

pdfcommand=['convert']
with open('varsFile.json') as var_json_file:
    myvar_dict = json.load(var_json_file)

nostat = False #no statistical error or not

varstr="nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1] MassAllj Mass0j Mass1j Mass2j Mass3j Mass34j Mass4j MassFull Mass0jFull Mass1jFull Mass2jFull Mass3jFull Mass34jFull Mass4jFull"
varlist = varstr.split(' ')
#varlist=["nJets"]
for var in varlist:
    prettyVar = myvar_dict[var]["prettyVars"]

    _yTitleTemp = '{prefix} \\frac{{d\\sigma_{{\\text{{fid}}}}}}{{d{xvar}}} {units}'

    if myvar_dict[var]["units"]:
        yt = _yTitleTemp.format(xvar=prettyVar,
                                    prefix='\\frac{1}{\\sigma_{\\text{fid}}}',
                                    units='\\, \\left( \\frac{{1}}{{\\text{{{unit}}}}} \\right)'.format(unit=myvar_dict[var]["units"].replace("[","").replace("]","")))
    else:
        yt = _yTitleTemp.format(prefix='\\frac{1}{\\sigma_{\\text{fid}}}',
                                    xvar=prettyVar, units='')

    fnames=[sys.argv[1]]
    labels=['Unfolded data/theo','RECO data/theo','MC shape RECO/Truth']
    hists=[]
    unfname='tot_%s_unf'%var
    truthname = 'tot_%s_true'%var
    dataname = 'tot_%s_data'%var
    MCname = 'tot_%s_SigMC'%var
    bkgname = 'tot_%s_bkg'%var
    #pdb.set_trace()

    f=r.TFile(fnames[0])
    hunf = f.Get(unfname).Clone()
    hdata = f.Get(dataname).Clone()
    hbkg = f.Get(bkgname).Clone()
    hMC = f.Get(MCname).Clone()
    htrue = f.Get(truthname).Clone()
    r.SetOwnership(f,False)
    checkZeroBin(htrue, 'Truth hist')
    checkZeroBin(hMC, 'RECO MC hist')
    hunf_d_t = hunf.Clone()
    hunf_d_t.Divide(htrue)
    hRECO_d_t = hdata.Clone()
    hRECO_d_t.Add(hbkg,-1)
    hRECO_d_t.Divide(hMC)
    
    hMC_R_t = hMC.Clone()
    areaMC = hMC.Integral(1,hMC.GetNbinsX())
    areatrue = htrue.Integral(1,htrue.GetNbinsX())
    hMC_R_t.Scale(1./areaMC)
    htrue_norm = htrue.Clone()
    htrue_norm.Scale(1./areatrue)
    hMC_R_t.Divide(htrue_norm)

    hists.append(hunf_d_t) #append in the orders of labels
    hists.append(hRECO_d_t)
    hists.append(hMC_R_t)
    
    colors = [3,2,4]
    markers = [1,2,3]
    maxs = []

    c1=r.TCanvas("canvas")
    c1.cd()
    #c1.SetLogy()
    if "Full" in var:
        c1.SetLogx()

    r.gStyle.SetOptDate(0)
    #pdb.set_trace()
    for i in range(len(hists)):
        hists[i].SetLineColor(colors[i])
        maxs.append(hists[i].GetMaximum())
        hists[i].SetTitle("")
        hists[i].GetYaxis().SetTitle("Ratio") #yt
        hists[i].GetYaxis().SetTitleOffset(1.5)
        hists[i].GetXaxis().SetTitle(prettyVar+" "+myvar_dict[var]["units"])
        hists[i].GetXaxis().SetTitleSize(0.04)
        hists[i].GetXaxis().SetTitleOffset(1.3)
        hists[i].SetStats(0)

    for i in range(len(hists)):
        hists[i].SetMaximum(max(maxs)*1.2)
        hists[i].SetMarkerStyle(1)
        if i == 0:   
            #hists[i].Draw("HIST P")
            if nostat:
                hists[i].Draw("HIST")
            else:
                hists[i].Draw()
            r.gStyle.SetLegendFont(42)
            r.gStyle.SetLegendTextSize(0.03)
            
            legend = r.TLegend (0.6 ,0.75 ,0.75 ,0.90)
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
    if 'nonreg' in sys.argv[1]:
        textbox_nonreg= getTextBox(0.7,0.9,"Matrix Inv",0.03)
    if "Full" in var:
        textbox= getTextBox(0.35,0.97,"Without 60 GeV < m_{Z1}, m_{Z2} < 120 GeV",0.03)
    #latex.DrawLatex(0.74,0.83 ,"59.7fb^{-1}")

    if not os.path.isdir('BottomLinePlots'):
        os.mkdir('BottomLinePlots')

    try:
        c1.SaveAs("BottomLinePlots/%s_BottomLineUnfolded.png"%(var))
    except:
        print("Problem saving plot.")
    c1.Clear()
    pdfcommand.append("BottomLinePlots/%s_BottomLineUnfolded.png"%(var))

pdfcommand.append('BottomLinePlots/BottomLine_plots.pdf')
subprocess.call(pdfcommand)
