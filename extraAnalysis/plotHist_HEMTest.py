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
    contents =[hist.GetBinContent(i) for i in range(1,hist.GetNbinsX()+1)]
    print(contents)
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
difflist = []
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

    fnames=[sys.argv[1],sys.argv[2]]
    labels=['Unfolded after/before']
    hists=[]
    unfname='tot_%s_unf'%var
    #pdb.set_trace()
    fa=r.TFile(fnames[0])
    fb=r.TFile(fnames[1])
    r.SetOwnership(fa,False)
    r.SetOwnership(fb,False)
    hunfa = fa.Get(unfname)#.Clone()
    hunfb = fb.Get(unfname)#+";1").Clone()
    
    checkZeroBin(hunfb, 'Unfolded hist before')
    checkZeroBin(hunfa, 'Unfolded hist after')
    
    hunf_a_b = hunfa.Clone()
    hunfaNorm = hunfa*(1./hunfa.Integral(1,hunfa.GetNbinsX()))
    hunfbNorm = hunfb*(1./hunfb.Integral(1,hunfa.GetNbinsX()))
    hunf_amb = hunfaNorm.Clone()
    hunf_a_b.Divide(hunfb)
    hunf_amb.Add(hunfbNorm,-1)
    diff = 0.
    for dbin in range(1,hunf_amb.GetNbinsX()+1):
        diff += abs(hunf_amb.GetBinContent(dbin)) #.Integral(1,hunf_amb.GetNbinsX())
    diff = round(diff*100,3)
    difflist.append(diff)
    print("Difference:{}%".format(diff))

    hists.append(hunf_a_b) #append in the orders of labels
    
    colors = [3]
    markers = [1]
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
    #if 'nonreg' in sys.argv[1]:
    textbox_diff= getTextBox(0.7,0.9,"diff {}%".format(diff),0.03)
    if "Full" in var:
        textbox= getTextBox(0.35,0.97,"Without 60 GeV < m_{Z1}, m_{Z2} < 120 GeV",0.03)
    #latex.DrawLatex(0.74,0.83 ,"59.7fb^{-1}")

    if not os.path.isdir('HEMPlots'):
        os.mkdir('HEMPlots')

    try:
        c1.SaveAs("HEMPlots/%s_HEMUnfolded.png"%(var))
    except:
        print("Problem saving plot.")
    c1.Clear()
    pdfcommand.append("HEMPlots/%s_HEMUnfolded.png"%(var))

pdfcommand.append('HEMPlots/HEM_plots.pdf')
subprocess.call(pdfcommand)
print(difflist,max(difflist),min(difflist))