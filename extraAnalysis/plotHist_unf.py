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

pdfcommand=['convert']
with open('varsFile.json') as var_json_file:
    myvar_dict = json.load(var_json_file)

nostat = False #no statistical error or not

varstr="nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1] MassAllj Mass0j Mass1j Mass2j Mass3j Mass34j Mass4j MassFull Mass0jFull Mass1jFull Mass2jFull Mass3jFull Mass34jFull Mass4jFull"
varlist = varstr.split(' ')

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

    fnames=['2016Full_MassFull_fixReplicas.root','2017Full_MassFull.root','2018Full_MassFull.root']
    labels=['2016','2017','2018']
    hists=[]
    histname='tot_%s_unf'%var
    #pdb.set_trace()
    for fn in fnames: 
        f=r.TFile(fn)
        htmp = f.Get(histname).Clone()
        htmp.Scale(1./htmp.Integral(1,htmp.GetNbinsX()))
        for i in range(1,htmp.GetNbinsX()+1):
            htmp.SetBinContent(i,htmp.GetBinContent(i)/htmp.GetBinWidth(i)) 
            htmp.SetBinError(i,htmp.GetBinError(i)/htmp.GetBinWidth(i))
        r.SetOwnership(f,False)
        hists.append(htmp)
    


    colors = [3,2,4]
    markers = [1,2,3]
    maxs = []

    c1=r.TCanvas("canvas")
    c1.cd()
    c1.SetLogy()
    if "Full" in var:
        c1.SetLogx()

    r.gStyle.SetOptDate(0)
    #pdb.set_trace()
    for i in range(len(hists)):
        hists[i].SetLineColor(colors[i])
        maxs.append(hists[i].GetMaximum())
        hists[i].SetTitle("")
        hists[i].GetYaxis().SetTitle(yt)
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
            
            legend = r.TLegend (0.8 ,0.75 ,0.95 ,0.90)
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
    if "Full" in var:
        textbox= getTextBox(0.35,0.9,"Without 60 GeV < m_{Z1}, m_{Z2} < 120 GeV",0.03)
    #latex.DrawLatex(0.74,0.83 ,"59.7fb^{-1}")

    if not os.path.isdir('OverlayPlots'):
        os.mkdir('OverlayPlots')

    try:
        c1.SaveAs("OverlayPlots/%s_overlayUnfolded.png"%(var))
    except:
        print("Problem saving plot.")
    c1.Clear()
    pdfcommand.append("OverlayPlots/%s_overlayUnfolded.png"%(var))

pdfcommand.append('OverlayPlots/overlay_plots.pdf')
subprocess.call(pdfcommand)
