import ROOT as r
import os,pdb,json,subprocess

def countXY(h,hrsp,xy): #value of xy: "x" or "y"
    h2 = h.Clone(h.GetName()+"Count"+xy)
    nbins = h.GetNbinsX()
    for i in range(1,nbins+1):
        if xy == "x":
            sumbin = sum([hrsp.GetBinContent(i,j) for j in range(1,nbins+1)])
            h2.SetBinContent(i,h.GetBinContent(i)-sumbin)
        
        if xy == "y":
            sumbin = sum([hrsp.GetBinContent(j,i) for j in range(1,nbins+1)])
            h2.SetBinContent(i,h.GetBinContent(i)-sumbin)

    return h2

def getLumiTextBox():
    texS = r.TLatex(0.65,0.96, str(59.7)+" fb^{-1} (13 TeV)")
    texS.SetNDC()
    texS.SetTextFont(42)
    texS.SetTextSize(0.045)
    texS.SetTextColor(r.kBlack)
    texS.Draw()
    texS1 = r.TLatex(0.15,0.96,"#bf{CMS}")
    texS1.SetNDC()
    texS1.SetTextFont(42)
    texS1.SetTextColor(r.kBlack)
    texS1.SetTextSize(0.045)
    texS1.Draw()

    texS2 = r.TLatex(0.25,0.96,"Preliminary")
    texS2.SetNDC()
    texS2.SetTextFont(52)
    texS2.SetTextColor(r.kBlack)
    texS2.SetTextSize(0.045)
    texS2.Draw()
    
    return texS,texS1,texS2

def redrawXaxis(h,varName):
    
    if "Full" in varName and "Mass" in varName:
            xaxis = r.TGaxis(h.GetXaxis().GetXmin(),h.GetMinimum(),h.GetXaxis().GetXmax(),h.GetMinimum(),h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax(),510,"G")
        
            xaxis.SetMoreLogLabels(True)
            xaxis.SetTickLength(0.03)
            #xaxis.SetLabelSize(0.025)
            xaxis.ChangeLabel(1,-1,0.,-1,-1,-1,"")
    elif varName == "nJets":
        xaxis = r.TGaxis(h.GetXaxis().GetXmin(),h.GetMinimum(),h.GetXaxis().GetXmax(),h.GetMinimum(),h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax(),505)
        xaxis.CenterLabels(True)
        xaxis.ChangeLabel(4,-1,-1,-1,-1,-1,"#geq 3")
    else:
        xaxis = r.TGaxis(h.GetXaxis().GetXmin(),h.GetMinimum(),h.GetXaxis().GetXmax(),h.GetMinimum(),h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax(),510)
    xaxis.SetTitle(prettyVars[varName]+''+units[varName])
    xaxis.SetLabelFont(42)
    xaxis.SetLabelOffset(0.01)
    #xaxis.SetTickLength(0.1)
    if "Full" in varName and "Mass" in varName:
        xaxis.SetLabelSize(0.04)
    else:
        xaxis.SetLabelSize(0.04)
    xaxis.SetTitleFont(42)
    xaxis.SetTitleSize(0.05)
    xaxis.SetTitleOffset(0.9)
    #if varName=="mass":
    #    xaxis.SetNoExponent(True)
    xaxis.Draw("SAME")
    #pdb.set_trace()
    return xaxis


def titleAndRatio(t,ra):
    tex = r.TLatex(0.43,0.9,t)
    tex.SetNDC()
    tex.SetTextFont(52)
    tex.SetTextColor(r.kBlue)
    tex.SetTextSize(0.035)
    tex.Draw()

    rtex = r.TLatex(0.52,0.73,ra)
    rtex.SetNDC()
    rtex.SetTextFont(52)
    rtex.SetTextColor(r.kBlue)
    rtex.SetTextSize(0.03)
    rtex.Draw()
    return tex,rtex

def plotHist(h1,h2,k1,k2):
    h1.GetYaxis().SetTitle("Events")
    h1.GetYaxis().SetTitleOffset(1.6)

    h1.Draw("HIST")
    h2.Draw("HIST SAME")
    legend = r.TLegend (0.5 ,0.75 ,0.9 ,0.90)
    #legend = r.TLegend (0.5 ,0.9 ,0.9 ,1.2)

    legend.SetFillColor(0)
    #legend.SetBorderSize(2)
    legend.AddEntry(h1,k1,"l")
    legend.AddEntry(h2,k2,"l")
    legend.SetTextSize(0.03)
    legend.SetLineWidth (0)
    legend.Draw("same")
    portion = h2.Integral(1,h2.GetNbinsX())/h1.Integral(1,h1.GetNbinsX())
    return portion,legend


varstr="nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1] MassAllj Mass0j Mass1j Mass2j Mass34j MassFull Mass0jFull Mass1jFull Mass2jFull Mass34jFull"
vars = varstr.split(" ")
#vars = ["MassFull"]
outdir = "2018FidPlots"
pdfcommand=['convert']
pdfcommand2=['convert']

if not os.path.isdir(outdir):
    os.mkdir(outdir)

with open('varsFile.json') as var_json_file:
    myvar_dict = json.load(var_json_file)
units = {}
prettyVars = {}
for key in myvar_dict.keys(): #key is the variable
    
    units[key] = myvar_dict[key]["units"]
    prettyVars[key] = myvar_dict[key]["prettyVars"]

fin = r.TFile("Fidcombined.root")
channels = ["eeee","eemm","mmmm"]
dict = {}
for var in vars:
    dict[var] = {}
    for chan in channels:
        hR = fin.Get("FidInfoRECO_%s_%s"%(var,chan))
        hRsp = fin.Get("FidInfoResp_%s_%s"%(var,chan))
        hT = fin.Get("FidInfoTruth_%s_%s"%(var,chan))
        if chan == channels[0]:
            dict[var]["Total"] = [hR.Clone("tot_RECO_%s"%var),hT.Clone("tot_Truth_%s"%var),hRsp.Clone("tot_Resp_%s"%var)]
        else:
            dict[var]["Total"][0].Add(hR)
            dict[var]["Total"][1].Add(hT)
            dict[var]["Total"][2].Add(hRsp)

        dict[var][chan] = [hR,hT,hRsp]

r.gStyle.SetOptDate(False)
r.gStyle.SetOptStat(0)
canvas_dimensions = [1000, 800]
c1 = r.TCanvas("c", "canvas",*canvas_dimensions)
c1.SetTopMargin(0.05)
c1.cd()
#
lineStyles = [1,2,1,2]
colors = [r.kOrange,2,r.kOrange,2]
for var in vars:
    for chan in channels + ["Total"]:
        hR,hT,hRsp = dict[var][chan]
        hdx = countXY(hR,hRsp,"x")
        hdy = countXY(hT,hRsp,"y")
        max1 = max(hdx.GetMaximum(),hR.GetMaximum())
        max2 = max(hdy.GetMaximum(),hT.GetMaximum())
        min1 = min(hdx.GetMinimum(),hR.GetMinimum(),0.)
        min2 = min(hdy.GetMinimum(),hT.GetMinimum(),0.)
        for i,h in enumerate([hdx,hR,hdy,hT]):
            if i<2:
                h.SetMaximum(1.2*max1)
                h.SetMinimum(min1)
            else:
                h.SetMaximum(1.2*max2)
                h.SetMinimum(min2)
            h.SetLineStyle(lineStyles[i])
            h.SetLineColor(colors[i])
            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTickLength(0)

        if "Full" in var:
            c1.SetLogx()
        else:
            c1.SetLogx(0)        
        c1.Divide(2,1)

        c1.cd(1)
        portionR,legend1 = plotHist(hR,hdx,"Total signal","Out of fiducial")
        t1,t2,t3 = getLumiTextBox()
        tex1,ratio1 = titleAndRatio("RECO Events","Out of fiducial portion %s"%(round(portionR,3)))
        xa1 = redrawXaxis(hdx,var)
        #pdb.set_trace()
        
        c1.cd(2)
        portionT,legend2 = plotHist(hT,hdy,"Total signal","Not reconstructed")
        t4,t5,t6 = getLumiTextBox()
        tex2,ratio2 = titleAndRatio("Truth Events","Non-reconstructed portion %s"%(round(portionT,3)))
        xa2 = redrawXaxis(hdy,var)
        
        c1.SaveAs(os.path.join(outdir,"%s_%s.png"%(var,chan)))

        if chan == "Total":
            pdfcommand.append(os.path.join(outdir,"%s_%s.png"%(var,chan)))
        else:
            pdfcommand2.append(os.path.join(outdir,"%s_%s.png"%(var,chan)))

        c1.Clear()

pdfcommand.append(os.path.join(outdir,"total.pdf"))      
pdfcommand2.append(os.path.join(outdir,"channels.pdf"))   
subprocess.call(pdfcommand)
subprocess.call(pdfcommand2)



