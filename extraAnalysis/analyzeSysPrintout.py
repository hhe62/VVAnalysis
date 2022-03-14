import ROOT as r
import pdb,subprocess,math,array
import sys,json,os
import numpy as np
from RooUnfoldBayes_reimplement import *

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

def makePlot(typesDiv,histsDiv,areasDiv,chan,c1,pdfcommand,plotlabel,folderName = 'SysDetailedPlots'):
    colors = [1,2,3,4,5,6]
    markers = [1,2,3,4,5,6]
    maxs = []
    mins = []
    for i in range(len(histsDiv)):
        histsDiv[i].SetLineColor(colors[i])
        maxs.append(histsDiv[i].GetMaximum())
        mins.append(histsDiv[i].GetMinimum())
        histsDiv[i].SetTitle("")
        histsDiv[i].GetYaxis().SetTitle('')
        histsDiv[i].GetYaxis().SetTitleOffset(1.5)
        histsDiv[i].GetXaxis().SetTitle(prettyVar+" "+ 'BinIndex')
        histsDiv[i].GetXaxis().SetTitleSize(0.04)
        histsDiv[i].GetXaxis().SetTitleOffset(1.3)
        histsDiv[i].SetStats(0)

    for i in range(len(histsDiv)):
        histsDiv[i].SetMaximum(max(maxs)*1.2)
        histsDiv[i].SetMinimum(min(mins))
        histsDiv[i].SetMarkerStyle(1)
        if i == 0:   
            #histsDiv[i].Draw("HIST P")
            if nostat:
                histsDiv[i].Draw("HIST")
            else:
                histsDiv[i].Draw()
            r.gStyle.SetLegendFont(42)
            r.gStyle.SetLegendTextSize(0.03)
            
            legend = r.TLegend (0.73 ,0.8 ,0.83 ,0.92)
        else:
            #histsDiv[i].SetMarkerStyle(markers[i])
            #histsDiv[i].Draw("HIST P SAME")
            if nostat:
                histsDiv[i].Draw("HIST SAME")
            else:
                histsDiv[i].Draw("SAME")
        legend.AddEntry(histsDiv[i],typesDiv[i]+' '+ str(round(areasDiv[i],2)))

    r.gStyle.SetLegendTextSize(0.018)
    r.gStyle.SetLegendFillColor(0)
    legend.SetLineWidth (0)
    legend.Draw("same")
    latex = r.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    #if "Full" in var:
    textbox= getTextBox(0.5,0.96,chan,0.03)
    textbox2= getTextBox(0.45,0.9,plotlabel,0.03)
    #latex.DrawLatex(0.74,0.83 ,"59.7fb^{-1}")
    
    if not os.path.isdir(folderName):
        os.mkdir(folderName)

    try:
        c1.SaveAs("%s/%s_%s_%s.png"%(folderName,var,chan,order))
    except:
        print("Problem saving plot.")
    c1.Clear()
    pdfcommand.append("%s/%s_%s_%s.png"%(folderName,var,chan,order))
    pdfcommand_global.append("%s/%s_%s_%s.png"%(folderName,var,chan,order))

def makeMatrixPlot(typesDiv,histsDiv,areasDiv,chan,c1,pdfcommand,plotlabel,folderName = 'SysDetailedPlots'):
    colors = [1,2,3,4,5,6]
    markers = [1,2,3,4,5,6]
    maxs = []
    mins = []
    for i in range(len(histsDiv)):
        histsDiv[i].SetLineColor(colors[i])
        maxs.append(histsDiv[i].GetMaximum())
        mins.append(histsDiv[i].GetMinimum())
        histsDiv[i].SetTitle("")
        histsDiv[i].GetYaxis().SetTitle('')
        histsDiv[i].GetYaxis().SetTitleOffset(1.5)
        histsDiv[i].GetXaxis().SetTitle(prettyVar+" "+ 'BinIndex')
        histsDiv[i].GetXaxis().SetTitleSize(0.04)
        histsDiv[i].GetXaxis().SetTitleOffset(1.3)
        histsDiv[i].SetStats(0)

    for i in range(len(histsDiv)):
        histsDiv[i].SetMaximum(max(maxs)*1.2)
        histsDiv[i].SetMinimum(min(mins))
        histsDiv[i].SetMarkerStyle(1)
        if i == 0:   
            #histsDiv[i].Draw("HIST P")
            if nostat:
                histsDiv[i].Draw("lego")
            else:
                histsDiv[i].Draw("lego")
            r.gStyle.SetLegendFont(42)
            r.gStyle.SetLegendTextSize(0.03)
            
            legend = r.TLegend (0.73 ,0.8 ,0.83 ,0.92)
        else:
            #histsDiv[i].SetMarkerStyle(markers[i])
            #histsDiv[i].Draw("HIST P SAME")
            if nostat:
                histsDiv[i].Draw("lego SAME")
            else:
                histsDiv[i].Draw("lego SAME")
        legend.AddEntry(histsDiv[i],typesDiv[i]+' '+ str(round(areasDiv[i],2)))

    r.gStyle.SetLegendTextSize(0.018)
    r.gStyle.SetLegendFillColor(0)
    legend.SetLineWidth (0)
    legend.Draw("same")
    latex = r.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    #if "Full" in var:
    textbox= getTextBox(0.5,0.96,chan,0.03)
    textbox2= getTextBox(0.45,0.9,plotlabel,0.03)
    #latex.DrawLatex(0.74,0.83 ,"59.7fb^{-1}")
    
    if not os.path.isdir(folderName):
        os.mkdir(folderName)

    try:
        c1.SaveAs("%s/%s_%s_%s.png"%(folderName,var,chan,order))
    except:
        print("Problem saving plot.")
    c1.Clear()
    pdfcommand.append("%s/%s_%s_%s.png"%(folderName,var,chan,order))
    pdfcommand_global.append("%s/%s_%s_%s.png"%(folderName,var,chan,order))

nostat = True
testUnf = True
pdfcommand=['convert']
pdfcommand_dmb=['convert']
pdfcommand_sig=['convert']
pdfcommand_truth=['convert']
pdfcommand_matrix = ['convert']
pdfcommand_global = ['convert']
commandlist = [pdfcommand,pdfcommand_dmb,pdfcommand_sig,pdfcommand_truth,pdfcommand_matrix,pdfcommand_global]

with open('varsFile.json') as var_json_file:
    myvar_dict = json.load(var_json_file)

varstr="nJets mjj dEtajj jetPt[0] jetPt[1] absjetEta[0] absjetEta[1] MassAllj Mass0j Mass1j Mass2j Mass3j Mass34j Mass4j MassFull Mass0jFull Mass1jFull Mass2jFull Mass3jFull Mass34jFull Mass4jFull"
varlist = varstr.split(' ')
var = sys.argv[2]#"Mass34jFull"

prettyVar = myvar_dict[var]["prettyVars"]

_yTitleTemp = '{prefix} \\frac{{d\\sigma_{{\\text{{fid}}}}}}{{d{xvar}}} {units}'

if myvar_dict[var]["units"]:
    yt = _yTitleTemp.format(xvar=prettyVar,
                                prefix='\\frac{1}{\\sigma_{\\text{fid}}}',
                                units='\\, \\left( \\frac{{1}}{{\\text{{{unit}}}}} \\right)'.format(unit=myvar_dict[var]["units"].replace("[","").replace("]","")))
else:
    yt = _yTitleTemp.format(prefix='\\frac{1}{\\sigma_{\\text{fid}}}',
                                xvar=prettyVar, units='')

labels=['2016','2017','2018']
channels =["eeee","eemm","mmmm"]
sysDict = {"eeee":{},"eemm":{},"mmmm":{}}
histDict = {"eeee":{},"eemm":{},"mmmm":{}}
for chan in sysDict.keys():
    sysDict[chan]['typeList'] = []
    sysDict[chan]['histList'] = []
    sysDict[chan]['areaList'] = []
for chan in histDict.keys():
    histDict[chan]['typeList'] = []
    histDict[chan]['dataList'] = []
    histDict[chan]['sigList'] = []
    histDict[chan]['bkgList'] = []
    histDict[chan]['truthList'] = []
    histDict[chan]['matrixList'] = []

unctype =''
channel = ''
recordbin = False

#first series of reading
fin = open(sys.argv[1])
for line in fin:
    if "Unc type before norm:" in line:
        unctype = line.strip().split('Unc type before norm:')[1]
        if not unctype:
            unctype = 'nominal'
        continue

    if 'Bin content before norm for channel ' in line:
        #pdb.set_trace()
        chanInd = line.strip().split('Bin content before norm for channel ')[1]
        if chanInd == '0':
            channel ='eeee'
        if chanInd == '1':
            channel ='eemm'
        if chanInd == '2':
            channel ='mmmm'
        recordbin = True
        continue

    if recordbin:
        bincontents=line.strip().split('[')[1].split(']')[0].split(', ')
        htmp = r.TH1F(unctype+channel,unctype+channel,len(bincontents),1,len(bincontents)+1)
        r.SetOwnership(htmp,False)
        for i in range(1,len(bincontents)+1):
            htmp.SetBinContent(i,float(bincontents[i-1]))
        if not 'PS' in unctype:
            sysDict[channel]['typeList'].append(unctype)
            sysDict[channel]['histList'].append(htmp)
            sysDict[channel]['areaList'].append(htmp.Integral(1,len(bincontents)))

        recordbin = False
fin.close()

#second series of reading
recordchannel = False
recordhist = False
recordmatrix = False
histtype = ''
fin = open(sys.argv[1])
for line in fin:
    if 'Printout record start here:' in line:
        recordchannel= True
        continue
    if recordchannel:
        channel=line.strip().split('channel:  ')[1]
        recordchannel=False
        continue
    if 'Position Indicator:' in line:
        unctype=line.strip().split('Position Indicator:')[1].replace(' ','')
        if not 'PS' in unctype:
            histDict[channel]['typeList'].append(unctype)
        continue
    if "Diagnostic bin contents of " in line:
        histtype = line.split('Diagnostic bin contents of ')[1].split(' ')[0]
        recordhist = True
        continue
    if recordhist:
        bincontents=line.strip().split('[')[1].split(']')[0].split(', ')
        htmp = r.TH1F(unctype+channel+histtype,unctype+channel+histtype,len(bincontents),1,len(bincontents)+1)
        #r.SetOwnership(htmp,False)
        for i in range(1,len(bincontents)+1):
            htmp.SetBinContent(i,float(bincontents[i-1]))
        if not 'PS' in unctype:
            histDict[channel][histtype+'List'].append(htmp)
        del htmp
        recordhist=False
        continue
    if 'Diagnostic matrix bin contents of response matrix' in line:
        histtype = 'matrix'
        recordmatrix = True
        continue
    if recordmatrix:
        bincontents=line.strip().split('[')[1].split(']')[0].split(', ')
        nbins = int((len(bincontents))**0.5)
        assert nbins == 10
        htmp = r.TH2F(unctype+channel+histtype,unctype+channel+histtype,nbins,1,nbins+1,nbins,1,nbins+1)
        indexmap = []
        for i in range(1,htmp.GetNbinsX()+1):
            for j in range(1,htmp.GetNbinsX()+1): 
                indexmap.append([i,j])
        for i in range(0,len(bincontents)):
            htmp.SetBinContent(indexmap[i][0],indexmap[i][1],float(bincontents[i]))
        if not 'PS' in unctype:
            histDict[channel][histtype+'List'].append(htmp)
        del htmp
        recordmatrix = False
        continue

#Making all the plots
c1=r.TCanvas("canvas")
c1.cd()
#c1.SetLogy()

r.gStyle.SetOptDate(0)
for chan in channels:
    types = sysDict[chan]['typeList']
    hists = sysDict[chan]['histList']
    areas = np.array(sysDict[chan]['areaList'])
    sortedInd = areas.argsort()

    types = [types[i] for i in sortedInd]
    hists = [hists[i] for i in sortedInd]
    areas = [areas[i] for i in sortedInd]
    
    assert len(types) == len(histDict[chan]['typeList'])
    sortedInd2 = [histDict[chan]['typeList'].index(types[i]) for i in range(0,len(types))]
    for prefix in ['type','data','sig','bkg','truth','matrix']:
        histDict[chan][prefix+'List'] = [histDict[chan][prefix+'List'][i] for i in sortedInd2]
    

    plotInd = [range(0,6), range(6,12),range(12,len(types))] #divide into groups
    #pdb.set_trace()
    for order,ranges in enumerate(plotInd):
        
        typesDiv = [types[i] for i in ranges]
        histsDiv = [hists[i] for i in ranges]
        areasDiv = [areas[i] for i in ranges]
        makePlot(typesDiv,histsDiv,areasDiv,chan,c1,pdfcommand,"Unfolded Results",folderName = 'SysDetailedPlots')

        divdict = {}
        for prefix in ['type','data','sig','bkg','truth','matrix']:
            divdict[prefix] = [histDict[chan][prefix+'List'][i] for i in ranges]
        divdict['dmb'] = [divdict['data'][i].Clone('dmb%s'%i) for i in range(0,len(divdict['data']))]
        for i in range(0,len(divdict['dmb'])):
            divdict['dmb'][i].Add(divdict['bkg'][i],-1)
        makePlot(divdict['type'],divdict['dmb'],areasDiv,chan,c1,pdfcommand_dmb,"Data minus bkg",folderName = 'DataMinusBkgPlots')    
        makePlot(divdict['type'],divdict['sig'],areasDiv,chan,c1,pdfcommand_sig,"Sig Plot",folderName = 'SignalPlots')
        makePlot(divdict['type'],divdict['truth'],areasDiv,chan,c1,pdfcommand_truth,"Truth Plot",folderName = 'TruthPlots')
        makeMatrixPlot(divdict['type'],divdict['matrix'],areasDiv,chan,c1,pdfcommand_matrix,"Resp Matrices",folderName = 'MatrixPlots')
        
        if testUnf and order == 0 and chan=='eeee':
            print('Testing '+ divdict['type'][0])
            #afile = r.TFile('ggZZxsec_down.root','RECREATE')
            #afile.cd()
            #divdict['sig'][0].Write()
            #divdict['truth'][0].Write()
            #divdict['matrix'][0].Write()
            #divdict['dmb'][0].Write()
            #sys.exit()
            testdmb = divdict['dmb'][0].Clone()
            dmb_printver = [-0.16745262345415665, -0.3142660111831682, -0.29620059239180124, 0.16902912781951573, -0.7315486655648783, 1.870785075540953, 2.6573732016988245, -0.020525411828250302, 0.014777588030264034, -0.013639418161901874]
            for i,x in enumerate(dmb_printver):
                testdmb.SetBinContent(i+1,x)
            denegative = False
            if denegative:
                for i in range(1,testdmb.GetNbinsX()+1):
                    if testdmb.GetBinContent(i) <0:
                        testdmb.SetBinContent(i,0.)

            respmatrix = RooUnfoldResponse(divdict['sig'][0],divdict['truth'][0],divdict['matrix'][0])
            unfold1 =RooUnfoldBayes(respmatrix, testdmb,4)
            hReco = unfold1.Hreco()
            #sys.exit()


        #pdb.set_trace()

        

pdfcommand.append('SysDetailedPlots'+"/SystematicPlots_%s.pdf"%var)
pdfcommand_dmb.append('DataMinusBkgPlots'+"/DataMinusBkgPlots_%s.pdf"%var)
pdfcommand_sig.append('SignalPlots'+"/SignalPlots_%s.pdf"%var)
pdfcommand_truth.append('TruthPlots'+"/TruthPlots_%s.pdf"%var)
pdfcommand_matrix.append('MatrixPlots'+"/MatrixPlots_%s.pdf"%var)
pdfcommand_global.append('AllPlots'+"/AllPlots_%s.pdf"%var)
if not os.path.isdir('AllPlots'):
        os.mkdir('AllPlots')
for comm in commandlist:
    subprocess.call(comm)

