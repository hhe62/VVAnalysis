#!/usr/bin/env python
import ROOT
import glob
import math
#from python import ConfigureJobs
from python import UserInput
import makeSimpleHtml
import os
import subprocess
import sys
import datetime
import array,json,pdb
from ROOT import vector as Vec


VFloat = Vec('float')
#from PlotTools import PlotStyle as Style, pdfViaTex
include_MiNNLO = True
EW_corr = True

#Divison and ticklength for redrawn y axes
yrdiv = 503
yrtl = 0.02

crossDrawOpt = "PESAME"
#style = Style()
#ROOT.gStyle.SetLineScalePS(1.8)
ROOT.gStyle.SetOptDate(False)
ROOT.gStyle.SetLineWidth(3)
#channels = ["eeee","eemm","mmmm"]
channels = []
def getComLineArgs():
    parser = UserInput.getDefaultParser()
    parser.add_argument("--lumi", "-l", type=float,
        default=41.5, help="luminosity value (in fb-1)")
    parser.add_argument('--thesis', action='store_true',
                        help='For making CMS thesis plots')
    parser.add_argument('--preliminary', action='store_true',
                        help='For making CMS Preliminary plots')
    parser.add_argument("--legend_left", action="store_true",
                        help="Put legend left or right")
    parser.add_argument("--titleOffset", type=float, default=1.0,
                        help="Scale default ymax by this amount")
    parser.add_argument("--scaleymax", type=float, default=1.0,
                        help="Scale default ymax by this amount")
    parser.add_argument("--scaleymin", type=float, default=1.0,
                        help="Scale default ymin by this amount")
    parser.add_argument("--scalelegy", type=float, default=1.0,
                        help="Scale default legend entry size by this amount")
    parser.add_argument("--scalelegx", type=float, default=1.0,
                        help="Scale default legend entry wdith by this amount")
    parser.add_argument("--output_file", "-o", type=str,
        default="", help="Output file name")
    parser.add_argument("--test", action='store_true',
        help="Run test job (no background estimate)")
    parser.add_argument("--variable", "-vr", type=str,
        default="all", help="variableName")
    parser.add_argument('--noNorm', action='store_true',
                        help='Leave differential cross sections in abolute normalization rather than normalizing to unity area.')
    parser.add_argument('--NormFb', action='store_true',
                        help='Normalize differential cross sections to the luminosity')
    parser.add_argument('--logy', '--logY', '--log', action='store_true',
                        help='Put vertical axis on a log scale.')
    parser.add_argument('--makeTotals', action='store_true',
                        help='plot total unfolded with uncertainities.')
    parser.add_argument('--unfoldDir', type=str, nargs='?',
                        default='/afs/cern.ch/user/h/hehe/www/ZZFullRun2/PlottingResults/ZZ4l2018/ZZSelectionsTightLeps/ANPlots/ZZ4l2018/DiffDist_FullRun2PreApproval_09Nov2019/',
                        help='Directory to put response and covariance plots in')
    return vars(parser.parse_args())

args = getComLineArgs()

today = datetime.date.today().strftime("%d%b%Y")

if "Full" in args['variable']:
    include_MiNNLO = False
#manager_path = ConfigureJobs.getManagerPath()
#Only MassAllj should plot EWK correction
EW_P4 = ("MassAllj" in args['variable'] and not "Full" in args['variable'] ) #or (args['variable'] == "nJets")

#Currently don't plot EWK for non-m4l jet variables
if not EW_P4:
    EW_corr = False
#Bottom margin for bottommost panel
if (EW_corr and EW_P4):
    bmg = 0.4
    bpwidth = 0.09
else:
    bmg = 0.45
    bpwidth = 0.11
lastbpw = round(bpwidth/(1.-bmg),2)
#y direction separation points of the panels
if not include_MiNNLO:
    #panel_breakpoints = [0.03,0.20,0.33]
    panel_breakpoints = [0.01,0.01+lastbpw,round(0.01+lastbpw+bpwidth,2)]
else:
    panel_breakpoints = [0.01,0.01+lastbpw] + [round(0.01+lastbpw+bpwidth*bpind,2) for bpind in range(1,3)]
    if EW_P4 and EW_corr:
        panel_breakpoints = [0.01,0.01+lastbpw] + [round(0.01+lastbpw+bpwidth*bpind,2) for bpind in range(1,4)]
        #panel_breakpoints = [0.01,0.15,0.227,0.304,0.381]
        #panel_breakpoints = [0.01,0.15,0.22,0.28,0.37]

pbpts = panel_breakpoints

analysis=args['analysis']
_binning = {
    'pt' : [25.*i for i in range(4)] + [100., 150., 200., 300.],
    #'mass' : [100.+100.*i for i in range(12)],
    'mass' : [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800.,1000.],
    'massFull' : [80.,100.,120.,130.,150.,180.,200.,240.,300.,400.,1000],
    'eta' : [6,0.,6.],
    'zmass' : [60., 80., 84., 86.] + [87.+i for i in range(10)] + [98., 102., 120.], #[12, 60., 120.],
    'z1mass' : [60., 80., 84., 86.] + [87.+i for i in range(10)] + [98., 102., 120.], #[12, 60., 120.],
    'z2mass' : [60., 75., 83.] + [84.+i for i in range(14)] + [105., 120.],#[12, 60., 120.],
    'z1pt' : [i * 25. for i in range(7)] + [200., 300.],
    'z2pt' : [i * 25. for i in range(7)] + [200., 300.],
    'zpt' : [i * 25. for i in range(7)] + [200., 300.],
    'zHigherPt' : [i * 25. for i in range(7)] + [200., 300.],
    'zLowerPt' : [i * 25. for i in range(7)] + [200., 300.],
    'leppt' : [i * 15. for i in range(11)],
    'l1Pt' : [0.,15.,30.,40.,50.]+[60.+15.*i for i in range(9)]+[195.,225.],#[14,0.,210.],#[15, 0., 150.],
    'dphiz1z2': [0.,1.5,2.0,2.25,2.5,2.75,3.0,3.25],
    'drz1z2': [0.,1.0,2.0,3.0,4.0,5.0,6.0]
    }
units = {
    'pt' : '[GeV]',
    'mass' : '[GeV]',
    'massFull' : '[GeV]',
    'eta' : '',
    'zmass' : '[GeV]',
    'z1mass' : '[GeV]',
    'z2mass' : '[GeV]',
    'zpt' : '[GeV]',
    'z1pt' : '[GeV]',
    'z2pt' : '[GeV]',
    'zHigherPt' : '[GeV]',
    'zLowerPt' : '[GeV]',
    'leppt' : '[GeV]',
    'l1Pt' : '[GeV]',
    'dphiz1z2': '',
    'drz1z2':'',
    }

yaxisunits = {
    'pt' : 'GeV',
    'mass' : 'GeV',
    'massFull' : 'GeV',
    'eta' : '',
    'zmass' : 'GeV',
    'z1mass' : 'GeV',
    'z2mass' : 'GeV',
    'zpt' : 'GeV',
    'z1pt' : 'GeV',
    'z2pt' : 'GeV',
    'zHigherPt' : 'GeV',
    'zLowerPt' : 'GeV',
    'leppt' : 'GeV',
    'l1Pt' : 'GeV',
    'dphiz1z2': '',
    'drz1z2':'',
    }

prettyVars = {
    'pt' : 'p_{T}^{4\\ell}',
    'mass' : 'm_{4\\ell}',
    'massFull' : 'm_{4\\ell}',
    'eta' : '\\eta_{4\\ell}',
    'zmass' : 'm_{Z}',
    'z1mass' : 'm_{Z_{1}}',
    'z2mass' : 'm_{Z_{2}}',
    'z1pt' : 'p_{T}^{Z_{1}}',
    'z2pt' : 'p_{T}^{Z_{2}}',
    'zpt' : 'p_{T}^{Z}',
    'zHigherPt' : 'p_\\text{T}^{\\text{Z}_{\\text{lead}}}',
    'zLowerPt' : 'p_\\text{T}^{\\text{Z}_{\\text{sublead}}}',
    'leppt' : 'p_{T}^{\\ell}',
    'l1Pt' : 'p_\\text{T}^{\\ell_1}', 
    'dphiz1z2': '\\Delta\\phi_{Z_{1},Z_{2}}',
    'drz1z2':'\\Delta\\text{R}_{Z_{1},Z_{2}}}',
    }

with open('varsFile.json') as var_json_file:
    myvar_dict = json.load(var_json_file)

ytilt_fac = myvar_dict[args['variable']]['ytilt_fac']
# list of variables not counting systematic shifts
varList=['Mass','ZZPt','ZPt','LepPt','dPhiZ1Z2','dRZ1Z2'] #With original list, histograms will be searched for all variables regardless of whether they are in runVariables
varNames={'mass': 'Mass','pt':'ZZPt','zpt':'ZPt','leppt':'LepPt','dphiz1z2':'dPhiZ1Z2','drz1z2':'dRZ1Z2'}

for key in myvar_dict.keys(): #key is the variable
    _binning[key] = myvar_dict[key]["_binning"]
    units[key] = myvar_dict[key]["units"]
    prettyVars[key] = myvar_dict[key]["prettyVars"]
    yaxisunits[key] = myvar_dict[key]["units"].replace("[","").replace("]","")
#    responseClassNames[key] = {c:myvar_dict[key]["responseClassNames"] for c in channels}
    if key == "MassAllj":
        #varNamesForResponseMaker[key] = {c:'Mass' for c in channels}
        varNames[key] = 'Mass'
    else:
        #varNamesForResponseMaker[key] = {c:str(key) for c in channels}
        varList.append(str(key))
        varNames[key] = str(key)


_xTitle = {}
_yTitle = {}
_yTitleNoNorm = {}

_yTitleTemp = '{prefix} \\frac{{d\\sigma_{{\\text{{fid}}}}}}{{d{xvar}}} {units}'
for var, prettyVar in prettyVars.iteritems():
    xt = prettyVar
    if yaxisunits[var]:
        xt += ' \\, \\left(\\text{{{}}}\\right)'.format(yaxisunits[var])
        yt = _yTitleTemp.format(xvar=prettyVar,
                                prefix='\\frac{1}{\\sigma_{\\text{fid}}}',
                                units='\\, \\left( \\frac{{1}}{{\\text{{{unit}}}}} \\right)'.format(unit=yaxisunits[var]))
        ytnn = _yTitleTemp.format(xvar=prettyVar, prefix='',
                                  units='\\, \\left( \\frac{{\\text{{fb}}}}{{\\text{{{unit}}}}} \\right)'.format(unit=yaxisunits[var]))
    else:
        yt = _yTitleTemp.format(prefix='\\frac{1}{\\sigma_{\\text{fid}}}',
                                xvar=prettyVar, units='')
        ytnn = _yTitleTemp.format(prefix='', xvar=prettyVar, units='\\left( \\text{fb} \\right)')

    _xTitle[var] = xt
    _yTitle[var] = yt
    if var == "dEtajj":
        _yTitle[var] = _yTitleTemp.format(xvar="\\vert \\Delta \\eta (j_1,j_2) \\vert",prefix='\\frac{1}{\\sigma_{\\text{fid}}}', units='')
    if var == "mjj":
        _yTitle[var] = _yTitleTemp.format(xvar="m_{jj}",prefix='\\frac{1}{\\sigma_{\\text{fid}}}', units='\\, \\left( \\frac{{1}}{{\\text{{{unit}}}}} \\right)'.format(unit='GeV'))
    _yTitleNoNorm[var] = ytnn

# list of variables not counting systematic shifts
#varList=['Mass','ZZPt','ZPt','LepPt','dPhiZ1Z2','dRZ1Z2']

# Sometimes need to more or resize legend
legDefaults = {
    'textsize' : 0.034, #.027,#2,
    'leftmargin' : 0.35,
    'entryheight' : 0.037,
    'rightmargin' : 0.03,
    }
legParams = {v.lower():legDefaults.copy() for v in varList}
legParams['z1mass'] = {
    'textsize' : .026,
    'leftmargin' : .03,
    'rightmargin' : .46,
    'entryheight' : .034,#23
    'entrysep' : .007,
    }
legParams['pt'] = legParams['zzpt'].copy()
legParams['zmass'] = legParams['z1mass'].copy()
legParams['z2mass'] = legParams['z1mass'].copy()
legParams['deltaEtajj'] = legParams['z1mass'].copy()
legParams['deltaEtajj']['leftmargin'] = .5
legParams['deltaEtajj']['rightmargin'] = .03
legParams['deltaEtajj']['topmargin'] = .05
legParams['eta'] = legParams['deltaEtajj'].copy()
#legParams['massFull']['leftmargin'] = 0.25

legParamsLogy = {v:p.copy() for v,p in legParams.iteritems()}
#legParamsLogy['l1Pt']['topmargin'] = 0.65
#legParamsLogy['l1Pt']['leftmargin'] = 0.2
#legParamsLogy['l1Pt']['rightmargin'] = 0.18
for key in legParamsLogy.keys():
    if "mass" in key.lower():
        print("legParamsLogy Key= %s=============================="%key)
        legParamsLogy[key]['topmargin'] = 0.075
        legParamsLogy[key]['leftmargin'] = 0.35
        legParamsLogy[key]['rightmargin'] = 0.025
        legParamsLogy[key]['textsize'] = 0.033
legParamsLogy['leppt']['topmargin'] = 0.05
#legParamsLogy['zHigherPt']['topmargin'] = 0.045
#legParamsLogy['massFull']['topmargin'] = 0.035

def normalizeBins(h):
    binUnit = 1 # min(h.GetBinWidth(b) for b in range(1,len(h)+1))
    for ib in range(1,h.GetNbinsX()+1):
        w = h.GetBinWidth(ib)
        h.SetBinContent(ib, h.GetBinContent(ib) * binUnit / w)
        h.SetBinError(ib, h.GetBinError(ib) * binUnit / w)
        if h.GetBinError(ib) > h.GetBinContent(ib):
            h.SetBinError(ib, h.GetBinContent(ib))
    h.Sumw2()

def unnormalizeBins(h):
    binUnit = 1 # min(h.GetBinWidth(b) for b in range(1,len(h)+1))
    for ib in range(1,h.GetNbinsX()+1):
        w = h.GetBinWidth(ib)
        h.SetBinContent(ib, h.GetBinContent(ib) * w / binUnit)
        h.SetBinError(ib, h.GetBinError(ib) * w / binUnit)
        if h.GetBinError(ib) > h.GetBinContent(ib):
            h.SetBinError(ib, h.GetBinContent(ib))
    h.Sumw2()

def createRatio(h1, h2):
    Nbins = h1.GetNbinsX()
    Ratio = h1.Clone("Ratio")
    hStackLast = h2.Clone("hStackLast")
    try:
        Ratio.Sumw2()
    except AttributeError:
        pass
    try:
        hStackLast.Sumw2()
    except AttributeError:
        pass
    for i in range(1,Nbins+1):
        stackcontent = hStackLast.GetBinContent(i)
        stackerror = hStackLast.GetBinError(i)
        datacontent = h1.GetBinContent(i)
        dataerror = h1.GetBinError(i)
        #print "bin: ",i
        #print "stackcontent: ",stackcontent," and data content: ",datacontent
        ratiocontent=0
        if(datacontent!=0):
            ratiocontent = datacontent/stackcontent
        if(datacontent!=0):
            error = ratiocontent*(math.sqrt(math.pow((dataerror/datacontent),2) + math.pow((stackerror/stackcontent),2)))
        else:
            error = 0. #set to 0 just in case, 0 bin cont shouldn't be drawn with PE1 option #2.07 #why 2.07?
        #print "ratio content: ",ratiocontent
        #print "stat error: ", error
        Ratio.SetBinContent(i,ratiocontent)
        Ratio.SetBinError(i,error)

    with open('varsFile.json') as var_json_file:
            myvar_dict = json.load(var_json_file)

    global my_varName
    Ratio.GetYaxis().SetRangeUser(myvar_dict[my_varName]['ratio_min'],myvar_dict[my_varName]['ratio_max'])
    #Ratio.GetYaxis().SetRangeUser(0.4,1.8)
    Ratio.SetStats(0)
    Ratio.GetYaxis().CenterTitle()
    Ratio.SetMarkerStyle(20)
    Ratio.SetMarkerSize(0.7)

    line = ROOT.TLine(h1.GetXaxis().GetXmin(), 1.,h1.GetXaxis().GetXmax(), 1.)
    line.SetLineStyle(7)

    Ratio.GetYaxis().SetLabelSize(0.14) #0.14
    Ratio.GetYaxis().SetTitleSize(0.16)
    Ratio.GetYaxis().SetLabelFont(42)
    Ratio.GetYaxis().SetTitleFont(42)
    Ratio.GetYaxis().SetTitleOffset(0.25)
    Ratio.GetYaxis().SetNdivisions(100)
    Ratio.GetYaxis().SetTickLength(0.0) #0.05

    Ratio.GetXaxis().SetLabelSize(0)
    Ratio.GetXaxis().SetTitleSize(0)

    return Ratio,line

def getPrettyLegend(hTrue, data_hist, hAltTrue, error_hist, coords,hTrueNNLO=None,hTrueEWC = None):
    tmpshift = 0.05
    if "MassAllj" in hTrue.GetName():
        tmpshift = 0.03
    elif "Mass" in hTrue.GetName():
        tmpshift = 0.04
    legend = ROOT.TLegend(coords[0]+tmpshift, coords[1]+0.01, coords[2]+tmpshift, coords[3]-0.01)
    ROOT.SetOwnership(legend, False)
    legend.SetName("legend")
    legend.SetFillStyle(0)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.047) #0.033 #0.025
    if "Full" in hTrue.GetName():
        legend.SetTextSize(0.042)    
    legend.SetTextColor(ROOT.kBlack)
    with open('listFile.json') as list_json_file:
        mylist_dict = json.load(list_json_file)
    sigLabel = mylist_dict["sigLabel"] #"POWHEG+MCFM+Pythia8" 
    sigLabelAlt = mylist_dict["sigLabelAlt"] #"MG5_aMC@NLO+MCFM+Pythia8"
    if data_hist:
        legend.AddEntry(data_hist, "Data + stat. unc.", "lep")
    legend.AddEntry(error_hist, "Stat. #oplus syst. unc.", "f")
    legend.AddEntry(hTrue, sigLabel,"lep")
    legend.AddEntry(hAltTrue, sigLabelAlt,"lep")
    if include_MiNNLO:
        legend.AddEntry(hTrueNNLO, "nNNLO+PS","lep")   
        if EW_corr:
            legend.AddEntry(hTrueEWC, "(nNNLO+PS)#times K_{EW}","lep")  
            #legend.AddEntry(hTrueEWC, "(nNNLO+PS) no GenWgt","lep")  

    #legend.AddEntry(hTrue, sigLabel,"lf")
    #legend.AddEntry(hAltTrue, sigLabelAlt,"l")
    return legend
def createCanvasPads(varName):
    #if matrix included
   # canvas_dimensions = [1000, 1400]
    canvas_dimensions = [1000, 1300]
    c = ROOT.TCanvas("c", "canvas",*canvas_dimensions)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetLegendBorderSize(0)
    # Upper histogram plot is pad1
    pad1 = ROOT.TPad("pad1", "pad1", 0.01, pbpts[-1], 0.99, 0.99)
    pad1.Draw()
    pad1.cd()
    if varName!="drz1z2":
        pad1.SetLogy()
    if "Full" in varName:
        pad1.SetLogx()
    pad1.SetFillColor(0)
    pad1.SetFrameLineWidth(3)
    pad1.SetFrameBorderMode(0)
    pad1.SetBorderMode(0)
    pad1.SetBottomMargin(0)  # joins upper and lower plot
    #pad1.SetGridx()
    #pad1.Draw()
    return c,pad1

def createPad2(canvas):
    # Lower ratio plot is pad2
    canvas.cd()  # returns to main canvas before defining pad2
    canvas.GetListOfPrimitives().SetOwner(True)
    
    pad2 = ROOT.TPad("pad2", "pad2", 0.01, pbpts[-2], 0.99, pbpts[-1])
   
    pad2.Draw()
    pad2.cd()
    pad2.SetFillColor(0)
    pad2.SetFrameLineWidth(3)
    pad2.SetFrameBorderMode(0)
    pad2.SetBorderMode(0)#bordermode = -1 box looks as it is behind the screen
   # bordermode = 0 no special effects
   # bordermode = 1 box looks as it is in front of the screen
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0)
    if "Full" in varName:
        pad2.SetLogx()

    #pad2.SetGridx()
    #pad2.Draw()
    return pad2

def createPad3(canvas):
    # Lower ratio plot is pad3
    canvas.cd()  # returns to main canvas before defining pad3
    
    pad3 = ROOT.TPad("pad3", "pad3", 0.01, pbpts[-3], 0.99, pbpts[-2])
    
    pad3.Draw()
    pad3.cd()
    pad3.SetFillColor(0)
    pad3.SetFrameLineWidth(3)
    pad3.SetFrameBorderMode(0)
    #pad3.SetFrameFillStyle(4000)
    if not include_MiNNLO:
        pad3.SetBorderMode(0)
        pad3.SetBottomMargin(bmg)
    else:
        pad3.SetBorderMode(0)
        pad3.SetBottomMargin(0)
    pad3.SetTopMargin(0)  # joins upper and lower plot
    
    #pad3.SetBottomMargin(0)
    if "Full" in varName:
        pad3.SetLogx()

    #pad3.SetGridx()
    #pad3.Draw()
    return pad3

def createPad4(canvas): #For MiNNLO
    # Lower ratio plot is pad4
    canvas.cd()  # returns to main canvas before defining pad4
    pad4 = ROOT.TPad("pad4", "pad4", 0.01, pbpts[-4], 0.99, pbpts[-3])
    pad4.Draw()
    pad4.cd()
    pad4.SetFillColor(0)
    pad4.SetFrameLineWidth(3)
    pad4.SetFrameBorderMode(0)
    #pad4.SetFrameFillStyle(4000)
    pad4.SetBorderMode(0)
    pad4.SetTopMargin(0)  # joins upper and lower plot
    if not (EW_P4 and EW_corr):
        pad4.SetBottomMargin(bmg)
    else:
        pad4.SetBottomMargin(0.)
    if "Full" in varName:
        pad4.SetLogx()

    #pad4.SetGridx()
    #pad4.Draw()
    return pad4

def createPad5(canvas): #For MiNNLO EWK m4l
    # Lower ratio plot is pad5
    canvas.cd()  # returns to main canvas before defining pad5
    pad5 = ROOT.TPad("pad5", "pad5", 0.01, pbpts[-5], 0.99, pbpts[-4])
    pad5.Draw()
    pad5.cd()
    pad5.SetFillColor(0)
    pad5.SetFrameLineWidth(3)
    pad5.SetFrameBorderMode(0)
    #pad5.SetFrameFillStyle(4000)
    pad5.SetBorderMode(0)
    pad5.SetTopMargin(0.)  # joins upper and lower plot
    
    pad5.SetBottomMargin(bmg) #This is currently the bottommost pad
  
    if "Full" in varName:
        pad5.SetLogx()

    #pad5.SetGridx()
    #pad5.Draw()
    return pad5

#rebin histos and take care of overflow bins
def rebin(hist,varName): #didn't handle error, but this function not actually used
    ROOT.SetOwnership(hist, False)
    #No need to rebin certain variables but still might need overflow check
    if varName not in ['eta']:
        bins=array.array('d',_binning[varName])
        Nbins=len(bins)-1 
        hist=hist.Rebin(Nbins,"",bins)
    else:
        Nbins = hist.GetSize() - 2
    add_overflow = hist.GetBinContent(Nbins) + hist.GetBinContent(Nbins + 1)
    hist.SetBinContent(Nbins, add_overflow)
    hist.SetBinContent(Nbins+1,0)
    return hist

#Draw Y axis with ticks
def getRYaxis(hUnf1,ratioErrorBand1,lastP):

    Ryaxis = ROOT.TGaxis(hUnf1.GetXaxis().GetXmax(),ratioErrorBand1.GetMinimum(),hUnf1.GetXaxis().GetXmax(),ratioErrorBand1.GetMaximum(),ratioErrorBand1.GetMinimum(),ratioErrorBand1.GetMaximum(),3,"+CS")
    Ryaxis.SetNdivisions(yrdiv)
    if lastP:
        Ryaxis.SetTickLength(yrtl/(1-bmg))
    else:
        Ryaxis.SetTickLength(yrtl)

    #axText3=getAxisTextBox(0.06,0.0,"Data/%s"%ratioName_alt,0.23,True)
    #Ryaxis.SetTitle("#scale[1.2]{Data/%s}"%ratioName_alt)
    #Ryaxis.SetTitle("Data/Theo.")
    Ryaxis.SetLabelFont(42)
    Ryaxis.SetLabelOffset(0.025) #0.01
    Ryaxis.SetLabelSize(0.) #0.1485
    Ryaxis.SetTitleFont(42)
    Ryaxis.SetTitleSize(0.)
    Ryaxis.SetTitleOffset(0.365)
    Ryaxis.Draw("SAME")

    return Ryaxis
    
def getLumiTextBox():
    texS = ROOT.TLatex(0.68,0.965, str(int(round(args['lumi'])))+" fb^{-1} (13 TeV)")
    texS.SetNDC()
    texS.SetTextFont(42)
    texS.SetTextSize(0.045)
    texS.SetTextColor(ROOT.kBlack)
    texS.Draw()
    texS1 = ROOT.TLatex(0.14,0.96,"#bf{CMS}")
    texS1.SetNDC()
    texS1.SetTextFont(42)
    texS1.SetTextColor(ROOT.kBlack)
    texS1.SetTextSize(0.1)
    texS1.Draw()

    texS2 = ROOT.TLatex(0.23,0.965,"Preliminary")
    texS2.SetNDC()
    texS2.SetTextFont(52)
    texS2.SetTextColor(ROOT.kBlack)
    texS2.SetTextSize(0.045)
    #texS2.Draw()
    return texS,texS1,texS2

def getSigTextBox(x,y,sigLabel,size): #check whether actually used
    if sigLabel=="POWHEG+MCFM+Pythia8":
        texS = ROOT.TLatex(x,y, "#bf{POWHEG+MCFM+Pythia8}")
    else:
        texS = ROOT.TLatex(x,y, "#bf{MG5_aMC@NLO+MCFM+Pythia8}")
    texS.SetNDC()
    texS.SetTextFont(42)
    texS.SetTextSize(size)
    texS.Draw()
    #return texS 
    #doesn't work without this last line. Uncomment it if want to use this function

def getAxisTextBox(x,y,axisLabel,size,rotated):
    texS = ROOT.TLatex(x,y,axisLabel)
    texS.SetNDC()
    #rotate for y-axis                                                                                                                                                                                             
    if rotated:
        texS.SetTextAngle(90)
    texS.SetTextFont(42)
    #texS.SetTextColor(ROOT.kBlack)                                                                                                                                                                                
    texS.SetTextSize(size)
    texS.Draw()
    return texS

ratioBand_count =0
def RatioErrorBand(Ratio,hUncUp,hUncDn,hTrueNoErrs,varName):
        global ratioBand_count 
        ratioBand_count+=1
        ratioGraph=ROOT.TGraphAsymmErrors(Ratio)
        ROOT.SetOwnership(ratioGraph,False)
        tmpData = Ratio.Clone("tmp")
        for i in range(1, tmpData.GetNbinsX()+1):
            if hTrueNoErrs.GetBinContent(i)==0:
                continue
            eUp=hUncUp.GetBinContent(i)
            eDn=hUncDn.GetBinContent(i)
            tru=hTrueNoErrs.GetBinContent(i)
            #print "eUp: ",eUp, "","eDn: ",eDn
            errorUp = tmpData.GetBinContent(i) + math.sqrt(math.pow(tmpData.GetBinError(i),2) + math.pow((eUp/tru),2))
            errorUp -= Ratio.GetBinContent(i) 
            errorDn = max(tmpData.GetBinContent(i) - math.sqrt(math.pow(tmpData.GetBinError(i),2) + math.pow((eDn/tru),2)),0)
            errorDn = Ratio.GetBinContent(i) - errorDn
            #print "Ratio (bin): ",i
            #print "stat. error: ",tmpData.GetBinError(i)
            #print "eUp/tru: ",eUp/tru
            #print "eDn/tru: ",eDn/tru
            #print "TotErrorUp: ",errorUp, "","TotErrorDn: ",errorDn
            
            #if "Mass" in varName and not "Full" in varName:
            #    if i==0:
            #        print("=======Test stat+sys vs stat========")
            #    print("Bin %s stat, sys up,dn"%i)
            #    print(tmpData.GetBinError(i),eUp/tru,eDn/tru)
            #    print("Ratio stat+sys/stat up/dn")
            #    print(errorUp/tmpData.GetBinError(i),errorDn/tmpData.GetBinError(i))

            #    errorUp = errorUp#*1.04
            #    errorDn = errorDn#*1.04


            ratioGraph.SetPointEYhigh(i-1, errorUp)
            ratioGraph.SetPointEYlow(i-1, errorDn)
        ratioGraph.SetFillColorAlpha(1,0.3)
        ratioGraph.SetFillStyle(3001)
        ratioGraph.GetXaxis().SetLabelSize(0)
        ratioGraph.GetXaxis().SetTitleSize(0)
        if ratioBand_count ==1:
            ratioGraph.GetYaxis().SetTickLength(0.0)
        else:
            ratioGraph.GetYaxis().SetTickLength(0.0)
        #ratioGraph.GetYaxis().SetLabelSize(0)
        #ratioGraph.GetYaxis().SetTitleSize(0)
        ratioGraph.GetXaxis().SetLimits(Ratio.GetXaxis().GetXmin(),Ratio.GetXaxis().GetXmax())
        with open('varsFile.json') as var_json_file:
            myvar_dict = json.load(var_json_file)
        
        if varName=="drz1z2":
            ratioGraph.SetMaximum(myvar_dict[my_varName]['ratio_max'])
            ratioGraph.SetMinimum(myvar_dict[my_varName]['ratio_min'])
        else:
            ratioGraph.SetMaximum(myvar_dict[my_varName]['ratio_max'])
            ratioGraph.SetMinimum(myvar_dict[my_varName]['ratio_min'])
        return ratioGraph

def MainErrorBand(hMain,hUncUp,hUncDn,varName,norm,normFb):
        with open('varsFile.json') as var_json_file:
            myvar_dict = json.load(var_json_file)
        MainGraph=ROOT.TGraphAsymmErrors(hMain)
        ROOT.SetOwnership(MainGraph,False)
        tmpData = hMain.Clone("tmp")
        for i in range(1, tmpData.GetNbinsX()+1):
            if hMain.GetBinContent(i)==0:
                continue
            eUp=hUncUp.GetBinContent(i)
            eDn=hUncDn.GetBinContent(i)
            #print "eUp: ",eUp, "","eDn: ",eDn
            errorUp = tmpData.GetBinContent(i) + math.sqrt(math.pow(tmpData.GetBinError(i),2) + math.pow(eUp,2))
            errorUp -= hMain.GetBinContent(i) 
            errorDn = max(tmpData.GetBinContent(i) - math.sqrt(math.pow(tmpData.GetBinError(i),2) + math.pow(eDn,2)),0)
            errorDn = hMain.GetBinContent(i) - errorDn
            #print "Main (bin): ",i
            #print "stat. error: ",tmpData.GetBinError(i)
            #print "eUp/tru: ",eUp/tru
            #print "eDn/tru: ",eDn/tru
            #print "TotErrorUp: ",errorUp, "","TotErrorDn: ",errorDn
            MainGraph.SetPointEYhigh(i-1, errorUp)
            MainGraph.SetPointEYlow(i-1, errorDn)
        MainGraph.SetFillColorAlpha(1,0.3)
#        MainGraph.SetFillColorAlpha(1,0.7)
        MainGraph.SetFillStyle(3001)
        if norm:
            drawyTitle = _yTitle[varName]
        elif normFb:
            drawyTitle = _yTitleNoNorm[varName]
        else:
            drawyTitle = "Events"
        MainGraph.GetYaxis().SetTitle(drawyTitle)
        #MainGraph.GetYaxis().CenterTitle()
        if "Allj" in hUncUp.GetName():
            tmpfac = 1.2
            tmpfac2 = 0.8
        elif "Full" in hUncUp.GetName():
            tmpfac = 1
            tmpfac2 = 0.9
        else:
            tmpfac = 1.2
            tmpfac2 = 0.8

        MainGraph.GetYaxis().SetTitleSize(tmpfac*ytilt_fac*hMain.GetYaxis().GetTitleSize())
        MainGraph.GetYaxis().SetLabelSize(1.1*hMain.GetYaxis().GetLabelSize()) #1.3
        if varName=="drz1z2":
            MainGraph.GetYaxis().SetTitleOffset(1.0)
        else:
            MainGraph.GetYaxis().SetLabelOffset(0.0)
            MainGraph.GetYaxis().SetTitleOffset(hMain.GetYaxis().GetTitleOffset()*args['titleOffset']*tmpfac2)
        #MainGraph.GetXaxis().SetLabelSize(0)
        #MainGraph.GetXaxis().SetTitleSize(0)
        #MainGraph.GetYaxis().SetLabelSize(0)
        #MainGraph.GetYaxis().SetTitleSize(0)
        if varName=="mass":
            MainGraph.GetYaxis().ChangeLabel(1,-1,0)
            MainGraph.GetYaxis().ChangeLabel(4,-1,0)
        MainGraph.GetXaxis().SetLimits(hMain.GetXaxis().GetXmin(),hMain.GetXaxis().GetXmax())
        #MainGraph.SetMaximum(1.5)

        #MainGraph.SetMaximum(1.2*(hMain.GetMaximum())*args["scaleymax"])
        MainGraph.SetMinimum(myvar_dict[varName]["ymin_fac"]*args['scaleymin']*(hMain.GetMinimum()))
        #if varName=="drz1z2":
        #    MainGraph.SetMinimum(0.0)
        #else:
            #MainGraph.SetMinimum(0.5*(hMain.GetMinimum()))
        return MainGraph

def generatePlots(hUnfolded,hUncUp,hUncDn,hTruth,hTruthAlt,varName,norm,normFb,lumi,unfoldDir,hTruthNNLO=None,hTruthEWC=None):
    global include_MiNNLO
    global EW_corr
    reset_include_MiNNLO = False
    reset_EW_corr = False
    if include_MiNNLO and not hTruthNNLO:
        include_MiNNLO = False
        reset_include_MiNNLO = True
    if EW_corr and not hTruthEWC:
        EW_corr = False
        reset_EW_corr = True
    with open('varsFile.json') as var_json_file:
        myvar_dict = json.load(var_json_file)
    
    top_xy = myvar_dict[varName]['top_xy'] #for MC labels positioning
    bottom_xy = myvar_dict[varName]['bottom_xy']
    xyP3 = myvar_dict[varName]['xyP3'] #for MC labels positioning
    xyP4 = myvar_dict[varName]['xyP4']
    
    top_fontsize=myvar_dict[varName]['top_size']
    bottom_fontsize=myvar_dict[varName]['bottom_size']
    fontsizeP3=myvar_dict[varName]['size_P3']  #Font sizes for 3rd and 4th ratio panel label
    fontsizeP4=myvar_dict[varName]['size_P4']
    ymax_fac=myvar_dict[varName]['ymax_fac']
    ymin_fac_extra=myvar_dict[varName]['ymin_fac_extra']

    UnfHists=[]
    TrueHists=[]
    # for normalization if needed
    nomArea = hUnfolded.Integral(1,hUnfolded.GetNbinsX())
    # Make uncertainties out of the unfolded histos
    ### plot
    hUnf = hUnfolded.Clone()
    hUnf.SetLineWidth(3)
    hTrue = hTruth.Clone()
    #Alt Signal 
    hTrueAlt = hTruthAlt.Clone()
    hTrueLeg = hTruthAlt.Clone() #doesn't seem to get used

    if include_MiNNLO:
        hTrueNNLO = hTruthNNLO.Clone()
        if EW_corr:
            hTrueEWC = hTruthEWC.Clone()    
    #lumi provided already in fb-1
    lumifb = lumi
    print("======================hUnf Integral before normalization: %s========================"%hUnf.Integral(1,hUnf.GetNbinsX()))
    if norm:
        hUnf.Scale(1.0/(hUnf.Integral(1,hUnf.GetNbinsX())))
    elif normFb:
        hUnf.Scale(1.0/lumifb)
    else:
        print "no special normalization"

    print ("hTrue histo here: ",hTrue)
    print ("unfoldDir: ",unfoldDir)
    xaxisSize = hUnf.GetXaxis().GetTitleSize()
    yaxisSize = hTrue.GetXaxis().GetTitleSize()
    if unfoldDir:
        #Create a ratio plot
        c,pad1 = createCanvasPads(varName)
        c.SetCanvasSize(1000, 1300);
        Unfmaximum = hUnf.GetMaximum()
        #hTrue.SetFillColor(ROOT.TColor.GetColor("#99ccff"))
        #hTrue.SetLineColor(ROOT.TColor.GetColor('#000099')) 
        hTrue.SetFillColor(ROOT.TColor.GetColor("#add8e6"))
        hTrue.SetLineColor(ROOT.TColor.GetColor('#377eb8'))
        hTrue.SetMarkerColor(ROOT.TColor.GetColor('#377eb8'))
        hTrue.SetLineStyle(1)
        hTrue.SetFillStyle(0)
        #AltSignal
        hTrueAlt.SetFillColor(2)
        hTrueAlt.SetLineStyle(10)#dashes
        hTrueAlt.SetFillStyle(0)#hollow
        hTrueAlt.SetLineColor(ROOT.kRed)
        hTrueAlt.SetMarkerColor(ROOT.kRed)
        if include_MiNNLO:
            hTrueNNLO.SetFillColor(8)
            hTrueNNLO.SetLineStyle(2)# special dashes
            hTrueNNLO.SetFillStyle(0)#hollow
            hTrueNNLO.SetLineColor(ROOT.kViolet)
            hTrueNNLO.SetMarkerColor(ROOT.kViolet)

            if EW_corr:
                hTrueEWC.SetFillColor(8)
                hTrueEWC.SetLineStyle(5)# special dashes
                hTrueEWC.SetFillStyle(0)#hollow
                hTrueEWC.SetLineColor(ROOT.kOrange)
                hTrueEWC.SetMarkerColor(ROOT.kOrange)
                
        print "Total Unf Data Integral",hUnf.Integral()
        Truthmaximum = hTrue.GetMaximum()
        Truthmaximum2 = hTrueAlt.GetMaximum()
        hTrue.SetLineWidth(4*hTrue.GetLineWidth())
        hTrueAlt.SetLineWidth(4*hTrueAlt.GetLineWidth())
        if include_MiNNLO:
            hTrueNNLO.SetLineWidth(4*hTrueNNLO.GetLineWidth())
            if EW_corr:
                hTrueEWC.SetLineWidth(4*hTrueEWC.GetLineWidth())

        if not norm and normFb:
            print "Inclusive fiducial cross section = {} fb".format(hUnf.Integral(1,hUnf.GetNbinsX()))
        if norm or normFb:
            normalizeBins(hUnf)

        if norm:
            hUncUp.Scale(1.0/hUnfolded.Integral(1,hUnfolded.GetNbinsX()))
            hUncDn.Scale(1.0/hUnfolded.Integral(1,hUnfolded.GetNbinsX()))
        elif normFb:
            hUncUp.Scale(1.0/lumifb)
            hUncDn.Scale(1.0/lumifb)
        else:
            print "no special normalization"

        if norm or normFb:
            normalizeBins(hUncUp)
            normalizeBins(hUncDn)

        #A good place to do extraction for HEPData, where hUnf and hUncUp and hUncDn are fully normalized
        quicksave = False
        if quicksave:
            if os.path.isfile("HEPData_extraction.root"): #redundant check
                extractionFile = ROOT.TFile("HEPData_extraction.root","UPDATE")
            else: 
                extractionFile = ROOT.TFile("HEPData_extraction.root","RECREATE")

            extractionFile.cd()
            hUnf.Write()
            hUncUp.Write()
            hUncDn.Write()
            print("Extracted hists for HEPData")
            sys.exit()

        if norm:
            trueInt = hTrue.Integral(1,hTrue.GetNbinsX())
            hTrue.Scale(1.0/trueInt)
            #hTrueUncUp /= trueInt # (trueInt + hTrueUncUp.Integral(0,hTrueUncUp.GetNbinsX()+1))
            #hTrueUncDn /= trueInt # (trueInt - hTrueUncDn.Integral(0,hTrueUncDn.GetNbinsX()+1))
            #Alt Signal
            AltTrueInt = hTrueAlt.Integral(1,hTrueAlt.GetNbinsX())
            hTrueAlt.Scale(1.0/AltTrueInt)

            if include_MiNNLO:
                NNLOTrueInt = hTrueNNLO.Integral(1,hTrueNNLO.GetNbinsX())
                hTrueNNLO.Scale(1.0/NNLOTrueInt)

                if EW_corr:
                    EWCTrueInt = hTrueEWC.Integral(1,hTrueEWC.GetNbinsX())
                    hTrueEWC.Scale(1.0/EWCTrueInt)

                    if varName == "MassAllj": #If m4l inclusive, after normalization, replace EWK with noEWK plot directly scaled by ratio [1.02244, etc.]
                        EWkfac = [1.02244,0.98414,0.97058,0.95705,0.95456,0.92758,0.91712,0.87614,0.81093]
                        
                        #for ifac in range(1,hTrueEWC.GetNbinsX()+1):
                        #    hTrueEWC.SetBinContent(ifac,hTrueNNLO.GetBinContent(ifac)*EWkfac[ifac-1])

        elif normFb:
            hTrue.Scale(1.0/lumifb)
            #hTrueUncUp /= lumifb
            #hTrueUncDn /= lumifb
            hTrueAlt.Scale(1.0/lumifb)
            if include_MiNNLO:
                hTrueNNLO.Scale(1.0/lumifb)
                if EW_corr:
                    hTrueEWC.Scale(1.0/lumifb)    
        else:
            print "no special normalization"

        print "Total Truth Integral",hTrue.Integral()
        print "Total Alt Truth Integral",hTrueAlt.Integral()
        if norm or normFb:
            normalizeBins(hTrue)
            #normalizeBins(hTrueUncUp)
            #normalizeBins(hTrueUncDn)
            normalizeBins(hTrueAlt)
            if include_MiNNLO:
                print("==================%s normalization check========================"%hTrueNNLO.GetName())
                print(hTrueNNLO.Integral(1,hTrueNNLO.GetNbinsX())) #print before normalization by bin width

                normalizeBins(hTrueNNLO)
                if EW_corr:
                    print("==================%s normalization check========================"%hTrueEWC.GetName())
                    print(hTrueEWC.Integral(1,hTrueEWC.GetNbinsX()))

                    normalizeBins(hTrueEWC)    

        #Don't know why draw twice. Commented the following two lines.
        #hTrue.Draw("HIST")
        #hTrueAlt.Draw("HIST")
        #pdb.set_trace()
        
        #if(Unfmaximum > Truthmaximum):
        #    hTrue.SetMaximum(Unfmaximum*args["scaleymax"]*ymax_fac)
        #else:
        #    hTrue.SetMaximum(Truthmaximum*args["scaleymax"]*ymax_fac)

        #Print bin content:
        xsecfac = 40.5
        print("==========Info: Unfolded bin content, total uncup,uncdn, and MiNNLO bin content and bin err=====================")
        print([xsecfac*hUnf.GetBinContent(infotmp) for infotmp in range(1,hUnf.GetNbinsX()+1)])
        print([xsecfac*math.sqrt((hUncUp.GetBinContent(infotmp))**2+(hUnf.GetBinError(infotmp))**2) for infotmp in range(1,hUnf.GetNbinsX()+1)])
        print([xsecfac*math.sqrt((hUncDn.GetBinContent(infotmp))**2+(hUnf.GetBinError(infotmp))**2) for infotmp in range(1,hUnf.GetNbinsX()+1)])
        if include_MiNNLO:
            print([xsecfac*hTrueNNLO.GetBinContent(infotmp) for infotmp in range(1,hUnf.GetNbinsX()+1)])
            print([xsecfac*hTrueNNLO.GetBinError(infotmp) for infotmp in range(1,hUnf.GetNbinsX()+1)])
        print("===================================================")

        hTrue.GetXaxis().SetTitle("")

        UnfErrBand = MainErrorBand(hUnf,hUncUp,hUncDn,varName,norm,normFb)
        if varName=="mass":
            UnfErrBand.SetMaximum(0.01*args['scaleymax']*ymax_fac)
        UnfErrBand.Draw("a2")
        hTrue.GetXaxis().SetLabelSize(0)
        hTrue.GetXaxis().SetTitleSize(0)
        #hTrue.GetYaxis().SetTitle("Events")
        #hTrue.GetYaxis().SetTitleOffset(1.0)
        hTrueAlt.GetXaxis().SetLabelSize(0)
        hTrueAlt.GetXaxis().SetTitleSize(0)
        hTrueAlt.Draw("E1 SAME") #drawing second time, for updating?
        hTrue.Draw("E1 SAME") #("PE1SAME")
        if include_MiNNLO:
            hTrueNNLO.GetXaxis().SetLabelSize(0)
            hTrueNNLO.GetXaxis().SetTitleSize(0)
            hTrueNNLO.Draw("E1 SAME") 
           
            if EW_corr:
                hTrueEWC.GetXaxis().SetLabelSize(0)
                hTrueEWC.GetXaxis().SetTitleSize(0)
                hTrueEWC.Draw("E1 SAME") 
                
#        hTrueAlt.Draw("HISTSAME") #drawing second time, for updating?
#        hTrue.Draw("HISTSAME")
        #hUnf.Sumw2(False)
        #hUnf.SetBinErrorOption(ROOT.TH1.kPoisson)
        hUnf.SetLineColor(ROOT.kBlack)
        hUnf.SetMarkerStyle(20)
        hUnf.SetMarkerSize(0.7)
        hUnf.GetXaxis().SetTitle("")
        hUnf.GetXaxis().SetLabelSize(0)
        hUnf.GetXaxis().SetTitleSize(0)
        hUnf.Draw(crossDrawOpt)

        axismaximum = max([hUnf.GetMaximum(),hTrue.GetMaximum(),hTrueAlt.GetMaximum(),UnfErrBand.GetMaximum()])
        if include_MiNNLO:
            axismaximum = max(axismaximum,hTrueNNLO.GetMaximum())
            if EW_corr:
                axismaximum = max(axismaximum,hTrueEWC.GetMaximum())    
            
        
        hTrue.SetMaximum(axismaximum*args["scaleymax"]*ymax_fac)
        hTrueAlt.SetMaximum(axismaximum*args["scaleymax"]*ymax_fac)
        if include_MiNNLO:
            hTrueNNLO.SetMaximum(axismaximum*args["scaleymax"]*ymax_fac)
            if EW_corr:
                hTrueEWC.SetMaximum(axismaximum*args["scaleymax"]*ymax_fac)
        hUnf.SetMaximum(axismaximum*args["scaleymax"]*ymax_fac)
        UnfErrBand.SetMaximum(axismaximum*args["scaleymax"]*ymax_fac)
        
        axisminimum = min([hUnf.GetMinimum(),hTrue.GetMinimum(),hTrueAlt.GetMinimum(),UnfErrBand.GetMinimum()])
        if include_MiNNLO:
            axisminimum = min(axisminimum,hTrueNNLO.GetMinimum())
            if EW_corr:
                axisminimum = min(axisminimum,hTrueEWC.GetMinimum())
        if not ymin_fac_extra==1.:
            hTrue.SetMinimum(axisminimum*ymin_fac_extra) #args["scaleymin"] is set to 0.3, not used here
            hTrueAlt.SetMinimum(axisminimum*ymin_fac_extra)
            if include_MiNNLO:
                hTrueNNLO.SetMinimum(axisminimum*ymin_fac_extra)
                if EW_corr:
                    hTrueEWC.SetMinimum(axisminimum*ymin_fac_extra)
            hUnf.SetMinimum(axisminimum*ymin_fac_extra)
            UnfErrBand.SetMinimum(axisminimum*ymin_fac_extra)
      
        #ROOT.dotrootImport('uhussain/CMSPlotDecorations')
        #scale_label = "Normalized to Unity" if args['lumi'] < 0 else \
        #    "%0.1f fb^{-1}" % args['lumi']
        #
        #lumi_text = ""
        #if args['thesis']:
        #    lumi_text = "Thesis" 
        #elif args['preliminary']:
        #    lumi_text = "Preliminary" 
        #
        #ROOT.CMSlumi(c, 0, 0, "%s (13 TeV)" % scale_label,lumi_text)
                #"Preliminary Simulation" if args.simulation else "Preliminary")
        
        offset = ROOT.gPad.GetLeftMargin() - 0.07 if args['legend_left'] else \
            ROOT.gPad.GetRightMargin() - 0.07 
        width = .33
        width *= args['scalelegx']
        xdist = 0.1 if args['legend_left'] else 0.91
        if varName=="leppt":
            xdist=0.15
        if varName=="drz1z2":
            xdist=0.07
        xcoords = [xdist+offset, xdist+width+offset] if args['legend_left'] \
            else [xdist-width-offset, xdist-offset]
        unique_entries = min(2, 8)
        ymax = 0.45 if varName=="leppt" else 0.915
        ycoords = [ymax, ymax - 0.08*unique_entries*args['scalelegy']]
        coordy_subtract = 0.05
        if varName == "MassAllj":
            coordy_subtract = 0.08
        if "Full" in varName:
            coordy_subtract = 0.0
        coords = [xcoords[0]-0.1, ycoords[0], xcoords[1], ycoords[1]-coordy_subtract] #extended legend frame
        if not include_MiNNLO:
            legend = getPrettyLegend(hTrue, hUnf, hTrueAlt, UnfErrBand, coords)
        else:
            if not EW_corr:
                legend = getPrettyLegend(hTrue, hUnf, hTrueAlt, UnfErrBand, coords,hTrueNNLO)
            else:
                legend = getPrettyLegend(hTrue, hUnf, hTrueAlt, UnfErrBand, coords,hTrueNNLO,hTrueEWC)
                
        legend.Draw()
        texS,texS1,texS2=getLumiTextBox()
        sigLabel = "POWHEG+MCFM+Pythia8" #used?
        sigLabelAlt = "MG5_aMC@NLO+MCFM+Pythia8"

        if varName in ["jetPt[0]","jetPt[1]","absjetEta[0]","absjetEta[1]","mjj","dEtajj"]:
            if varName in ["jetPt[0]","absjetEta[0]"]:
                nJetsText=getAxisTextBox(0.17,0.1,"Events with #geq 1 jet",0.06,False)
            if varName in ["jetPt[1]","absjetEta[1]","mjj","dEtajj"]:
                nJetsText=getAxisTextBox(0.17,0.1,"Events with #geq 2 jets",0.06,False)
        
        if "Mass" in varName:
            if "0" in varName:
                nJetsText=getAxisTextBox(0.17,0.1,"Events with 0 jet",0.06,False)
            if "1" in varName:
                nJetsText=getAxisTextBox(0.17,0.1,"Events with 1 jet",0.06,False)
            if "2" in varName:
                nJetsText=getAxisTextBox(0.17,0.1,"Events with 2 jets",0.06,False)
            if "34" in varName:
                nJetsText=getAxisTextBox(0.17,0.1,"Events with #geq 3 jets",0.06,False)

        #if varName=="dphiz1z2" or varName=="drz1z2":
        #    leg = ROOT.TLegend(0.15,0.60,0.15+0.015*len(sigLabelAlt),0.90,"")
        #elif varName=="leppt":
        #    leg = ROOT.TLegend(0.20,0.18,0.20+0.015*len(sigLabelAlt),0.48,"")
        #else:
        #    leg = ROOT.TLegend(0.55,0.60,0.55+0.015*len(sigLabelAlt),0.90,"")
        #leg.AddEntry(hUnf,"Data + stat. unc.","lep")
        #leg.AddEntry(UnfErrBand, "Stat. #oplus syst. unc.","f")
        #leg.AddEntry(hTrue, sigLabel,"lf")

        #hTrueLeg.SetFillColor(2)
        #hTrueLeg.SetLineStyle(10)#dashes
        #hTrueLeg.SetFillColorAlpha(2,0.4)
        #hTrueLeg.SetFillStyle(3001)#solid
        #hTrueLeg.SetLineColor(ROOT.kRed)
        #hTrueLeg.SetLineWidth(4*hTrueLeg.GetLineWidth())
        #leg.AddEntry(hTrueLeg, sigLabelAlt,"l")
        #leg.SetFillColor(ROOT.kWhite)
        #leg.SetBorderSize(1)
        #leg.SetFillStyle(1001)
        #leg.SetTextSize(0.025)
        #leg.Draw()

        #SecondPad
        #nominal sample
        with open('listFile.json') as list_json_file:
            mylist_dict = json.load(list_json_file)
        ratioName_nom =mylist_dict["sigLabel"]
        ratioName_alt =mylist_dict["sigLabelAlt"]

        pad2 = createPad2(c)

        hTrueNoErrs = hTrue.Clone() # need central value only to keep ratio uncertainties consistent
        nbins=hTrueNoErrs.GetNbinsX()
        print("trueNbins: ",nbins)

        Unfbins=hUnf.GetNbinsX()
        print("UnfNbins: ",Unfbins)
        #hTrueNoErrs.SetError(array.array('d',[0.]*nbins))
        #Starting the ratio proceedure
        Ratio,line = createRatio(hUnf, hTrueNoErrs)
        ratioErrorBand = RatioErrorBand(Ratio,hUncUp,hUncDn,hTrueNoErrs,varName)
        ratioErrorBand.GetYaxis().SetLabelSize(0)
        ratioErrorBand.GetYaxis().SetTitleSize(0)
        Ratio.GetYaxis().SetLabelSize(0)
        Ratio.GetYaxis().SetTitleSize(0)
        ratioErrorBand.Draw("a2")
        
        sigTex = getSigTextBox(0.15,0.8,sigLabel,0.14) #used?
        Ratio.Draw(crossDrawOpt)
        line.SetLineColor(ROOT.kBlack)
        #line.SetLineColor(ROOT.TColor.GetColor('#377eb8'))
        line.Draw("same")

        Altyaxis = ROOT.TGaxis(hUnf.GetXaxis().GetXmin(),ratioErrorBand.GetMinimum(),hUnf.GetXaxis().GetXmin(),ratioErrorBand.GetMaximum(),ratioErrorBand.GetMinimum(),ratioErrorBand.GetMaximum(),3,"CS")
        Altyaxis.SetNdivisions(yrdiv)
        dataTheoSize = 0.223
        tmpy = -0.01
        if varName == "MassAllj":
            dataTheoSize *= 1.
        elif "Full" in varName:
            dataTheoSize *= 0.9
            tmpy = 0.0
        else:
            dataTheoSize *= 0.97
            tmpy = 0.0
        axText2=getAxisTextBox(0.06,tmpy,"Data/Pred.",dataTheoSize,True)
        MCTextNom=getAxisTextBox(top_xy[0],top_xy[1],ratioName_nom,top_fontsize,False)
        

        #Altyaxis.SetTitle("#scale[1.2]{Data/%s}"%ratioName_nom)
        Altyaxis.SetTickLength(yrtl)
        Altyaxis.SetLabelFont(42)
        Altyaxis.SetLabelOffset(0.025) #0.01
        if not include_MiNNLO:
            Altyaxis.SetLabelSize(0.189)
        else:
            Altyaxis.SetLabelSize(0.167)
        Altyaxis.SetTitleFont(42)
        Altyaxis.SetTitleSize(0.16) #0.16
        Altyaxis.SetTitleOffset(0.29) #0.29
        #Altyaxis.ChangeLabel(2,-1,0.189,-1,-1,-1,"2")
        Altyaxis.Draw("SAME")
        AltyaxisR = getRYaxis(hUnf,ratioErrorBand,False)
        
        #ThirdPad
        pad3 = createPad3(c)
        


        hTrueAltNoErrs = hTrueAlt.Clone() # need central value only to keep ratio uncertainties consistent
        #nbins=hTrueNoErrs.GetNbinsX()
        #print("trueNbins: ",nbins)

        #Unfbins=hUnf.GetNbinsX()
        #print("UnfNbins: ",Unfbins)

        #hTrueNoErrs.SetError(array.array('d',[0.]*nbins))
        #Starting the ratio proceedure
        AltRatio,Altline = createRatio(hUnf, hTrueAltNoErrs)
        AltRatioErrorBand = RatioErrorBand(AltRatio,hUncUp,hUncDn,hTrueAltNoErrs,varName) 
        AltRatioErrorBand.GetYaxis().SetLabelSize(0)
        AltRatioErrorBand.GetYaxis().SetTitleSize(0)

        #Currently this is the line that will change bottom tick length
        #AltRatioErrorBand.GetXaxis().SetTickLength(0.1)
        
        AltRatioErrorBand.Draw("a2")
        #if varName == "nJets":
        #    AltRatio.GetXaxis().SetNdivisions(505)
        #    AltRatio.GetXaxis().CenterLabels(True)
        
        AltRatio.Draw(crossDrawOpt)
        #ratioErrorBand.Draw("p")
        Altline.SetLineColor(ROOT.kBlack)
        Altline.Draw("same")
        
        if include_MiNNLO:
            MCTextAlt=getAxisTextBox(bottom_xy[0],bottom_xy[1],ratioName_alt,bottom_fontsize,False)
        else:
            MCTextAlt=getAxisTextBox(bottom_xy[0],bottom_xy[1],ratioName_alt,bottom_fontsize*(1-bmg),False)

        AltTex = getSigTextBox(0.15,0.85,sigLabelAlt,0.11)

        yaxis = ROOT.TGaxis(hUnf.GetXaxis().GetXmin(),ratioErrorBand.GetMinimum(),hUnf.GetXaxis().GetXmin(),ratioErrorBand.GetMaximum(),ratioErrorBand.GetMinimum(),ratioErrorBand.GetMaximum(),3,"CS")
        yaxis.SetNdivisions(yrdiv)
        yaxis.SetTickLength(yrtl)
        #axText3=getAxisTextBox(0.06,0.0,"Data/%s"%ratioName_alt,0.23,True)
        #yaxis.SetTitle("#scale[1.2]{Data/%s}"%ratioName_alt)
        #yaxis.SetTitle("Data/Theo.")
        yaxis.SetLabelFont(42)
        yaxis.SetLabelOffset(0.025) #0.01
        if not include_MiNNLO:
            yaxis.SetLabelSize(0.108)
        else:
            yaxis.SetLabelSize(0.163)
        yaxis.SetTitleFont(42)
        yaxis.SetTitleSize(0.12)
        yaxis.SetTitleOffset(0.365)
        yaxis.Draw("SAME")
        
        yaxisR = getRYaxis(hUnf,ratioErrorBand,not include_MiNNLO)

        if include_MiNNLO:
            #Fourth pad
            pad4 = createPad4(c)
            
            hTrueNNLONoErrs = hTrueNNLO.Clone() # need central value only to keep ratio uncertainties consistent
            if EW_corr:
                hTrueEWCNoErrs = hTrueEWC.Clone() # need central value only to keep ratio uncertainties consistent

            #nbins=hTrueNoErrs.GetNbinsX()
            #print("trueNbins: ",nbins)

            #Unfbins=hUnf.GetNbinsX()
            #print("UnfNbins: ",Unfbins)

            #hTrueNoErrs.SetError(array.array('d',[0.]*nbins))
            #Starting the ratio proceedure
            NNLORatio,NNLOline = createRatio(hUnf, hTrueNNLONoErrs)
            if EW_corr:
                EWCRatio,EWCline = createRatio(hUnf,hTrueEWCNoErrs)
                #EWCRatio,EWCline = createRatio(hTrueEWCNoErrs, hTrueNNLONoErrs)

            NNLORatioErrorBand = RatioErrorBand(NNLORatio,hUncUp,hUncDn,hTrueNNLONoErrs,varName) 
            NNLORatioErrorBand.GetYaxis().SetLabelSize(0)
            NNLORatioErrorBand.GetYaxis().SetTitleSize(0)
            NNLORatioErrorBand.Draw("a2")
            NNLORatio.Draw(crossDrawOpt)

            if EW_P4 and EW_corr:
                EWCRatioErrorBand = RatioErrorBand(EWCRatio,hUncUp,hUncDn,hTrueEWCNoErrs,varName) 
                EWCRatioErrorBand.GetYaxis().SetLabelSize(0)
                EWCRatioErrorBand.GetYaxis().SetTitleSize(0)

            if EW_corr and not EW_P4:
                EWCRatio.SetLineColor(ROOT.kOrange)
                EWCRatio.Draw(crossDrawOpt)

            #ratioErrorBand.Draw("p")
            NNLOline.SetLineColor(ROOT.kBlack)
            NNLOline.Draw("same")
            
            ratioName_NNLO = "nNNLO+PS"
            sigLabelNNLO = "nNNLO+PS"
            MCTextNNLO=getAxisTextBox(xyP3[0],xyP3[1],ratioName_NNLO,fontsizeP3,False)
            NNLOTex = getSigTextBox(0.15,0.85,sigLabelNNLO,0.11)

            NNLOyaxis = ROOT.TGaxis(hUnf.GetXaxis().GetXmin(),ratioErrorBand.GetMinimum(),hUnf.GetXaxis().GetXmin(),ratioErrorBand.GetMaximum(),ratioErrorBand.GetMinimum(),ratioErrorBand.GetMaximum(),3,"CS")
            NNLOyaxis.SetNdivisions(yrdiv)
            NNLOyaxis.SetTickLength(yrtl)
            if not (EW_corr and EW_P4):
                NNLOyaxis.SetTickLength(yrtl/(1.-bmg))    
            #axText3=getAxisTextBox(0.06,0.0,"Data/%s"%ratioName_alt,0.23,True)
            #NNLOyaxis.SetTitle("#scale[1.2]{Data/%s}"%ratioName_alt)
            #NNLOyaxis.SetTitle("Data/Theo.")
            NNLOyaxis.SetLabelFont(42)
            NNLOyaxis.SetLabelOffset(0.025) #0.01
            if not (EW_P4 and EW_corr):
                NNLOyaxis.SetLabelSize(0.1) #0.1485
            else:
                NNLOyaxis.SetLabelSize(0.165) #0.1485

            NNLOyaxis.SetTitleFont(42)
            NNLOyaxis.SetTitleSize(0.12)
            NNLOyaxis.SetTitleOffset(0.365)
            NNLOyaxis.Draw("SAME")

            NNLOyaxisR = getRYaxis(hUnf,ratioErrorBand, not (EW_corr and EW_P4))

            if EW_P4 and EW_corr:
                pad5 = createPad5(c)

                EWCRatioErrorBand.Draw("a2")
                EWCRatio.Draw(crossDrawOpt)
                EWCline.SetLineColor(ROOT.kBlack)
                EWCline.Draw("same")

                ratioName_EWC = "(nNNLO+PS)#times K_{EW}"
                sigLabelEWC = "(nNNLO+PS)#times K_{EW}"
                MCTextEWC=getAxisTextBox(xyP4[0],xyP4[1],ratioName_EWC,fontsizeP4,False)
                EWCTex = getSigTextBox(0.15,0.85,sigLabelEWC,0.11)

                EWCyaxis = ROOT.TGaxis(hUnf.GetXaxis().GetXmin(),ratioErrorBand.GetMinimum(),hUnf.GetXaxis().GetXmin(),ratioErrorBand.GetMaximum(),ratioErrorBand.GetMinimum(),ratioErrorBand.GetMaximum(),3,"CS")
                EWCyaxis.SetNdivisions(yrdiv)
                EWCyaxis.SetTickLength(yrtl/(1-bmg))
                #axText3=getAxisTextBox(0.06,0.0,"Data/%s"%ratioName_alt,0.23,True)
                #EWCyaxis.SetTitle("#scale[1.2]{Data/%s}"%ratioName_alt)
                #EWCyaxis.SetTitle("Data/Theo.")
                EWCyaxis.SetLabelFont(42)
                EWCyaxis.SetLabelOffset(0.025) #0.01
                EWCyaxis.SetLabelSize(0.1) #0.1485
                EWCyaxis.SetTitleFont(42)
                EWCyaxis.SetTitleSize(0.12)
                EWCyaxis.SetTitleOffset(0.365)
                EWCyaxis.Draw("SAME")
                EWCyaxisR = getRYaxis(hUnf,ratioErrorBand, True)


        #redraw axis
        if "Full" in varName and "Mass" in varName:
            xaxis = ROOT.TGaxis(hUnf.GetXaxis().GetXmin(),ratioErrorBand.GetMinimum(),hUnf.GetXaxis().GetXmax(),ratioErrorBand.GetMinimum(),hUnf.GetXaxis().GetXmin(),hUnf.GetXaxis().GetXmax(),510,"G")
        
            xaxis.SetMoreLogLabels(True)
            xaxis.SetTickLength(0.07)
            #xaxis.SetLabelSize(0.025)
            xaxis.ChangeLabel(1,-1,0.,-1,-1,-1,"")
        elif varName == "nJets":
            xaxis = ROOT.TGaxis(hUnf.GetXaxis().GetXmin(),ratioErrorBand.GetMinimum(),hUnf.GetXaxis().GetXmax(),ratioErrorBand.GetMinimum(),hUnf.GetXaxis().GetXmin(),hUnf.GetXaxis().GetXmax(),505)
            xaxis.CenterLabels(True)
            xaxis.ChangeLabel(4,-1,-1,-1,-1,-1,"#geq 3")
        else:
            xaxis = ROOT.TGaxis(hUnf.GetXaxis().GetXmin(),ratioErrorBand.GetMinimum(),hUnf.GetXaxis().GetXmax(),ratioErrorBand.GetMinimum(),hUnf.GetXaxis().GetXmin(),hUnf.GetXaxis().GetXmax(),510)
        xaxis.SetTitle(prettyVars[varName]+' '+units[varName])
        #labelTex = getSigTextBox(0.9,0.8,prettyVars[varName]+''+units[varName])
        #if varName=="mass":
            #xaxis.ChangeLabel(1,-1,0.1)
            #xaxis.ChangeLabel(2,0.2)
            #xaxis.ChangeLabel(3,0.3)
            #xaxis.ChangeLabel(4,0.4)
            #xaxis.ChangeLabel(5,0.5)
            #xaxis.ChangeLabel(6,0.6)
            #xaxis.ChangeLabel(7,0.7)
            #xaxis.ChangeLabel(8,0.8)
            #xaxis.ChangeLabel(9,0.9)
            #xaxis.ChangeLabel(10,1.0)
        xaxis.SetLabelFont(42)
        #xaxis.SetLabelOffset(0.03)
        xaxis.SetLabelOffset(0.03)
        #xaxis.SetTickLength(0.1)
        if "Full" in varName and "Mass" in varName:
            xaxis.SetLabelSize(0.12)
        elif "Mass" in varName and not "All" in varName:
            xaxis.SetLabelSize(0.15)
        elif "mjj" in varName:
            xaxis.SetLabelSize(0.15)
        else:
            xaxis.SetLabelSize(0.162)
        xaxis.SetTitleFont(42)
        xaxis.SetTitleSize(0.18)
        if EW_corr and EW_P4:
            xaxis.SetTitleOffset(1.)
        else:
            xaxis.SetTitleOffset(1.2)

        if varName=="mass":
            xaxis.SetNoExponent(True)
        xaxis.Draw("SAME")

        
        c.Update()
        print("CanvasWidth: ", c.GetWw())
        print("CanvasHeight: ", c.GetWh())
        plotName="Ratio_%s" % (varName)
        output_name="/".join([unfoldDir,plotName])
        #c.Print("%s/Ratio_%s.png" % (unfoldDir,varName))
        #c.Print("%s/Ratio_%s.eps" % (unfoldDir,varName))
        c.Print(output_name+".eps")
        c.Print(output_name+".png")
        subprocess.call(["epstopdf", "--outfile=%s" % output_name+".pdf", output_name+".eps"],env={})
        os.remove(output_name+".eps")
        del c

    if reset_include_MiNNLO:
        include_MiNNLO = True
        
    if reset_EW_corr:
        EW_corr = True
def mkdir(plotDir):
    for outdir in [plotDir]:
        try:
            os.makedirs(os.path.expanduser(outdir))
        except OSError as e:
            print e
            pass

#varNames={'mass': 'Mass','pt':'ZZPt','zpt':'ZPt','leppt':'LepPt','dphiz1z2':'dPhiZ1Z2','drz1z2':'dRZ1Z2'}
UnfoldDir=args['unfoldDir']
UnfoldOutDirs={}
#Make differential cross sections normalizing to unity area.')
norm = not args['noNorm']
#normalize with the luminosity instead of area under the curve
normFb = args['NormFb']
#def main():
runVariables=[]
runVariables.append(args['variable'])
print "runVariables: ",runVariables
#Plot histograms from these respective root files generated wiht saveUnfolded.py
with open('listFile.json') as list_json_file:
    mylist_dict = json.load(list_json_file)

if analysis=="ZZ4l2016":
    fUse = ROOT.TFile(mylist_dict['f16'],"read")
    #fUse = ROOT.TFile("UnfHistsFinal-18Jun2021-ZZ4l2016.root","read")
    #fUse = ROOT.TFile("UnfHistsFinal-18Apr2020-ZZ4l2016.root","read")
elif analysis=="ZZ4l2017":
    fUse = ROOT.TFile(mylist_dict['f17'],"read")
    #fUse = ROOT.TFile("UnfHistsFinal-22May2021-ZZ4l2017.root","read")
    #fUse = ROOT.TFile("UnfHistsFinal-18Apr2020-ZZ4l2017.root","read")
elif analysis=="ZZ4l2018":
    if args['lumi'] < 100. :
        fUse = ROOT.TFile(mylist_dict['f18'],"read")
        #fUse = ROOT.TFile("UnfHistsFinal-2018NewqqZZMC-14Jul2021-ZZ4l2018.root","read")
        #fUse = ROOT.TFile("UnfHistsFinal-22May2021-ZZ4l2018.root","read")
        #fUse = ROOT.TFile("UnfHistsFinal-18Apr2020-ZZ4l2018.root","read")
    else:
        fUse = ROOT.TFile(mylist_dict['fFull'],"read")
        #fUse = ROOT.TFile("allyear_UnfHist.root","read") 
        #fUse = ROOT.TFile("UnfHistsFinal-18Apr2020-ZZ4lFullRun2.root","read")
#fUse = ROOT.TFile.Open("UnfHistsFull09Nov2019-ZZ4l2018.root","update")

#keep a "global" copy of variable name
my_varName = ''
for varName in runVariables:

    my_varName = varName 
    print "varName:", varNames[varName]
    # save unfolded distributions by channel, then systematic
    hUnfolded = {}
    hTrue = {}
    hTrueAlt = {}
    hErr = {}
    hErrTrue = {}
    for chan in channels:
        print "channel: ",chan
        #print "hUnfolded: ",hUnfolded
        #print "hTrue: ",hTrue
        UnfoldOutDir=UnfoldDir+"/"+chan+"/plots"
        if chan not in UnfoldOutDirs:
            UnfoldOutDirs[chan]=UnfoldOutDir
        if not os.path.exists(UnfoldOutDir):
            mkdir(UnfoldOutDir)
        hUnfolded[chan] = fUse.Get(chan+"_"+varName+"_unf")
        hTrue[chan] = fUse.Get(chan+"_"+varName+"_true")
        hTrueAlt[chan] = fUse.Get(chan+"_"+varName+"_trueAlt")
        #print("returning unfolded? ",hUnfolded[chan])
        #print("returning truth? ",hTrue[chan])
        #print("returning Alt truth? ",hTrueAlt[chan])
        #Get the total UncUp and total UncDown histograms from the file as well
        hUncUp = fUse.Get(chan+"_"+varName+"_totUncUp") 
        hUncDn = fUse.Get(chan+"_"+varName+"_totUncDown")
        #print "UnfoldOutDir:",UnfoldOutDir
        print("Somehow channel plots get called")
        sys.exit()
        generatePlots(hUnfolded[chan],hUncUp,hUncDn,hTrue[chan],hTrueAlt[chan],varName,norm,normFb,args['lumi'],UnfoldOutDir)
    
    if args['makeTotals']:
        #Now access the histograms for all channels combined together
        #While saving in makeResponseMatrix.py, make a "total" category as wel
        hUnfTot = fUse.Get("tot_"+varName+"_unf") 
        hTrueTot = fUse.Get("tot_"+varName+"_true")
        hTrueAltTot = fUse.Get("tot_"+varName+"_trueAlt") 
        if include_MiNNLO:
            hTrueNNLOTot = fUse.Get("tot_"+varName+"_trueNNLO")    
            if EW_corr:
                hTrueEWCTot = fUse.Get("tot_"+varName+"_trueEWC")    
                #hTrueEWCTot = fUse.Get("tot_"+varName+"_trueNNLONoGenW") 
        hTotUncUp = fUse.Get("tot_"+varName+"_totUncUp")
        hTotUncDn = fUse.Get("tot_"+varName+"_totUncDown") 
        UnfoldOutDir=UnfoldDir+"/"+"tot"+"/plots"
        if "tot" not in UnfoldOutDirs:
            UnfoldOutDirs["tot"]=UnfoldOutDir
        if not os.path.exists(UnfoldOutDir):
            mkdir(UnfoldOutDir)
        if not include_MiNNLO:
            generatePlots(hUnfTot,hTotUncUp,hTotUncDn,hTrueTot,hTrueAltTot,varName,norm,normFb,args['lumi'],UnfoldOutDir)
        else:
            if not EW_corr:
                generatePlots(hUnfTot,hTotUncUp,hTotUncDn,hTrueTot,hTrueAltTot,varName,norm,normFb,args['lumi'],UnfoldOutDir,hTrueNNLOTot)
            else:
                generatePlots(hUnfTot,hTotUncUp,hTotUncDn,hTrueTot,hTrueAltTot,varName,norm,normFb,args['lumi'],UnfoldOutDir,hTrueNNLOTot,hTrueEWCTot)

#Show plots nicely on my webpages
for cat in ["tot"]:   
#for cat in ["eeee","eemm","mmmm","tot"]:   
    #This is where we put all the plots in html format for quick access/debugging
    makeSimpleHtml.writeHTML(os.path.expanduser(UnfoldOutDirs[cat].replace("/plots", "")), "Unfolded Distributions (from MC)")
fUse.Close()
#exit(0)

#if __name__ == "__main__":
#    main()
