#!/usr/bin/env python
import ROOT
import glob
import math
from python import SelectorTools
from python import UserInput
from python import OutputTools
from python import ConfigureJobs
from python import HistTools
import makeSimpleHtml
import os,pdb,copy
import subprocess
import sys,json
import datetime
import array
import time
from ROOT import vector as Vec
VFloat = Vec('float')
from PlotTools import PlotStyle as Style, pdfViaTex
style = Style()
ROOT.gStyle.SetLineScalePS(1.8)

with open('listFile.json') as list_json_file:
    mylist_dict = json.load(list_json_file)

channels = ["eeee","eemm","mmmm"]
mynominalName=mylist_dict['nomname']
myaltname= mylist_dict['altname']
applyreg=mylist_dict["reg"]

#channels = ["eeee"]
def getComLineArgs():
    parser = UserInput.getDefaultParser()
    parser.add_argument("--lumi", "-l", type=float,
        default=41.5, help="luminosity value (in fb-1)")
    parser.add_argument("--output_file", "-o", type=str,
        default="", help="Output file name")
    parser.add_argument("--test", action='store_true',
        help="Run test job (no background estimate)")
    parser.add_argument("--uwvv", action='store_true',
        help="Use UWVV format ntuples in stead of NanoAOD")
    parser.add_argument("--with_background", action='store_true',
        help="Don't run background selector")
    parser.add_argument("--noHistConfig", action='store_true',
        help="Don't rely on config file to specify hist info")
    parser.add_argument("-j", "--numCores", type=int, default=1,
        help="Number of cores to use (parallelize by dataset)")
    parser.add_argument("--input_tier", type=str,
        default="", help="Selection stage of input files")
    parser.add_argument("--year", type=str,
        default="default", help="Year of Analysis")
    parser.add_argument("-c", "--channels", 
                        type=lambda x : [i.strip() for i in x.split(',')],
                        default=["eee","eem","emm","mmm"], help="List of channels"
                        "separated by commas. NOTE: set to Inclusive for NanoAOD")
    parser.add_argument("--scalefactors_file", "-sf", type=str,
        default="", help="ScaleFactors file name")
    parser.add_argument("--leptonSelections", "-ls", type=str,
        default="TightLeptons", help="Either All Loose or Tight")
    parser.add_argument("--output_selection", type=str,
        default="", help="Selection stage of output file "
        "(Same as input if not give)")
    parser.add_argument("-b", "--hist_names", 
                        type=lambda x : [i.strip() for i in x.split(',')],
                        default=["all"], help="List of histograms, "
                        "as defined in ZZ4lRun2DatasetManager, separated "
                        "by commas")
    parser.add_argument("--variable", "-vr", type=str,
        default="all", help="variableName")
    parser.add_argument('--noNorm', action='store_true',
                        help='Leave differential cross sections in abolute normalization rather than normalizing to unity area.')
    parser.add_argument('--plotResponse', action='store_true',
                        help='plot Response Matrices and covariance matrices.')
    parser.add_argument('--plotSystResponse', action='store_true',
                        help='plotResponse Matrices varied up or down with systematics.')
    parser.add_argument('--makeTotals', action='store_true',
                        help='plot total unfolded with uncertainities.')
    parser.add_argument('--noSyst', action='store_true',
                        help='No Systematics calculations.')
    parser.add_argument('--diagnostic', action='store_true',
                        help='Save diagnostic histograms and exit in _generateUncertainties function')
    parser.add_argument('--plotDir', type=str, nargs='?',
                        default='/afs/cern.ch/user/h/hehe/www/ZZFullRun2/PlottingResults/ZZ4l2017/ZZSelectionsTightLeps/ANPlots/ZZ4l2017/RespMat_Moriond2019IDMuSF',
                        help='Directory to put response and covariance plots in')
    parser.add_argument('--unfoldDir', type=str, nargs='?',
                        default='/afs/cern.ch/user/h/hehe/ZZFullRun2Unfolding',
                        help='Directory to put unfolded histo files in')
    parser.add_argument('--nIter', type=int, nargs='?', default=8,
                        help='Number of iterations for D\'Agostini method')
    return vars(parser.parse_args())

args = getComLineArgs()
manager_path = ConfigureJobs.getManagerPath()
selection = args['selection']
if selection == "":
    selection = "LooseLeptons"
    #commented_print "Info: Using BasicZZSelections for hist defintions"
#analysis = "/".join([args['analysis'], selection])
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
# Names of compiled C++ classes to make response matrices fast
# (this is extremely slow in Python because it requires a combination of
# information from multiple trees, which can't be done with TTree::Draw())
responseClassNames = {
    'mass' : {c:'FloatBranchResponseMatrixMaker' for c in channels},
    #'massFull' : {c:'FullSpectrumFloatResponseMatrixMaker' for c in channels},
    'pt' : {c:'FloatBranchResponseMatrixMaker' for c in channels},
    'eta' : {c:'AbsFloatBranchResponseMatrixMaker' for c in channels},
    'z1mass' : {'eeee':'FloatBranchResponseMatrixMaker',
                'mmmm':'FloatBranchResponseMatrixMaker',
                'eemm':'Z1ByMassResponseMatrixMaker',},
    'z2mass' : {'eeee':'FloatBranchResponseMatrixMaker',
                'mmmm':'FloatBranchResponseMatrixMaker',
                'eemm':'Z2ByMassResponseMatrixMaker',},
    'z1pt' : {'eeee':'FloatBranchResponseMatrixMaker',
              'mmmm':'FloatBranchResponseMatrixMaker',
              'eemm':'Z1ByMassResponseMatrixMaker',},
    'z2pt' : {'eeee':'FloatBranchResponseMatrixMaker',
              'mmmm':'FloatBranchResponseMatrixMaker',
              'eemm':'Z2ByMassResponseMatrixMaker',},
    #'zHigherPt' : {c:'Z1ByPtResponseMatrixMaker' for c in channels},
    #'zLowerPt' : {c:'Z2ByPtResponseMatrixMaker' for c in channels},
    'leppt' : {c:'AllLeptonBranchResponseMatrixMaker' for c in channels},
    #'l1Pt' : {c:'LeptonMaxBranchResponseMatrixMaker' for c in channels},
    'zpt' : {c:'BothZsBranchResponseMatrixMaker' for c in channels},
    'zmass' : {c:'BothZsBranchResponseMatrixMaker' for c in channels},
    'dphiz1z2':{c:'ZZAbsDeltaPhiResponseMatrixMaker' for c in channels},
    'drz1z2':{c:'ZZDeltaRResponseMatrixMaker' for c in channels},
    }

# Variable names usable by response maker classes
varNamesForResponseMaker = {
    'mass' : {c:'Mass' for c in channels},
    #'massFull' : {c:'Mass' for c in channels},
    'pt' : {c:'Pt' for c in channels},
    'eta' : {c:'Eta' for c in channels},
    'z1mass' : {'eeee':'e1_e2_Mass','mmmm':'m1_m2_Mass','eemm':'Mass'}, # 4e/4mu just use 1 variable because that's easy
    'z2mass' : {'eeee':'e3_e4_Mass','mmmm':'m3_m4_Mass','eemm':'Mass'}, # for 2e2mu, the response maker class will figure it out
    'z1pt' : {'eeee':'e1_e2_Pt','mmmm':'m1_m2_Pt','eemm':'Pt'}, # 4e/4mu just use 1 variable because that's easy
    'z2pt' : {'eeee':'e3_e4_Pt','mmmm':'m3_m4_Pt','eemm':'Pt'}, # for 2e2mu, the response maker class will figure it out
    'zpt' : {c:'Pt' for c in channels},
    'zmass' : {c:'Mass' for c in channels},
    #'zHigherPt' : {c:'Pt' for c in channels},
    #'zLowerPt' : {c:'Pt' for c in channels},
    'leppt' : {c:'Pt' for c in channels},
    #'l1Pt' : {c:'Pt' for c in channels},
    #'zPt' : {c:'Pt' for c in channels},
    'dphiz1z2': {c:'' for c in channels}, #variable names already set in ResponseMatrix.cxx
    'drz1z2': {c:'' for c in channels},
}

with open('varsFile.json') as var_json_file:
    myvar_dict = json.load(var_json_file)

# list of variables not counting systematic shifts
varList=['Mass'] #['Mass','ZZPt','ZPt','LepPt','dPhiZ1Z2','dRZ1Z2'] #With original list, histograms will be searched for all variables regardless of whether they are in runVariables
varNames={'mass': 'Mass','pt':'ZZPt','zpt':'ZPt','leppt':'LepPt','dphiz1z2':'dPhiZ1Z2','drz1z2':'dRZ1Z2'}

for key in myvar_dict.keys(): #key is the variable
    _binning[key] = myvar_dict[key]["_binning"]
    units[key] = myvar_dict[key]["units"]
    prettyVars[key] = myvar_dict[key]["prettyVars"]
    responseClassNames[key] = {c:myvar_dict[key]["responseClassNames"] for c in channels}
    if key == "MassAllj":
        varNamesForResponseMaker[key] = {c:'Mass' for c in channels}
        varNames[key] = 'Mass'
    else:
        varNamesForResponseMaker[key] = {c:str(key) for c in channels}
        varList.append(str(key))
        varNames[key] = str(key)

varNamesCopy = copy.deepcopy(varNames) #varNames is also defined and used in a function. Just in case, make a copy and put this in the function

def getLumiTextBox():
    texS = ROOT.TLatex(0.76,0.955, str(args['lumi'])+" fb^{-1} (13 TeV)")
    texS.SetNDC()
    texS.SetTextFont(42)
    texS.SetTextSize(0.030)
    texS.Draw()
    texS1 = ROOT.TLatex(0.11,0.955,"#bf{CMS} #it{Preliminary}")
    texS1.SetNDC()
    texS1.SetTextFont(42)
    texS1.SetTextSize(0.030)
    texS1.Draw()
    return texS,texS1


def generateAnalysisInputs():    
    #dictionary of SF histograms
    hSF = {}
    year = analysis[4:]
    eLowRecoFile = ROOT.TFile.Open('data/Ele_Reco_LowEt_%s.root' % (year)) 
    hSF['eLowReco'] = eLowRecoFile.Get('EGamma_SF2D').Clone()
    hSF['eLowReco'].SetDirectory(0)
    eLowRecoFile.Close()
    
    eRecoFile = ROOT.TFile.Open('data/Ele_Reco_%s.root' % (year))
    hSF['eReco'] = eRecoFile.Get('EGamma_SF2D').Clone()
    hSF['eReco'].SetDirectory(0)
    eRecoFile.Close()
    
    eIdFile = ROOT.TFile.Open('data/ElectronSF_Legacy_%s_NoGap.root' % (year))
    hSF['eSel'] = eIdFile.Get('EGamma_SF2D').Clone() 
    hSF['eSel'].SetDirectory(0)
    eIdFile.Close()

    eIdGapFile = ROOT.TFile.Open('data/ElectronSF_Legacy_%s_Gap.root' % (year))
    hSF['eSelGap'] = eIdGapFile.Get('EGamma_SF2D').Clone() 
    hSF['eSelGap'].SetDirectory(0)
    eIdGapFile.Close()
    
    if year=="2016":
        mIdFile = ROOT.TFile.Open('data/final_HZZ_SF_2016_legacy_mupogsysts_newLoose_noTracking_1610.root')
    elif year=="2017":
        mIdFile = ROOT.TFile.Open('data/ScaleFactors_mu_Moriond2018_final.root')
    elif year=="2018":
        mIdFile = ROOT.TFile.Open('data/final_HZZ_muon_SF_2018_IsBDT_0610.root')
    hSF['m'] = mIdFile.Get('FINAL').Clone()
    hSF['m'].SetDirectory(0)

    hSF['mErr'] = mIdFile.Get('ERROR').Clone()
    hSF['mErr'].SetDirectory(0)
    mIdFile.Close()

    #dictionary of PU weights
    hPU={}
    pileupFile = ROOT.TFile.Open('data/PileupWeights%s/PU_Central.root' % (year))
    hPU[''] = pileupFile.Get('pileup')
    hPU[''].SetDirectory(0)
    pileupFile.Close()
    
    pileupFileUp = ROOT.TFile.Open('data/PileupWeights%s/PU_minBiasUP.root' % (year))
    hPU['Up'] = pileupFileUp.Get('pileup')
    hPU['Up'].SetDirectory(0)
    pileupFileUp.Close()

    pileupFileDown = ROOT.TFile.Open('data/PileupWeights%s/PU_minBiasDOWN.root' % (year))
    hPU['Down'] = pileupFileDown.Get('pileup')
    hPU['Down'].SetDirectory(0)
    pileupFileDown.Close()

    return hSF,hPU
#ROOT.gSystem.Load('Utilities/scripts/ResponseMatrixMaker_cxx')
#sigSamples is a dictionary containing sample names and kfactor*cross-section
#sumW is a dictionary with sigsample:sumweights stored
ROOT.gSystem.Load('Utilities/scripts/ResponseMatrixMaker_cxx')
def generateResponseClass(varName, channel,sigSamples,sigSamplesPath,sumW,hPUWt,hSF={}):
    
    className = responseClassNames[varName][channel]

    for h in hSF.values()+hPUWt.values():
        ROOT.SetOwnership(h,False)

    if hSF:
        className = 'SFHist'+className

    #if not hasattr(ROOT,className):
    #        ROOT.gSystem.Load('Utilities/scripts/ResponseMatrixMaker_cxx','kTRUE')
    
    #for example C=<class 'ROOT.BranchValueResponseMatrixMaker<float>'>     
    C = getattr(ROOT, className)
    #commented_print("className:",C)
    
    #filelist=["zz4l-powheg"]
    filelist=[str(i) for i in sigSamples.keys()] 
    #improve this by getting this info from ZZDatasetManager just like its done in makeCompositeHists
    #sigConstWeights = {sample : (1.256*35900*1.0835)/sumW
    #                   for sample in ConfigureJobs.getListOfFiles(filelist, selection)}
   
    sigConstWeights = {sample : (sigSamples[sample.split("__")[0]]*1000*args['lumi'])/sumW[sample]
                       for sample in [str(i) for i in sigSamples.keys()] }
    #print "sigConstWeights: ",sigConstWeights
    #print("_binning: ",_binning)
    binning = _binning[varName]
    vBinning = VFloat()
    ROOT.SetOwnership(vBinning,False)
    #print("vBinning: ",vBinning)
    #print("Content of the ROOT vector object: {}".format([x for x in vBinning]))
    #print("binning: ",binning)
    if len(binning) == 3:
        binningTemp = [binning[1] + i * (binning[2] - binning[1])/float(binning[0]) for i in xrange(binning[0]+1)]
        #print("binningTemp: ",binningTemp)
        for b in binningTemp:
            #print("b: ",b)
            vBinning.push_back(b)
    else:
        for b in binning:
            vBinning.push_back(b)

    #commented_print("Content of the ROOT vector object: {}".format([x for x in vBinning]))
    #print("vBinning: ",vBinning)
    responseMakers = {}
    #for sample, file_path in sigFileNames.items():
    #for sample in ConfigureJobs.getListOfFiles(filelist,selection):
    year = int(analysis[4:])
    for sample in sigSamplesPath.keys():
        if sample==myaltname:
            continue
        #print "sample:", sample #expect zz4l-powheg
        #file_path = ConfigureJobs.getInputFilesPath(sample,selection,analysis)
        file_path=sigSamplesPath[sample]
        #print("where are the histos leaking")
        resp = C(channel, varNamesForResponseMaker[varName][channel], vBinning)
        #ROOT.SetOwnership(resp,False)
        file_path=file_path.encode("utf-8")
        #print "file_path: ",file_path
        fileList=glob.glob(file_path)
        #print "type in fileList should be str: ",type(fileList[0])
        for fName in fileList:
            resp.registerFile(fName)
        resp.registerPUWeights(hPUWt[''])
        resp.registerPUWeights(hPUWt['Up'], 'Up')
        resp.registerPUWeights(hPUWt['Down'], 'Down')
        resp.setConstantScale(sigConstWeights[sample])

        resp.setYear(year)
        if hSF:
            resp.registerElectronSelectionSFHist(hSF['eSel'])
            resp.registerElectronSelectionGapSFHist(hSF['eSelGap'])
            resp.registerElectronLowRecoSFHist(hSF['eLowReco'])
            resp.registerElectronRecoSFHist(hSF['eReco'])
            resp.registerMuonSFHist(hSF['m'])
            resp.registerMuonSFErrorHist(hSF['mErr'])
            #print("scale factors are being added")

        responseMakers[sample] = resp

    altResponseMakers = {}
    for sample in [myaltname]:#["zz4l-amcatnlo"]:
        #we only need to make new responseMatrix for zz4l-amcatnlo, ggZZ responseMatrices are already done above. 
        file_path=sigSamplesPath[sample]
        #print("where are the histos leaking")
        resp = C(channel, varNamesForResponseMaker[varName][channel], vBinning)
        #ROOT.SetOwnership(resp,False)
        file_path=file_path.encode("utf-8")
        #print "file_path: ",file_path
        fileList=glob.glob(file_path)
        #print "type in fileList should be str: ",type(fileList[0])
        for fName in fileList:
            resp.registerFile(fName)
        resp.registerPUWeights(hPUWt[''])
        resp.registerPUWeights(hPUWt['Up'], 'Up')
        resp.registerPUWeights(hPUWt['Down'], 'Down')
        resp.setConstantScale(sigConstWeights[sample])
        resp.setYear(year)
        if hSF:
            resp.registerElectronSelectionSFHist(hSF['eSel'])
            resp.registerElectronSelectionGapSFHist(hSF['eSelGap'])
            resp.registerElectronLowRecoSFHist(hSF['eLowReco'])
            resp.registerElectronRecoSFHist(hSF['eReco'])
            resp.registerMuonSFHist(hSF['m'])
            resp.registerMuonSFErrorHist(hSF['mErr'])
            #print("scale factors are being added")

        altResponseMakers[sample] = resp
        #print "resp: ",resp
        #del resp

    for sample in responseMakers.keys():
        print "sigSamples: " ,sample
    
    for sample in altResponseMakers.keys():
        print "altsigSamples: " ,sample

    for Resp in responseMakers.values()+altResponseMakers.values():
        ROOT.SetOwnership(Resp,False)
    
    return responseMakers,altResponseMakers

_printCounter = 0
#Load the RooUnfold library into ROOT
ROOT.gSystem.Load("RooUnfold/libRooUnfold")
def unfold(varName,chan,responseMakers,altResponseMakers,hSigDic,hAltSigDic,hSigSystDic,hTrueDic,hTrueSystDic_qqZZonly,hAltTrueDic,hDataDic,hbkgDic,hbkgMCDic,hbkgMCSystDic,nIter,plotDir=''):
    global _printCounter
    #get responseMakers from the function above- this is the whole game.
    #responseMakers = generateResponseClass(varName, chan,sigSamples,sumW,hSF)
    #if chan == 'mmmm':
    #    pdb.set_trace()
    # outputs
    hUnfolded = {}
    hTruth={}
    hTrueAlt = {}
    hResponseNominal={}
    print "responseMakers: ",responseMakers
    hResponseNominal = {s:resp for s,resp in responseMakers.items()}
    print "hResponseNominal:",hResponseNominal
    
    #Setup() is called here for all signals?
    hResponseSig1 = hResponseNominal["ggZZ4e"].getResponse("pu_Up")
    hResponseSig2 = hResponseNominal["ggZZ4m"].getResponse("pu_Up")
    hResponseSig3 = hResponseNominal["ggZZ4t"].getResponse("pu_Up")
    hResponseSig4 = hResponseNominal["ggZZ2e2tau"].getResponse("pu_Up")
    hResponseSig5 = hResponseNominal["ggZZ2e2mu"].getResponse("pu_Up")
    hResponseSig6 = hResponseNominal[mynominalName].getResponse("pu_Up")
    #This will pop the powheg response matrix from the hResponseNominal Dictionary
    hResponseNominalTotal = hResponseNominal.pop(mynominalName)
    #print "hRespNominalTotal: ",hResponseNominalTotal
    #This gets us the response matrix as a TH2D for "zz4l-powheg"
    #commented_print "This hResponse is full of leaks here"
    hResponse = hResponseNominalTotal.getResponse('nominal').Clone()
    hResponse.SetDirectory(0)
    #ROOT.SetOwnership(hResponse,False)
    #print("where are all the leaks") 
    #print "type of hResp: " ,hResponse
    #Now we need to add the rest of the response matrices (MCFMs) to this POHWEG matrix
    #Looping over the values of the dictionary (it doesn't have powheg anymore)
    #print "hResponseNominal after zz-pohwheg:",hResponseNominal
    for response in hResponseNominal.values():
        #commented_print "Is the leak here"
        #commented_print "response: ",response
        respMat = response.getResponse('nominal').Clone()
        #ROOT.SetOwnership(respMat,False)
        hResponse.Add(respMat)
        respMat.SetDirectory(0)
        #commented_print "Is the leak where"
        #ROOT.SetOwnership(respMat,False)  
        del respMat
        #respMat.Delete()

    #print ("The leaks happen in this for loop")
    #hResponseNominalTotal = sum(resp for resp in hResponseNominal.values())
    

    #print "type of Total hResp: " ,hResponse
    # But its better to use RooUnfoldResponse here
    #RooUnfoldResponse constructor - create from already-filled histograms
    # "response" gives the response matrix, measured X truth.
    #  "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
    #   but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
    #    in "truth" for unmeasured events (inefficiency).
    #     "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
    #      to indicate, respectively, no fakes and/or no inefficiency.


    ## Give hSig and hTrue in the form of histograms
    global varNamesCopy    #use a copy from global variable instead
    varNames= varNamesCopy #{'mass': 'Mass','pt':'ZZPt','zpt':'ZPt','leppt':'LepPt','dphiz1z2':'dPhiZ1Z2','drz1z2':'dRZ1Z2'}
    #varNames={'zmass':'ZMass','mass': 'Mass','pt':'ZZPt','eta':'ZZEta','z1mass':'Z1Mass','z1pt':'Z1Pt','z2mass':'Z2Mass','z2pt':'Z2Pt','zpt':'ZPt','leppt':'LepPt'}

    hSigNominal = hSigDic[chan][varNames[varName]]
    #print "sigHist: ", hSigNominal,", ",hSigNominal.Integral()
    #pdb.set_trace()
    hTrue = hTrueDic[chan]["Gen"+varNames[varName]]
    #histTrue.Scale((1.256*35900*1.0835)/zzSumWeights) 
    hData = hDataDic[chan][varNames[varName]]
    #print "dataHist: ",hData,", ",hData.Integral()
    #Get the background hists - #Get the histName_Fakes_chan histos
    hBkgNominal = hbkgDic[chan][varNames[varName]+"_Fakes"].Clone()
    #print "NonPromptHist: ",hBkgNominal,", ",hBkgNominal.Integral()
    hBkgMCNominal = hbkgMCDic[chan][varNames[varName]].Clone()
    #print "VVVHist: ",hBkgMCNominal,", ",hBkgMCNominal.Integral()
    #Add the two backgrounds

    hBkgNominal = rebin(hBkgNominal,varName)
    truncateTH1(hBkgNominal)
    hBkgMCNominal = rebin(hBkgMCNominal,varName)
    hBkgTotal=hBkgNominal.Clone()
    hBkgTotal.Add(hBkgMCNominal)
    #print "TotBkgHist: ",hBkgTotal,", ",hBkgTotal.Integral()

    #No need to rebin certain variables
    #if varNames[varName] not in ['ZZEta']:
        #bins=array.array('d',_binning[varName])
        #Nbins=len(bins)-1 
        #hSigNominal=hSigNominal.Rebin(Nbins,"",bins)
        #hTrue=hTrue.Rebin(Nbins,"",bins)
        #hData=hData.Rebin(Nbins,"",bins)
        #hBkgTotal=hBkgTotal.Rebin(Nbins,"",bins)
    #using the rebin function which takes care of overflow bins
    hSigNominal=rebin(hSigNominal,varName)
    hTrue=rebin(hTrue,varName)
    hData=rebin(hData,varName)
    #hBkgTotal=rebin(hBkgTotal,varName)
    
    xaxisSize = hSigNominal.GetXaxis().GetTitleSize()
    yaxisSize = hTrue.GetXaxis().GetTitleSize()
    #print "xaxisSize: ",xaxisSize
    #commented_print "trueHist: ",hTrue,", ",hTrue.Integral()
    #print "TotBkgHist after rebinning: ",hBkgTotal,", ",hBkgTotal.Integral()
    hTruth['']=hTrue

    hUnfolded[''], hCov, hResp = getUnfolded(hSigNominal.Clone(),hBkgTotal,hTruth[''],hResponse,hData, nIter,withRespAndCov=True,isNom=True)
    print("Position Indicator: nominal")
    printTH1(hData,'data nominal')
    printTH1(hSigNominal,'sig nominal')
    printTH1(hBkgTotal,'bkg nominal')
    printTH1(hTruth[''],'truth nominal')
    printTH2(hResponse,'response matrix nominal')
    printTH1(hUnfolded[''],'unfolded results nominal')
    printTH1Error(hUnfolded[''],'unfolded errors nominal')
    #print "hUnfolded['']: ",hUnfolded[''].Integral()
    
    #print("hResp: ",hResp) 
    #del hResponse
    # plot covariance and response
    if plotDir and args['plotResponse']:
        cRes = ROOT.TCanvas("c","canvas")
        if varName == 'massFull':
            cRes.SetLogx()
            cRes.SetLogy()
        draw_opt = "colz text45"
        hResp.GetXaxis().SetTitle('Reco '+prettyVars[varName]+''+units[varName])
        hResp.GetXaxis().SetTitleOffset(1.2)
        hResp.GetYaxis().SetTitle('True '+prettyVars[varName]+''+units[varName])
        hResp.GetXaxis().SetTitleSize(0.75*xaxisSize)
        hResp.GetYaxis().SetTitleSize(0.75*yaxisSize)
        hResp.Draw(draw_opt)
        texS,texS1=getLumiTextBox()
        #style.setCMSStyle(cRes, '', dataType='  Preliminary', intLumi=35900.)
        #print "plotDir: ",plotDir 
        plotName="response_%s" % (varName)
        output_name="/".join([plotDir,plotName])
        #cRes.Print("%s/response_%s.png" % (plotDir,varName))
        cRes.Print(output_name+".eps")
        cRes.Print(output_name+".png") 
        subprocess.call(["epstopdf", "--outfile=%s" % output_name+".pdf", output_name+".eps"],env={})
        os.remove(output_name+".eps")

        del cRes

        cCov = ROOT.TCanvas("c","canvas")
        if varName == 'massFull':
            cCov.SetLogx()
            cCov.SetLogy()
        draw_opt = "colz text45"
        hCov.Draw(draw_opt)
        texS,texS1=getLumiTextBox()
        covName="covariance_%s" % (varName)
        cov_name="/".join([plotDir,covName])
        cCov.Print(cov_name+".eps")
        cCov.Print(cov_name+".png") 
        subprocess.call(["epstopdf", "--outfile=%s" % cov_name+".pdf", cov_name+".eps"],env={})
        os.remove(cov_name+".eps")
        #style.setCMSStyle(cCov, '', dataType='Preliminary', intLumi=35900.)
        #cCov.Print("%s/covariance_%s.pdf" % (plotDir,varName))
        del cCov
    if not args['noSyst']: 
        #ggZZ xsec
        global hTrueDic_ggZZup
        global hTrueDic_ggZZdn
        global hSigDic_ggZZup
        global hSigDic_ggZZdn

        xsecScale = {'Up':1.+0.18,'Down':1.-0.14}
        for sys, scale in xsecScale.iteritems():
            #print "lumi uncert.",sys
            #print "scale: ",scale
            
            if "Up" in sys: #no need to clone, only used once
                hSigGX = hSigDic_ggZZup[chan][varNames[varName]]
                hTrueGX= hTrueDic_ggZZup[chan]["Gen"+varNames[varName]]
            else:
                hSigGX = hSigDic_ggZZdn[chan][varNames[varName]]
                hTrueGX= hTrueDic_ggZZdn[chan]["Gen"+varNames[varName]]
            hSigGX.SetDirectory(0)
            hTrueGX.SetDirectory(0)

            hBkgGX = hbkgDic[chan][varNames[varName]+"_Fakes"].Clone() 
            hBkgGX.SetDirectory(0)
            hBkgMCGX = hbkgMCDic[chan][varNames[varName]].Clone()
            hBkgMCGX.SetDirectory(0)
            hBkgGX = rebin(hBkgGX,varName)
            truncateTH1(hBkgGX)
            hBkgMCGX = rebin(hBkgMCGX,varName)
            hBkgTotalGX=hBkgGX.Clone() #clone in case
            hBkgTotalGX.Add(hBkgMCGX)

            hResponseGX = {s:resp for s,resp in responseMakers.items()}
            hRespGXTot = hResponseGX.pop(mynominalName)
            
            hRespGX = hRespGXTot.getResponse('nominal').Clone()
            hRespGX.SetDirectory(0)
            for resp in hResponseGX.values(): #ggZZ content
                respMatGX = resp.getResponse('nominal').Clone()
                respMatGX.Scale(scale)
                hRespGX.Add(respMatGX)
                respMatGX.SetDirectory(0)
                del respMatGX

            
            hSigGX=rebin(hSigGX,varName)
            hTrueGX=rebin(hTrueGX,varName)
            #hBkgTotalGX=rebin(hBkgTotalGX,varName)
            

            #print "trueHist: ",hTrueLumiShift,", ",hTrueLumiShift.Integral()
            #print "TotBkgHistLumi after rebinning: ",hBkgTotalLumi,", ",hBkgTotalLumi.Integral()
            #pdb.set_trace()
            hUnfolded['ggZZxsec_'+sys] = getUnfolded(hSigGX,hBkgTotalGX,hTrueGX,hRespGX,hData, nIter)
            print("Position Indicator:"+'ggZZxsec_'+sys)            
            printTH1(hData,'data '+'ggZZxsec_'+sys)
            printTH1(hSigGX,'sig '+'ggZZxsec_'+sys)
            printTH1Error(hSigGX,'sig Error '+'ggZZxsec_'+sys)
            printTH1(hBkgTotalGX,'bkg '+'ggZZxsec_'+sys)
            printTH1(hBkgGX,'bkg nonprompt only '+'ggZZxsec_'+sys)
            printTH1(hTrueGX,'truth '+'ggZZxsec_'+sys)
            printTH2(hRespGX,'response matrix '+'ggZZxsec_'+sys)
            del hSigGX
            del hBkgMCGX
            del hTrueGX
            del hBkgGX
            del hRespGX

# fake rate
        fakeUnc = 0.4
        fakeScale = {'Up':1.+fakeUnc,'Down':1.-fakeUnc}
        for sys, scale in fakeScale.iteritems():
            #print "fake uncert.",sys
            #print "scale: ",scale
            hSigFake = hSigNominal.Clone()
            hSigFake.SetDirectory(0)

            hBkgFake = hbkgDic[chan][varNames[varName]+"_Fakes"].Clone()
            hBkgFake.Scale(scale)
            hBkgFake.SetDirectory(0)
            hBkgFake=rebin(hBkgFake,varName)
            truncateTH1(hBkgFake)
            hBkgMCFake = hbkgMCDic[chan][varNames[varName]].Clone()
            hBkgMCFake.SetDirectory(0)
            hBkgMCFake=rebin(hBkgMCFake,varName)
            hBkgTotalFake=hBkgFake.Clone()
            hBkgTotalFake.Add(hBkgMCFake)

            hTrueFakeShift = hTruth[''].Clone()
            hTrueFakeShift.SetDirectory(0)
            hResponseFake = hResponse.Clone()
            hResponseFake.SetDirectory(0)
            #print "SigFakeHist: ",hSigFake,", ",hSigFake.Integral()
            #print "VVVHist: ",hBkgMCFake,", ",hBkgMCFake.Integral()
            #Add the two backgrounds

            #hBkgTotalFake=rebin(hBkgTotalFake,varName)

            #print "trueHist: ",hTrueFakeShift,", ",hTrueFakeShift.Integral()
            #print "TotBkgHistFake after rebinning: ",hBkgTotalFake,", ",hBkgTotalFake.Integral()

            hUnfolded['fake_'+sys],hCovFake,hRespFake = getUnfolded(hSigFake,hBkgTotalFake,hTrueFakeShift,hResponseFake,hData, nIter,withRespAndCov=True)
            #hUnfolded['fake_'+sys] = getUnfolded(hSigFake,hBkgTotalFake,hTrueFakeShift,hResponseFake,hData, nIter)
            print("Position Indicator:"+'fake_'+sys)
            printTH1(hData,'data '+'fake_'+sys)
            printTH1(hSigFake,'sig '+'fake_'+sys)
            printTH1(hBkgTotalFake,'bkg '+'fake_'+sys)
            printTH1(hTrueFakeShift,'truth '+'fake_'+sys)
            printTH2(hResponseFake,'response matrix '+'fake_'+sys)
            del hSigFake
            del hBkgMCFake
            del hBkgFake
            del hTrueFakeShift

        # luminosity
        lumiUnc = 0.023
        lumiScale = {'Up':1.+lumiUnc,'Down':1.-lumiUnc}
        for sys, scale in lumiScale.iteritems():
            #print "lumi uncert.",sys
            #print "scale: ",scale
            hSigLumi = hSigNominal * scale
            hSigLumi.SetDirectory(0)

            hBkgLumi = hbkgDic[chan][varNames[varName]+"_Fakes"].Clone()
            hBkgLumi.SetDirectory(0)
            hBkgLumi=rebin(hBkgLumi,varName)
            truncateTH1(hBkgLumi)
            hBkgMCLumi = hbkgMCDic[chan][varNames[varName]].Clone()
            hBkgMCLumi.SetDirectory(0)
            hBkgMCLumi.Scale(scale)
            hBkgMCLumi=rebin(hBkgMCLumi,varName)
            hBkgTotalLumi=hBkgLumi.Clone()
            hBkgTotalLumi.Add(hBkgMCLumi)

            hTrueLumiShift = hTruth[''] * scale
            hTrueLumiShift.SetDirectory(0)
            hResponseLumi = hResponse.Clone()
            hResponseLumi.Scale(scale)
            hResponseLumi.SetDirectory(0)
            #print "SigLumiHist: ",hSigLumi,", ",hSigLumi.Integral()
            #print "VVVHist: ",hBkgMCLumi,", ",hBkgMCLumi.Integral()
            #Add the two backgrounds

            #hBkgTotalLumi=rebin(hBkgTotalLumi,varName)

            #print "trueHist: ",hTrueLumiShift,", ",hTrueLumiShift.Integral()
            #print "TotBkgHistLumi after rebinning: ",hBkgTotalLumi,", ",hBkgTotalLumi.Integral()

            hUnfolded['lumi_'+sys],hCovLumi,hRespLumi = getUnfolded(hSigLumi,hBkgTotalLumi,hTrueLumiShift,hResponseLumi,hData, nIter,withRespAndCov=True)
            print("Position Indicator:"+'lumi_'+sys)
            printTH1(hData,'data '+'lumi_'+sys)
            printTH1(hSigLumi,'sig '+'lumi_'+sys)
            printTH1(hBkgTotalLumi,'bkg '+'lumi_'+sys)
            printTH1(hTrueLumiShift,'truth '+'lumi_'+sys)
            printTH2(hResponseLumi,'response matrix '+'lumi_'+sys)
            del hSigLumi
            del hBkgMCLumi
            del hBkgLumi
            del hTrueLumiShift

            #print "hUnfolded['lumi_']: ",hUnfolded['lumi_'+sys].Integral()

            if plotDir and args['plotSystResponse']:
                cResLumi = ROOT.TCanvas("c","canvas",800,800)
                if varName == 'massFull':
                    cRes.SetLogx()
                    cRes.SetLogy()
                draw_opt = "colz text45"
                hRespLumi.GetXaxis().SetTitle('Reco '+prettyVars[varName]+''+units[varName])
                hRespLumi.GetYaxis().SetTitle('True '+prettyVars[varName]+''+units[varName])
                hRespLumi.GetXaxis().SetTitleSize(0.75*xaxisSize)
                hRespLumi.GetYaxis().SetTitleSize(0.75*yaxisSize)
                hRespLumi.Draw(draw_opt)
                style.setCMSStyle(cResLumi, '', dataType='  Preliminary', intLumi=35900.)
                #print "plotDir: ",plotDir
                cResLumi.Print("%s/response_%s_%s.png" % (plotDir,varName,'lumi'+sys))
                cResLumi.Print("%s/response_%s_%s.pdf" % (plotDir,varName,'lumi'+sys))
            
                del cResLumi
        
        hResponsePU = {s:resp for s,resp in responseMakers.items()}
        hRespPUTot = hResponsePU.pop(mynominalName)
        #commented_print "No errors in PU chain?"
        # PU reweight uncertainty
        for sys in ['Up','Down']:
            hRespPU = hRespPUTot.getResponse('pu_'+sys)
            hRespPU.SetDirectory(0)
            for resp in hResponsePU.values():
                respMatPU = resp.getResponse('pu_'+sys)
                hRespPU.Add(respMatPU)
                respMatPU.SetDirectory(0)
                del respMatPU
            #print "hSigSystDic: ",hSigSystDic
            hSigPU = hSigSystDic[chan][varNames[varName]+"_CMS_pileup"+sys]
            hSigPU.SetDirectory(0)
            #print 'pu_'+sys 
            #print "sigHist: ", hSigPU,", ",hSigPU.Integral()
            hBkgPU = hbkgDic[chan][varNames[varName]+"_Fakes"].Clone()
            hBkgPU.SetDirectory(0)
            hBkgPU=rebin(hBkgPU,varName)
            truncateTH1(hBkgPU)
            #print "NonPromptHist: ",hBkgPU,", ",hBkgPU.Integral()
            hBkgMCPU = hbkgMCSystDic[chan][varNames[varName]+"_CMS_pileup"+sys].Clone()
            hBkgMCPU.SetDirectory(0)
            hBkgMCPU=rebin(hBkgMCPU,varName)
            #print "VVVHist: ",hBkgMCPU,", ",hBkgMCPU.Integral()
            hBkgPUTotal=hBkgPU.Clone()
            hBkgPUTotal.Add(hBkgMCPU) #don't think hBkgMCPU needs to be cloned
            #print "TotBkgPUHist: ",hBkgPUTotal,", ",hBkgPUTotal.Integral()
            hSigPU=rebin(hSigPU,varName)
            #hBkgPUTotal=rebin(hBkgPUTotal,varName)
            #print "TotBkgPUHist after Rebinning: ",hBkgPUTotal,", ",hBkgPUTotal.Integral()
            
            hUnfolded['pu_'+sys] = getUnfolded(hSigPU,
                                                     hBkgPUTotal,
                                                     hTruth[''],
                                                     hRespPU,
                                                     hData, nIter)
            print("Position Indicator:"+'pu_'+sys)
            printTH1(hData,'data '+'pu_'+sys)
            printTH1(hSigPU,'sig '+'pu_'+sys)
            printTH1(hBkgPUTotal,'bkg '+'pu_'+sys)
            printTH1(hTruth[''],'truth '+'pu_'+sys)
            printTH2(hRespPU,'response matrix '+'pu_'+sys)
            del hSigPU
            del hBkgMCPU
            del hBkgPU
            del hRespPU
            
#Add systematics for JES and JER
        hResponseJET = {s:resp for s,resp in responseMakers.items()}
        hRespJETTot = hResponseJET.pop(mynominalName)
        
        hSigJET = hSigSystDic[chan][varNames[varName]+"_jetsysts"] #TH2
        hSigJET.SetDirectory(0)
            
        hBkgJET = hbkgDic[chan][varNames[varName]+"_Fakes"].Clone()
        hBkgJET.SetDirectory(0)
        hBkgJET=rebin(hBkgJET,varName)
        truncateTH1(hBkgJET)
            
        hBkgMCJET = hbkgMCSystDic[chan][varNames[varName]+"_jetsysts"].Clone() #TH2
        hBkgMCJET.SetDirectory(0)

        for i,sys in enumerate(["jes_up","jes_dn","jer_up","jer_dn"]):
            hRespJET = hRespJETTot.getResponse(sys)
            hRespJET.SetDirectory(0)
            for resp in hResponseJET.values():
                respMatJET = resp.getResponse(sys)
                hRespJET.Add(respMatJET)
                respMatJET.SetDirectory(0)
                del respMatJET
            
                    
            hSigJETsub = hSigJET.ProjectionX("JET_%s"%i,i+1,i+1,"e") #histogram bin count starts with 1
            hSigJETsub.SetDirectory(0)
            
            hBkgMCJETsub = hBkgMCJET.ProjectionX("JETBkg_%s"%i,i+1,i+1,"e")
            hBkgMCJETsub=rebin(hBkgMCJETsub,varName)
            hBkgJETTotal=hBkgJET.Clone()
            hBkgJETTotal.Add(hBkgMCJETsub)
            
            hSigJETsub=rebin(hSigJETsub,varName)
            #hBkgJETTotal=rebin(hBkgJETTotal,varName)
           
            
            hUnfolded[sys] = getUnfolded(hSigJETsub,
                                                     hBkgJETTotal,
                                                     hTruth[''],
                                                     hRespJET,
                                                     hData, nIter)
            print("Position Indicator:"+sys)
            printTH1(hData,'data '+sys)
            printTH1(hSigJETsub,'sig '+sys)
            printTH1(hBkgJETTotal,'bkg '+sys)
            printTH1(hTruth[''],'truth '+sys)
            printTH2(hRespJET,'response matrix '+sys)
            del hSigJETsub
            del hBkgMCJETsub
            del hRespJET
        del hBkgJET
        

#Add systematics for scales and PDF 
        #pdb.set_trace()
        hResponsePS = {s:resp for s,resp in responseMakers.items()} #sample and response classes
        hRespPSTot = hResponsePS.pop(mynominalName)
        # PDF and scale uncertainty
        nscales = 9 # six scale variation indices 1,2,3,4,6,8
        npdf = 103
        #indices run from 0 to 111 for nominal, scale, pdf+alpha_s for 2017,18. In 2016 no nominal PDF and last relevant index is 110.
        PSlist = []
        hRespPS = hRespPSTot.getScaleResponses() #vec<TH2D>
        for histPS in hRespPS:
            histPS.SetDirectory(0)
            PSlist.append(histPS)
        for resp in hResponsePS.values():
            respMatPS = resp.getResponse('nominal').Clone()
            for i, item in enumerate(PSlist):
                PSlist[i].Add(respMatPS)
            respMatPS.SetDirectory(0)
            del respMatPS
        
        #scale_pdf_syslist = ['scale%s'%i for i in range(1,nscales)] + ['pdf%s'%i for i in range(10,109)] + ['alphas_up','alphas_dn']
        hSigPSt = hSigSystDic_qqZZonly[chan][varNames[varName]+"_lheWeights"]
        hSigPSt.SetDirectory(0)
        hSig_ggZZonly = hSigDic_ggZZonly[chan][varNames[varName]]
        hSig_ggZZonly.SetDirectory(0)

        hTrueLHEt = hTrueSystDic_qqZZonly[chan]["Gen"+varNames[varName]+"_lheWeights"] #TH2
        hTrueLHEt.SetDirectory(0)
        hTrue_ggZZonly = hTrueDic_ggZZonly[chan]["Gen"+varNames[varName]]
        hTrue_ggZZonly.SetDirectory(0)
        #print 'pu_'+sys 
        #print "sigHist: ", hSigPS,", ",hSigPS.Integral()
        hBkgPSt = hbkgDic[chan][varNames[varName]+"_Fakes"].Clone()
        hBkgPSt.SetDirectory(0)
        hBkgPSt=rebin(hBkgPSt,varName)
        truncateTH1(hBkgPSt)
        #print "NonPromptHist: ",hBkgPS,", ",hBkgPS.Integral()
        hBkgMCPSt = hbkgMCDic[chan][varNames[varName]].Clone() #hbkgMCSystDic[chan][varNames[varName]+"_lheWeights"]
        hBkgMCPSt.SetDirectory(0)
        hBkgMCPSt=rebin(hBkgMCPSt,varName)
        hBkgPSTotal=hBkgPSt.Clone()
        hBkgPSTotal.Add(hBkgMCPSt)
        #hBkgPSTotal=rebin(hBkgPSTotal,varName)
        #print "VVVHist: ",hBkgMCPS,", ",hBkgMCPS.Integral()
        #hBkgPSTotalt=hBkgPSt.Clone()
        #hBkgPSTotalt.Add(hBkgMCPSt)
        #print "TotBkgPSHist: ",hBkgPSTotal,", ",hBkgPSTotal.Integral()
        #hSigPSt=rebin(hSigPSt,varName)
        #hBkgPSTotalt=rebin(hBkgPSTotalt,varName)
        year = analysis[4:]
        if year=='2018' or year=='2017':
            ind_last = 111
        
        if year=='2016':
            ind_last = 110
        

        for i in range(0,ind_last+1): #pick all indices other than nominal
            if i == 0: #or i == 9: 
                continue #only skip nominal 0 since 9 is not nominal PDF for 2016
            #print "hSigSystDic: ",hSigSystDic

            binnum = i + 1 # 1st bin corresponds to 0
            hSigPS = hSigPSt.ProjectionX("PS%s"%i,i+1,i+1,"e") #expect e option to instruct computing the error
            hSigPS.SetDirectory(0)
            hSigPS.Add(hSig_ggZZonly)
            #print 'pu_'+sys 
            #print "sigHist: ", hSigPS,", ",hSigPS.Integral()
            hTrueLHE = hTrueLHEt.ProjectionX("PS%s"%i,i+1,i+1,"e")
            hTrueLHE.SetDirectory(0)
            hTrueLHE.Add(hTrue_ggZZonly)
            #print "NonPromptHist: ",hBkgPS,", ",hBkgPS.Integral()
            
            #hBkgMCPS = hBkgMCPSt.ProjectionX("PSBkg%s"%i,i+1,i+1,"e")
            #hBkgMCPS.SetDirectory(0)
            
            #print "VVVHist: ",hBkgMCPS,", ",hBkgMCPS.Integral()
            
            #hBkgPSTotal=hBkgPSt.Clone()
            #hBkgPSTotal.Add(hBkgMCPS)
            
            #print "TotBkgPSHist: ",hBkgPSTotal,", ",hBkgPSTotal.Integral()
            hSigPS=rebin(hSigPS,varName)
            
            hTrueLHE = rebin(hTrueLHE,varName)
            #print "TotBkgPSHist after Rebinning: ",hBkgPSTotal,", ",hBkgPSTotal.Integral()
            
            hUnfolded['PS_'+str(i)] = getUnfolded(hSigPS,
                                                     hBkgPSTotal,
                                                     hTrueLHE,
                                                     PSlist[i],
                                                     hData, nIter)
            del hSigPS # This just deletes the variable not the histogram? What is the purpose?
            #del hBkgMCPS
            del hTrueLHE
            print("Position Indicator:"+'PS_'+str(i))
            #del hBkgPS
            #del hRespPS 
            #don't delete hBkgTotal in the original codes?
        del hBkgPSt
        for item in PSlist:
            del item
        

        # lepton efficiency uncertainty
        for lep in set(chan):
            hResponseSyst = {s:resp for s,resp in responseMakers.items()}
            hRespSystTot = hResponseSyst.pop(mynominalName)
            #commented_print "No errors in systematics chain?"
            for sys in ['Up','Down']:
                hRespSyst = hRespSystTot.getResponse(lep+'Eff_'+sys)
                hRespSyst.SetDirectory(0)
                for response in hResponseSyst.values():
                    respMatSyst = response.getResponse(lep+'Eff_'+sys)
                    hRespSyst.Add(respMatSyst)
                    respMatSyst.SetDirectory(0)
                    del respMatSyst
                #print "hSigSystDic: ",hSigSystDic
                hSigSyst = hSigSystDic[chan][varNames[varName]+"_CMS_eff_"+lep+sys]
                hSigSyst.SetDirectory(0)
                #print lep+'Eff_'+sys 
                #print "sigHist: ", hSigSyst,", ",hSigSyst.Integral()
                hBkgSyst = hbkgDic[chan][varNames[varName]+"_Fakes"].Clone()
                hBkgSyst.SetDirectory(0)
                hBkgSyst=rebin(hBkgSyst,varName)
                truncateTH1(hBkgSyst)
                #print "NonPromptHist: ",hBkgSyst,", ",hBkgSyst.Integral()
                hBkgMCSyst = hbkgMCSystDic[chan][varNames[varName]+"_CMS_eff_"+lep+sys].Clone()
                hBkgMCSyst.SetDirectory(0)
                hBkgMCSyst=rebin(hBkgMCSyst,varName)
                #print "VVVHist: ",hBkgMCSyst,", ",hBkgMCSyst.Integral()
                hBkgSystTotal=hBkgSyst.Clone()
                hBkgSystTotal.Add(hBkgMCSyst)
                #print "TotBkgSystHist: ",hBkgSystTotal,", ",hBkgSystTotal.Integral()
                hSigSyst=rebin(hSigSyst,varName)
                #hBkgSystTotal=rebin(hBkgSystTotal,varName)
                #print "TotBkgSystHist after Rebinning: ",hBkgSystTotal,", ",hBkgSystTotal.Integral()
                #if lep+'Eff_'+sys=="eEff_Up":
                #    pdb.set_trace()
                hUnfolded[lep+'Eff_'+sys] ,hCovLep,hRespLep = getUnfolded(hSigSyst,
                                                         hBkgSystTotal,
                                                         hTruth[''],
                                                         hRespSyst,
                                                         hData, nIter,withRespAndCov=True)
                print("Position Indicator:"+lep+'Eff_'+sys)
                printTH1(hData,'data '+lep+'Eff_'+sys)
                printTH1(hSigSyst,'sig '+lep+'Eff_'+sys)
                printTH1(hBkgSystTotal,'bkg '+lep+'Eff_'+sys)
                printTH1(hTruth[''],'truth '+lep+'Eff_'+sys)
                printTH2(hRespSyst,'response matrix '+lep+'Eff_'+sys)
                del hSigSyst
                del hBkgMCSyst
                del hBkgSyst
                del hRespSyst

                if plotDir and args['plotSystResponse']:
                    cResSyst = ROOT.TCanvas("c","canvas",1200,1200)
                    if varName == 'massFull':
                        cResSyst.SetLogx()
                        cResSyst.SetLogy()
                    draw_opt = "colz text45"
                    hRespLep.GetXaxis().SetTitle('Reco '+prettyVars[varName]+sys+''+units[varName])
                    hRespLep.GetYaxis().SetTitle('True '+prettyVars[varName]+sys+''+units[varName])
                    hRespLep.GetXaxis().SetTitleSize(0.75*xaxisSize)
                    hRespLep.GetYaxis().SetTitleSize(0.75*yaxisSize)
                    hRespLep.Draw(draw_opt)
                    style.setCMSStyle(cResSyst, '', dataType='  Preliminary', intLumi=35900.)
                    #print "plotDir: ",plotDir
                    cResSyst.Print("%s/response_%s_%s.png" % (plotDir,varName,lep+'Eff_'+sys))
                    cResSyst.Print("%s/response_%s_%s.pdf" % (plotDir,varName,lep+'Eff_'+sys))
                    del cResSyst
    
    hResponseDebug = hResponse.Clone()
    del hResponse 
    
    #Alternative signal zz4l-amcatnlo
    hResponseAltNominal={} 
    print "AltresponseMakers: ",altResponseMakers
    hResponseAltNominal = {s:resp for s,resp in altResponseMakers.items()}
    print "hResponseNominal:",hResponseNominal
    print "hResponseAltNominal:",hResponseAltNominal

    hResponseSig7 = hResponseAltNominal[myaltname].getResponse("pu_Up")
    #This will pop the amcnlo response matrix from the hResponseAltNominal Dictionary
    hResponseAltNominalTotal = hResponseAltNominal.pop(myaltname)

    hAltResponse = hResponseAltNominalTotal.getResponse('nominal').Clone()
    hAltResponse.SetDirectory(0)

    #Looping over the values of the dictionary (it doesn't have powheg anymore)
    #it only has the MFCFM signals in it so add them to the amcatnlo now!
    for response in hResponseNominal.values():
        print "response: ",response
        respMat = response.getResponse('nominal')
        hAltResponse.Add(respMat)
        respMat.SetDirectory(0)
        del respMat

    hAltSigNominal = hAltSigDic[chan][varNames[varName]]

    hAltTrue = hAltTrueDic[chan]["Gen"+varNames[varName]]

    hAltSigNominal=rebin(hAltSigNominal,varName)
    hAltTrue=rebin(hAltTrue,varName)

    #commented_print "AltTrueHist: ",hAltTrue,", ",hAltTrue.Integral()
    
    hTrueAlt['']=hAltTrue
    #if chan=='mmmm':
    #    pdb.set_trace()
    
    hUnfolded['generator']  = getUnfolded(hAltSigNominal,hBkgTotal,hTrueAlt[''],hAltResponse,hData, nIter) 
    print("Position Indicator:generator")
    printTH1(hData,'data '+'generator')
    printTH1(hAltSigNominal,'sig '+'generator')
    printTH1(hBkgTotal,'bkg '+'generator')
    printTH1(hTrueAlt[''],'truth '+'generator')
    printTH2(hAltResponse,'response matrix '+'generator')
    del hResponseDebug
    # make everything local (we'll cache copies)
    
    for h in hUnfolded.values()+hTruth.values()+hTrueAlt.values():
        ROOT.SetOwnership(h,False)
        #commented_print("histos: ",h)
        #commented_print ("hTruthOut of Unfold: ",h.Integral())
        #h.SetDirectory(0)
    #commented_print "hUnfolded (inside unfold): ",hUnfolded
    #commented_print "hTrueAlt (inside unfold): ",hTrueAlt
    return hUnfolded,hTruth,hTrueAlt,hData,hSigNominal,hBkgTotal

def printTH1(hist,label):
    contents = []
    for i in range(1,hist.GetNbinsX()+1):
        contents.append(hist.GetBinContent(i))
    time.sleep(0.2) #prevent mixing with potential warning message
    print('Diagnostic bin contents of %s:'%label)
    print(contents)

def truncateTH1(hist):
    
    for i in range(1,hist.GetNbinsX()+1):
        if hist.GetBinContent(i)<0.:
            tmperror = abs(hist.GetBinError(i)) #should be positive error, just in case
            hist.SetBinContent(i,0.)
            hist.SetBinError(i,tmperror)

def printTH1Error(hist,label):
    contents = []
    for i in range(1,hist.GetNbinsX()+1):
        contents.append(hist.GetBinError(i))
    time.sleep(0.2) #prevent mixing with potential warning message
    print('Diagnostic bin error of %s:'%label)
    print(contents)

def printTH2(hist,label):
    contents = []
    for i in range(1,hist.GetNbinsX()+1):
        for j in range(1,hist.GetNbinsX()+1): #same xy dimension
            contents.append(hist.GetBinContent(i,j))
    time.sleep(0.2) #prevent mixing with potential warning message
    print('Diagnostic matrix bin contents of %s:'%label)
    print(contents)


#rebin histos and take care of overflow bins
def rebin(hist,varName):
    ROOT.SetOwnership(hist, False)
    #No need to rebin certain variables but still might need overflow check
    if varName not in ['eta']:
        bins=array.array('d',_binning[varName])
        Nbins=len(bins)-1 
        hist=hist.Rebin(Nbins,"",bins)
    else:
        Nbins = hist.GetSize() - 2
    add_overflow = hist.GetBinContent(Nbins) + hist.GetBinContent(Nbins + 1)
    lastbin_error = math.sqrt(math.pow(hist.GetBinError(Nbins),2)+math.pow(hist.GetBinError(Nbins+1),2))
    hist.SetBinContent(Nbins, add_overflow)
    hist.SetBinError(Nbins, lastbin_error)
    hist.SetBinContent(Nbins+1,0)
    hist.SetBinError(Nbins+1,0)
    #if not hist.GetSumw2(): hist.Sumw2()
    return hist

def getUnfolded(hSig, hBkg, hTrue, hResponse, hData, nIter,withRespAndCov=False,isNom=False):
    Response = getattr(ROOT,"RooUnfoldResponse")
    specBkg = True #a separate bkg treatment according to D'Agostini
    clean0 = False #clean bins with negative value
    if not applyreg:
        specBkg = False
    #commented_print "TrueBeforeResponse: ", hTrue,", ",hTrue.Integral()
    #commented_print "SigBeforeResponse: ", hSig,", ",hSig.Integral()
    hSig = hSig.Clone() #reassign local variable
    if specBkg:
        hSig.Add(hBkg)
    response = Response(hSig, hTrue.Clone(), hResponse.Clone()) 
    ROOT.SetOwnership(response,False)
    ROOT.SetOwnership(hData,False)
    #Response matrix as a 2D-histogram: (x,y)=(measured,truth)
    hResp = response.Hresponse()
    #hResp = hResponse
    #print "hResp out of response: ",hResp

    #RooUnfoldIter = getattr(ROOT,"RooUnfoldBayes")

    RooUnfoldInv = getattr(ROOT,"RooUnfoldInvert")
    RooUnfoldBayes = getattr(ROOT,"RooUnfoldBayes")
    #RooUnfoldBinbyBin = getattr(ROOT,"RooUnfoldBinByBin")
    try:
        svd = ROOT.TDecompSVD(response.Mresponse())
        sig = svd.GetSig()
        try:
            condition = sig.Max() / max(0., sig.Min())
        except ZeroDivisionError:
            condition = float('inf')
            raise

        print("Printout record start here:")
        print "channel: ",chan
        print "variable: ",varNames[varName]
        print "hResp out of response: ",hResp
        print ''
        print 'condition: {}'.format(condition)
        print ''

    except:
        #commented_print "It broke! #commented_printing debug info"
        #commented_print "Sig: {}, bkg: {}, true: {}, response: {}".format(hSig.Integral(), hBkg.Integral(), hTrue.Integral(), hResponse.Integral())
        c = ROOT.TCanvas("c1","canvas",800,800)
        hSig.Draw()
        style.setCMSStyle(c, '', dataType='Debug', intLumi=35900.)
        c.Print("DebugPlots/sig{}.png".format(_printCounter))
        hBkg.draw()
        _style.setCMSStyle(c, '', dataType='Debug', intLumi=35900.)
        c.Print("bkg{}.png".format(_printCounter))
        hTrue.Draw()
        style.setCMSStyle(c, '', dataType='Debug', intLumi=35900.)
        c.Print("DebugPlots/true{}.png".format(_printCounter))
        hData.Draw()
        style.setCMSStyle(c, '', dataType='Debug', intLumi=35900.)
        c.Print("DebugPlots/data{}.png".format(_printCounter))
        draw_opt = "colz text45"
        hResponse.Draw(draw_opt)
        style.setCMSStyle(c, '', dataType='Debug', intLumi=35900.)
        c.Print("DebugPlots/resp{}.png".format(_printCounter))
        c.Print("DebugPlots/resp{}.root".format(_printCounter))
        _printCounter += 1

    #calculate chi-2 in smeared space for bottom-line test
    chi2_smear = 0.
    chi2_goodness = 0.
    mResp = response.Mresponse()
    exclude_bt_bin = []
    if isNom:
        
        for i in range(1,hData.GetNbinsX()+1):
            #if not hData.GetBinError(i) == 0.:
            error_denominator = hData.GetBinError(i) #don't fix zero division for now
            #elif not hSig.GetBinError(i) == 0.:
            #    error_denominator = hSig.GetBinError(i)
            #else:
            #    exclude_bt_bin.append(i)
            #    continue
            try:
                if specBkg:
                    chi2_smear += (hData.GetBinContent(i)-hSig.GetBinContent(i))**2/(error_denominator**2)
                else:
                    chi2_smear += (hData.GetBinContent(i)-hBkg.GetBinContent(i)-hSig.GetBinContent(i))**2/(error_denominator**2)
            except ZeroDivisionError:
                print("Zero division, chi2_smear set to inf")
                chi2_smear = float('inf')
                break

    #commented_print "hData: ", hData.Integral()
    hDataMinusBkg = hData.Clone()
    hDataMinusBkg.Reset()
    #commented_print "hBkg: ", hBkg.Integral()
    hDataMinusBkg.Add(hData,1)
    if not specBkg:
        hDataMinusBkg.Add(hBkg,-1)
    if clean0:
        for ind in range(1,hDataMinusBkg.GetNbinsX()+1):
            if hDataMinusBkg.GetBinContent(ind) <0.:
                hDataMinusBkg.SetBinContent(ind,0.)
    #hDataMinusBkg.Add(hBkg,-1)
    #HistTools.zeroNegativeBins(hDataMinusBkg)
    #commented_print "DataMinusbkgIntegral: ",hDataMinusBkg, ", ",hDataMinusBkg.Integral()
    #Unfolding using 4 iterations and then stopping
    #if varNames[varName] not in ["Z1Mass","Z2Mass"]:
    #    nIter=8
    #print "No.of iterations: ",nIter
    #commented_print "response: ",response

    #Simply inverting the matrix
    if not applyreg:
        unf = RooUnfoldInv(response, hDataMinusBkg)
    else:
        unf3 = RooUnfoldBayes(response, hDataMinusBkg,4)
    #unf3 = RooUnfoldBayes(response, hSig,4)
    #unf = RooUnfoldIter(response, hDataMinusBkg, nIter)
    #commented_print "unf: ",unf 

    #Unfolds using the method of correction factors
    #unf = RooUnfoldBinbyBin(response, hSig)

    #del hDataMinusBkg
    #This is the unfolded "data" distribution
    if not applyreg:
        hOut = unf.Hreco()
    else:
        hOut3=unf3.Hreco()
    #ROOT.SetOwnership(hOut,False)
    #if not hOut:
        #commented_print hOut
    #    raise ValueError("The unfolded histogram got screwed up somehow!")
    #commented_print("hOut: ",hOut,"",hOut.Integral())
    #commented_print("Check first bin simple:",hOut.GetBinContent(1))
    #commented_print("Check first bin reg:",hOut.GetBinContent(1))
    #Returns covariance matrices for error calculation of type withError
    #0: Errors are the square root of the bin content
    #1: Errors from the diagonals of the covariance matrix given by the unfolding
    #2: Errors from the covariance matrix given by the unfolding => We use this one for now
    #3: Errors from the covariance matrix from the variation of the results in toy MC tests
    if not applyreg:
        hCov = unf.Ereco(2)
    else:
        hCov3=unf3.Ereco(2)

    #Now calculate chi-2 at unfolded space for bottom-line test
    if applyreg:
        hCovInv = hCov3.Clone()
        hresult = hOut3.Clone()
    else:
        hCovInv = hCov.Clone()
        hresult = hOut.Clone()
    
    chi2_unf = 0.
    try:
        hCovInv.Invert()
    except:
        print("Covariance matrix not invertible,chi2_unf set to -1")
        chi2_unf = -1.

    if isNom:
        for nrow in range(1,hresult.GetNbinsX()+1):
            tmpSmear = 0.
            for ncol in range(1,hresult.GetNbinsX()+1):
                tmpSmear += mResp(nrow-1,ncol-1)* hresult.GetBinContent(ncol)
            try:
                chi2_goodness += (hData.GetBinContent(nrow)-hBkg.GetBinContent(nrow)-tmpSmear)**2/tmpSmear
            except ZeroDivisionError:
                print("Zero division, chi2_goodness set to inf")
                chi2_goodness = float('inf')
                break
        print("Channel %s for variable %s, chi2_goodness is %s"%(chan,varNames[varName],chi2_goodness))
        

    if not chi2_unf ==-1.:
        if isNom:
            for nrow in range(1,hresult.GetNbinsX()+1):
                for ncol in range(1,hresult.GetNbinsX()+1):
                    if nrow in exclude_bt_bin or ncol in exclude_bt_bin:
                        continue
                    try:
                        chi2_unf+= (hresult.GetBinContent(nrow)-hTrue.GetBinContent(nrow))*(hresult.GetBinContent(ncol)-hTrue.GetBinContent(ncol))*hCovInv(nrow-1,ncol-1) #different indexing
                    except ZeroDivisionError:
                        print("Zero division, chi2_unf set to inf")
                        chi2_unf = float('inf')
                        break
    if isNom:
        print("Channel %s for variable %s, chi2_smear is %s, chi2_unf is %s"%(chan,varNames[varName],chi2_smear,chi2_unf))
        #print("chi2_smear - chi2_unf is %s"%(chi2_smear-chi2_unf))


    #hOut.SetDirectory(0)
    #hResp.SetDirectory(0)
    #ROOT.SetOwnership(hCov,False)
    #print("hCov: ",hCov) 
    #print("where is the crash happening?")
    #return hCov.Clone(),hResp.Clone()
    #return hOut
    if withRespAndCov:
        if applyreg:
            return hOut3,hCov3.Clone(),hResp.Clone()
        else:
            return hOut,hCov.Clone(),hResp.Clone()
    
    #del hDataMinusBkg
    #print "DataMinusbkgIntegral: ",hDataMinusBkg, ", ",hDataMinusBkg.Integral()
    if applyreg:
        return hOut3
    else:
        return hOut

_generateUncertainties_count=0 #global counter for distinguishing histogram to suppress warning message in the printout

def _generateUncertainties(hDict,varName,norm): #hDict is hUnfolded dict
    #if called after rebin, the overflow bin is 0 already, and do we include underflow bin?
    #rebin was not called for histograms in hUnfolded, so use nbins+1 -> but there shouldn't be overflow bin content for unfolded hist.
    global _generateUncertainties_count
    _generateUncertainties_count+=1
    nominalArea = hDict[''].Integral(1,hDict[''].GetNbinsX()+1) #original codes start with 0th bin.
    hn = hDict[''].Clone()
    if norm:
        hn.Scale(1.0/nominalArea)
    hErr = {'Up':{},'Down':{}}

    #Note on indexing: 0 to 8 always for scale, the 9th index is pdf nominal for 2018, but the pdf nominal is not included in 2016  
    #so in 2016 the 9th index is the first pdf variation, and the other pdf indices shift by 1

    year = analysis[4:]
    if year=='2018' or year=='2017':
        ind_as1 = 110
        ind_as2 = 111
        ind_pdf1 = 10
    if year=='2016':
        ind_as1 = 109
        ind_as2 = 110
        ind_pdf1 = 9
    
    pdf_strs = ['PS_'+ str(x) for x in range(ind_pdf1,ind_as1)]
        


    pdflist={}
    scaleindlist = [str(i) for i in range(1,9)] # 1 to 8 always correspond to scales
    scalenamelist = ['PS_'+str(i) for i in [1,2,3,4,6,8]]
    if norm:
        scalehists = []
        for sys in scalenamelist:
            hsc = hDict[sys].Clone()
            hsc.Scale(1.0/(hsc.Integral(1,hsc.GetNbinsX()+1)))
            scalehists.append(hsc)
    else:
        scalehists = [hDict[sys].Clone() for sys in scalenamelist]

    if varName == "eta":
        histbins=array.array('d',[0.,1.0,2.0,3.0,4.0,5.0,6.0])
    else:
        histbins=array.array('d',_binning[varName])
    
    hUncScaleUp=ROOT.TH1D("hUncScaleUp"+str(_generateUncertainties_count),"QCD scale uncertainty up.",len(histbins)-1,histbins)
    hUncScaleDn=ROOT.TH1D("hUncScaleDn"+str(_generateUncertainties_count),"QCD scale uncertainty down.",len(histbins)-1,histbins)


    avghist = None
    firstadd = True
    
    #calculate avg hist for pdf variations MC replica case, used for 2016 replicas error
    for sys, h in hDict.iteritems():
        if not sys:
            continue
        he = h.Clone()
        if norm: #divide by nominal area for pdf variations (not QCD scale or alpha_s)
            he.Scale(1.0/(he.Integral(1,he.GetNbinsX()+1)))
            #he.Scale(nominalArea/(he.Integral(0,he.GetNbinsX()+1)))

        if sys in pdf_strs:  #not ('111' in sys or '110' in sys or sys.replace('PS_','') in scaleindlist):
            if firstadd:
                avghist = he.Clone()
                firstadd = False
            else:
                avghist.Add(he)
    
    avghist.Scale(0.01) #divided by 100

    #Output some diagnostics histogram in _generateUncertainties function
    global diagnostic_count
    diagnostic_count +=1
    if diagnostic_count==4 and args['diagnostic']: #4 corresponds to total case
        
        #hUncPDFsum=ROOT.TH1D("hUncPDFsum","PDF Uncertainty",len(histbins)-1,histbins)
        #tmpPDFSum = range(0,hn.GetNbinsX()+1)
        #for i in range(0,hn.GetNbinsX()+1):
        #    tmpPDFSum[i] = 0.

        myoutputFile=ROOT.TFile("diagnostics_%s.root"%varName,"RECREATE")
        myoutputFile.cd()
        hnomtmp = hn.Clone("Nominal")
        hnomtmp.Write()
        for sys, h in hDict.iteritems():
            if not sys:
                continue
            he = h.Clone(sys)

            if norm: #divide by nominal area for pdf variations (not QCD scale or alpha_s)
                if 'eEff' in sys or 'mEff' in sys or 'generator' in sys:
                    he.Scale(1.0/(he.Integral(1,he.GetNbinsX()+1)))
                    he.Write()
                    print(he)
                #if 'PS' in sys:
                #    if not ('111' in sys or '110' in sys or sys.replace('PS_','') in scaleindlist):
                #        he.Scale(1.0/(he.Integral(1,he.GetNbinsX()+1)))
                #        for i in range(1,hn.GetNbinsX()+1):
                #            tmpPDFSum[i] += (hn.GetBinContent(i)-he.GetBinContent(i))**2
                #        he.Write()
                #        print(he)

        #for i in range(1,hn.GetNbinsX()+1):
        #    hUncPDFsum.SetBinContent(i,math.sqrt(tmpPDFSum[i]))

        #hUncPDFsum.Write()
        myoutputFile.Close()
        print("File written")
        del sys
        import sys
        sys.exit()

#main systematics calculation
    for sys, h in hDict.iteritems():
        if not sys:
            continue
        he = h.Clone()

        if norm: #divide by nominal area for pdf variations (not QCD scale or alpha_s)
            if 'PS' in sys:
                if sys in pdf_strs:#not ('111' in sys or '110' in sys or sys.replace('PS_','') in scaleindlist):
                    he.Scale(1.0/(he.Integral(1,he.GetNbinsX()+1)))
            else:
                he.Scale(1.0/(he.Integral(1,he.GetNbinsX()+1)))
            #he.Scale(nominalArea/(he.Integral(0,he.GetNbinsX()+1)))

        #if sys.replace('PS_','') in scaleindlist:
        #    pdb.set_trace()
        #print "systematic: ",sys
        #if sys=="generator":
        #print "Unc: ",he.Integral()
        #print "nominalArea:",nominalArea
        #Subtract the nominal histogram from the SysUp or Down and pdf/scale histograms other than alpha_s variations

        if not 'PS' in sys:
            if not norm:
                he.Add(hDict[''],-1)
            else:
                he.Add(hn,-1)

        elif sys in pdf_strs:#not ('111' in sys or '110' in sys or sys.replace('PS_','') in scaleindlist): #last two correspond to alpha_s variation
            if not norm:
                if year=='2016':
                    he.Add(avghist,-1) #avghist is already normalized or not normalized according to norm
                else:
                    he.Add(hDict[''],-1)
            else:
                if year=='2016':
                    he.Add(avghist,-1)
                else:
                    he.Add(hn,-1)

        sysName = sys.replace('_Up','').replace('_Down','') #used for systematics other than the added scale/pdf and jet systs
        sysName2 = sys.replace('_up','').replace('_dn','')#should only work for "jes_up","jes_dn","jer_up","jer_dn"
        #if sys=="generator":
        #commented_print "Unc after nominal subtraction: ",he.Integral()
        if '_Up' in sys:
            hErr['Up'][sysName] = he
        elif '_Down' in sys:
            hErr['Down'][sysName] = he
        elif '_up' in sys:  #should only work for "jes_up","jes_dn","jer_up","jer_dn"
            hErr['Up'][sysName2] = he
        elif '_dn' in sys:
            hErr['Down'][sysName2] = he
        elif 'PS_' in sys: 
            if sys in pdf_strs: #not ('111' in sys or '110' in sys or sys.replace('PS_','') in scaleindlist):
                hErr['Up'][sys] = he
                hErr['Down'][sys] = he.Clone() #these pdf variations get into the square sum 
        else:
            hErr['Up'][sysName] = he
            he2 = he.Clone()
            hErr['Down'][sysName] = he2
    #pdb.set_trace()
    h_alphas_up = hDict['PS_'+str(ind_as2)].Clone()
    h_alphas_down = hDict['PS_'+str(ind_as1)].Clone()
    if norm:
        h_alphas_up.Scale(1.0/(h_alphas_up.Integral(1,h_alphas_up.GetNbinsX()+1)))
        h_alphas_down.Scale(1.0/(h_alphas_down.Integral(1,h_alphas_down.GetNbinsX()+1)))
    h_alphas = h_alphas_up.Clone()
    h_alphas.Add(h_alphas_down,-1)
    h_alphas.Scale(0.5)
    
    hErr['Up']['alpha_s'] = h_alphas
    hErr['Down']['alpha_s'] = h_alphas.Clone()

    for i in range(1,hUncScaleUp.GetNbinsX()+1):
        tmpcontent = [h.GetBinContent(i)-hn.GetBinContent(i) for h in scalehists]
        maxi=max(tmpcontent)
        mini=min(tmpcontent)
        hUncScaleUp.SetBinContent(i,maxi)
        hUncScaleDn.SetBinContent(i,mini)

    #if norm:
    #    hUncScaleUp.Add(hn,-1)
    #    hUncScaleDn.Add(hn,-1)
    #else:
    #    hUncScaleUp.Add(hDict[''],-1)
    #    hUncScaleDn.Add(hDict[''],-1)
    hErr['Up']['QCD_scales'] = hUncScaleUp
    hErr['Down']['QCD_scales'] = hUncScaleDn
    return hErr

def _sumUncertainties(errDict,varName):
    if varName == "eta":
        histbins=array.array('d',[0.,1.0,2.0,3.0,4.0,5.0,6.0])
    else:
        histbins=array.array('d',_binning[varName])
    #print "histbins: ",histbins
    hUncUp=ROOT.TH1D("hUncUp","Total Up Uncert.",len(histbins)-1,histbins)
    hUncDn=ROOT.TH1D("hUncDn","Total Dn Uncert.",len(histbins)-1,histbins)
    sysList = errDict['Up'].keys()
    #pdb.set_trace()
    print "sysList: ",sysList
    #print "hUncUp: ",hUncUp,"",hUncUp.Integral()
    #print "hUncDown: ",hUncDn,"",hUncDn.Integral()
    totUncUp=totUncDn=0.
    UncUpHistos= [errDict['Up'][sys] for sys in sysList]
    UncDnHistos= [errDict['Down'][sys] for sys in sysList]
    #commented_print "UncUpHistos: ",UncUpHistos
    #commented_print "UncDnHistos: ",UncDnHistos
    for i,sys in enumerate(sysList):
        print "systematic: ",sys
        print "UncUp: ",UncUpHistos[i].Integral()
        print "UncDn: ",UncDnHistos[i].Integral()
    LumiUp = errDict['Up']['generator'] #lumi ->generator??
    LumiDn = errDict['Down']['generator']
    print "GeneratorUp: ",LumiUp.Integral()
    print "GeneratorDn: ",LumiDn.Integral()
    for i in range(1,hUncUp.GetNbinsX()+1):
        for h1, h2 in zip(UncUpHistos,UncDnHistos):
            print "histUp: ",h1,"",h1.GetBinContent(i)
            print "histDn: ",h2,"",h2.GetBinContent(i)
            totUncUp += max(h1.GetBinContent(i),h2.GetBinContent(i))**2
            totUncDn += min(h1.GetBinContent(i),h2.GetBinContent(i))**2

        totUncUp = math.sqrt(totUncUp)
        totUncDn = math.sqrt(totUncDn)
        #print "totUncUp: ",totUncUp
        #print "totUncDn: ",totUncDn
        hUncUp.SetBinContent(i,totUncUp)
        hUncDn.SetBinContent(i,totUncDn)
    #commented_print("hUncUp: ",hUncUp,"",hUncUp.Integral()) 
    #commented_print("hUncDown: ",hUncDn,"",hUncDn.Integral())

    return hUncUp, hUncDn

_errfolder = "ErrorLog_%s"%(str(datetime.datetime.now()).split(" ")[0]) #replace(" ", '-').replace(".","-").replace(":","-")
if not os.path.isdir(_errfolder):
    os.mkdir(_errfolder)
_sumUncertainties_info_count=0 # counter used to distinguish histograms to suppress warning message
def _sumUncertainties_info(norm,errDict,varName,hUnf,chan=''): #same as above but used to printout info
    global _sumUncertainties_info_count
    _sumUncertainties_info_count +=1
    year = analysis[4:]
    systSum = {}
    
    
    ferrinfo=open(_errfolder+"/"+"ErrorInfo_%s%s.log"%(varName,chan),'w')
    ferrinfo.write("Var: %s \n"%varName)
    tmparea0= hUnf.Integral(1,hUnf.GetNbinsX()) #used for shorthand
    ferrinfo.write("Area: %s \n"%tmparea0)
    tmparea = tmparea0 #used to multiply the unc to cancel the division in plotUnfolded script to reduce code modification
    tmparea2 = tmparea0 # This is used to calculate the portion of error
    tmparea3 = tmparea0 # used to normalize stat unc for printout only
    totalbw = 0.
    for i in range(1,hUnf.GetNbinsX()+1):
        totalbw += hUnf.GetBinWidth(i)

    if not norm:
        tmparea = 1.0
        tmparea3 = 1.0
    if norm:
        tmparea2 = 1.0


    if varName == "eta":
        histbins=array.array('d',[0.,1.0,2.0,3.0,4.0,5.0,6.0])
    else:
        histbins=array.array('d',_binning[varName])
    #print "histbins: ",histbins
    hUncUp=ROOT.TH1D("hUncUp"+str(_sumUncertainties_info_count),"Total Up Uncert.",len(histbins)-1,histbins)
    hUncDn=ROOT.TH1D("hUncDn"+str(_sumUncertainties_info_count),"Total Dn Uncert.",len(histbins)-1,histbins)
    sysList = errDict['Up'].keys()

    extraSystSumList = ['pdf','stat','total']
    for key in extraSystSumList:
        systSum[key] = {}
        systSum[key]['Up'] = systSum[key]['Down'] = 0.
    
    for key in sysList:
        if not 'PS' in key:
            systSum[key] ={}
            systSum[key]['Up'] = 0.
            systSum[key]['Down'] = 0.

    #pdb.set_trace()
    print("Current channel:%s"%chan)
    print("Current variable:%s"%varName)
    print("Current sysList: ",sysList)
    #print "hUncUp: ",hUncUp,"",hUncUp.Integral()
    #print "hUncDown: ",hUncDn,"",hUncDn.Integral()
    totUncUp=totUncDn=0.
    UncUpHistos= [errDict['Up'][sys] for sys in sysList]
    UncDnHistos= [errDict['Down'][sys] for sys in sysList]
    #print "UncUpHistos: ",UncUpHistos
    #print "UncDnHistos: ",UncDnHistos
    print("Current Unf: ",[hUnf.GetBinContent(x) for x in range(1,hUnf.GetNbinsX()+1)])
    print("Current stat error: ",[hUnf.GetBinError(x) for x in range(1,hUnf.GetNbinsX()+1)])
    for i,sys in enumerate(sysList):
        if not 'PS' in sys:
            ferrinfo.write("Source Up %s:%s\n"%(sys,[UncUpHistos[i].GetBinContent(x) for x in range(1,UncUpHistos[i].GetNbinsX()+1)]))
            ferrinfo.write("Source Dn %s:%s\n"%(sys,[UncDnHistos[i].GetBinContent(x) for x in range(1,UncDnHistos[i].GetNbinsX()+1)]))
        print("Current systematic: ",sys)
        print("Current UncUp: ",[UncUpHistos[i].GetBinContent(x) for x in range(1,UncUpHistos[i].GetNbinsX()+1)])
        print("Current UncDn: ",[UncDnHistos[i].GetBinContent(x) for x in range(1,UncDnHistos[i].GetNbinsX()+1)])
    ferrinfo.write("Source Stat unc:%s\n"%([hUnf.GetBinError(x)/tmparea3 for x in range(1,hUnf.GetNbinsX()+1)]))
    LumiUp = errDict['Up']['generator'] #lumi ->generator??
    LumiDn = errDict['Down']['generator']
    #print "GeneratorUp: ",LumiUp.Integral()
    #print "GeneratorDn: ",LumiDn.Integral()
    
    pdfbincontent = []
    for i in range(1,hUncUp.GetNbinsX()+1):
        if norm:# divided by normalized height in that bin to get the portion of shape change, then weighted by bin width for averaging
            #extraNfac = 1.0/(hUnf.GetBinContent(i)/tmparea0)*hUnf.GetBinWidth(i)/totalbw
            extraNfac = 1.0 #currently don't use average weighted by bin width
        ferrinfo.write("Bin %s\n"%i)
        totUncUp=totUncDn=0. #should be reset each time but not done in original codes?
        PS_sum = 0.
        for j,(h1, h2) in enumerate(zip(UncUpHistos,UncDnHistos)):
            #if chan == 'mmmm' and sysList[j]=='lumi':
            #    pdb.set_trace()
            if not 'PS' in sysList[j]: #only 10 to 109 in the list for 17, 18 and 9 to 108 for 2016
                ferrinfo.write("Syst: %s \n"%sysList[j])
                ferrinfo.write("histUp: %s \n"%(h1.GetBinContent(i)))
                ferrinfo.write("histDn: %s \n"%(h2.GetBinContent(i)))
                #ferrinfo.write("Chosen max/min %s/%s \n"%(max(h1.GetBinContent(i),h2.GetBinContent(i)),min(h1.GetBinContent(i),h2.GetBinContent(i))))
                if h1.GetBinContent(i)*h2.GetBinContent(i)<=0: #opposite sign then determine up/down envelop by +up/-down
                    totUncUp += max(h1.GetBinContent(i),h2.GetBinContent(i))**2
                    totUncDn += min(h1.GetBinContent(i),h2.GetBinContent(i))**2
                    systSum[sysList[j]]['Up'] +=abs(max(h1.GetBinContent(i),h2.GetBinContent(i)))*extraNfac
                    systSum[sysList[j]]['Down'] += abs(min(h1.GetBinContent(i),h2.GetBinContent(i)))*extraNfac
                else: #same sign then follow original up/down order
                    totUncUp +=(h1.GetBinContent(i))**2
                    totUncDn +=(h2.GetBinContent(i))**2
                    systSum[sysList[j]]['Up'] +=abs(h1.GetBinContent(i))*extraNfac
                    systSum[sysList[j]]['Down'] += abs(h2.GetBinContent(i))*extraNfac
            else: #already satify in generateUncertainties#if not ('111' in sys or '110' in sys or sys.replace('PS_','') in scaleindlist):
                PS_sum += h1.GetBinContent(i)**2 
        
        if year=='2016':
            PS_sum = PS_sum/(100.-1.) # for MC replicas
        totUncUp += PS_sum
        totUncDn += PS_sum
        totUncUp = math.sqrt(totUncUp)
        totUncDn = math.sqrt(totUncDn)
        PS_sum = math.sqrt(PS_sum)
        pdfbincontent.append(PS_sum)
        systSum['pdf']['Up'] += PS_sum*extraNfac
        systSum['pdf']['Down'] += PS_sum*extraNfac
        ferrinfo.write("Syst: PDF (not alpha_s) Value:%s \n"%PS_sum)
        
        ferrinfo.write("totSysUncUp: %s \n"%totUncUp)
        ferrinfo.write("totSysUncDn: %s \n"%totUncDn)

        ferrinfo.write("Normalized Stat Unc: %s \n"%(hUnf.GetBinError(i)/tmparea3))
        systSum['stat']['Up'] += abs(hUnf.GetBinError(i))/tmparea3   #divide by area to show the portion in the normalized case for stat error 
        systSum['stat']['Down'] += abs(hUnf.GetBinError(i))/tmparea3 #abs just in case. Shouldn't need it.

        finalup = (totUncUp**2+(hUnf.GetBinError(i)/tmparea3)**2)**0.5
        finaldn = (totUncDn**2+(hUnf.GetBinError(i)/tmparea3)**2)**0.5
        systSum['total']['Up'] += finalup*extraNfac
        systSum['total']['Down'] += finaldn*extraNfac

        #normup = finalup/tmparea/hUnf.GetBinWidth(i)
        #normdn = finaldn/tmparea/hUnf.GetBinWidth(i)
        ferrinfo.write("Sys+stat Combined up: %s\n"%finalup)
        ferrinfo.write("Sys+stat Combined dn: %s\n"%finaldn)
        #ferrinfo.write("norm Combined up: %s\n"%(format(normup,".2e")))
        #ferrinfo.write("norm Combined dn: %s\n"%(format(normdn,".2e")))
        
        hUncUp.SetBinContent(i,totUncUp*tmparea)
        hUncDn.SetBinContent(i,totUncDn*tmparea)
    #commented_print("hUncUp: ",hUncUp,"",hUncUp.Integral()) 
    #commented_print("hUncDown: ",hUncDn,"",hUncDn.Integral())
    ferrinfo.write("Source pdf unc:%s\n"%(pdfbincontent))
    ferrinfo.write("======================\n")
    ferrinfo.write("Error Summary \n")
    sumportup = 0.
    sumportdn = 0.
    for key in systSum.keys():
        tmpup = systSum[key]['Up']
        tmpdn = systSum[key]['Down']
        portup = tmpup/tmparea2
        portdn = tmpdn/tmparea2
        #portup = tmpup/systSum['total']['Up']
        #portdn = tmpdn/systSum['total']['Down']
        if not key=='total':
            sumportup += portup
            sumportdn += portdn
        tup = round(tmpup,3)
        tdn = round(tmpdn,3)
        pup = round(portup,4)
        pdn = round(portdn,4)
        ferrinfo.write("%s: PortionUp %s PortionDn %s \n"%(key,pup,pdn))
    ferrinfo.write("Sum portion up and down: %s %s \n"%(sumportup, sumportdn))
    return hUncUp, hUncDn

def sumUnfoldDict(*hUnfoldeds):
    #pdb.set_trace()
    hUnfoldedtot = {}
    uncList = []
    for dic in hUnfoldeds: #one dict per channel
        uncList += dic.keys()
    keys = set(uncList)
    
    for key in keys:
        if varName == "eta":
            histbins=array.array('d',[0.,1.0,2.0,3.0,4.0,5.0,6.0])
        else:
            histbins=array.array('d',_binning[varName])
            #print "histTot histbins: ",histbins
        #histTot=ROOT.TH1D("histTotUnfolded_tmp","Tot hist",len(histbins)-1,histbins)
        #ROOT.SetOwnership(histTot,False)
        hUnfoldedtot[key] = ROOT.TH1D("histTotUnfolded_tmp","Tot hist "+key,len(histbins)-1,histbins)
        hUnfoldedtot[key].SetName(key+'tot_Sum') 
        
        print("Unc type before norm:%s"%key)
        
    
        for i,dic in enumerate(hUnfoldeds):
            try:
                hUnfoldedtot[key].Add(dic[key])
                print("Bin content before norm for channel %s"%i)
                print([dic[key].GetBinContent(j) for j in range(1,dic[key].GetNbinsX()+1)])
            except KeyError:
                hUnfoldedtot[key].Add(dic['']) #For eEff and mEff, add nominal if the variation doesn't exist in that channel
    
    return hUnfoldedtot


def _combineChannelUncertainties(*errDicts):
    hUncTot = {}
    uncList = []
    for errDict in errDicts:
        for sys in ['Up','Down']:
            uncList += errDict[sys].keys()
    uncList = set(uncList)
    #commented_print "uncList:",uncList
    for sys in ['Up','Down']:
        hUncTot[sys] = {}
        for unc in uncList:
            if varName == "eta":
                histbins=array.array('d',[0.,1.0,2.0,3.0,4.0,5.0,6.0])
            else:
                histbins=array.array('d',_binning[varName])
            #print "histTot histbins: ",histbins
            histTot=ROOT.TH1D("histTot","Tot Uncert.",len(histbins)-1,histbins)
            ROOT.SetOwnership(histTot,False)
            hUncTot[sys][unc] = histTot.Clone() #clone just in case in new codes
            for errDict in errDicts:
                try:
                    hUncTot[sys][unc].Add(errDict[sys][unc])
                except KeyError:
                    continue

    return hUncTot

def mkdir(plotDir):
    for outdir in [plotDir]:
        try:
            os.makedirs(os.path.expanduser(outdir))
        except OSError as e:
            print e
            pass

plotDir=args['plotDir']
nIterations=args['nIter']



#pdb.set_trace()
#Dictionary where signal samples are keys with cross-section*kfactors as values
sigSampleDic=ConfigureJobs.getListOfFilesWithXSec(ConfigureJobs.getListOfEWK())
sigSampleList=[str(i) for i in sigSampleDic.keys()]
print "sigSamples: ",sigSampleList

AltsigSampleDic=ConfigureJobs.getListOfFilesWithXSec([myaltname,])
AltsigSampleList=[str(i) for i in AltsigSampleDic.keys()]
print "AltsigSamples: ",AltsigSampleList

#Combine sigSamples
TotSigSampleList = list(set(sigSampleList) | set(AltsigSampleList))
sigSampleDic.update(AltsigSampleDic)
#Replace fOut with fUse once you have run all the data samples and the backgrounds - right now unfolded data looking really big- subtract backgrounds
if args['test']:
    sigSamplesPath={}
    if analysis=="ZZ4l2016":
        fUse = ROOT.TFile("SystGenFiles/Fullsys_fullrange16.root","update") 
        #fUse = ROOT.TFile("SystGenFiles/For_unfolding_Hists02Feb2022_ZZ4l2016_Moriond_fullSyst.root","update") 
        #fUse = ROOT.TFile("SystGenFiles/For_unfolding_Hists15Dec2021_ZZ4l2016_Moriond_fullSyst.root","update")
        #fUse = ROOT.TFile("SystGenFiles/Hists25Jun2020-ZZ4l2016_Moriond.root","update") #Note one more recent file is not in the commented list 
        #fUse = ROOT.TFile("SystGenFiles/Hists31Mar2020-ZZ4l2016_Moriond.root","update")
    elif analysis=="ZZ4l2017":
        fUse = ROOT.TFile("SystGenFiles/Fullsys_fullrange17_full.root","update")
        #fUse = ROOT.TFile("SystGenFiles/For_unfolding_Hists02Feb2022_ZZ4l2017_Moriond_fullSyst.root","update")
        #fUse = ROOT.TFile("SystGenFiles/For_unfolding_Hists20Dec2021_ZZ4l2017_Moriond_fullSyst.root","update")
        #fUse = ROOT.TFile("SystGenFiles/For_unfolding_Hists17May2021_ZZ4l2017_Moriond_fullSyst.root","update")
        #fUse = ROOT.TFile("SystGenFiles/Hists07Jun2020-ZZ4l2017_Moriond.root","update")
    elif analysis=="ZZ4l2018":
        #fUse = ROOT.TFile("SystGenFiles/HEM1516.root","update") 
        fUse = ROOT.TFile("SystGenFiles/Fullsys_fullrange18_full.root")#,"update")
        #fUse = ROOT.TFile("SystGenFiles/For_unfolding_Hists02Feb2022_ZZ4l2018_MVA_fullSyst.root","update")
        #fUse = ROOT.TFile("SystGenFiles/Syst_qqZZNewMCadded_Hists18Oct2021-ZZ4l2018_MVA.root","update")
        #fUse = ROOT.TFile("SystGenFiles/Syst_qqZZNewMCadded_Hists30Aug2021-ZZ4l2018_MVA.root","update") #Most recent before jet syst and pdf/scale syst
        #fUse = ROOT.TFile("SystGenFiles/Syst_qqZZNewMCadded_Hists14Jul2021-ZZ4l2018_MVA.root","update") #Recent results after adding new MC qqZZ first Round
        #fUse = ROOT.TFile("SystGenFiles/For_unfolding_Hists17May2021_ZZ4l2018_MVA_fullSyst.root","update") #Recent results before adding new MC qqZZ
        #fUse = ROOT.TFile("SystGenFiles/Hists08Jun2020-ZZ4l2018_MVA.root","update")
    fOut=fUse
    #pdb.set_trace()
    for dataset in TotSigSampleList:
        file_path = ConfigureJobs.getInputFilesPath(dataset,selection, analysis)
        #print "file_path:",file_path
        sigSamplesPath[dataset]=file_path

#Sum all data and return a TList of all histograms that are booked. And an empty datSumW dictionary as there are no sumWeights
alldata,dataSumW = HistTools.makeCompositeHists(fOut,"AllData", 
    ConfigureJobs.getListOfFilesWithXSec([args['analysis']+"data"],manager_path), args['lumi'],
    underflow=False, overflow=False)


#all ewkmc/this is also allSignal histos, scaled properly, kind of a repeat of above but with ggZZ added
ewkmc,ewkSumW = HistTools.makeCompositeHists(fOut,"AllEWK", ConfigureJobs.getListOfFilesWithXSec(
    ConfigureJobs.getListOfEWK(), manager_path), args['lumi'],
    underflow=False, overflow=False)

ewkmc_qqZZonly,ewkSumW_qqZZonly = HistTools.makeCompositeHists(fOut,"AllEWKqqZZonly", ConfigureJobs.getListOfFilesWithXSec(
    ConfigureJobs.getListOfEWK()[:1], manager_path), args['lumi'],
    underflow=False, overflow=False)

ewkmc_ggZZonly,ewkSumW_ggZZonly = HistTools.makeCompositeHists(fOut,"AllEWKggZZonly", ConfigureJobs.getListOfFilesWithXSec(
    ConfigureJobs.getListOfEWK()[1:], manager_path), args['lumi'],
    underflow=False, overflow=False)


ewkmc_ggZZup,ewkSumW_ggZZup = HistTools.makeCompositeHists_scaling(fOut,"AllEWKggZZup", ConfigureJobs.getListOfFilesWithXSec(
    ConfigureJobs.getListOfEWK(), manager_path), args['lumi'],
    underflow=False, overflow=False,scale_fac=1.+0.18)

ewkmc_ggZZdn,ewkSumW_ggZZdn = HistTools.makeCompositeHists_scaling(fOut,"AllEWKggZZdn", ConfigureJobs.getListOfFilesWithXSec(
    ConfigureJobs.getListOfEWK(), manager_path), args['lumi'],
    underflow=False, overflow=False,scale_fac=1.-0.14)

altSigmc,altSigSumW = HistTools.makeCompositeHists(fOut,"AltSig", ConfigureJobs.getListOfFilesWithXSec(
    ConfigureJobs.getListOfaltSig(), manager_path), args['lumi'],
    underflow=False, overflow=False)

#Update ewkSumW dictionary with sumWeights value of zz4l-amcatnlo from altSigSumW, the common keys should not be duplicated
ewkSumW.update(altSigSumW)

#all mcbkg that needs to be subtracted
allVVVmc,VVVSumW = HistTools.makeCompositeHists(fOut,"AllVVV", ConfigureJobs.getListOfFilesWithXSec(
    ConfigureJobs.getListOfVVV(), manager_path), args['lumi'],
    underflow=False, overflow=False)

#This is the non-prompt background
ewkcorr = HistTools.getDifferenceDirect(fOut, "DataEWKCorrected", alldata, ewkmc)

print "Signals: ",ewkSumW
#print the sum for a sample (zz4l-powheg)
zzSumWeights = ewkSumW[mynominalName]  
#print "sumW (zz4l-powheg): ",zzSumWeights

#getHistInDic function also takes care of adding the histograms in eemm+mmee, hence the input here is channels=[eeee,eemm,mmmm]
#dataHists dictionary
hDataDic=OutputTools.getHistsInDic(alldata,varList,channels)

hSigDic=OutputTools.getHistsInDic(ewkmc,varList,channels)
hSigDic_ggZZonly=OutputTools.getHistsInDic(ewkmc_ggZZonly,varList,channels)
hSigDic_ggZZup=OutputTools.getHistsInDic(ewkmc_ggZZup,varList,channels)
hSigDic_ggZZdn=OutputTools.getHistsInDic(ewkmc_ggZZdn,varList,channels)

#Alt signals containing zzl4-amcatnlo instead of zz4l-powheg
hAltSigDic=OutputTools.getHistsInDic(altSigmc,varList,channels)

#TrueHists dictionary
#hTrueDic=OutputTools.getHistsInDic(allzzPowheg,["Gen"+s for s in varList],channels)
#pdb.set_trace()
hTrueDic=OutputTools.getHistsInDic(ewkmc,["Gen"+s for s in varList],channels)
hTrueDic_ggZZonly=OutputTools.getHistsInDic(ewkmc_ggZZonly,["Gen"+s for s in varList],channels)
hTrueDic_ggZZup=OutputTools.getHistsInDic(ewkmc_ggZZup,["Gen"+s for s in varList],channels)
hTrueDic_ggZZdn=OutputTools.getHistsInDic(ewkmc_ggZZdn,["Gen"+s for s in varList],channels)

#Alt signals containing zzl4-amcatnlo instead of zz4l-powheg
hAltTrueDic=OutputTools.getHistsInDic(altSigmc,["Gen"+s for s in varList],channels)

#Non-prompt background dictionary
hbkgDic=OutputTools.getHistsInDic(ewkcorr,[s+"_Fakes" for s in varList],channels)
#strange python debug
print "channels: ",channels
if "mmee" in channels:
    channels.remove("mmee")
#VVV background dictionary
hbkgMCDic=OutputTools.getHistsInDic(allVVVmc,varList,channels)
print "hbkgMCDic: ",hbkgMCDic

runVariables=[]
runVariables.append(args['variable'])
print "runVariables: ",runVariables

##Systematic histos
systList=[]
gensystList=[]
for chan in channels:
    for s in runVariables:
        systList.append(varNames[s]+"_lheWeights")
        gensystList.append("Gen"+varNames[s]+"_lheWeights")
        systList.append(varNames[s]+"_jetsysts")
    for sys in ["Up","Down"]: 
        for s in runVariables:
            systList.append(varNames[s]+"_CMS_pileup"+sys)
            for lep in set(chan):         
                systList.append(varNames[s]+"_CMS_eff_"+lep+sys)

print systList
if not args['noSyst']:  #systList has repeated variables, but shouldn't matter as it will just reassigin same value in the dictionary
    hSigSystDic=OutputTools.getHistsInDic(ewkmc,systList,channels)
    hSigSystDic_qqZZonly=OutputTools.getHistsInDic(ewkmc_qqZZonly,systList,channels)
    hTrueSystDic_qqZZonly=OutputTools.getHistsInDic(ewkmc_qqZZonly,gensystList,channels)
    hbkgMCSystDic=OutputTools.getHistsInDic(allVVVmc,systList,channels)
else:
    hSigSystDic = None #Since they are still put into function arguments
    hbkgMCSystDic = None

OutputDirs={}
diagnostic_count=0
SFhistos,PUhistos = generateAnalysisInputs()
#So even when selector runs, run the unfolding procedure only on the variables provided.

#Make differential cross sections normalizing to unity area.')
norm = not args['noNorm']
savehists=[]
for varName in runVariables:
    print "varName:", varNames[varName]
    # save unfolded distributions by channel, then systematic
    hDataDict = {}
    hMCSigDict = {}
    hBkgTotDict={}
    hUnfolded = {"eeee":{},"eemm":{},"mmmm":{}}
    hTrue = {"eeee":{},"eemm":{},"mmmm":{}}
    hTrueAlt = {"eeee":{},"eemm":{},"mmmm":{}}
    hErr = {}
    hErrTrue = {}
    for chan in channels:
        print "channel: ",chan
        print "hUnfolded: ",hUnfolded
        print "hTrue: ",hTrue
        OutputDir=plotDir+"/"+chan+"/plots"
        if chan not in OutputDirs:
            OutputDirs[chan]=OutputDir
        if not os.path.exists(OutputDir):
            mkdir(OutputDir)
            OutputDirs[chan]=OutputDir
         
        responseMakers,altResponseMakers = generateResponseClass(varName, chan,sigSampleDic,sigSamplesPath,ewkSumW,PUhistos,SFhistos)
        print "hUnfolded in main: ", hUnfolded
        print "hTrue in main: ", hTrue
        print "hTrue in main: ", hTrueAlt
        hUnfolded[chan], hTrue[chan],hTrueAlt[chan],hDataDict[chan],hMCSigDict[chan],hBkgTotDict[chan] = unfold(varName,chan,responseMakers,altResponseMakers,hSigDic,hAltSigDic,hSigSystDic,hTrueDic,hTrueSystDic_qqZZonly,hAltTrueDic,hDataDic,hbkgDic,hbkgMCDic,hbkgMCSystDic,nIterations,OutputDir)
        print("returning unfolded? ",hUnfolded[chan])
        #print("returning truth? ",hTrue[chan])

        if not args['noSyst']: 
            hErr[chan]= _generateUncertainties(hUnfolded[chan],varName,norm)
            print "hErr[",chan,"]: ",hErr[chan]
            (hUncUp, hUncDn) = _sumUncertainties_info(norm,hErr[chan],varName,hUnfolded[chan][''],chan)
        #hErrTrue[chan] = _generateUncertainties(hTrue[chan],norm)
        #(hTrueUncUp, hTrueUncDn) = _sumUncertainties(hErrTrue[chan],varName)
            hDataSave = hDataDict[chan].Clone()
            dataName = chan+"_"+varName+"_data"
            hDataSave.SetName(dataName)
            savehists.append(hDataSave)
            #save histograms instead of plotting them
            hUnf = hUnfolded[chan][''].Clone()
            unfName = chan+"_"+varName+"_unf"
            hUnf.SetName(unfName)
            savehists.append(hUnf)
            #hTrue
            hTheo = hTrue[chan][''].Clone()
            truName = chan+"_"+varName+"_true"
            hTheo.SetName(truName)
            savehists.append(hTheo)
            #hTrueAlt
            hTheoAlt = hTrueAlt[chan][''].Clone()
            AltTruName = chan+"_"+varName+"_trueAlt"
            hTheoAlt.SetName(AltTruName)
            savehists.append(hTheoAlt)
            #hUncUp
            UncUp = hUncUp.Clone()
            UncUpName = chan+"_"+varName+"_totUncUp"
            UncUp.SetName(UncUpName)
            savehists.append(UncUp)
            #hUncDn
            UncDn = hUncDn.Clone()
            UncDnName = chan+"_"+varName+"_totUncDown"
            UncDn.SetName(UncDnName)
            savehists.append(UncDn)
    print "savehists: ",savehists
    if args['makeTotals']:
        if not args['noSyst']:
            
            hUnfoldedtotDic = sumUnfoldDict(*[hUnfolded['eeee'],hUnfolded['eemm'],hUnfolded['mmmm']])
            hErr_tot= _generateUncertainties(hUnfoldedtotDic,varName,norm)
            hErrTot = hErr_tot#_combineChannelUncertainties(*hErr.values())
            hTotUncUp, hTotUncDn = _sumUncertainties_info(norm,hErrTot,varName,hUnfoldedtotDic[''].Clone()) #hTot.Clone())

        if "eeee" in channels:
            hTotData = hDataDict["eeee"]
            hTotBkg = hBkgTotDict["eeee"]
            hTotMCSig = hMCSigDict["eeee"]
            hTot = hUnfolded["eeee"]['']
            hTrueTot = hTrue["eeee"]['']
            hTrueAltTot = hTrueAlt["eeee"]['']
            #channels.remove("eeee")
        print("channels before adding histos: ",channels)
        for c in ["eemm","mmmm"]:
            hTotData.Add(hDataDict[c])
            hTotBkg.Add(hBkgTotDict[c])
            hTotMCSig.Add(hMCSigDict[c])
            hTot.Add(hUnfolded[c][''])
            hTrueTot.Add(hTrue[c][''])
            hTrueAltTot.Add(hTrueAlt[c][''])
        print "hErr.values(): ",hErr.values()
        #Saving Total histograms
        hTotalData = hTotData.Clone()
        TotDatName = "tot_"+varName+"_data"
        hTotalData.SetName(TotDatName)
        savehists.append(hTotalData)

        hTotalBkg = hTotBkg.Clone()
        TotBkgName = "tot_"+varName+"_bkg"
        hTotalBkg.SetName(TotBkgName)
        savehists.append(hTotalBkg)

        hTotalMCSig = hTotMCSig.Clone()
        TotSigMCName = "tot_"+varName+"_SigMC"
        hTotalMCSig.SetName(TotSigMCName)
        savehists.append(hTotalMCSig)

        hTotUnf = hTot.Clone()
        TotName = "tot_"+varName+"_unf"
        hTotUnf.SetName(TotName)
        savehists.append(hTotUnf)
        #hTrue
        hTotTrue = hTrueTot.Clone()
        TotTruName = "tot_"+varName+"_true"
        hTotTrue.SetName(TotTruName)
        savehists.append(hTotTrue)
        #hTrueAlt
        hTotTrueAlt = hTrueAltTot.Clone()
        TotAltTruName = "tot_"+varName+"_trueAlt"
        hTotTrueAlt.SetName(TotAltTruName)
        savehists.append(hTotTrueAlt)
        #hUncUp
        if not args['noSyst']:
            TotUncUp = hTotUncUp.Clone()
            TotUncUpName = "tot_"+varName+"_totUncUp"
            TotUncUp.SetName(TotUncUpName)
            savehists.append(TotUncUp)
            #hUncDn
            TotUncDn = hTotUncDn.Clone()
            TotUncDnName = "tot_"+varName+"_totUncDown"
            TotUncDn.SetName(TotUncDnName)
            savehists.append(TotUncDn)

if args['plotResponse']:       
    for cat in ["eeee","eemm","mmmm"]:  
    #This is where we put all the response plots in html format for quick access/debugging
        makeSimpleHtml.writeHTML(os.path.expanduser(OutputDirs[cat].replace("/plots", "")), "2D ResponseMatrices (from MC)")

today = datetime.date.today().strftime("%d%b%Y")
tmpFileName = "UnfHistsFinal-NewMCRedo-26Jul2021-%s.root" % (analysis) 
#tmpFileName = "UnfHistsFinal-%s-%s.root" % (today, analysis) 
#tmpFileName = "UnfHistsFinal-18Apr2020-%s.root" % (analysis) 
fHistOut = ROOT.TFile.Open(tmpFileName, "update")
#fOut = ROOT.TFile.Open("/".join([outputFolder, outputFile]), "RECREATE")
fHistOut.cd()
for newhists in savehists:
    newhists.Write()
fHistOut.Close()

