#!/usr/bin/env python
import ROOT
from python import SelectorTools
from python import UserInput
from python import OutputTools
from python import ConfigureJobs
from python import HistTools
import os
import logging
import sys
import datetime
import subprocess,pdb

#logging.basicConfig(level=logging.DEBUG)

def getComLineArgs():
    parser = UserInput.getDefaultParser()
    parser.add_argument("--lumi", "-l", type=float,
        default=35.87, help="luminosity value (in fb-1)")
    parser.add_argument("--output_file", "-o", type=str,
        default="test.root", help="Output file name")
    parser.add_argument("--muIDType", type=str,
        default="", help="Moriond Muon ID or MVA ID")
    parser.add_argument("--test", action='store_true',
        help="Run test job (no background estimate)")
    parser.add_argument("--uwvv", action='store_true',
        help="Use UWVV format ntuples in stead of NanoAOD")
    parser.add_argument("--with_background", action='store_true',
        help="Don't run background selector")
    parser.add_argument("--with_Gen", action='store_true',
        help="Don't run ZZGen selector")
    parser.add_argument("--noHistConfig", action='store_true',
        help="Don't rely on config file to specify hist info")
    parser.add_argument("-j", "--numCores", type=int, default=1,
        help="Number of cores to use (parallelize by dataset)")
    parser.add_argument("--input_tier", type=str,
        default="", help="Selection stage of input files")
    parser.add_argument("--year", type=str,
        default="default", help="Year of Analysis")
    parser.add_argument("--scalefactors_file", "-sf", type=str,
        default="", help="ScaleFactors file name")
    parser.add_argument("-c", "--channels", 
                        type=lambda x : [i.strip() for i in x.split(',')],
                        default=["eee","eem","emm","mmm"], help="List of channels"
                        "separated by commas. NOTE: set to Inclusive for NanoAOD")
    parser.add_argument("-b", "--hist_names", 
                        type=lambda x : [i.strip() for i in x.split(',')],
                        default=["all"], help="List of histograms, "
                        "as defined in ZZ4lRun2DatasetManager, separated "
                        "by commas")
    return vars(parser.parse_args())

def makeHistFile(args):
    ROOT.gROOT.SetBatch(True)

    manager_path = ConfigureJobs.getManagerPath()
    if manager_path not in sys.path:
        sys.path.insert(0, "/".join([manager_path, 
            "ZZ4lRun2DatasetManager", "Utilities/python"]))

    today = datetime.date.today().strftime("%d%b%Y")
    
    if args['test']:
        tmpFileName = "Hists%s-%s.root" % (today, args['output_file']) 
    else:
        tmpFileName = "Hists%s-%s.root" % (today, args['output_file']) if args['selection'] == "SignalSync" \
            else "Hists%s-%s_%s.root" % (today, args['analysis'],args['muIDType'])
    toCombine = args['with_background'] or args['with_Gen']
    fOut = ROOT.TFile(tmpFileName if not toCombine else tmpFileName.replace(".root", "sel.root"), "recreate")
    combinedNames = [fOut.GetName()]

    addScaleFacs = False
    if args['analysis'] == "WZxsec2016" or args['analysis'] == 'Zstudy_2016' or args['scalefactors_file']:
        addScaleFacs = True
    fjetPUSF = ROOT.TFile("data/jetSF/scalefactorsPUID_81Xtraining.root")
    fjetPUeff= ROOT.TFile("data/jetSF/effcyPUID_81Xtraining.root")
    sf_inputs = [ROOT.TParameter(bool)("applyScaleFacs", False)]
    fr_inputs = []
    if addScaleFacs:
        fScales = ROOT.TFile(args['scalefactors_file'])
        if "ZZ4l" in args['analysis']:
            mZZTightFakeRate = fScales.Get("mZZTightFakeRate")
            eZZTightFakeRate = fScales.Get("eZZTightFakeRate")
            if mZZTightFakeRate:
                mZZTightFakeRate.SetName("fakeRate_allMu")
            if eZZTightFakeRate:
                eZZTightFakeRate.SetName("fakeRate_allE")
            if "2018" in args['scalefactors_file']:
                muonRunSF= fScales.Get('muonRun18SF')
                muonRunSF.SetName("muonRunSF")
                electronLowRecoSF = fScales.Get('electronLowReco18SF')
                electronLowRecoSF.SetName("electronLowRecoSF")
                electronRecoSF = fScales.Get('electronReco18SF')
                electronRecoSF.SetName("electronRecoSF")
                electronRunSF = fScales.Get('electronRun18SF')
                electronRunSF.SetName("electronRunSF")
                electronRunGapSF = fScales.Get('electronRun18GapSF')
                electronRunGapSF.SetName("electronRunGapSF")
                jetPUSF = fjetPUSF.Get("h2_eff_sf2018_T")
                jetPUSF.SetName("jetPUSF")
                jetPUeff = fjetPUeff.Get("h2_eff_mc2018_T")
                jetPUeff.SetName("jetPUeff")
            elif "2016" in args['scalefactors_file']:
                muonRunSF= fScales.Get('muonRun16SF')
                muonRunSF.SetName("muonRunSF")
                electronLowRecoSF = fScales.Get('electronLowReco16SF')
                electronLowRecoSF.SetName("electronLowRecoSF")
                electronRecoSF = fScales.Get('electronReco16SF')
                electronRecoSF.SetName("electronRecoSF")
                electronRunSF = fScales.Get('electronRun16SF')
                electronRunSF.SetName("electronRunSF")
                electronRunGapSF = fScales.Get('electronRun16GapSF')
                electronRunGapSF.SetName("electronRunGapSF")
                jetPUSF = fjetPUSF.Get("h2_eff_sf2016_T")
                jetPUSF.SetName("jetPUSF")
                jetPUeff = fjetPUeff.Get("h2_eff_mc2016_T")
                jetPUeff.SetName("jetPUeff")
            elif "2017" in args['scalefactors_file']:
                muonRunSF= fScales.Get('muonRun17SF')
                muonRunSF.SetName("muonRunSF")
                electronLowRecoSF = fScales.Get('electronLowReco17SF')
                electronLowRecoSF.SetName("electronLowRecoSF")
                electronRecoSF = fScales.Get('electronReco17SF')
                electronRecoSF.SetName("electronRecoSF")
                electronRunSF = fScales.Get('electronRun17SF')
                electronRunSF.SetName("electronRunSF")
                electronRunGapSF = fScales.Get('electronRun17GapSF')
                electronRunGapSF.SetName("electronRunGapSF")
                jetPUSF = fjetPUSF.Get("h2_eff_sf2017_T")
                jetPUSF.SetName("jetPUSF")
                jetPUeff = fjetPUeff.Get("h2_eff_mc2017_T")
                jetPUeff.SetName("jetPUeff")
            else: 
                print "what scale factors you want?"
                sys.exit()
            pileupSF = fScales.Get('pileupSF')
            jetPUSF.SetDirectory(0)
            jetPUeff.SetDirectory(0)
            ROOT.SetOwnership(jetPUSF,False)
            ROOT.SetOwnership(jetPUeff,False)
            fr_inputs = [eZZTightFakeRate, mZZTightFakeRate,]
            sf_inputs = [electronLowRecoSF,electronRecoSF,electronRunSF, electronRunGapSF,muonRunSF,pileupSF,jetPUSF,jetPUeff]
        else:
            fScales = ROOT.TFile('data/scaleFactors.root')
            mCBTightFakeRate = fScales.Get("mCBTightFakeRate")
            eCBTightFakeRate = fScales.Get("eCBTightFakeRate")
            useSvenjasFRs = False
            useJakobsFRs = False
            if useSvenjasFRs:
                mCBTightFakeRate = fScales.Get("mCBTightFakeRate_Svenja")
                eCBTightFakeRate = fScales.Get("eCBTightFakeRate_Svenja")
            elif useJakobsFRs:
                mCBTightFakeRate = fScales.Get("mCBTightFakeRate_Jakob")
                eCBTightFakeRate = fScales.Get("eCBTightFakeRate_Jakob")
            # For medium muons
            #mCBMedFakeRate.SetName("fakeRate_allMu")
            if mCBTightFakeRate:
                mCBTightFakeRate.SetName("fakeRate_allMu")
            if eCBTightFakeRate:
                eCBTightFakeRate.SetName("fakeRate_allE")

            muonIsoSF = fScales.Get('muonIsoSF')
            muonIdSF = fScales.Get('muonTightIdSF')
            electronTightIdSF = fScales.Get('electronTightIdSF')
            electronGsfSF = fScales.Get('electronGsfSF')
            pileupSF = fScales.Get('pileupSF')

            #fPrefireEfficiency = ROOT.TFile('data/Map_Jet_L1FinOReff_bxm1_looseJet_JetHT_Run2016B-H.root')
            fPrefireEfficiency = ROOT.TFile('data/Map_Jet_L1FinOReff_bxm1_looseJet_SingleMuon_Run2016B-H.root')
            prefireEff = fPrefireEfficiency.Get('prefireEfficiencyMap')

            fr_inputs = [eCBTightFakeRate, mCBTightFakeRate,]
            sf_inputs = [electronTightIdSF, electronGsfSF, muonIsoSF, muonIdSF, pileupSF, prefireEff]

        sf_inputs.append(ROOT.TParameter(bool)("applyScaleFacs", True))

    if args['input_tier'] == '':
        args['input_tier'] = args['selection']
    selection = args['selection'].split("_")[0]

    if selection == "Inclusive2Jet":
        selection = "Wselection"
        print "Info: Using Wselection for hist defintions"
    analysis = "/".join([args['analysis'], selection])
    hists, hist_inputs = UserInput.getHistInfo(analysis, args['hist_names'], args['noHistConfig'])

    selector = SelectorTools.SelectorDriver(args['analysis'], args['selection'], args['input_tier'], args['year'])
    selector.setOutputfile(fOut.GetName())
    selector.setInputs(sf_inputs+hist_inputs)

    selector.setNtupeType("UWVV" if args['uwvv'] else "NanoAOD")
    if args['uwvv']:
        logging.debug("Processing channels " % args['channels'])
        selector.setChannels(args['channels'])
    selector.setNumCores(args['numCores'])

    if args['filenames']:
        selector.setDatasets(args['filenames'])
    else:
        selector.setFileList(*args['inputs_from_file'])
    mc = selector.applySelector()

    print("Pause here")
    #sys.exit()

    if args['with_background']:
        selector.isBackground()
        selector.setInputs(sf_inputs+hist_inputs+fr_inputs)
        output_name = tmpFileName.replace(".root", "bkgd.root")
        selector.setOutputfile(output_name)
        bkgd = selector.applySelector()
        combinedNames.append(output_name)
    #pdb.set_trace()
    if args['with_Gen']:
        selector.isGen()
        selector.setChannels([c+"Gen" for c in args['channels']])
        # Make sure to remove data from the dataset lists
        selector.setInputs(hist_inputs)
        output_name = tmpFileName.replace(".root", "gen.root")
        selector.setOutputfile(output_name)
        combinedNames.append(output_name)
        if args['filenames']:
            #selector.setDatasets(args['filenames'])
            selector.setDatasets(ConfigureJobs.getListOfGenFilenames(args['analysis']))
        else:
            selector.setFileList(*args['inputs_from_file'])
        gen = selector.applySelector()
        selector.setChannels(args['channels'])
        selector.outputFile().Close()

    if len(combinedNames) > 1:
        rval = subprocess.call(["hadd", "-f", tmpFileName] + combinedNames)
        if rval == 0:
            map(os.remove, combinedNames)

    if args['test']:
        fOut.Close()
        sys.exit(0)

    fOut.Close()
    
    #! Can exit here if don't want to do final nonprompt background calc. (e.g. pure MCs)
    #sys.exit()
    
    fOut = ROOT.TFile.Open(tmpFileName, "update")
    #pdb.set_trace()
    alldata = HistTools.makeCompositeHists(fOut,"AllData", 
        ConfigureJobs.getListOfFilesWithXSec([args['analysis']+"data"], manager_path), args['lumi'],
        underflow=False, overflow=False)
    OutputTools.writeOutputListItem(alldata, fOut)
    alldata.Delete()

    if "ZZ4l" not in args['analysis']:
        nonpromptmc = HistTools.makeCompositeHists(fOut, "NonpromptMC", ConfigureJobs.getListOfFilesWithXSec( 
            ConfigureJobs.getListOfNonpromptFilenames(), manager_path), args['lumi'],
            underflow=False, overflow=False)
        nonpromptmc.Delete()

        OutputTools.writeOutputListItem(nonpromptmc, fOut)

    ewkmc = HistTools.makeCompositeHists(fOut,"AllEWK", ConfigureJobs.getListOfFilesWithXSec(
        ConfigureJobs.getListOfEWKFilenames(args['analysis']), manager_path), args['lumi'],
        underflow=False, overflow=False)
    OutputTools.writeOutputListItem(ewkmc, fOut)
    ewkmc.Delete()

    ewkcorr = HistTools.getDifference(fOut, "DataEWKCorrected", "AllData", "AllEWK")
    OutputTools.writeOutputListItem(ewkcorr, fOut)
    ewkcorr.Delete()

def main():
    args = getComLineArgs()
    makeHistFile(args)
    exit(0)

if __name__ == "__main__":
    main()
