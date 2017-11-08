#!/usr/bin/env python
import ROOT
import argparse
import os
import sys
from Utilities.python import ApplySelection
from Utilities.python.prettytable import PrettyTable

def getComLineArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--selections", type=str,
                        required=True, help="Name of selections to apply, "
                        "separated by commas. They must be"
                        " mapped to a cuts json via Cuts/definitions.json")
    parser.add_argument("-t", "--trigger", type=str, default="",
                        choices=["DoubleEG", "DoubleMuon", "MuonEG", 
                            "SingleMuon", "SingleElectron", "MonteCarlo", ""],
                        help="Name of trigger to select in data")
    parser.add_argument("-f", "--filelist", type=str,
                        required=True, help="List of input file names "
                        "to be processed (separated by commas)")
    parser.add_argument("-a", "--analysis", type=str,
                        required=True, help="Analysis name, used"
                        " in selection the cut json")
    parser.add_argument("-o", "--output_file_name", type=str,
                        required=True, help="Name of output file")
    parser.add_argument("-d", "--no_deduplicate", action='store_true',
                        help="Don't remove duplicated events from ntuple")
    return vars(parser.parse_args())
def writeNtupleToFile(output_file, tree, state, cut_string, deduplicate):
    state_dir = output_file.mkdir(state)
    state_dir.cd()
    save_tree = tree.CopyTree(cut_string) if not deduplicate else \
        getDeduplicatedTree(tree, state, cut_string).CopyTree("")
    save_tree.Write()
    # Remove AutoSaved trees
    output_file.Purge()
    ROOT.gROOT.cd()
    return save_tree.GetEntries()
def getDeduplicatedTree(tree, state, cut_string):
    selector = ROOT.disambiguateFinalStates()
    if state.count('e') > 2:
        l1_l2_cand_mass = "e1_e2_Mass"
        l1_l2_cand_pt = "e1_e2_Pt"
        l3_l4_cand_mass = "e3_e4_Mass"
        l3_l4_cand_pt = "e3_e4_Pt"
    elif state.count('m') > 2:
        l1_l2_cand_mass = "m1_m2_Mass"
        l1_l2_cand_pt = "m1_m2_Pt"
        l3_l4_cand_mass = "m3_m4_Mass"
        l3_l4_cand_pt = "m3_m4_Pt"
    else: 
        l1_l2_cand_mass = "e1_e2_Mass"
        l1_l2_cand_pt = "e1_e2_Pt"
        l3_l4_cand_mass = "m1_m2_Mass"
        l3_l4_cand_pt = "m1_m2_Pt"
    #zcand_name = "e1_e2_Mass" if state.count('e') >= 2 else "m1_m2_Mass"
    selector.setZCandidateBranchName(l1_l2_cand_mass,l1_l2_cand_pt,l3_l4_cand_mass,l3_l4_cand_pt)
    new_tree = tree.CopyTree(cut_string)
    new_tree.Process(selector)#, cut_string)
    entryList = selector.GetOutputList().FindObject('bestCandidates')
    new_tree.SetEntryList(entryList)
    return new_tree
def writeMetaTreeToFile(output_file, metaTree):
    output_file.cd()
    meta_dir = output_file.mkdir("metaInfo")
    meta_dir.cd()
    save_mt = metaTree.CopyTree("")
    save_mt.Write()
def skimNtuple(selections, analysis, trigger, filelist, output_file_name, deduplicate):
    current_path = os.getcwd()
    os.chdir(sys.path[0])
    ROOT.gROOT.SetBatch(True)
    output_file = ROOT.TFile(output_file_name, "RECREATE")
    ROOT.gROOT.cd()
    with open(filelist) as input_file:
        input_files = [('root://cmsxrootd.hep.wisc.edu/' + i.strip()) \
            if "store" in i[:6] else i.strip() for i in input_file.readlines()]
    metaTree = ROOT.TChain("metaInfo/metaInfo")
    for file_path in input_files:
        metaTree.Add(file_path)
    event_counts = {"initial" : {}, "selected" : {}}
    for state in ["eeee", "eemm", "mmmm"]:
        tree = ROOT.TChain("%s/ntuple" % state)
        for file_path in input_files:
            tree.Add(file_path)
        event_counts["initial"][state] = tree.GetEntries()
        cuts = ApplySelection.CutString()
        cuts.append(ApplySelection.buildCutString(state, 
            selections.split(","), analysis, trigger).getString())
        cut_string = cuts.getString()
        print "Cut string is %s " % cut_string
        #ApplySelection.setAliases(tree, state, "Cuts/aliases.json")
        event_counts["selected"][state] = writeNtupleToFile(output_file, tree, state, cut_string,
            deduplicate)
    writeMetaTreeToFile(output_file, metaTree)
    event_info = PrettyTable(["Selection", "eeee", "eemm","mmmm"])
    for selection, events in event_counts.iteritems():
        event_info.add_row([selection, events["eeee"], events["eemm"],events["mmmm"]])
    print "\nResults for selection: %s" % selections
    if deduplicate:
        print "NOTE: Events deduplicated by choosing the ordering with m_l1_l2 " \
            "closest to m_{Z}^{PDG} AFTER making full selection\n"
    else:
        print "NOTE: Events NOT deduplicated! Event may appear in multiple rows of ntuple!\n"
    print event_info.get_string()
    
    os.chdir(current_path)
def main():
    args = getComLineArgs()
    skimNtuple(args['selections'], args['analysis'], args['trigger'], args['filelist'], 
        args['output_file_name'], not args['no_deduplicate'])
    exit(0)
if __name__ == "__main__":
    main()
