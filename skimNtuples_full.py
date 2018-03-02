#!/usr/bin/env python
import ROOT
import argparse
import os
import json
import sys
from pprint import pprint
from Utilities.python import ApplySelection
from Utilities.python import ApplyAllSelections
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
def calculateTreeEntries(tree, state, cut_string):
    #state_dir = output_file.mkdir(state)
    #state_dir.cd()
    save_tree = tree.CopyTree(cut_string) 
    #if cut_string != "deduplicate" else \
    #    getDeduplicatedTree(tree, state, cut_string).CopyTree("")
    #save_tree.Write()
    # Remove AutoSaved trees
    #output_file.Purge()
    #ROOT.gROOT.cd()
    return save_tree.GetEntries()
def getDeduplicatedTree(tree, state, cut_string):
    selector = ROOT.disambiguateFinalStates()
    if state.count('e') > 2:
        l1_l2_cand_mass = "e1_e2_Mass"
        l1_cand_pt = "e1Pt"
        l2_cand_pt = "e2Pt"
        l3_l4_cand_mass = "e3_e4_Mass"
        l3_cand_pt = "e3Pt"
        l4_cand_pt = "e4Pt"
    elif state.count('m') > 2:
        l1_l2_cand_mass = "m1_m2_Mass"
        l1_cand_pt = "m1Pt"
        l2_cand_pt = "m2Pt"
        l3_l4_cand_mass = "m3_m4_Mass"
        l3_cand_pt = "m3Pt"
        l4_cand_pt = "m4Pt"
    else: 
        l1_l2_cand_mass = "e1_e2_Mass"
        l1_cand_pt = "e1Pt"
        l2_cand_pt = "e2Pt"
        l3_l4_cand_mass = "m1_m2_Mass"
        l3_cand_pt = "m1Pt"
        l4_cand_pt = "m2Pt"
    #zcand_name = "e1_e2_Mass" if state.count('e') >= 2 else "m1_m2_Mass"
    selector.setZCandidateBranchName(l1_l2_cand_mass,l1_cand_pt,l2_cand_pt,l3_l4_cand_mass,l3_cand_pt,l4_cand_pt)
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
        print input_files
    metaTree = ROOT.TChain("metaInfo/metaInfo")
    for file_path in input_files:
        metaTree.Add(file_path)
    #event_counts for writing the tree to the File with all selections together
    event_counts = {"initial" : {}, "selected" : {}}
    #This dictionary is for saving entries for each selection
    EventCounts={}
    #creating different event_counts dictionaries for each selection
    EventCounts["NoSelection"] = {"eeee" : {}, "eemm" : {},"mmmm" : {}}
    MakeSelection=''
    for printS in selections.split(","):
        MakeSelection+=printS
        EventCounts[MakeSelection] = {"eeee" : {}, "eemm" : {},"mmmm" : {}}
    #for selection in selections.split(","):
    #    EventCounts[selection] = {"eeee" : {}, "eemm" : {},"mmmm" : {}}
    #print EventCounts
    for state in ["eeee", "eemm", "mmmm"]:
        tree = ROOT.TChain("%s/ntuple" % state)
        #Counttree is not written to the Ntuple but just used for calculating entries
        Counttree = ROOT.TChain("%s/ntuple" % state)
        for file_path in input_files:
            tree.Add(file_path)
            Counttree.Add(file_path)
        event_counts["initial"][state] = tree.GetEntries()
        cuts = ApplyAllSelections.CutString()
        cuts.append(ApplyAllSelections.buildCutString(state, 
           selections.split(","), analysis, trigger).getString())
        cuts_string = cuts.getString()
        #print "Cut string is %s " % cuts_string
        event_counts["selected"][state] = writeNtupleToFile(output_file, tree, state, cuts_string,deduplicate)
        selection=[]
        selectionstring=''
        #EventCounts will take the Counttree
        EventCounts["NoSelection"][state] = Counttree.GetEntries()
        #count=0
        for s in selections.split(","):
            selection.append(s)
            selectionstring+=s
            #print selection
            cut = ApplySelection.CutString()
            cut.append(ApplySelection.buildCutString(state, 
            selection, analysis, trigger).getString())
            cut_separate = cut.getString()
            #print "Cut string is %s " % cut_separate
            #if('empty' in selection and count==0):
            #    count+=1
            #    EventCounts[selectionstring][state]["PassTrigger"] = calculateTreeEntries(Counttree, state, cut_separate)
            if('deduplicated' in selection):
                print selectionstring
                #print cut_separate
                specialtree = getDeduplicatedTree(Counttree, state, cut_separate).CopyTree("") 
                EventCounts[selectionstring][state] = specialtree.GetEntries()
                #deduplication=selection
                #EventCounts[selection][state]["%s"%selection] = calculateTreeEntries(tree, state, cut_separate, deduplication)
            else:
                #print selection
                #print cut_separate
                EventCounts[selectionstring][state] = calculateTreeEntries(Counttree, state, cut_separate)
                #event_counts[state]["%s"%selection] = writeNtupleToFile(output_file, tree, state, cut_separate,deduplicate)
    #print selectionstring
    pprint(EventCounts)
    with open('CutFlow.json', 'w') as fp:
        json.dump(EventCounts,fp)
    #print event_counts
    #This loop over selections again to produce .json files for each selection containing events selected
    #Copy the logic above to produce json files with subsequent cuts
    #printSelection=''
    #Count=0
    #for printS in selections.split(","):
    #    printSelection+=printS
    #    print printSelection
        #if('empty' in printSelection and Count==0):
        #    Count+=1
            #with open('PassTrigger.json', 'w') as fp:

            #event_info = PrettyTable(["Channels", "Initial", "PassTrigger"])
            ##print EventCounts[selection]
            #for channel, events in EventCounts[printSelection].iteritems():
            #    event_info.add_row([channel, Events["NoSelection"][channel]["initial"], events["PassTrigger"]])
    
            #print "\nResults for selection: PassTrigger"

            #print event_info.get_string()
        #else:
        #    #with open('%s.json' %printSelection, 'w') as fp:
        #    with open('CutFlow.json', 'a') as fp:
        #        json.dump(EventCounts[printSelection],fp)
        
           # event_info = PrettyTable(["Channels", "Initial", "%s"%printSelection])

           # for channel, events in EventCounts[printSelection].iteritems():
           #     event_info.add_row([channel, events["initial"], events["%s"%printSelection]])
    
           # print "\nResults for selection: %s" % printSelection

           # print event_info.get_string()

    #for states,status in event_counts.iteritems():
    #        #print states,status
    #        print json.dumps({states:{'processed': status['initial'], 'selected': status['selected']},})
    writeMetaTreeToFile(output_file, metaTree)
    with open('AllSelections.json', 'w') as fp:
        json.dump(event_counts,fp)
    Allevent_info = PrettyTable(["Selection", "eeee", "eemm","mmmm"]) 
    for allselections, events in event_counts.iteritems():
        Allevent_info.add_row([allselections, events["eeee"], events["eemm"],events["mmmm"]]) 
    print "\nResults for selection: %s" % selections
    if deduplicate:
        print "NOTE: Events deduplicated by choosing the ordering with m_l1_l2 " \
            "closest to m_{Z}^{PDG} AFTER making full selection\n"
    else:
        print "NOTE: Events NOT deduplicated! Event may appear in multiple rows of ntuple!\n"
    print Allevent_info.get_string()
    
    os.chdir(current_path)
def main():
    args = getComLineArgs()
    print args['filelist']
    skimNtuple(args['selections'], args['analysis'], args['trigger'], args['filelist'], 
        args['output_file_name'], not args['no_deduplicate'])
    exit(0)
if __name__ == "__main__":
    main()