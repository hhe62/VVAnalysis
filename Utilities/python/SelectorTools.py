#!/usr/bin/env python
import ROOT
import glob
import datetime
import ConfigureJobs, OutputTools
import sys
import os
import multiprocessing
import subprocess
import logging
import UserInput

class SelectorDriver(object):
    def __init__(self, analysis, selection, input_tier, year):
        # TODO: Make this a configurable argument
        #logging.basicConfig(level=logging.DEBUG)
        logging.basicConfig(level=logging.WARNING)

        selector_map = {
            "WZxsec2016" : "WZSelector",
            "Zstudy" : "ZSelector",
            "Zstudy_2016" : "ZSelector",
            "Zstudy_2017" : "ZSelector",
            "ZZ4l2016" : "ZZSelector",
            "ZZ4l2017" : "ZZSelector",
            "ZZ4l2018" : "ZZSelector",
            "WGen" : "WGenSelector",
            "ZGen" : "ZGenSelector",
            "ThreeLep" : "ThreeLepSelector",
        }

        self.analysis = analysis
        self.selection = selection
        self.input_tier = input_tier
        if analysis not in selector_map.keys():
            raise ValueError("Analysis does not point to " \
                "a defined selector. Please edit "
                "Utilities/python/SelectorTools.py to add it.")
        self.selector_name = selector_map[analysis]
        self.addSumweights = True
        self.ntupleType = "NanoAOD"
        self.year = year
        self.numCores = 1
        self.channels = ["Inclusive"]
        self.outfile_name = "temp.root"
        self.outfile = 0
        self.datasets = {}

    # Needed to parallelize class member function, see
    # https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map
    def __call__(self, args):
        self.processDatasetHelper(args)

    def tempfileName(self, dataset):
        return "temp_%s_%s" % (dataset, self.outfile_name)

    def setChannels(self, channels):
        self.channels = channels

    def isBackground(self):
        self.selector_name = self.selector_name.replace("Selector", "BackgroundSelector")

    def isGen(self):
        self.selector_name = self.selector_name.replace("BackgroundSelector", "GenSelector")
    
    def isFake(self):
        self.selector_name = self.selector_name.replace("ZZ", "FakeRate")

    def outputFile(self):
        return self.outfile
    
    def setOutputfile(self, outfile_name):
        if self.outfile:
            self.outfile.Close()
        self.outfile_name = outfile_name
        self.outfile = ROOT.gROOT.FindObject(outfile_name)
        if not self.outfile:
            self.outfile = ROOT.TFile.Open(outfile_name, "recreate")
        self.current_file = self.outfile

    def addTNamed(self, name, title):
        if not self.inputs:
            raise ValueError("Input list is not defined!")
        named = self.inputs.FindObject(name)
        if named:
            named.SetTitle(title)
        else:
            named = ROOT.TNamed(name, title)
            self.inputs.Add(named)

    def setInputs(self, inputs):
        self.inputs = ROOT.TList()
        for inp in inputs:
            self.inputs.Add(inp)
        self.addTNamed("ntupleType", self.ntupleType)
        self.addTNamed("selection", self.selection)
        self.addTNamed("year", self.year)
        
    def setNtupeType(self, ntupleType):
        self.ntupleType = ntupleType
        self.addTNamed("ntupleType", self.ntupleType)

    def setNumCores(self, numCores):
        self.numCores = numCores

    def setSumWeights(self, addSumweights):
        self.addSumweights = addSumweights

    def setFileList(self, list_of_files, nPerJob, jobNum):
        if not os.path.isfile(list_of_files):
            raise ValueError("%s is not a valid file." % list_of_files)
        filelist = [f.split("#")[0].strip() for f in open(list_of_files).readlines()]
        # Remove empty/commented lines
        filelist = filter(lambda  x: len(x) > 2, filelist)
        nPerJob = int(nPerJob)
        if nPerJob < 1:
            raise ValueError("Number of files per job must be >= 1.")
        jobNum = int(jobNum)
        maxNum = len(filelist)
        firstEntry = nPerJob*jobNum
        if firstEntry > maxNum:
            raise ValueError("The first file to process (nPerJob*jobNum) = (%i*%i)" % (nPerJob, jobNum) \
                    + " is greater than the number of entries in file %s (%s)." % (list_of_files, maxNum))
        lastEntry = min(nPerJob*(jobNum+1), maxNum)
        
        for line in filelist[firstEntry:lastEntry]:
            if "@" not in line:
                dataset = "Unknown"
                file_path = line
            else:
                # Intended for running specified files, use the format name:file
                dataset, file_path = line.split("@")
            if dataset not in self.datasets.keys():
                self.datasets[dataset] = [file_path]
            else:
                self.datasets[dataset].append(file_path)

    def setDatasets(self, datalist):
        self.clearDatasets()
        datasets = ConfigureJobs.getListOfFiles(datalist, self.input_tier)
        for dataset in datasets:
            if "@" in dataset:
                dataset, file_path = [f.strip() for f in dataset.split("@")]
            else:
                try:
                    file_path = ConfigureJobs.getInputFilesPath(dataset, 
                        self.input_tier, self.analysis)
                except ValueError as e:
                    logging.warning(e)
                    continue
            self.datasets[dataset] = [file_path]

    def clearDatasets(self):
        self.datasets.clear()

    def applySelector(self):
        for chan in self.channels:
            self.addTNamed("channel", chan)
            logging.info("Processing channel %s" % chan)
            if self.numCores > 1:
                self.processParallelByDataset(self.datasets, chan)
            else: 
                for dataset, file_path in self.datasets.iteritems():
                    self.processDataset(dataset, file_path, chan)

                #Multiple files in 1-thread case also, and do the combination here instead of inside processParallel
                tempTreeFiles = ["FilledNtuples/TreeFile_"+self.tempfileName(dataset) for dataset in self.datasets]
                
                self.combineParallelTreeFiles(tempTreeFiles, chan)
        if len(self.channels) > 1 and self.numCores > 1:
            tempfiles = [self.outfile_name.replace(".root", "_%s.root" % c) for c in self.channels]
            self.combineParallelFiles(tempfiles, "Inclusive")

        #Same tree files regardless of single of multi-thread 
        tempTreefiles = ["FilledNtuples/TreeFile_"+self.outfile_name.replace(".root", "_%s.root" % c) for c in self.channels]
        self.combineParallelTreeFiles(tempTreefiles, "Inclusive")

    def processDataset(self, dataset, file_path, chan):
        logging.info("Processing dataset %s" % dataset)
        select = getattr(ROOT, self.selector_name)() #! selector initialized here
        select.SetInputList(self.inputs)
        self.addTNamed("name", dataset)
        # Only add for one channel
        addSumweights = self.addSumweights and self.channels.index(chan) == 0 and "data" not in dataset and "Background" not in self.selector_name and "Gen" not in self.selector_name
        if addSumweights:
            #print "Should only go here once!"
            sumweights_hist = ROOT.gROOT.FindObject("sumweights")
            # Avoid accidentally combining sumweights across datasets
            if sumweights_hist:
                del sumweights_hist
            if not self.outfile:
                self.outfile = ROOT.TFile.Open(self.outfile_name)
            sumweights_hist = self.outfile.Get("%s/sumweights" % dataset)
            if not sumweights_hist:
                sumweights_hist = ROOT.TH1D("sumweights", "sumweights", 1000, 0, 1000)
            sumweights_hist.SetDirectory(ROOT.gROOT)
        self.processLocalFiles(select, file_path, addSumweights, chan)
        output_list = select.GetOutputList()
        name = self.inputs.FindObject("name").GetTitle()
        dataset_list = output_list.FindObject(name)
        ntuple_written = output_list.FindObject(name + "_fTreeNtuple_"+ chan)
        if not dataset_list or dataset_list.ClassName() != "TList":
            logging.warning("No output found for dataset %s" % dataset)
            dataset_list = output_list.FindObject("Unknown")
            if dataset_list and dataset_list.ClassName() == "TList":
                logging.warning('Falling back to dataset "Unknown"')
            else:
                logging.warning('Skipping dataset %s' % dataset)
                return False

        if not ntuple_written or ntuple_written.ClassName() != "TTree":
            logging.warning("No filled ntuple found for dataset %s" % dataset)
            
            #logging.warning('For filled ntuple skipping dataset %s' % dataset)
            

        if addSumweights:
            dataset_list.Add(ROOT.gROOT.FindObject("sumweights"))
        if self.numCores > 1:
            self.outfile.Close()
            chanNum = self.channels.index(chan)
            #Shouldn't matter whether recreate or update for hists, since the temp files will be deleted after hadd
            self.current_file = ROOT.TFile.Open(self.tempfileName(dataset), "recreate" if chanNum == 0 else "update")
        
        #Always create temp file for ntuple case regardless of whether multithread
        if not os.path.isdir("FilledNtuples"):
            try:
                os.mkdir("FilledNtuples")
            except:
                print("Failed to create FilledNtuples folder. Could be created by other thread at the same time.")
        
        OutputTools.writeOutputListItem(dataset_list, self.current_file)
        
        if self.current_file != self.outfile:
            self.current_file.Close()

        self.current_Treefile = ROOT.TFile.Open("FilledNtuples/TreeFile_"+self.tempfileName(dataset), "recreate")
        self.current_Treefile.cd()
        ntuple_written.Write() #write TTree
        self.current_Treefile.Close()
        
        #dataset_list.Delete()
        #ntuple_written.Delete()
        output_list.Delete()
        return True

    def getFileNames(self, file_path):

        #! Don't know why it was implemented this way, but xrootd and xrootd_user here are expected to be false even if path is correct
        #! Major part that is working is the glob.glob(file_path) part

        xrootd = "/store" in file_path.split("/hdfs/")[0][:7]
        xrootd_user = "/store/user" in file_path.split("/hdfs/")[0][:12]
        if not (xrootd or os.path.isfile(file_path) or os.path.isdir(file_path.rsplit("/", 1)[0].rstrip("/*"))):
            raise ValueError("Invalid path! Skipping dataset. Path was %s" 
                % file_path)

        # Assuming these are user files on HDFS, otherwise it won't work
        if (xrootd and not xrootd_user):
            xrd = 'root://%s/' % ConfigureJobs.getXrdRedirector()
            filenames = [xrd + file_path]
            return filenames
        filenames =  glob.glob(file_path) if not xrootd_user else \
                ConfigureJobs.getListOfHDFSFiles(file_path)
        filenames = ['root://cmsxrootd.hep.wisc.edu/' + f if "/store/user" in f[0:12] else f for f in filenames]
        return filenames

    def getTreeName(self, chan):
        # TODO: Fix this! This is an extremely ineffient way to separate the eemm and mmee
        # since it involves reading the file an extra time
        channel = chan if "mmee" not in chan else chan.replace("mmee", "eemm")
        return "Events" if self.ntupleType == "NanoAOD" else ("%s/ntuple" % channel)

    def combineParallelFiles(self, tempfiles, chan):
        tempfiles = filter(os.path.isfile, tempfiles)
        outfile = self.outfile_name
        if chan != "Inclusive" and len(self.channels) != 1:
            outfile = self.outfile_name.replace(".root", "_%s.root" % chan)
        rval = subprocess.call(["hadd", "-f", outfile] + tempfiles)
        if rval == 0:
            map(os.remove, tempfiles)
        else:
            raise RuntimeError("Failed to collect data from parallel run")

    def combineParallelTreeFiles(self, tempfiles, chan):
        tempfiles = filter(os.path.isfile, tempfiles)
        
        if chan != "Inclusive":
            outfile = "FilledNtuples/TreeFile_"+self.outfile_name.replace(".root", "_%s.root" % chan)
        else:
            outfile = "FilledNtuples/TreeFile_%s_"%self.selector_name + self.outfile_name.replace(".root", "_%s.root" % chan)
            
        outfileMerged = ROOT.TFile.Open(outfile,"recreate")
        outfileMerged.Close() 
        for f in tempfiles:
            os.system("rootcp %s:* %s"%(f,outfile))
        
       
        map(os.remove, tempfiles)
       

    def processParallelByDataset(self, datasets, chan):
        numCores = min(self.numCores, len(datasets))
        p = multiprocessing.Pool(processes=self.numCores)
        p.map(self, [[dataset, f, chan] for dataset, f in datasets.iteritems()])
        # Store arrays in temp files, since it can get way too big to keep around in memory
        tempfiles = [self.tempfileName(d) for d in datasets] 
        tempTreeFiles = ["FilledNtuples/TreeFile_"+self.tempfileName(d) for d in datasets]
        self.combineParallelFiles(tempfiles, chan)
        self.combineParallelTreeFiles(tempTreeFiles, chan)

    # Pool.map can only take in one argument, so expand the array
    def processDatasetHelper(self, args):
        self.processDataset(*args)

    def processLocalFiles(self, selector, file_path, addSumweights, chan,):
        filenames = []
        for entry in file_path:
            filenames.extend(self.getFileNames(entry))
        for i, filename in enumerate(filenames):
            self.processFile(selector, filename, addSumweights, chan, i+1)

    def processFile(self, selector, filename, addSumweights, chan, filenum=1):
        logging.debug("Processing file: %s" % filename)
        rtfile = ROOT.TFile.Open(filename)
        if not rtfile or not rtfile.IsOpen() or rtfile.IsZombie():
            raise IOError("Failed to open file %s!" % filename)
        tree_name = self.getTreeName(chan)
        tree = rtfile.Get(tree_name)
        if not tree:
            raise ValueError(("tree %s not found for file %s. " \
                    "Either the file is corrupted or the ntupleType (%s) is wrong.") 
                % (tree_name, filename, self.ntupleType)
            )
        logging.debug("Processing tree %s for file %s." % (tree.GetName(), rtfile.GetName()))
        
        #Codes for trigger test
        #aliases = UserInput.readJson("Cuts/%s/aliases.json" % self.analysis)
        #for nameAlias, valueAlias in aliases["Event"].iteritems():
        #    tree.SetAlias(nameAlias, valueAlias)
        
        #TriggerStr = "(singleIsoMuPass || doubleMuDZPass || tripleMuPass)"
        #TriggerStr = "((singleIsoMuPass || doubleMuDZPass || tripleMuPass) && (singleEPass || doubleEPass || tripleEPass))"
        #TriggerStr = "(!(singleIsoMuPass || doubleMuDZPass || tripleMuPass) && (singleEPass || doubleEPass || tripleEPass))"
        TriggerStr = "((singleIsoMuPass || doubleMuDZPass || tripleMuPass) && (singleEPass))"
        #TriggerStr = "((singleIsoMuPass || doubleMuDZPass || tripleMuPass) && (doubleEPass))"
        #TriggerStr = "((singleIsoMuPass || doubleMuDZPass || tripleMuPass) && (singleEPass || doubleEPass))"
        #TriggerStr = "((singleIsoMuPass || doubleMuDZPass || tripleMuPass) && (tripleEPass))"
        tree.Process(selector, "")
        #tree.Process(selector, TriggerStr)
        logging.debug("Processed with selector %s." % selector.GetName())
        if addSumweights:
            self.fillSumweightsHist(rtfile, filenum)
        logging.debug("Added sumweights hist.")
        rtfile.Close()

    # You can use filenum to index the files and sum separately, but it's not necessary
    def fillSumweightsHist(self, rtfile, filenum=1):
        if self.ntupleType == "NanoAOD":
            sumweights_branch = "genEventSumw"
            meta_tree_name = "Runs"
        else:
            sumweights_branch = "summedWeights"
            meta_tree_name = "metaInfo/metaInfo"
        meta_tree = rtfile.Get(meta_tree_name)
        ROOT.gROOT.cd()
        sumweights_hist = ROOT.gROOT.FindObject("sumweights")
        tmplabel = sumweights_hist.GetName()+"_i"
        tmpweights_hist = sumweights_hist.Clone(tmplabel)
        meta_tree.Draw("%i>>%s" % (filenum, tmplabel), sumweights_branch)
        sumweights_hist.Add(tmpweights_hist)
