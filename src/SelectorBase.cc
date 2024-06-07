#include "Analysis/VVAnalysis/interface/SelectorBase.h"
#include <boost/algorithm/string.hpp>
#include <TStyle.h>
#include <regex>
#include "TParameter.h"

void SelectorBase::Begin(TTree * /*tree*/)
{
    TString option = GetOption();
}

void SelectorBase::SlaveBegin(TTree * /*tree*/)
{
    if (GetInputList() != nullptr)
    {
        TParameter<bool> *applyScaleFactors = (TParameter<bool> *)GetInputList()->FindObject("applyScaleFacs");
        if (applyScaleFactors != nullptr && applyScaleFactors->GetVal())
        {
            SetScaleFactors();
        }
    }
}

void SelectorBase::Init(TTree *tree)
{
    if (!tree)
        return;
    fChain = tree;

    TString option = GetOption();

    if (GetInputList() != nullptr)
    {
        TNamed *ntupleType = (TNamed *)GetInputList()->FindObject("ntupleType");
        TNamed *name = (TNamed *)GetInputList()->FindObject("name");
        TNamed *chan = (TNamed *)GetInputList()->FindObject("channel");
        TNamed *selection = (TNamed *)GetInputList()->FindObject("selection");
        TNamed *year = (TNamed *)GetInputList()->FindObject("year");

        if (ntupleType != nullptr)
        {
            std::string ntupleName = ntupleType->GetTitle();
            if (ntupleName == "NanoAOD")
                ntupleType_ = NanoAOD;
            else if (ntupleName == "UWVV")
                ntupleType_ = UWVV;
            else
                throw std::invalid_argument("Unsupported ntuple type!");
        }
        else
        {
            std::cerr << "INFO: Assuming NanoAOD ntuples" << std::endl;
            ntupleType_ = NanoAOD;
        }

        if (name != nullptr)
        {
            name_ = name->GetTitle();
        }
        else
        {
            name_ = GetNameFromFile();
        }
        if (name_ == "")
        {
            std::cerr << "INFO: Using default name \"Unknown\" for file" << std::endl;
            name_ = "Unknown";
        }
        if (year != nullptr)
        {
            year_ = yearMap_[year->GetTitle()];
        }

        if (chan != nullptr)
        {
            channelName_ = chan->GetTitle();
        }
        else if (ntupleType_ == UWVV)
            channelName_ = fChain->GetTree()->GetDirectory()->GetName();
        if (selection != nullptr)
        {
            selectionName_ = selection->GetTitle();
        }
    }

    if (selectionMap_.find(selectionName_) != selectionMap_.end())
    {
        selection_ = selectionMap_[selectionName_];
    }
    else
        throw std::invalid_argument("Invalid selection!");

    isMC_ = false;
    if (name_.find("data") == std::string::npos)
    {
        isMC_ = true;
    }
    if (doSystematics_ && isMC_ && !isNonPrompt_) // isNonpromptEstimate?
        variations_.insert(systematics_.begin(), systematics_.end());

    currentHistDir_ = dynamic_cast<TList *>(fOutput->FindObject(name_.c_str()));

    if (channelMap_.find(channelName_) != channelMap_.end())
        channel_ = channelMap_[channelName_];
    else
    {
        std::string message = "Invalid channel choice! ";
        message += "Choice was " + channelName_ + "\n";
        message += "Valid choices are: ";
        for (const auto &chan : channelMap_)
            message += chan.first + ", ";
        throw std::invalid_argument(message);
    }

    ftntpName_ = name_ + "_fTreeNtuple_"+ channelName_;
    ftntp_ = dynamic_cast<TTree *>(fOutput->FindObject(ftntpName_.c_str()));
    if (ftntp_ == nullptr)
    {
        ftntp_ = new TTree(ftntpName_.c_str(),"Ntuple of selected events");
        fOutput->Add(ftntp_);
    }

    if (currentHistDir_ == nullptr)
    {
        currentHistDir_ = new TList();
        currentHistDir_->SetName(name_.c_str());
        fOutput->Add(currentHistDir_);
        // Watch for something that I hope never happens,
        size_t existingObjectPtrsSize = allObjects_.size();
        SetupNewDirectory();
        if (existingObjectPtrsSize > 0 && allObjects_.size() != existingObjectPtrsSize)
        {
            std::invalid_argument(Form("SelectorBase: Size of allObjects has changed!: %lu to %lu", existingObjectPtrsSize, allObjects_.size()));
        }
    }
    UpdateDirectory();

    SetBranches();
}

void SelectorBase::SetScaleFactors()
{
    std::invalid_argument("No scale factors defined for selector!");
}

void SelectorBase::SetBranches()
{
    if (ntupleType_ == UWVV)
        SetBranchesUWVV();
    else if (ntupleType_ == NanoAOD)
        SetBranchesNanoAOD();
}

void SelectorBase::LoadBranches(Long64_t entry, std::pair<Systematic, std::string> variation)
{
    if (ntupleType_ == UWVV)
    {
        LoadBranchesUWVV(entry, variation);
    }
    else if (ntupleType_ == NanoAOD)
        LoadBranchesNanoAOD(entry, variation);
}

Bool_t SelectorBase::Process(Long64_t entry)
{
    // std::cout<<"Does it enter Process"<<std::endl;
    for (const auto &variation : variations_)
    {
        // std::cout<<"Does it enter variations"<<std::endl;
        LoadBranches(entry, variation);
        FillHistograms(entry, variation);
    }
    return kTRUE;
}

Bool_t SelectorBase::Notify()
{
    return kTRUE;
}

float SelectorBase::GetPrefiringEfficiencyWeight(
    std::vector<float> *jetPt, std::vector<float> *jetEta)
{
    float prefire_weight = 1;
    for (size_t i = 0; i < jetPt->size(); i++)
    {
        float jPt = jetPt->at(i);
        float jEta = std::abs(jetEta->at(i));
        prefire_weight *= (1 - prefireEff_->GetEfficiency(prefireEff_->FindFixBin(jEta, jPt)));
    }
    return prefire_weight;
}

void SelectorBase::Terminate()
{
}

void SelectorBase::SlaveTerminate()
{
}
void SelectorBase::UpdateDirectory()
{
    for (TNamed **objPtrPtr : allObjects_)
    {
        if (*objPtrPtr == nullptr)
            std::invalid_argument("SelectorBase: Call to UpdateObject but existing pointer is null");
        *objPtrPtr = (TNamed *)currentHistDir_->FindObject((*objPtrPtr)->GetName());
        if (*objPtrPtr == nullptr)
            std::invalid_argument("SelectorBase: Call to UpdateObject but current directory has no instance");
    }
}

template <typename T>
void SelectorBase::InitializeHistMap(std::vector<std::string> &labels, std::map<std::string, T *> &histMap)
{
    for (auto &label : labels)
    {
        if (channel_ != Inclusive)
        {
            auto histName = getHistName(label, "", channelName_);
            histMap[histName] = {};
        }
        else
        {
            for (auto &chan : allChannels_)
            {
                auto histName = getHistName(label, "", chan);
                histMap[histName] = {};
            }
        }
    }
}

void SelectorBase::InitializeHistogramsFromConfig()
{
    TList *histInfo = (TList *)GetInputList()->FindObject("histinfo");
    if (histInfo == nullptr)
        throw std::domain_error("Can't initialize histograms without passing histogram information to TSelector");

    InitializeHistMap(hists1D_, histMap1D_);
    InitializeHistMap(weighthists1D_, weighthistMap1D_);
    InitializeHistMap(jethists1D_, jethistMap1D_);
    InitializeHistMap(jetTest2D_, jetTestMap2D_);

    for (auto &&entry : *histInfo)
    {
        TNamed *currentHistInfo = dynamic_cast<TNamed *>(entry);
        std::string name = currentHistInfo->GetName();
        std::vector<std::string> histData = ReadHistDataFromConfig(currentHistInfo->GetTitle());

        std::vector<std::string> channels = {channelName_};
        if (channel_ == Inclusive)
        {
            channels = allChannels_;
        }

        for (auto &chan : channels)
        {
            auto histName = getHistName(name, "", chan);
            if (hists2D_.find(histName) != hists2D_.end() || histMap1D_.find(histName) != histMap1D_.end())
            {
                InitializeHistogramFromConfig(name, chan, histData);
            }
            // No need to print warning for every channel
            else if (chan == channels.front())
                std::cerr << "Skipping invalid histogram " << name << std::endl;
        }
    }
}

void SelectorBase::InitializeHistogramFromConfig(std::string name, std::string channel, std::vector<std::string> histData)
{
    if (histData.size() != 4 && histData.size() != 7)
    {
        std::cerr << "Malformed data string for histogram '" << name
                  << ".' Must have form: 'Title; (optional info) $ nbins, xmin, xmax'"
                  << "\n   OR form: 'Title; (optional info) $ nbins, xmin, xmax nbinsy ymin ymax'"
                  << std::endl;
        exit(1);
    }
    std::string histName = getHistName(name, "", channel);

    int nbins = std::stoi(histData[1]);
    float xmin = std::stof(histData[2]);
    float xmax = std::stof(histData[3]);

    if (histData.size() == 4)
    {
        AddObject<TH1D>(histMap1D_[histName], histName.c_str(), histData[0].c_str(), nbins, xmin, xmax);
        // std::cout<<"doSystematics: "<<doSystematics_<<std::endl;
        // std::cout<<systHists_.size()<<std::endl;
        // std::cout<<"histName: "<<histName<<std::endl;
        if (doSystematics_ && !isNonPrompt_ && std::find(systHists_.begin(), systHists_.end(), name) != systHists_.end())
        {
            // std::cout<<"are there systHists_"<<std::endl;
            for (auto &syst : systematics_)
            {
                // std::cout<<"systHists getting filled?"<<std::endl;
                std::string syst_histName = getHistName(name, syst.second, channel);
                histMap1D_[syst_histName] = {};
                AddObject<TH1D>(histMap1D_[syst_histName], syst_histName.c_str(),
                                histData[0].c_str(), nbins, xmin, xmax);
                // TODO: Cleaner way to determine if you want to store systematics for weighted entries
                // if (isaQGC_ && doaQGC_ && (weighthistMap1D_.find(name) != weighthistMap1D_.end())) {
                //    std::string weightsyst_histName = name+"_lheWeights_"+syst.second;
                //    AddObject<TH2D>(weighthistMap1D_[syst_histName],
                //        (weightsyst_histName+"_"+channel).c_str(), histData[0].c_str(),
                //        nbins, xmin, xmax, 1000, 0, 1000);
                //}
            }
        }
        // Weight hists must be subset of 1D hists!
        // std::cout<<"size of weighthistMap1D_: "<<weighthistMap1D_.size()<<std::endl;
        if (isMC_ && !isNonPrompt_ && (weighthistMap1D_.find(histName) != weighthistMap1D_.end()))
        {
            // std::cout<<"Is weightHists getting filled?"<<std::endl;
            AddObject<TH2D>(weighthistMap1D_[histName],
                            (name + "_lheWeights_" + channel).c_str(), histData[0].c_str(),
                            nbins, xmin, xmax, 1000, 0, 1000);
        }

        if (isMC_ && !isNonPrompt_ && (jethistMap1D_.find(histName) != jethistMap1D_.end()))
        {
            // std::cout<<"Is weightHists getting filled?"<<std::endl;
            AddObject<TH2D>(jethistMap1D_[histName],
                            (name + "_jetsysts_" + channel).c_str(), histData[0].c_str(),
                            nbins, xmin, xmax, 4, 0, 4); // 4 systs for JES and JER
        }

        // jetPt vs jetEta for single nJets, also define hist in hists1D_ to pass checks in InitializeHistogramsFromConfig()
        // Now also include jet Eta vs jet Phi for HEM15/16 studies
        float jEta_binning[] = {0.0, 1.5, 2.4, 3.2, 4.7};
        float jPt_binning[] = {30., 50.,70., 100., 150.,200., 300., 500.};
        float jEta_binning2[] = {0.0, 1.5, 3.0, 4.7};
        float jPt_binning2[] = {30., 50., 100., 170., 300.};
        float jPhi_binning[] = {-3.15, -2.5, -1.57, -0.87, 0.0, 0.87, 1.57, 2.5, 3.15};
        float jEta_binning3[] = {-4.7, -3.0, -2.5, -1.3, 0.0, 1.3, 2.5, 3.0, 4.7};
        if (jetTestMap2D_.find(histName) != jetTestMap2D_.end())
        {
            if (histName.find("HEM") != std::string::npos)
            {
                if (histName.find("HEM2") != std::string::npos)
                {
                    AddObject<TH2D>(jetTestMap2D_[histName],
                                    (name + "_JetEtaVsPhiFillPt_" + channel).c_str(), histData[0].c_str(),
                                    8, jPhi_binning, 8, jEta_binning3);
                }
                else
                {
                    AddObject<TH2D>(jetTestMap2D_[histName],
                                    (name + "_JetEtaVsPhi_" + channel).c_str(), histData[0].c_str(),
                                    8, jPhi_binning, 8, jEta_binning3);
                }
            }

            else if (histName.find("N1") != std::string::npos)
            {
                AddObject<TH2D>(jetTestMap2D_[histName],
                                (name + "_vsJetEta_" + channel).c_str(), histData[0].c_str(),
                                5, jPt_binning, 4, jEta_binning);
            }
            else
            {
                AddObject<TH2D>(jetTestMap2D_[histName],
                                (name + "_vsJetEta_" + channel).c_str(), histData[0].c_str(),
                                4, jPt_binning2, 3, jEta_binning2);
            }
            // std::cout<<"Is weightHists getting filled?"<<std::endl;
            // AddObject<TH2D>(jetTestMap2D_[histName],
            //     (name+"_vsJetEta_"+channel).c_str(), histData[0].c_str(),
            //     nbins, xmin, xmax, 50, 0, 5);
        }
    }
    else
    {
        int nbinsy = std::stoi(histData[4]);
        float ymin = std::stof(histData[5]);
        float ymax = std::stof(histData[6]);
        AddObject<TH2D>(hists2D_[histName], histName.c_str(), histData[0].c_str(), nbins, xmin, xmax,
                        nbinsy, ymin, ymax);
        if (doSystematics_ && !isNonPrompt_ && std::find(systHists2D_.begin(), systHists2D_.end(), histName) != systHists2D_.end())
        {
            for (auto &syst : systematics_)
            {
                std::string syst_hist_name = name + "_" + syst.second + "_" + channel;
                hists2D_[syst_hist_name] = {};
                AddObject<TH2D>(hists2D_[syst_hist_name], syst_hist_name.c_str(),
                                histData[0].c_str(), nbins, xmin, xmax, nbinsy, ymin, ymax);
            }
        }
        // 3D weight hists must be subset of 2D hists!
        if (isMC_ && !isNonPrompt_ && (weighthistMap2D_.find(histName) != weighthistMap2D_.end()))
        {
            AddObject<TH3D>(weighthistMap2D_[histName],
                            (name + "_lheWeights_" + channel).c_str(), histData[0].c_str(),
                            nbins, xmin, xmax, nbinsy, ymin, ymax, 1000, 0, 1000);
        }
    }
}

std::vector<std::string> SelectorBase::ReadHistDataFromConfig(std::string histDataString)
{
    std::vector<std::string> histData;
    boost::split(histData, histDataString, boost::is_any_of("$"));
    std::vector<std::string> binInfo;
    if (histData.size() != 2)
        return {};

    boost::split(binInfo, histData[1], boost::is_any_of(","));

    histData.pop_back();
    for (const auto &x : binInfo)
    {
        histData.push_back(x);
    }

    return histData;
}

void SelectorBase::SetupNewDirectory()
{
    if (addSumweights_)
        AddObject<TH1D>(sumWeightsHist_, "sumweights", "sumweights", 1, 0, 10);
}

std::string SelectorBase::getHistName(std::string histName, std::string variationName)
{
    return getHistName(histName, variationName, "");
}

std::string SelectorBase::getHistName(std::string histName, std::string variationName, std::string channel)
{
    if (channel == "")
        channel = channelName_;
    if (variationName != "")
        return histName + "_" + variationName + "_" + channel;
    return histName + "_" + channel;
}

//Copy of getHistName, used for naming ntuple TTree branch
std::string SelectorBase::getBranchName(std::string bName, std::string variationName)
{
    return getBranchName(bName, variationName, "");
}

std::string SelectorBase::getBranchName(std::string bName, std::string variationName, std::string channel)
{
    //if (channel == "")
    //    channel = channelName_;
    if (channel != ""){
        channel = "_" + channel;  
    }

    //*Don't add variation name to branch name except for weight, but use weight branches to distinguish syst from central *
    if (variationName != "" && bName == "weight")
        return bName + "_" + variationName + channel;

    return bName + channel;
}
