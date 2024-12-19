#ifndef ZZSelector_h
#define ZZSelector_h

#include <vector>
#include "Analysis/VVAnalysis/interface/ZZSelectorBase.h"
#include <TH3.h>
#include <TTreeFormula.h>

class ZZSelector : public ZZSelectorBase
{
public:
    bool applyFullSelection_ = true;
    bool isaTGC_ = false;
    bool doaTGC_ = false;
    bool applyPUSF_ = false;
    bool applyPUSFNtp_ = true;

    std::vector<float> *scaleWeights = NULL;
    std::vector<float> *pdfWeights = NULL;
    std::vector<float> lheWeights;
    unsigned int weight_info_;
    TTreeFormula   *fCutFormula;
    bool passCurrentTrig;
    float dEtajj;
    float dEtajj_jesUp;
    float dEtajj_jesDown;
    float dEtajj_jerUp;
    float dEtajj_jerDown;
    float dRjj;
    float dRjj_jesUp;
    float dRjj_jesDown;
    float dRjj_jerUp;
    float dRjj_jerDown;
    float zep3l;
    float zep3l_jesUp;
    float zep3l_jesDown;
    float zep3l_jerUp;
    float zep3l_jerDown;
    Float_t mjj;
    Float_t jetPUSFmulfac;
    Float_t mjj_jesUp;
    Float_t mjj_jesDown;
    Float_t mjj_jerUp;
    Float_t mjj_jerDown;
    unsigned int nJets;
    unsigned int nJets_jesUp;
    unsigned int nJets_jesDown;
    unsigned int nJets_jerUp;
    unsigned int nJets_jerDown;
    std::vector<float> *jetEta = NULL;
    std::vector<float> *jetEta_jesUp = NULL;
    std::vector<float> *jetEta_jesDown = NULL;
    std::vector<float> *jetEta_jerUp = NULL;
    std::vector<float> *jetEta_jerDown = NULL;
    std::vector<float> *jetPhi = NULL;
    std::vector<float> *jetPhi_jesUp = NULL;
    std::vector<float> *jetPhi_jesDown = NULL;
    std::vector<float> *jetPhi_jerUp = NULL;
    std::vector<float> *jetPhi_jerDown = NULL;
    std::vector<float> *jetPt = NULL;
    std::vector<float> *jetPt_jesUp = NULL;
    std::vector<float> *jetPt_jesDown = NULL;
    std::vector<float> *jetPt_jerUp = NULL;
    std::vector<float> *jetPt_jerDown = NULL;
    std::vector<float> *jetCSVv2 = NULL;
    UInt_t nvtx;
    Float_t Mass;
    Float_t Pt;
    Float_t Eta;
    float dPhiZZ; //DeltaPhi between Z1 and Z2
    float dRZZ;   //DeltaR between Z1 and Z2

    std::vector<int> *jetPUID = NULL;
    std::vector<int> *isGenJetMatched = NULL;

    TBranch *b_jetPUID;
    TBranch *b_isGenJetMatched;
    TBranch *b_jetPUSFmulfac;
    TBranch *b_mjj;
    TBranch *b_mjj_jesUp;
    TBranch *b_mjj_jesDown;
    TBranch *b_mjj_jerUp;
    TBranch *b_mjj_jerDown;

    TBranch *b_nJets;
    TBranch *b_nJets_jesUp;
    TBranch *b_nJets_jesDown;
    TBranch *b_nJets_jerUp;
    TBranch *b_nJets_jerDown;

    TBranch *b_jetPt;
    TBranch *b_jetPt_jesUp;
    TBranch *b_jetPt_jesDown;
    TBranch *b_jetPt_jerUp;
    TBranch *b_jetPt_jerDown;

    TBranch *b_jetEta;
    TBranch *b_jetEta_jesUp;
    TBranch *b_jetEta_jesDown;
    TBranch *b_jetEta_jerUp;
    TBranch *b_jetEta_jerDown;

    TBranch *b_nvtx;
    TBranch *b_Mass;
    TBranch *b_Pt;
    TBranch *b_Eta;
    //TBranch* b_jetPt;
    //TBranch* b_jetEta;
    //TBranch* b_mjj;
    TBranch *b_jetPhi;
    int cen_count = 0;
    float count = 0;

    // Readers to access the data (delete the ones you do not need).
    virtual void Init(TTree *tree) override;
    virtual void SetupNewDirectory() override;

    ClassDefOverride(ZZSelector, 0);

protected:
    virtual void SetBranchesUWVV() override;
    void LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void ApplyScaleFactors();
    bool PassesZZSelection(bool nonPrompt);
    bool PassesZZSelectionLoose(bool nonPrompt);
    bool Passes4eExtraCut();
    bool Passes2e2mExtraCut(Long64_t entry);
    bool PassesZZjjSelection();
    bool PassesZZjjBaselineSelection();
    bool PassesHZZSelection(bool nonPrompt);
    unsigned int GetLheWeightInfo();
    void SetVariables(Long64_t entry);
    //bool DuplicatedEvent();
    bool TightZZLeptons();
    bool ZZSelection();
    bool Z4lSelection();
    bool ZSelection();
    bool HZZSIPSelection();
    //bool SmartCut();
    bool HZZLowMass();
    bool HZZMediumMass();
    bool TestMuons();
    void ShiftEfficiencies(Systematic variation);
};

#endif
