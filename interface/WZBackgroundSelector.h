#ifndef WZBackgroundSelector_h
#define WZBackgroundSelector_h

#include <vector>
#include "Analysis/WZAnalysis/interface/ScaleFactor.h"
#include "Analysis/WZAnalysis/interface/WZSelector.h"

class WZBackgroundSelector : public WZSelector {
public :
    ScaleFactor * fakeRate_allE_;
    ScaleFactor * fakeRate_allMu_;
    
    TH1D* mjjHistPPF_;
    TH1D* mjjHistPFP_;
    TH1D* mjjHistFPP_;
    TH1D* mjjHistFFP_;
    TH1D* mjjHistFPF_;
    TH1D* mjjHistPFF_;
    TH1D* mjjHistFFF_;
    TH1D* ZMassHistPPF_;
    TH1D* ZMassHistPFP_;
    TH1D* ZMassHistFPP_;
    TH1D* ZMassHistFFP_;
    TH1D* ZMassHistFPF_;
    TH1D* ZMassHistPFF_;
    TH1D* ZMassHistFFF_;
    
    WZBackgroundSelector(TTree * /*tree*/ =0) { }
    virtual ~WZBackgroundSelector() { }
    virtual Bool_t  Process(Long64_t entry) override;
    virtual void    SetupNewDirectory() override;
    virtual void    SlaveBegin(TTree *tree) override;

    ClassDefOverride(WZBackgroundSelector,0);
private:
    float getl1FakeRate();
    float getl2FakeRate();
    float getl3FakeRate();
    bool IsFPPRegion();
    bool IsPFPRegion();
    bool IsPPFRegion();
    bool IsFPFRegion();
    bool IsFFPRegion();
    bool IsPFFRegion();
    bool IsFFFRegion();
    float getEventWeight();
};

#endif
