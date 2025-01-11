#include "Analysis/VVAnalysis/interface/ZZSelector.h"
#include "TLorentzVector.h"
#include <boost/algorithm/string.hpp>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Important global variables:
// doSystematics_ :turn on syst calc., much more time consuming
// applyPUSFNtp_  :apply PU SF from ntuple calculated in UWVV; default to true
// writeNtp_      :also write selected events in ntuple files
//
// Local variables:
// processFullRange  : if false then processing for full m4l range in FillHistograms function is not done.
// writeNtpFullRange : if true, only write ntuple for full mass range; false -- only for on-shell region
// fillLHEwgt : if false then don't fill hists with alternative LHE weights for both full range and on-shell region
// tempNA  : used to set weight to N/A value. Set to float -99999.
// tempNAv : used to set variable values to N/A. Set to float -99999.
//
// Additional cut:
// PassesZZjjBaselineSelection() : For EWK ZZjj analsysis, should be turned off (commented out) for other analyses
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ZZSelector::Init(TTree *tree)
{
  systematics_ = {
      {electronRecoEffUp, "CMS_RecoEff_eUp"},
      {electronRecoEffDown, "CMS_RecoEff_eDown"},
      {electronEfficiencyUp, "CMS_eff_eUp"},
      {electronEfficiencyDown, "CMS_eff_eDown"},
      {muonEfficiencyUp, "CMS_eff_mUp"},
      {muonEfficiencyDown, "CMS_eff_mDown"},
      {pileupUp, "CMS_pileupUp"},
      {pileupDown, "CMS_pileupDown"},
  };
  doSystematics_ = true; // false;//true;

  // This would be set true inside ZZBackground Selector
  // isNonPrompt_ = false;

  systHists_ = {"yield", "Mass", "MassFull", "nJets", "jetPt[1]", "jetPt[0]", "jetEta[0]", "jetEta[1]", "absjetEta[0]", "absjetEta[1]", "mjj", "dEtajj", "Mass0j", "Mass1j", "Mass2j", "Mass3j", "Mass34j", "Mass4j", "Mass0jFull", "Mass1jFull", "Mass2jFull", "Mass3jFull", "Mass34jFull", "Mass4jFull"};
  // hists1D_ = {
  //      "yield", "backgroundControlYield","nTruePU","nvtx","ZMass","Z1Mass","Z2Mass","ZZPt",
  //      "Z1Pt","Z2Pt","Z1Phi","Z2Phi","dPhiZ1Z2","ZPt","LepPt","LepEta","Lep12Pt","Lep12Eta",
  //      "Lep34Pt","Lep34Eta","Z1lep1_Eta","Z1lep1_Phi","Z1lep1_Pt","Z1lep1_PdgId","Z1lep2_Eta",
  //      "Z1lep2_Phi","Z1lep2_Pt","Z1lep2_PdgId","Z2lep1_Eta","Z2lep1_Phi","Z2lep1_Pt","Z2lep1_PdgId",
  //      "Z2lep2_Eta","Z2lep2_Phi","Z2lep2_Pt","Z2lep2_PdgId","Mass","nJets",
  // };

  hists1D_ = {"yield", "Mass", "MassFull", "nJets", "jetPt[1]", "jetPt[0]", "jetEta[0]", "jetEta[1]", "absjetEta[0]", "absjetEta[1]", "mjj", "dEtajj", "Mass0j", "Mass1j", "Mass2j", "Mass3j", "Mass34j", "Mass4j", "Mass0jFull", "Mass1jFull", "Mass2jFull", "Mass3jFull", "Mass34jFull", "Mass4jFull", "nJets_central"};

  jetTest2D_ = {}; // also defined in hists1D_ to pass checks in InitializeHistogramsFromConfig()
  jethists1D_ = {"Mass", "MassFull", "nJets", "jetPt[1]", "jetPt[0]", "jetEta[0]", "jetEta[1]", "absjetEta[0]", "absjetEta[1]", "mjj", "dEtajj", "Mass0j", "Mass1j", "Mass2j", "Mass3j", "Mass34j", "Mass4j", "Mass0jFull", "Mass1jFull", "Mass2jFull", "Mass3jFull", "Mass34jFull", "Mass4jFull"};

  weighthists1D_ = {"yield", "Mass", "MassFull", "nJets", "jetPt[1]", "jetPt[0]", "jetEta[0]", "jetEta[1]", "absjetEta[0]", "absjetEta[1]", "mjj", "dEtajj", "Mass0j", "Mass1j", "Mass2j", "Mass3j", "Mass34j", "Mass4j", "Mass0jFull", "Mass1jFull", "Mass2jFull", "Mass3jFull", "Mass34jFull", "Mass4jFull"};
  ZZSelectorBase::Init(tree);
  // fCutFormula = new TTreeFormula("CutFormula", fOption, fChain);
  // fCutFormula->SetQuickLoad(kTRUE);
  // if (!fCutFormula->GetNdim())
  //{
  //   delete fCutFormula;
  //   fCutFormula = 0;
  // }
  fCutFormula = 0; // turn off trigger cut test
}

void ZZSelector::SetBranchesUWVV()
{
  ZZSelectorBase::SetBranchesUWVV();
  if (isMC_)
  {
    weight_info_ = GetLheWeightInfo();
    if (weight_info_ > 0)
      fChain->SetBranchAddress("scaleWeights", &scaleWeights, &b_scaleWeights);
    if ((weight_info_ == 2 || weight_info_ == 3) && doSystematics_ && !isNonPrompt_)
      fChain->SetBranchAddress("pdfWeights", &pdfWeights, &b_pdfWeights);
  }
  fChain->SetBranchAddress("Mass", &Mass, &b_Mass);
  fChain->SetBranchAddress("Pt", &Pt, &b_Pt);
  fChain->SetBranchAddress("Eta", &Eta, &b_Eta);

  fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);

  fChain->SetBranchAddress("jetPUID", &jetPUID, &b_jetPUID);
  if (isMC_)
  {

    applyPUSF_ = false;
    if (applyPUSF_ || applyPUSFNtp_)
    {
      try
      { //! catching doesn't seem to work. Crash directly in either here or getentry()
        fChain->SetBranchAddress("isGenJetMatched", &isGenJetMatched, &b_isGenJetMatched);
        fChain->SetBranchAddress("jetPUSFmulfac", &jetPUSFmulfac, &b_jetPUSFmulfac);
      }
      catch (...)
      {
        applyPUSF_ = false;
        std::cout << "Cannot set address; PU SF is not applied for this sample" << std::endl;
      }
    }
    fChain->SetBranchAddress("jetPt_jesUp", &jetPt_jesUp, &b_jetPt_jesUp);
    fChain->SetBranchAddress("jetPt_jesDown", &jetPt_jesDown, &b_jetPt_jesDown);
    fChain->SetBranchAddress("jetPt_jerUp", &jetPt_jerUp, &b_jetPt_jerUp);
    fChain->SetBranchAddress("jetPt_jerDown", &jetPt_jerDown, &b_jetPt_jerDown);
  }

  fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);

  fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
  if (isMC_)
  {
    fChain->SetBranchAddress("jetEta_jesUp", &jetEta_jesUp, &b_jetEta_jesUp);
    fChain->SetBranchAddress("jetEta_jesDown", &jetEta_jesDown, &b_jetEta_jesDown);
    fChain->SetBranchAddress("jetEta_jerUp", &jetEta_jerUp, &b_jetEta_jerUp);
    fChain->SetBranchAddress("jetEta_jerDown", &jetEta_jerDown, &b_jetEta_jerDown);
  }

  fChain->SetBranchAddress("mjj", &mjj, &b_mjj);
  if (isMC_)
  {
    fChain->SetBranchAddress("mjj_jesUp", &mjj_jesUp, &b_mjj_jesUp);
    fChain->SetBranchAddress("mjj_jesDown", &mjj_jesDown, &b_mjj_jesDown);
    fChain->SetBranchAddress("mjj_jerUp", &mjj_jerUp, &b_mjj_jerUp);
    fChain->SetBranchAddress("mjj_jerDown", &mjj_jerDown, &b_mjj_jerDown);
  }

  fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
  if (isMC_)
  {
    fChain->SetBranchAddress("nJets_jesUp", &nJets_jesUp, &b_nJets_jesUp);
    fChain->SetBranchAddress("nJets_jesDown", &nJets_jesDown, &b_nJets_jesDown);
    fChain->SetBranchAddress("nJets_jerUp", &nJets_jerUp, &b_nJets_jerUp);
    fChain->SetBranchAddress("nJets_jerDown", &nJets_jerDown, &b_nJets_jerDown);
  }
}

unsigned int ZZSelector::GetLheWeightInfo()
{
  std::vector<std::string> noLheWeights = {
      "ggZZ2e2mu", "ggZZ4e", "ggZZ4m", "ggZZ4t", "ggZZ2e2tau", "ggZZ2mu2tau", "zz4l-sherpa", "ZZJJTo2e2mu-EWK-phantom", "ZZJJTo4e-EWK-phantom", "ZZJJTo4mu-EWK-phantom"};
  std::vector<std::string> scaleAndPdfWeights = {
      "wz3lnu-powheg", "wz3lnu-mg5amcnlo",
      "ZZZ", "WZZ", "WWZ",
      "zz4l-powheg", "zz4l-amcatnlo",
      "ZZJJTo4L-EWK"};
  std::vector<std::string> allLheWeights = {
      //"wzjj-aqgcft", "wzjj-aqgcfm", "wzjj-aqgcfs",
      //"wz-atgc_pt0-200", "wz-atgc_pt200-300",
      //"wz-atgc_pt300"
  };

  if ((std::find(noLheWeights.begin(), noLheWeights.end(), name_) != noLheWeights.end()) || (isaTGC_))
    return 0;
  if (std::find(scaleAndPdfWeights.begin(), scaleAndPdfWeights.end(), name_) != scaleAndPdfWeights.end())
    return 2;
  if (std::find(allLheWeights.begin(), allLheWeights.end(), name_) != allLheWeights.end())
    return 3;

  if (isUL_L1check)
  {
    return 0;
  }
  return 1;
}
void ZZSelector::LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation)
{

  ZZSelectorBase::LoadBranchesUWVV(entry, variation);

  if (fCutFormula && fCutFormula->EvalInstance() > 0.)
  {
    passCurrentTrig = true;
  }
  else
  {
    passCurrentTrig = false;
  }
  // b_MtToMET->GetEntry(entry);
  // b_l1Pt->GetEntry(entry);
  // b_l2Pt->GetEntry(entry);
  // b_l3Pt->GetEntry(entry);
  b_Mass->GetEntry(entry);
  b_Pt->GetEntry(entry);
  b_Eta->GetEntry(entry);
  // b_jetPt->GetEntry(entry);
  b_jetPhi->GetEntry(entry);
  // b_jetEta->GetEntry(entry);

  b_mjj->GetEntry(entry);
  if (isMC_)
  {
    b_mjj_jesUp->GetEntry(entry);
    b_mjj_jesDown->GetEntry(entry);
    b_mjj_jerUp->GetEntry(entry);
    b_mjj_jerDown->GetEntry(entry);
  }

  b_jetPt->GetEntry(entry);
  b_jetPUID->GetEntry(entry);
  if (isMC_)
  {
    if (applyPUSF_ || applyPUSFNtp_)
    {
      try
      { //! catching doesn't seem to work. Crash directly somewhere for sample without isGenJetMatched.
        b_isGenJetMatched->GetEntry(entry);
        b_jetPUSFmulfac->GetEntry(entry);
      }
      catch (...)
      {
        applyPUSF_ = false;
        std::cout << "PU SF is not applied for this sample" << std::endl;
      }
    }
    b_jetPt_jesUp->GetEntry(entry);
    b_jetPt_jesDown->GetEntry(entry);
    b_jetPt_jerUp->GetEntry(entry);
    b_jetPt_jerDown->GetEntry(entry);
  }

  b_jetEta->GetEntry(entry);
  if (isMC_)
  {
    b_jetEta_jesUp->GetEntry(entry);
    b_jetEta_jesDown->GetEntry(entry);
    b_jetEta_jerUp->GetEntry(entry);
    b_jetEta_jerDown->GetEntry(entry);
  }

  // erase jet with pt<50GeV for testing
  bool raisejPt = false;
  if (raisejPt)
  {
    if (jetPt->size() == jetEta->size())
    {
      auto jetit = jetPt->begin();
      auto jetait = jetEta->begin();
      while (jetit != jetPt->end())
      {
        if (*jetit < 50)
        {
          jetit = jetPt->erase(jetit);
          jetait = jetEta->erase(jetait);
        }
        else
        {
          ++jetit;
          ++jetait;
        }
      }
    }
    else
    {
      std::cout << "Something Wrong jetPt vs jetEta size" << jetPt->size() << " " << jetEta->size() << std::endl;
    }
  } // raise jet pt

  b_nJets->GetEntry(entry);
  if (isMC_)
  {
    b_nJets_jesUp->GetEntry(entry);
    b_nJets_jesDown->GetEntry(entry);
    b_nJets_jerUp->GetEntry(entry);
    b_nJets_jerDown->GetEntry(entry);
  }

  // std::cout<<"channel in LoadBranches function: "<<channel_<<std::endl;
  if (channel_ == eemm || channel_ == mmee)
  {
    // if(TightZZLeptons()){//i don't think this condition is needed even though it might save time but it messes up sf application for CRs in eemm,mmee states
    if (TightZZLeptons())
    {
      SetVariables(entry);
    }
  }

  if (isMC_)
  {
    ApplyScaleFactors();
  }

  if (isMC_)
  {
    bool applyjetMatch = false;
    if (applyjetMatch)
    {
      // apply gen jet matching for remaining processing if MC. For data, should apply PU id at ntuplization step but can also apply here
      if (jetPt->size() == jetEta->size() && jetPt->size() == jetPUID->size() && jetPUID->size() == isGenJetMatched->size())
      {
        auto jetit = jetPt->begin();
        auto jetait = jetEta->begin();
        auto jPUIDit = jetPUID->begin();
        auto jmatchit = isGenJetMatched->begin();

        while (jetit != jetPt->end())
        {

          if (*jmatchit < 1)
          {
            jetit = jetPt->erase(jetit);
            jetait = jetEta->erase(jetait);
            jPUIDit = jetPUID->erase(jPUIDit);
            jmatchit = isGenJetMatched->erase(jmatchit);
            // std::cout<<"A jet is removed by PU id"<<std::endl;
          }
          else
          {
            ++jetit;
            ++jetait;
            ++jPUIDit;
            ++jmatchit;
          }
        }
      }
      else
      {
        std::cout << "Something Wrong jetPt vs jetEta, jetPUID, isGenJetMatched size" << jetPt->size() << " " << jetEta->size() << " " << jetPUID->size() << std::endl; //<< " " << isGenJetMatched->size() << std::endl;
      }
    }
  }

  bool applyjetPUID = false;
  if (applyjetPUID)
  {
    // apply jet PU id for remaining processing if MC. For data, should apply PU id at ntuplization step but can also apply here
    if (jetPt->size() == jetEta->size() && jetPt->size() == jetPUID->size()) //&& jetPUID->size() == isGenJetMatched->size())
    {
      auto jetit = jetPt->begin();
      auto jetait = jetEta->begin();
      auto jPUIDit = jetPUID->begin();

      while (jetit != jetPt->end())
      {

        if (*jetit < 50 && *jPUIDit < 7)
        {
          jetit = jetPt->erase(jetit);
          jetait = jetEta->erase(jetait);
          jPUIDit = jetPUID->erase(jPUIDit);
          // std::cout<<"A jet is removed by PU id"<<std::endl;
        }
        else
        {
          ++jetit;
          ++jetait;
          ++jPUIDit;
        }
      }
    }
    else
    {
      std::cout << "Something Wrong jetPt vs jetEta, jetPUID, isGenJetMatched size" << jetPt->size() << " " << jetEta->size() << " " << jetPUID->size() << std::endl; //<< " " << isGenJetMatched->size() << std::endl;
    }
  }

  if (variation.first == Central)
  {
    if (isMC_ && doSystematics_ && !isNonPrompt_)
    {
      if (isMC_ && weight_info_ > 0)
      {
        b_scaleWeights->GetEntry(entry);
        lheWeights = *scaleWeights;
        if (weight_info_ == 2)
        {
          b_pdfWeights->GetEntry(entry);
          // Only keep NNPDF weights
          lheWeights.insert(lheWeights.end(), pdfWeights->begin(),
                            pdfWeights->begin() + std::min(static_cast<size_t>(103), pdfWeights->size()));
        }
        else if (weight_info_ == 3)
        {
          b_pdfWeights->GetEntry(entry);
          lheWeights.insert(lheWeights.end(), pdfWeights->begin(), pdfWeights->end());
        }
      }
    }
  }
  else if (isMC_)
  {
    // Systematic uncertainties and creating shiftUp and shiftDown histograms
    // Starting with lepton Efficiencies
    if (variation.first == electronEfficiencyUp || variation.first == electronEfficiencyDown || variation.first == electronRecoEffUp || variation.first == electronRecoEffDown ||
        variation.first == muonEfficiencyUp || variation.first == muonEfficiencyDown)
    {
      ShiftEfficiencies(variation.first);
    }
    else if (variation.first == pileupUp)
    {
      weight *= pileupSF_->Evaluate1D(nTruePU, ScaleFactor::ShiftUp) / pileupSF_->Evaluate1D(nTruePU);
    }
    else if (variation.first == pileupDown)
    {
      weight *= pileupSF_->Evaluate1D(nTruePU, ScaleFactor::ShiftDown) / pileupSF_->Evaluate1D(nTruePU);
    }
  }

  auto deltaPhiZZ = [](float phi1, float phi2)
  {
    float pi = TMath::Pi();
    float dphi = fabs(phi1 - phi2);
    if (dphi > pi)
      dphi = 2.0 * pi - dphi;
    return dphi;
  };
  auto deltaEtajj = [](std::vector<float> *jEta)
  {
    if (jEta->size() < 2)
      return -1.;
    double etaDiff = jEta->at(0) - jEta->at(1);
    return std::abs(etaDiff);
  };

  auto deltaRZZ = [](float eta1, float eta2, float dPhi)
  {
    float dEta = eta1 - eta2;
    return std::sqrt(dPhi * dPhi + dEta * dEta);
  };

  dEtajj = deltaEtajj(jetEta);
  dPhiZZ = deltaPhiZZ(Z1Phi, Z2Phi);
  dRZZ = deltaRZZ(Z1Eta, Z2Eta, dPhiZZ);
}

void ZZSelector::ApplyScaleFactors()
{

  if (applyPUSFNtp_)
  {

    weight *= jetPUSFmulfac;
  }
  // bool applyPUSF = true;
  if (applyPUSF_)
  {

    for (std::size_t ind = 0; ind < jetPt->size(); ind++)
    {

      if (isGenJetMatched->at(ind) < 1)
      { // only used real jets for applying SF
        continue;
      }

      if (jetPt->at(ind) < 50)
      {
        auto jetbin = jetPUSF_->FindBin(jetPt->at(ind), jetEta->at(ind));
        float jetPUSF = (float)jetPUSF_->GetBinContent(jetbin);
        auto jeffbin = jetPUeff_->FindBin(jetPt->at(ind), jetEta->at(ind));
        float jeffPU = (float)jetPUeff_->GetBinContent(jeffbin);
        float mulfac = 1.;

        if (jetPUID->at(ind) < 7)
        { // doesn't pass PU id
          mulfac = (1. - jetPUSF * jeffPU) / (1. - jeffPU);
        }
        else
        {
          mulfac = jetPUSF;
        }
        weight *= mulfac;
      }
    }
  }

  // In order to get around the Overflow issue, set the Pt, not ideal.
  // std::cout<<"weight before SF: "<<weight<<std::endl;
  if (channel_ == eeee)
  {
    float pt_e1 = l1Pt < EleSF_MAX_PT_ ? l1Pt : EleSF_MAX_PT_ - 0.01;
    float pt_e2 = l2Pt < EleSF_MAX_PT_ ? l2Pt : EleSF_MAX_PT_ - 0.01;
    float pt_e3 = l3Pt < EleSF_MAX_PT_ ? l3Pt : EleSF_MAX_PT_ - 0.01;
    float pt_e4 = l4Pt < EleSF_MAX_PT_ ? l4Pt : EleSF_MAX_PT_ - 0.01;
    if (eIdSF_ != nullptr and eGapIdSF_ != nullptr)
    {
      // In order to get around the Overflow issue, set the Pt, not ideal.
      if (l1IsGap)
      {
        weight *= eGapIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
        // std::cout<<"weight for GapE: "<<weight<<std::endl;
      }
      else
      {
        weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
      }
      if (l2IsGap)
      {
        weight *= eGapIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
      }
      else
      {
        weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
      }
      if (l3IsGap)
      {
        weight *= eGapIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
      }
      else
      {
        weight *= eIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
      }
      if (l4IsGap)
      {
        weight *= eGapIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
      }
      else
      {
        weight *= eIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
      }
    }
    if (eRecoSF_ != nullptr and eLowRecoSF_ != nullptr)
    {
      // Applying Electron Reco SFs
      if (l1Pt < 20)
      {
        weight *= eLowRecoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
      }
      else
      {
        weight *= eRecoSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
      }
      if (l2Pt < 20)
      {
        weight *= eLowRecoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
      }
      else
      {
        weight *= eRecoSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
      }
      if (l3Pt < 20)
      {
        weight *= eLowRecoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
      }
      else
      {
        weight *= eRecoSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
      }
      if (l4Pt < 20)
      {
        weight *= eLowRecoSF_->Evaluate2D(std::abs(l4Eta), l4Pt);
      }
      else
      {
        weight *= eRecoSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
      }
    }
  }
  else if (channel_ == eemm)
  {
    float pt_e1 = l1Pt < EleSF_MAX_PT_ ? l1Pt : EleSF_MAX_PT_ - 0.01;
    float pt_e2 = l2Pt < EleSF_MAX_PT_ ? l2Pt : EleSF_MAX_PT_ - 0.01;
    float pt_m3 = l3Pt < MuSF_MAX_PT_ ? l3Pt : MuSF_MAX_PT_ - 0.01;
    float pt_m4 = l4Pt < MuSF_MAX_PT_ ? l4Pt : MuSF_MAX_PT_ - 0.01;
    if (eIdSF_ != nullptr and eGapIdSF_ != nullptr)
    {
      // In order to get around the Overflow issue, set the Pt, not ideal.
      if (l1IsGap)
      {
        weight *= eGapIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
      }
      else
      {
        weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
      }
      if (l2IsGap)
      {
        weight *= eGapIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
      }
      else
      {
        weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
      }
    }
    // Applying Electron Reco SFs
    if (eRecoSF_ != nullptr and eLowRecoSF_ != nullptr)
    {
      if (l1Pt < 20)
      {
        weight *= eLowRecoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
      }
      else
      {
        weight *= eRecoSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
      }
      if (l2Pt < 20)
      {
        weight *= eLowRecoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
      }
      else
      {
        weight *= eRecoSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
      }
    }
    if (mIdSF_ != nullptr)
    {
      weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3);
      weight *= mIdSF_->Evaluate2D(std::abs(l4Eta), pt_m4);
    }
  }
  else if (channel_ == mmee)
  {
    float pt_m1 = l1Pt < MuSF_MAX_PT_ ? l1Pt : MuSF_MAX_PT_ - 0.01;
    float pt_m2 = l2Pt < MuSF_MAX_PT_ ? l2Pt : MuSF_MAX_PT_ - 0.01;
    float pt_e3 = l3Pt < EleSF_MAX_PT_ ? l3Pt : EleSF_MAX_PT_ - 0.01;
    float pt_e4 = l4Pt < EleSF_MAX_PT_ ? l4Pt : EleSF_MAX_PT_ - 0.01;
    if (mIdSF_ != nullptr)
    {
      weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1);
      weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), pt_m2);
    }
    if (eIdSF_ != nullptr and eGapIdSF_ != nullptr)
    {
      if (l3IsGap)
      {
        weight *= eGapIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
      }
      else
      {
        weight *= eIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
      }
      if (l4IsGap)
      {
        weight *= eGapIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
      }
      else
      {
        weight *= eIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
      }
    }
    // Applying Electron Reco SFs
    if (eRecoSF_ != nullptr and eLowRecoSF_ != nullptr)
    {
      if (l3Pt < 20)
      {
        weight *= eLowRecoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
      }
      else
      {
        weight *= eRecoSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
      }
      if (l4Pt < 20)
      {
        weight *= eLowRecoSF_->Evaluate2D(std::abs(l4Eta), l4Pt);
      }
      else
      {
        weight *= eRecoSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
      }
    }
  }
  else
  {
    float pt_m1 = l1Pt < MuSF_MAX_PT_ ? l1Pt : MuSF_MAX_PT_ - 0.01;
    float pt_m2 = l2Pt < MuSF_MAX_PT_ ? l2Pt : MuSF_MAX_PT_ - 0.01;
    float pt_m3 = l3Pt < MuSF_MAX_PT_ ? l3Pt : MuSF_MAX_PT_ - 0.01;
    float pt_m4 = l4Pt < MuSF_MAX_PT_ ? l4Pt : MuSF_MAX_PT_ - 0.01;
    if (mIdSF_ != nullptr)
    {
      weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1);
      weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), pt_m2);
      weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3);
      weight *= mIdSF_->Evaluate2D(std::abs(l4Eta), pt_m4);
    }
  }
  if (pileupSF_ != nullptr)
  {
    weight *= pileupSF_->Evaluate1D(nTruePU);
  }
}
// Similar to Kenneth's SetShiftedMasses function which i will need later as well
void ZZSelector::SetVariables(Long64_t entry)
{
  if (!(e1e2IsZ1(entry)))
  {
    // std::cout<<"e1e2IsZ1 is working"<<std::endl;
    //! Should switch ZEta and ZPhi etc. as well, but doesn't affect results for the purpose so far
    float tempMass = Z1mass;
    Z1mass = Z2mass;
    Z2mass = tempMass;
    float tempPt = Z1pt;
    Z1pt = Z2pt;
    Z2pt = tempPt;
    bool templ1IsTight = l1IsTight;
    l1IsTight = l3IsTight;
    l3IsTight = templ1IsTight;
    bool templ2IsTight = l2IsTight;
    l2IsTight = l4IsTight;
    l4IsTight = templ2IsTight;
    float templ1IsIso = l1IsIso;
    l1IsIso = l3IsIso;
    l3IsIso = templ1IsIso;
    float templ2IsIso = l2IsIso;
    l2IsIso = l4IsIso;
    l4IsIso = templ2IsIso;
    bool templ1IsGap = l1IsGap;
    l1IsGap = l3IsGap;
    l3IsGap = templ1IsGap;
    bool templ2IsGap = l2IsGap;
    l2IsGap = l4IsGap;
    l4IsGap = templ2IsGap;
    float templ1Pt = l1Pt;
    l1Pt = l3Pt;
    l3Pt = templ1Pt;
    float templ2Pt = l2Pt;
    l2Pt = l4Pt;
    l4Pt = templ2Pt;
    float templ1PVDZ = l1PVDZ;
    l1PVDZ = l3PVDZ;
    l3PVDZ = templ1PVDZ;
    float templ2PVDZ = l2PVDZ;
    l2PVDZ = l4PVDZ;
    l4PVDZ = templ2PVDZ;
    float templ1Eta = l1Eta;
    l1Eta = l3Eta;
    l3Eta = templ1Eta;
    float templ2Eta = l2Eta;
    l2Eta = l4Eta;
    l4Eta = templ2Eta;
    float templ1Phi = l1Phi;
    l1Phi = l3Phi;
    l3Phi = templ1Phi;
    float templ2Phi = l2Phi;
    l2Phi = l4Phi;
    l4Phi = templ2Phi;
    float templ1SIP3D = l1SIP3D;
    l1SIP3D = l3SIP3D;
    l3SIP3D = templ1SIP3D;
    float templ2SIP3D = l2SIP3D;
    l2SIP3D = l4SIP3D;
    l4SIP3D = templ2SIP3D;
    int templ1PdgId = l1PdgId;
    l1PdgId = l3PdgId;
    l3PdgId = templ1PdgId;
    int templ2PdgId = l2PdgId;
    l2PdgId = l4PdgId;
    l4PdgId = templ2PdgId;
    float templ1Mass = l1Mass;
    l1Mass = l3Mass;
    l3Mass = templ1Mass;
    float templ2Mass = l2Mass;
    l2Mass = l4Mass;
    l4Mass = templ2Mass;
  }
}
void ZZSelector::ShiftEfficiencies(Systematic variation)
{
  ScaleFactor::Variation shift = ScaleFactor::Variation::ShiftUp;
  if (variation == electronEfficiencyDown || variation == electronRecoEffDown || variation == muonEfficiencyDown)
    shift = ScaleFactor::Variation::ShiftDown;

  if (channel_ == eeee)
  {
    float pt_e1 = l1Pt < EleSF_MAX_PT_ ? l1Pt : EleSF_MAX_PT_ - 0.01;
    float pt_e2 = l2Pt < EleSF_MAX_PT_ ? l2Pt : EleSF_MAX_PT_ - 0.01;
    float pt_e3 = l3Pt < EleSF_MAX_PT_ ? l3Pt : EleSF_MAX_PT_ - 0.01;
    float pt_e4 = l4Pt < EleSF_MAX_PT_ ? l4Pt : EleSF_MAX_PT_ - 0.01;
    // std::cout<<eRecoSF_<<std::endl;
    if (variation == electronRecoEffUp || variation == electronRecoEffDown)
    {
      if (eRecoSF_ != nullptr and eLowRecoSF_ != nullptr)
      {
        // std::cout<<"Is it a null ptr"<<std::endl;
        // Applying Electron Reco SFs Up/Down for ElectronRecoSyst
        if (l1Pt < 20)
        {
          weight *= eLowRecoSF_->Evaluate2D(std::abs(l1Eta), l1Pt, shift) / eLowRecoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
        }
        else
        {
          weight *= eRecoSF_->Evaluate2D(std::abs(l1Eta), pt_e1, shift) / eRecoSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
        }
        if (l2Pt < 20)
        {
          weight *= eLowRecoSF_->Evaluate2D(std::abs(l2Eta), l2Pt, shift) / eLowRecoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        }
        else
        {
          weight *= eRecoSF_->Evaluate2D(std::abs(l2Eta), pt_e2, shift) / eRecoSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
        }
        if (l3Pt < 20)
        {
          weight *= eLowRecoSF_->Evaluate2D(std::abs(l3Eta), l3Pt, shift) / eLowRecoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
        else
        {
          weight *= eRecoSF_->Evaluate2D(std::abs(l3Eta), pt_e3, shift) / eRecoSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
        }
        if (l4Pt < 20)
        {
          weight *= eLowRecoSF_->Evaluate2D(std::abs(l4Eta), l4Pt, shift) / eLowRecoSF_->Evaluate2D(std::abs(l4Eta), l4Pt);
        }
        else
        {
          weight *= eRecoSF_->Evaluate2D(std::abs(l4Eta), pt_e4, shift) / eRecoSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
        }
      }
    }
    // Applying Electron ID SFs Up/Down for ElectronIDEffSyst
    else if (variation == electronEfficiencyUp || variation == electronEfficiencyDown)
    {
      if (eIdSF_ != nullptr and eGapIdSF_ != nullptr)
      {
        if (l1IsGap)
        {
          weight *= eGapIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1, shift) / eGapIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
        }
        else
        {
          weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1, shift) / eIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
        }
        if (l2IsGap)
        {
          weight *= eGapIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2, shift) / eGapIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
        }
        else
        {
          weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2, shift) / eIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
        }
        if (l3IsGap)
        {
          weight *= eGapIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3, shift) / eGapIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
        }
        else
        {
          weight *= eIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3, shift) / eIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
        }
        if (l4IsGap)
        {
          weight *= eGapIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4, shift) / eGapIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
        }
        else
        {
          weight *= eIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4, shift) / eIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
        }
      }
    }
  } // channel eeee
  else if (channel_ == eemm)
  {
    // In order to get around the Overflow issue, set the Pt, not ideal.
    float pt_e1 = l1Pt < EleSF_MAX_PT_ ? l1Pt : EleSF_MAX_PT_ - 0.01;
    float pt_e2 = l2Pt < EleSF_MAX_PT_ ? l2Pt : EleSF_MAX_PT_ - 0.01;
    float pt_m3 = l3Pt < MuSF_MAX_PT_ ? l3Pt : MuSF_MAX_PT_ - 0.01;
    float pt_m4 = l4Pt < MuSF_MAX_PT_ ? l4Pt : MuSF_MAX_PT_ - 0.01;
    if (variation == electronRecoEffUp || variation == electronRecoEffDown)
    {
      // Applying Electron Reco SFs Up/Down for ElectronRecoEffSyst
      if (l1Pt < 20)
      {
        weight *= eLowRecoSF_->Evaluate2D(std::abs(l1Eta), l1Pt, shift) / eLowRecoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
      }
      else
      {
        weight *= eRecoSF_->Evaluate2D(std::abs(l1Eta), pt_e1, shift) / eRecoSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
      }
      if (l2Pt < 20)
      {
        weight *= eLowRecoSF_->Evaluate2D(std::abs(l2Eta), l2Pt, shift) / eLowRecoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
      }
      else
      {
        weight *= eRecoSF_->Evaluate2D(std::abs(l2Eta), pt_e2, shift) / eRecoSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
      }
    }
    // Applying Electron ID SFs Up/Down for ElectronIDEffSyst
    else if (variation == electronEfficiencyUp || variation == electronEfficiencyDown)
    {
      if (l1IsGap)
      {
        weight *= eGapIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1, shift) / eGapIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
      }
      else
      {
        weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1, shift) / eIdSF_->Evaluate2D(std::abs(l1Eta), pt_e1);
      }
      if (l2IsGap)
      {
        weight *= eGapIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2, shift) / eGapIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
      }
      else
      {
        weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2, shift) / eIdSF_->Evaluate2D(std::abs(l2Eta), pt_e2);
      }
    }
    else if (variation == muonEfficiencyUp || variation == muonEfficiencyDown)
    {
      weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3, shift) / mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3);
      weight *= mIdSF_->Evaluate2D(std::abs(l4Eta), pt_m4, shift) / mIdSF_->Evaluate2D(std::abs(l4Eta), pt_m4);
    }
  }
  else if (channel_ == mmee)
  {
    float pt_m1 = l1Pt < MuSF_MAX_PT_ ? l1Pt : MuSF_MAX_PT_ - 0.01;
    float pt_m2 = l2Pt < MuSF_MAX_PT_ ? l2Pt : MuSF_MAX_PT_ - 0.01;
    float pt_e3 = l3Pt < EleSF_MAX_PT_ ? l3Pt : EleSF_MAX_PT_ - 0.01;
    float pt_e4 = l4Pt < EleSF_MAX_PT_ ? l4Pt : EleSF_MAX_PT_ - 0.01;
    if (variation == muonEfficiencyUp || variation == muonEfficiencyDown)
    {
      weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1, shift) / mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1);
      weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), pt_m2, shift) / mIdSF_->Evaluate2D(std::abs(l2Eta), pt_m2);
    }
    else if (variation == electronRecoEffUp || variation == electronRecoEffDown)
    {
      // Applying Electron Reco SFs Up/Down for ElectronRecoEffSyst
      if (l3Pt < 20)
      {
        weight *= eLowRecoSF_->Evaluate2D(std::abs(l3Eta), l3Pt, shift) / eLowRecoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
      }
      else
      {
        weight *= eRecoSF_->Evaluate2D(std::abs(l3Eta), pt_e3, shift) / eRecoSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
      }
      if (l4Pt < 20)
      {
        weight *= eLowRecoSF_->Evaluate2D(std::abs(l4Eta), l4Pt, shift) / eLowRecoSF_->Evaluate2D(std::abs(l4Eta), l4Pt);
      }
      else
      {
        weight *= eRecoSF_->Evaluate2D(std::abs(l4Eta), pt_e4, shift) / eRecoSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
      }
    }
    // Applying Electron ID SFs Up/Down for ElectronIDEffSyst
    else if (variation == electronEfficiencyUp || variation == electronEfficiencyDown)
    {
      if (l3IsGap)
      {
        weight *= eGapIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3, shift) / eGapIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
      }
      else
      {
        weight *= eIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3, shift) / eIdSF_->Evaluate2D(std::abs(l3Eta), pt_e3);
      }
      if (l4IsGap)
      {
        weight *= eGapIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4, shift) / eGapIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
      }
      else
      {
        weight *= eIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4, shift) / eIdSF_->Evaluate2D(std::abs(l4Eta), pt_e4);
      }
    }
  }
  else if (channel_ == mmmm && (variation == muonEfficiencyUp || variation == muonEfficiencyDown))
  {
    float pt_m1 = l1Pt < MuSF_MAX_PT_ ? l1Pt : MuSF_MAX_PT_ - 0.01;
    float pt_m2 = l2Pt < MuSF_MAX_PT_ ? l2Pt : MuSF_MAX_PT_ - 0.01;
    float pt_m3 = l3Pt < MuSF_MAX_PT_ ? l3Pt : MuSF_MAX_PT_ - 0.01;
    float pt_m4 = l4Pt < MuSF_MAX_PT_ ? l4Pt : MuSF_MAX_PT_ - 0.01;
    weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1, shift) / mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1);
    // std::cout<<"l1 Shift: "<< mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1, shift)<<std::endl;
    // std::cout<<"l1 sf: "<< mIdSF_->Evaluate2D(std::abs(l1Eta), pt_m1)<<std::endl;
    weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), pt_m2, shift) / mIdSF_->Evaluate2D(std::abs(l2Eta), pt_m2);
    weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3, shift) / mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3);
    // std::cout<<"l3 Shift: "<< mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3, shift)<<std::endl;
    // std::cout<<"l3 sf: "<< mIdSF_->Evaluate2D(std::abs(l3Eta), pt_m3)<<std::endl;
    weight *= mIdSF_->Evaluate2D(std::abs(l4Eta), pt_m4, shift) / mIdSF_->Evaluate2D(std::abs(l4Eta), pt_m4);
  }
}

bool ZZSelector::PassesZZjjSelection()
{
  if ((jetPt->size() != jetEta->size() || jetPt->size() < 2) || (mjj < 100))
    return false;
  else
    return true;
}

bool ZZSelector::PassesZZjjBaselineSelection() // Define a separate selection in updated codes
{
  if ((jetPt->size() != jetEta->size() || jetPt->size() < 2) || (mjj <= 100) || (Mass <= 180.0))
    return false;
  else
    return true;
}

bool ZZSelector::Passes4eExtraCut()
{

  float lpt_arraySort[] = {l1Pt, l2Pt, l3Pt, l4Pt};
  std::sort(lpt_arraySort, lpt_arraySort + 4, std::greater<float>());

  if (channel_ == eeee)
  {

    if (lpt_arraySort[0] > 23 && lpt_arraySort[1] > 12)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  else
  {
    return true;
  }
}

bool ZZSelector::Passes2e2mExtraCut(Long64_t entry)
{

  float lpt1Tmp = 0.;
  float lpt2Tmp = 0.;

  if (channel_ == eemm)
  {
    float lpt_arraySort[] = {l1Pt, l2Pt}; // By default l1pt,l2pt set to e1pt,e2pt for eemm and mmee
    std::sort(lpt_arraySort, lpt_arraySort + 2, std::greater<float>());
    lpt1Tmp = lpt_arraySort[0];
    lpt2Tmp = lpt_arraySort[1];
  }
  else if (channel_ == mmee)
  {
    float lpt_arraySort[] = {l3Pt, l4Pt}; // pt switched by SetVariables function
    std::sort(lpt_arraySort, lpt_arraySort + 2, std::greater<float>());
    lpt1Tmp = lpt_arraySort[0];
    lpt2Tmp = lpt_arraySort[1];
  }

  if (channel_ == eemm || channel_ == mmee)
  {

    if (lpt1Tmp > 23 && lpt2Tmp > 12)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  else
  {
    return true;
  }
}

bool ZZSelector::PassesZZSelection(bool nonPrompt)
{
  // This nonPrompt boolean is for ZZBackgroundSelector
  // When running ZZBackgroundSelector, FillHistograms should run just with ZZSelection, we cannot require TightZZLeptons by definition
  if (nonPrompt)
  {
    // Because we are using Cut Based Moriond ID for Muons, we need SIP cut for 2016,2017
    if (year_ == yr2016 || year_ == yr2017)
    {
      if (ZZSelection() && HZZSIPSelection())
        return true;
      else
        return false;
    }
    else
    {
      if (ZZSelection())
        return true;
      else
        return false;
    }
  }
  else
  {
    // std::cout<<"nonPrompt inside function: "<<nonPrompt<<std::endl;
    if (year_ == yr2016 || year_ == yr2017)
    {
      if (ZZSelection() && TightZZLeptons() && HZZSIPSelection())
        return true;
      else
        return false;
    }
    else
    {
      if (ZZSelection() && TightZZLeptons())
        return true;
      else
        return false;
    }
  }
}

bool ZZSelector::PassesZZSelectionLoose(bool nonPrompt)
{
  // This nonPrompt boolean is for ZZBackgroundSelector
  // When running ZZBackgroundSelector, FillHistograms should run just with ZZSelection, we cannot require TightZZLeptons by definition
  if (nonPrompt)
  {
    // Because we are using Cut Based Moriond ID for Muons, we need SIP cut for 2016,2017
    if (year_ == yr2016 || year_ == yr2017)
    {
      if (HZZSIPSelection())
        return true;
      else
        return false;
    }
    else
    {
      if (true)
        return true;
      else
        return false;
    }
  }
  else
  {
    // std::cout<<"nonPrompt inside function: "<<nonPrompt<<std::endl;
    if (year_ == yr2016 || year_ == yr2017)
    {
      if (TightZZLeptons() && HZZSIPSelection())
        return true;
      else
        return false;
    }
    else
    {
      if (TightZZLeptons())
        return true;
      else
        return false;
    }
  }
}

bool ZZSelector::PassesHZZSelection(bool nonPrompt)
{
  if (nonPrompt)
  {
    if (ZSelection())
      return true;
    else
      return false;
  }
  else
  {
    if (ZSelection() && TightZZLeptons())
      return true;
    else
      return false;
  }
}
bool ZZSelector::TightZZLeptons()
{
  if (tightZ1Leptons() && tightZ2Leptons())
    return true;
  else
    return false;
}
bool ZZSelector::ZZSelection()
{
  if ((Z1mass > 60.0 && Z1mass < 120.0) && (Z2mass > 60.0 && Z2mass < 120.0))
  {
    return true;
  }
  else
  {
    return false;
  }
}

// We already require 4 < Z1,Z2 < 120  in the "Loose Skim"
bool ZZSelector::ZSelection()
{
  if (Z1mass > 40.0 && Z2mass > 12.0)
    return true;
  else
    return false;
}
bool ZZSelector::Z4lSelection()
{
  if (Mass > 80.0 && Mass < 100.0)
    return true;
  else
    return false;
}
// This is no longer needed as I apply SIP < 4.0 on electrons while skimming and we don't need SIP cut on muons with latest HZZID
bool ZZSelector::HZZSIPSelection()
{
  if ((l1SIP3D < 4.0 && l2SIP3D < 4.0 && l3SIP3D < 4.0 && l4SIP3D < 4.0))
    return true;
  else
    return false;
}
bool ZZSelector::HZZLowMass()
{
  if (Mass > 130.0 && Mass < 170.0)
    return true;
  else
    return false;
}
bool ZZSelector::HZZMediumMass()
{
  if (Mass > 138.0 && Mass < 300.0)
    return true;
  else
    return false;
}
bool ZZSelector::TestMuons()
{
  if ((Z1mass > 82.0 && Z1mass < 102.0) && (Z2mass < 40.0))
    return true;
  else
    return false;
}

void ZZSelector::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation)
{

  if (!PassesZZSelectionLoose(isNonPrompt_))
  {
    return;
  }

  if (!Passes4eExtraCut()) // Apply extra 23/12 GeV cut to 4e channel
  {
    return;
  }

  //! Specific to EWK ZZjj analysis and should be turned off otherwise
  if (!PassesZZjjBaselineSelection())
  {
    return;
  }

  //! Fill variables for full mass range and on-shell ZZ ntuple filling to avoid accidental modification etc.

  float tempNAv = -99999.;
  int nJets_tmp = jetPt->size();
  float l1pt_tmp = l1Pt;
  float l2pt_tmp = l2Pt;
  float l3pt_tmp = l3Pt;
  float l4pt_tmp = l4Pt;

  float l1eta_tmp = l1Eta;
  float l2eta_tmp = l2Eta;
  float l3eta_tmp = l3Eta;
  float l4eta_tmp = l4Eta;
  
  float jpt0_tmp;
  float jeta0_tmp;

  if (nJets_tmp >= 1)
  {
    jpt0_tmp = jetPt->at(0);
    jeta0_tmp = jetEta->at(0);
  }
  else
  {
    jpt0_tmp = tempNAv;
    jeta0_tmp = tempNAv;
  }

  float jpt1_tmp;
  float jeta1_tmp;
  float jetaproduct_tmp; 
  float zep1_tmp;
  float zep2_tmp;

  //! One value for each shift (JES up, JER dn etc.)
  std::vector<float> jetPt0_syst_tmp = {tempNAv, tempNAv, tempNAv, tempNAv};
  std::vector<float> jetEta0_syst_tmp = {tempNAv, tempNAv, tempNAv, tempNAv};
  std::vector<float> jetPt1_syst_tmp = {tempNAv, tempNAv, tempNAv, tempNAv};
  std::vector<float> jetEta1_syst_tmp = {tempNAv, tempNAv, tempNAv, tempNAv};
  std::vector<float> mjj_syst_tmp = {tempNAv, tempNAv, tempNAv, tempNAv};
  std::vector<float> dEtajj_syst_tmp = {tempNAv, tempNAv, tempNAv, tempNAv};
  std::vector<float> jetaproduct_syst_tmp = {tempNAv, tempNAv, tempNAv, tempNAv};
  std::vector<float> zep1_syst_tmp = {tempNAv, tempNAv, tempNAv, tempNAv};
  std::vector<float> zep2_syst_tmp = {tempNAv, tempNAv, tempNAv, tempNAv};

  //! mjj and dEtajj has default values already (but dEtajj=-1 is small although unphysical... unlike -9999), but need to assign values for temporary jetPt[1] and jetEta[1]
  //! => mjj and dEtajj NaN values are set here as well for consistency; dPhiZZ is always defined with the two Z's present.
  if (nJets_tmp >= 2)
  {
    jpt1_tmp = jetPt->at(1);
    jeta1_tmp = jetEta->at(1);
    jetaproduct_tmp = jetEta->at(0)*jetEta->at(1);
    zep1_tmp = Z1Eta - 0.5*(jetEta->at(0)+jetEta->at(1));
    zep2_tmp = Z2Eta - 0.5*(jetEta->at(0)+jetEta->at(1));
  }
  else
  {
    jpt1_tmp = tempNAv;
    jeta1_tmp = tempNAv;
    dEtajj = tempNAv;
    mjj = tempNAv;
    jetaproduct_tmp = tempNAv;
    zep1_tmp = tempNAv;
    zep2_tmp = tempNAv;
  }

  std::vector<std::vector<float> *> vjetEta = {jetEta_jesUp, jetEta_jesDown, jetEta_jerUp, jetEta_jerDown};
  std::vector<std::vector<float> *> vjetPt = {jetPt_jesUp, jetPt_jesDown, jetPt_jerUp, jetPt_jerDown};
  std::vector<unsigned int> vnJets = {nJets_jesUp, nJets_jesDown, nJets_jerUp, nJets_jerDown};
  std::vector<float> vmjj = {mjj_jesUp, mjj_jesDown, mjj_jerUp, mjj_jerDown};
  std::vector<std::string> jetSystNames = {"jesUp", "jesDown", "jerUp", "jerDown"};

  bool processFullRange = false;
  bool writeNtpFullRange = false;
  bool fillLHEwgt = false;

  //==============================================================================
  //   __       _ _             _  _   _
  //  / _|_   _| | |  _ __ ___ | || | | |  _ __ __ _ _ __   __ _  ___
  // | |_| | | | | | | '_ ` _ \| || |_| | | '__/ _` | '_ \ / _` |/ _ \
// |  _| |_| | | | | | | | | |__   _| | | | | (_| | | | | (_| |  __/
  // |_|  \__,_|_|_| |_| |_| |_|  |_| |_| |_|  \__,_|_| |_|\__, |\___|
  //                                                       |___/
  //
  //
  //===============================================================================

  if (processFullRange)
  {
    // int nJets_tmp2 = jetPt->size();
    if (writeNtp_ && writeNtpFullRange)
    {
      SafeSetBranch(ftntp_, getBranchName("weight", variation.second), &weight);
      SafeSetBranch(ftntp_, getBranchName("Mass", variation.second), &Mass);
      SafeSetBranch(ftntp_, getBranchName("nJets", variation.second), &nJets_tmp);
      SafeSetBranch(ftntp_, getBranchName("jetPt0", variation.second), &jpt0_tmp);
      SafeSetBranch(ftntp_, getBranchName("jetEta0", variation.second), &jeta0_tmp);

      SafeSetBranch(ftntp_, getBranchName("run", variation.second), &run);
      SafeSetBranch(ftntp_, getBranchName("lumi", variation.second), &lumi);
      SafeSetBranch(ftntp_, getBranchName("evt", variation.second), &evt);

      SafeSetBranch(ftntp_, getBranchName("jetPt1", variation.second), &jpt1_tmp);
      SafeSetBranch(ftntp_, getBranchName("jetEta1", variation.second), &jeta1_tmp);

      SafeSetBranch(ftntp_, getBranchName("mjj", variation.second), &mjj);
      SafeSetBranch(ftntp_, getBranchName("dEtajj", variation.second), &dEtajj);

      if (isMC_)
      {
        SafeSetBranch(ftntp_, getBranchName("genWeight", variation.second), &genWeight);
        SafeSetBranch(ftntp_, getBranchName("L1prefiringWeight", variation.second), &L1prefiringWeight);
        SafeSetBranch(ftntp_, getBranchName("L1prefiringWeightUp", variation.second), &L1prefiringWeightUp);
        SafeSetBranch(ftntp_, getBranchName("L1prefiringWeightDn", variation.second), &L1prefiringWeightDn);
      }
      ftntp_->Fill();
    }
    // if (!Passes2e2mExtraCut(entry)) //Apply extra 23/12 GeV cut to electrons in 2e2m channel
    //   {
    //  return;
    //   }

    // if (!passCurrentTrig)
    //{
    //   return;
    // }

    if ((variation.first == Central || (doaTGC_ && isaTGC_)) && isMC_)
    {
      //! Do jet systematics JES and JER for full m4l range
      if (isMC_) //! duplicated condition
      {
        for (size_t i = 0; i < vjetEta.size(); i++)
        { //! No actual syst for full m4l but just for consistency
          SafeHistFill(jethistMap1D_, getHistName("MassFull", variation.second), Mass, i, weight);

          if (vnJets[i] == 0)
          {
            SafeHistFill(jethistMap1D_, getHistName("Mass0jFull", variation.second), Mass, i, weight);
          }
          else if (vnJets[i] == 1)
          {
            SafeHistFill(jethistMap1D_, getHistName("Mass1jFull", variation.second), Mass, i, weight);
          }
          else if (vnJets[i] == 2)
          {
            SafeHistFill(jethistMap1D_, getHistName("Mass2jFull", variation.second), Mass, i, weight);
          }
          else
          {
            if (vnJets[i] == 3)
            {
              SafeHistFill(jethistMap1D_, getHistName("Mass3jFull", variation.second), Mass, i, weight);
            }

            if (vnJets[i] >= 3)
            {
              SafeHistFill(jethistMap1D_, getHistName("Mass34jFull", variation.second), Mass, i, weight);
            }

            if (vnJets[i] >= 4)
            {
              SafeHistFill(jethistMap1D_, getHistName("Mass4jFull", variation.second), Mass, i, weight);
            }
          }
        }
      }

      //! LHE reweighting, related to syst. unc. with PDF and QCD scales
      if (fillLHEwgt)
      {
        for (size_t i = 0; i < lheWeights.size(); i++) // expect 0 to 111 currently
        {
          SafeHistFill(weighthistMap1D_, getHistName("MassFull", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
          if (jetPt->size() == 0 && jetPt->size() == jetEta->size())
          {
            SafeHistFill(weighthistMap1D_, getHistName("Mass0jFull", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
          }
          else if (jetPt->size() == 1 && jetPt->size() == jetEta->size())
          {
            SafeHistFill(weighthistMap1D_, getHistName("Mass1jFull", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
          }
          else if (jetPt->size() == 2 && jetPt->size() == jetEta->size())
          {
            SafeHistFill(weighthistMap1D_, getHistName("Mass2jFull", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
          }
          else
          {

            if (jetPt->size() == 3 && jetPt->size() == jetEta->size())
            {
              SafeHistFill(weighthistMap1D_, getHistName("Mass3jFull", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
            }

            if (jetPt->size() >= 3 && jetPt->size() == jetEta->size())
            {
              SafeHistFill(weighthistMap1D_, getHistName("Mass34jFull", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
            }

            if (jetPt->size() >= 4 && jetPt->size() == jetEta->size())
            {
              SafeHistFill(weighthistMap1D_, getHistName("Mass4jFull", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
            }
          }
        } //! For-loop for filling LHE weights
      } //! If fillLHEwgt
    }

    //! Nominal full mass range histograms being filled
    SafeHistFill(histMap1D_, getHistName("MassFull", variation.second), Mass, weight);
    if (jetPt->size() == 0 && jetPt->size() == jetEta->size())
    {
      SafeHistFill(histMap1D_, getHistName("Mass0jFull", variation.second), Mass, weight);
    }
    else if (jetPt->size() == 1 && jetPt->size() == jetEta->size())
    {
      SafeHistFill(histMap1D_, getHistName("Mass1jFull", variation.second), Mass, weight);
    }
    else if (jetPt->size() == 2 && jetPt->size() == jetEta->size())
    {
      SafeHistFill(histMap1D_, getHistName("Mass2jFull", variation.second), Mass, weight);
    }
    else
    {
      if (jetPt->size() == 3 && jetPt->size() == jetEta->size())
      {
        SafeHistFill(histMap1D_, getHistName("Mass3jFull", variation.second), Mass, weight);
      }

      if (jetPt->size() >= 3 && jetPt->size() == jetEta->size())
      {
        SafeHistFill(histMap1D_, getHistName("Mass34jFull", variation.second), Mass, weight);
      }

      if (jetPt->size() >= 4 && jetPt->size() == jetEta->size())
      {
        SafeHistFill(histMap1D_, getHistName("Mass4jFull", variation.second), Mass, weight);
      }
    }

    /*  // sort lepton pt
      float lpt_array[] = {l1Pt, l2Pt, l3Pt, l4Pt};
      std::sort(lpt_array, lpt_array + 4, std::greater<float>());
      float l1PtTmp, l2PtTmp, l3PtTmp, l4PtTmp;
      l1PtTmp = lpt_array[0];
      l2PtTmp = lpt_array[1];
      l3PtTmp = lpt_array[2];
      l4PtTmp = lpt_array[3];

      float e1PtTmp = 0.;
      float e2PtTmp = 0.;
      if (channel_ == eemm)
      {
        float lpt_arraySort[] = {l1Pt, l2Pt}; // By default l1pt,l2pt set to e1pt,e2pt for eemm and mmee
        std::sort(lpt_arraySort, lpt_arraySort + 2, std::greater<float>());
        e1PtTmp = lpt_arraySort[0];
        e2PtTmp = lpt_arraySort[1];
      }
      else if (channel_ == mmee)
      {
        float lpt_arraySort[] = {l3Pt, l4Pt}; // pt switched by SetVariables function
        std::sort(lpt_arraySort, lpt_arraySort + 2, std::greater<float>());
        e1PtTmp = lpt_arraySort[0];
        e2PtTmp = lpt_arraySort[1];
      }

      //finish sorting lepton pt
    */

    if (80 < Mass && Mass < 110)
    {

      //    SafeHistFill(histMap1D_, getHistName("PassTriggerFull", variation.second), 1, weight);

      //    SafeHistFill(histMap1D_, getHistName("LepPtFull", variation.second), l1PtTmp, weight);
      //    SafeHistFill(histMap1D_, getHistName("LepPtFull", variation.second), l2PtTmp, weight);
      //    SafeHistFill(histMap1D_, getHistName("LepPtFull", variation.second), l3PtTmp, weight);
      //    SafeHistFill(histMap1D_, getHistName("LepPtFull", variation.second), l4PtTmp, weight);
      //    SafeHistFill(histMap1D_, getHistName("LepPt1Full", variation.second), l1PtTmp, weight);
      //    SafeHistFill(histMap1D_, getHistName("LepPt2Full", variation.second), l2PtTmp, weight);
      //    SafeHistFill(histMap1D_, getHistName("LepPt3Full", variation.second), l3PtTmp, weight);
      //    SafeHistFill(histMap1D_, getHistName("LepPt4Full", variation.second), l4PtTmp, weight);
      //    SafeHistFill(histMap1D_, getHistName("e1PtSortedFull", variation.second), e1PtTmp, weight);
      //    SafeHistFill(histMap1D_, getHistName("e2PtSortedFull", variation.second), e2PtTmp, weight);
    }

    // eta for all jets in full mass range
    for (unsigned int ind = 0; ind < jetPt->size(); ind++)
    {
      //    SafeHistFill(histMap1D_, getHistName("jetEtaAllj", variation.second), jetEta->at(ind), weight);
      //    SafeHistFill(histMap1D_, getHistName("absjetEtaAllj", variation.second), std::abs(jetEta->at(ind)), weight);

      if (jetPt->at(ind) > 50)
      {
        //      SafeHistFill(histMap1D_, getHistName("jetEtaAllj50", variation.second), jetEta->at(ind), weight);
        //      SafeHistFill(histMap1D_, getHistName("absjetEtaAllj50", variation.second), std::abs(jetEta->at(ind)), weight);
      }
    }

    if (Mass < 180)
    {
      for (unsigned int ind = 0; ind < jetPt->size(); ind++)
      {
        //      SafeHistFill(histMap1D_, getHistName("jetEtaAllj_180", variation.second), jetEta->at(ind), weight);
        //      SafeHistFill(histMap1D_, getHistName("absjetEtaAllj_180", variation.second), std::abs(jetEta->at(ind)), weight);

        if (jetPt->at(ind) > 50)
        {
          //        SafeHistFill(histMap1D_, getHistName("jetEtaAllj50_180", variation.second), jetEta->at(ind), weight);
          //        SafeHistFill(histMap1D_, getHistName("absjetEtaAllj50_180", variation.second), std::abs(jetEta->at(ind)), weight);
        }
      }
    }
    // bool noBlind = true;
    // Applying the ZZ Selection here
    // std::cout<<"Is fillHistograms working?"<<std::endl;
    // std::cout<<"isNonPrompt_ in FillHistograms:"<<isNonPrompt_<<std::endl;

    //=========================================================
  } //! Turn full m4l range processing on/off
  //=========================================================

  if (!PassesZZSelection(isNonPrompt_))
  {
    return;
  }

  //===========================================================================================
  //
  //                       _          _ _                  _
  //   ___  _ __        ___| |__   ___| | |  _ __ ___  __ _(_) ___  _ __
  //  / _ \| '_ \ _____/ __| '_ \ / _ \ | | | '__/ _ \/ _` | |/ _ \| '_ \ 
  // | (_) | | | |_____\__ \ | | |  __/ | | | | |  __/ (_| | | (_) | | | |
  //  \___/|_| |_|     |___/_| |_|\___|_|_| |_|  \___|\__, |_|\___/|_| |_|
  //                                                  |___/
  //
  // On-shell ZZ regions starts here after the "PassesZZSelection" requirement above
  //
  //
  //
  //===========================================================================================

  for (unsigned int ind = 0; ind < jetPt->size(); ind++)
  {
    //    SafeHistFill(histMap1D_, getHistName("jetEtaAllj_120", variation.second), jetEta->at(ind), weight);
    //    SafeHistFill(histMap1D_, getHistName("absjetEtaAllj_120", variation.second), std::abs(jetEta->at(ind)), weight);

    if (jetPt->at(ind) > 50)
    {
      //      SafeHistFill(histMap1D_, getHistName("jetEtaAllj50_120", variation.second), jetEta->at(ind), weight);
      //      SafeHistFill(histMap1D_, getHistName("absjetEtaAllj50_120", variation.second), std::abs(jetEta->at(ind)), weight);
    }
  }

  // std::cout<<"eventWeight in ZZSelector: "<<weight<<std::endl;
  if ((variation.first == Central || (doaTGC_ && isaTGC_)) && isMC_)
  {
    //! Do jet systematics JES and JER in on-shell region
    if (isMC_) //! duplicated condition
    {

      for (size_t i = 0; i < vjetEta.size(); i++)
      {
        SafeHistFill(jethistMap1D_, getHistName("nJets", variation.second), vnJets[i], i, weight);

        if (vnJets[i] > 0)
        {
          SafeHistFill(jethistMap1D_, getHistName("jetPt[0]", variation.second), vjetPt[i]->at(0), i, weight);

          SafeHistFill(jethistMap1D_, getHistName("jetEta[0]", variation.second), vjetEta[i]->at(0), i, weight);
          SafeHistFill(jethistMap1D_, getHistName("absjetEta[0]", variation.second), std::abs(vjetEta[i]->at(0)), i, weight);
        }
        if (vnJets[i] > 1)
        {
          SafeHistFill(jethistMap1D_, getHistName("jetPt[1]", variation.second), vjetPt[i]->at(1), i, weight);

          SafeHistFill(jethistMap1D_, getHistName("jetEta[1]", variation.second), vjetEta[i]->at(1), i, weight); //}
          SafeHistFill(jethistMap1D_, getHistName("absjetEta[1]", variation.second), std::abs(vjetEta[i]->at(1)), i, weight);

          SafeHistFill(jethistMap1D_, getHistName("dEtajj", variation.second), std::abs(vjetEta[i]->at(0) - vjetEta[i]->at(1)), i, weight);
          SafeHistFill(jethistMap1D_, getHistName("mjj", variation.second), vmjj[i], i, weight);
        }

        //! No actual syst for m4l but just for consistency
        SafeHistFill(jethistMap1D_, getHistName("Mass", variation.second), Mass, i, weight);

        if (vnJets[i] == 0)
        {
          SafeHistFill(jethistMap1D_, getHistName("Mass0j", variation.second), Mass, i, weight);
        }
        else if (vnJets[i] == 1)
        {
          SafeHistFill(jethistMap1D_, getHistName("Mass1j", variation.second), Mass, i, weight);
        }
        else if (vnJets[i] == 2)
        {
          SafeHistFill(jethistMap1D_, getHistName("Mass2j", variation.second), Mass, i, weight);
        }
        else
        {
          if (vnJets[i] == 3)
          {
            SafeHistFill(jethistMap1D_, getHistName("Mass3j", variation.second), Mass, i, weight);
          }

          if (vnJets[i] >= 3)
          {
            SafeHistFill(jethistMap1D_, getHistName("Mass34j", variation.second), Mass, i, weight);
          }

          if (vnJets[i] >= 4)
          {
            SafeHistFill(jethistMap1D_, getHistName("Mass4j", variation.second), Mass, i, weight);
          }
        }
      } // loop over syst indices

    } // if (isMC_) for jet systs

    // std::cout<<"does it go into lheWeights"<<std::endl;
    // std::cout << "lheWeights.size() " << lheWeights.size() << std::endl;

    //! LHE reweighting, related to syst. unc. with PDF and QCD scales
    if (fillLHEwgt)
    {
      for (size_t i = 0; i < lheWeights.size(); i++) // expect 0 to 111 currently
      {
        SafeHistFill(weighthistMap1D_, getHistName("yield", variation.second), 1, i, lheWeights[i] / lheWeights[0] * weight);
        SafeHistFill(weighthistMap1D_, getHistName("Mass", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
        //      SafeHistFill(weighthistMap1D_, getHistName("ZZPt", variation.second), Pt, i, lheWeights[i] / lheWeights[0] * weight);
        //      SafeHistFill(weighthistMap1D_, getHistName("ZPt", variation.second), Z1pt, i, lheWeights[i] / lheWeights[0] * weight);
        //      SafeHistFill(weighthistMap1D_, getHistName("ZPt", variation.second), Z2pt, i, lheWeights[i] / lheWeights[0] * weight);
        //      SafeHistFill(weighthistMap1D_, getHistName("LepPt", variation.second), l1Pt, i, lheWeights[i] / lheWeights[0] * weight);
        //      SafeHistFill(weighthistMap1D_, getHistName("LepPt", variation.second), l2Pt, i, lheWeights[i] / lheWeights[0] * weight);
        //      SafeHistFill(weighthistMap1D_, getHistName("LepPt", variation.second), l3Pt, i, lheWeights[i] / lheWeights[0] * weight);
        //      SafeHistFill(weighthistMap1D_, getHistName("LepPt", variation.second), l4Pt, i, lheWeights[i] / lheWeights[0] * weight);
        //      SafeHistFill(weighthistMap1D_, getHistName("dPhiZ1Z2", variation.second), dPhiZZ, i, lheWeights[i] / lheWeights[0] * weight);
        //      SafeHistFill(weighthistMap1D_, getHistName("dRZ1Z2", variation.second), dRZZ, i, lheWeights[i] / lheWeights[0] * weight);

        SafeHistFill(weighthistMap1D_, getHistName("nJets", variation.second), jetPt->size(), i, lheWeights[i] / lheWeights[0] * weight);

        if (jetPt->size() > 0 && jetPt->size() == jetEta->size())
        {
          SafeHistFill(weighthistMap1D_, getHistName("jetPt[0]", variation.second), jetPt->at(0), i, lheWeights[i] / lheWeights[0] * weight);

          SafeHistFill(weighthistMap1D_, getHistName("jetEta[0]", variation.second), jetEta->at(0), i, lheWeights[i] / lheWeights[0] * weight);
          SafeHistFill(weighthistMap1D_, getHistName("absjetEta[0]", variation.second), std::abs(jetEta->at(0)), i, lheWeights[i] / lheWeights[0] * weight);
        }
        if (jetPt->size() > 1 && jetPt->size() == jetEta->size())
        {
          SafeHistFill(weighthistMap1D_, getHistName("jetPt[1]", variation.second), jetPt->at(1), i, lheWeights[i] / lheWeights[0] * weight);

          SafeHistFill(weighthistMap1D_, getHistName("jetEta[1]", variation.second), jetEta->at(1), i, lheWeights[i] / lheWeights[0] * weight); //}
          SafeHistFill(weighthistMap1D_, getHistName("absjetEta[1]", variation.second), std::abs(jetEta->at(1)), i, lheWeights[i] / lheWeights[0] * weight);

          SafeHistFill(weighthistMap1D_, getHistName("dEtajj", variation.second), dEtajj, i, lheWeights[i] / lheWeights[0] * weight);
          SafeHistFill(weighthistMap1D_, getHistName("mjj", variation.second), mjj, i, lheWeights[i] / lheWeights[0] * weight);
        }

        if (jetPt->size() == 0 && jetPt->size() == jetEta->size())
        {
          SafeHistFill(weighthistMap1D_, getHistName("Mass0j", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
        }
        else if (jetPt->size() == 1 && jetPt->size() == jetEta->size())
        {
          SafeHistFill(weighthistMap1D_, getHistName("Mass1j", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
        }
        else if (jetPt->size() == 2 && jetPt->size() == jetEta->size())
        {
          SafeHistFill(weighthistMap1D_, getHistName("Mass2j", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
        }
        else
        {

          if (jetPt->size() == 3 && jetPt->size() == jetEta->size())
          {
            SafeHistFill(weighthistMap1D_, getHistName("Mass3j", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
          }

          if (jetPt->size() >= 3 && jetPt->size() == jetEta->size())
          {
            SafeHistFill(weighthistMap1D_, getHistName("Mass34j", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
          }

          if (jetPt->size() >= 4 && jetPt->size() == jetEta->size())
          {
            SafeHistFill(weighthistMap1D_, getHistName("Mass4j", variation.second), Mass, i, lheWeights[i] / lheWeights[0] * weight);
          }
        }
      } //! For-loop for filling LHE weight
    } //! if fillLHEwgt
  }
  // std::cout<<"isNonPrompt_ in FillHistograms after ZZSelection:"<<isNonPrompt_<<std::endl;
  // std::cout<<run<<":"<<lumi<<":"<<evt<<std::endl;
  // std::cout << "variation.second: "<<variation.second;

  //!=====================A place where the on-shell selections have been applied and we fill the ntuple ====================================================

  float tempNA = -99999.;

  if (writeNtp_ && !writeNtpFullRange)
  {

    //*Set weight branch of central and all systs to a N/A variable, then only fill nonzero weight for the current syst branch*
    for (const auto &syst_temp : variations_)
    {
      SafeSetBranch(ftntp_, getBranchName("weight", syst_temp.second), &tempNA);
    }

    SafeSetBranch(ftntp_, getBranchName("weight", variation.second), &weight);

    // getBranchName function only changes the branch name for weight
    SafeSetBranch(ftntp_, getBranchName("Mass", variation.second), &Mass);
    SafeSetBranch(ftntp_, getBranchName("nJets", variation.second), &nJets_tmp);
    SafeSetBranch(ftntp_, getBranchName("jetPt0", variation.second), &jpt0_tmp);
    SafeSetBranch(ftntp_, getBranchName("jetEta0", variation.second), &jeta0_tmp);

    SafeSetBranch(ftntp_, getBranchName("run", variation.second), &run);
    SafeSetBranch(ftntp_, getBranchName("lumi", variation.second), &lumi);
    SafeSetBranch(ftntp_, getBranchName("evt", variation.second), &evt);

    SafeSetBranch(ftntp_, getBranchName("l1pt", variation.second), &l1pt_tmp);
    SafeSetBranch(ftntp_, getBranchName("l2pt", variation.second), &l2pt_tmp);
    SafeSetBranch(ftntp_, getBranchName("l3pt", variation.second), &l3pt_tmp);
    SafeSetBranch(ftntp_, getBranchName("l4pt", variation.second), &l4pt_tmp);

    SafeSetBranch(ftntp_, getBranchName("l1eta", variation.second), &l1eta_tmp);
    SafeSetBranch(ftntp_, getBranchName("l2eta", variation.second), &l2eta_tmp);
    SafeSetBranch(ftntp_, getBranchName("l3eta", variation.second), &l3eta_tmp);
    SafeSetBranch(ftntp_, getBranchName("l4eta", variation.second), &l4eta_tmp);
    SafeSetBranch(ftntp_, getBranchName("jEtaProduct", variation.second), &jetaproduct_tmp);
    SafeSetBranch(ftntp_, getBranchName("Zep1", variation.second), &zep1_tmp);
    SafeSetBranch(ftntp_, getBranchName("Zep2", variation.second), &zep2_tmp);

    SafeSetBranch(ftntp_, getBranchName("jetPt1", variation.second), &jpt1_tmp);
    SafeSetBranch(ftntp_, getBranchName("jetEta1", variation.second), &jeta1_tmp);

    SafeSetBranch(ftntp_, getBranchName("mjj", variation.second), &mjj);
    SafeSetBranch(ftntp_, getBranchName("dEtajj", variation.second), &dEtajj);
    SafeSetBranch(ftntp_, getBranchName("dPhiZZ", variation.second), &dPhiZZ);

    if (isMC_)
    {
      SafeSetBranch(ftntp_, getBranchName("genWeight", variation.second), &genWeight);
      SafeSetBranch(ftntp_, getBranchName("L1prefiringWeight", variation.second), &L1prefiringWeight);
      SafeSetBranch(ftntp_, getBranchName("L1prefiringWeightUp", variation.second), &L1prefiringWeightUp);
      SafeSetBranch(ftntp_, getBranchName("L1prefiringWeightDn", variation.second), &L1prefiringWeightDn);

      //! fill jet syst variables
      for (size_t i = 0; i < vjetEta.size(); i++)
      {

        SafeSetBranch(ftntp_, getBranchName("nJets" + jetSystNames[i], variation.second), &vnJets[i]);

        if (vnJets[i] > 0)
        {
          jetPt0_syst_tmp[i] = vjetPt[i]->at(0);
          jetEta0_syst_tmp[i] = vjetEta[i]->at(0);
        }

        if (vnJets[i] > 1)
        {
          jetPt1_syst_tmp[i] = vjetPt[i]->at(1);
          jetEta1_syst_tmp[i] = vjetEta[i]->at(1);
          mjj_syst_tmp[i] = vmjj[i];
          dEtajj_syst_tmp[i] = std::abs(vjetEta[i]->at(0) - vjetEta[i]->at(1));
          jetaproduct_syst_tmp[i] = vjetEta[i]->at(0) * vjetEta[i]->at(1);
          zep1_syst_tmp[i] = Z1Eta- (vjetEta[i]->at(0) + vjetEta[i]->at(1));
          zep2_syst_tmp[i] = Z2Eta- (vjetEta[i]->at(0) + vjetEta[i]->at(1));
        }

        SafeSetBranch(ftntp_, getBranchName("jetPt0" + jetSystNames[i], variation.second), &jetPt0_syst_tmp[i]);
        SafeSetBranch(ftntp_, getBranchName("jetEta0" + jetSystNames[i], variation.second), &jetEta0_syst_tmp[i]);
        SafeSetBranch(ftntp_, getBranchName("jetPt1" + jetSystNames[i], variation.second), &jetPt1_syst_tmp[i]);
        SafeSetBranch(ftntp_, getBranchName("jetEta1" + jetSystNames[i], variation.second), &jetEta1_syst_tmp[i]);
        SafeSetBranch(ftntp_, getBranchName("mjj" + jetSystNames[i], variation.second), &mjj_syst_tmp[i]);
        SafeSetBranch(ftntp_, getBranchName("dEtajj" + jetSystNames[i], variation.second), &dEtajj_syst_tmp[i]);
        SafeSetBranch(ftntp_, getBranchName("jEtaProduct" + jetSystNames[i], variation.second), &jetaproduct_syst_tmp[i]);
        SafeSetBranch(ftntp_, getBranchName("Zep1" + jetSystNames[i], variation.second), &zep1_syst_tmp[i]);
        SafeSetBranch(ftntp_, getBranchName("Zep2" + jetSystNames[i], variation.second), &zep2_syst_tmp[i]);
      }
    }

    ftntp_->Fill();
  }

  SafeHistFill(histMap1D_, getHistName("Mass", variation.second), Mass, weight);
  //  SafeHistFill(histMap1D_, getHistName("PVDZ", variation.second), l1PVDZ, weight);
  //  SafeHistFill(histMap1D_, getHistName("PVDZ", variation.second), l2PVDZ, weight);
  //  SafeHistFill(histMap1D_, getHistName("PVDZ", variation.second), l3PVDZ, weight);
  //  SafeHistFill(histMap1D_, getHistName("PVDZ", variation.second), l4PVDZ, weight);
  //  SafeHistFill(histMap1D_, getHistName("deltaPVDZ_sameZ", variation.second), std::abs(l1PVDZ - l2PVDZ), weight);
  //  SafeHistFill(histMap1D_, getHistName("deltaPVDZ_sameZ", variation.second), std::abs(l3PVDZ - l4PVDZ), weight);
  //  SafeHistFill(histMap1D_, getHistName("deltaPVDZ_diffZ", variation.second), std::abs(l1PVDZ - l3PVDZ), weight);
  //  SafeHistFill(histMap1D_, getHistName("deltaPVDZ_diffZ", variation.second), std::abs(l2PVDZ - l4PVDZ), weight);
  //  SafeHistFill(histMap1D_, getHistName("deltaPVDZ_diffZ", variation.second), std::abs(l1PVDZ - l4PVDZ), weight);
  //  SafeHistFill(histMap1D_, getHistName("deltaPVDZ_diffZ", variation.second), std::abs(l2PVDZ - l3PVDZ), weight);

  // if (!PassesZZjjSelection()){
  //   return;
  // }

  // if (Mass <= 180){
  //  return;
  // }

  // if (Mass <= 70.0 || Mass>=1000.0){
  //   return;
  // }

  SafeHistFill(histMap1D_, getHistName("yield", variation.second), 1, weight);
  // SafeHistFill(histMap1D_, getHistName("Mass", variation.second), Mass,weight);
  //  SafeHistFill(histMap1D_, getHistName("Z1Mass", variation.second), Z1mass, weight);
  //  SafeHistFill(histMap1D_, getHistName("Z2Mass", variation.second), Z2mass, weight);
  //  SafeHistFill(histMap1D_, getHistName("ZMass", variation.second), Z1mass, weight);
  //  SafeHistFill(histMap1D_, getHistName("ZMass", variation.second), Z2mass, weight);
  //  SafeHistFill(histMap1D_, getHistName("ZPt", variation.second), Z1pt, weight);
  //  SafeHistFill(histMap1D_, getHistName("ZPt", variation.second), Z2pt, weight);
  //  SafeHistFill(histMap1D_, getHistName("dPhiZ1Z2", variation.second), dPhiZZ, weight);
  //  SafeHistFill(histMap1D_, getHistName("dRZ1Z2", variation.second), dRZZ, weight);
  //  SafeHistFill(histMap1D_, getHistName("ZZPt", variation.second), Pt, weight);
  //  SafeHistFill(histMap1D_, getHistName("ZZEta", variation.second), Eta, weight);

  //  SafeHistFill(histMap1D_, getHistName("SIP3D", variation.second), l1SIP3D, weight);
  //  SafeHistFill(histMap1D_, getHistName("SIP3D", variation.second), l2SIP3D, weight);
  //  SafeHistFill(histMap1D_, getHistName("SIP3D", variation.second), l3SIP3D, weight);
  //  SafeHistFill(histMap1D_, getHistName("SIP3D", variation.second), l4SIP3D, weight);
  // Making LeptonPt and Eta plots
  //  SafeHistFill(histMap1D_, getHistName("LepPt", variation.second), l1PtTmp, weight);
  //  SafeHistFill(histMap1D_, getHistName("LepPt", variation.second), l2PtTmp, weight);
  //  SafeHistFill(histMap1D_, getHistName("LepPt", variation.second), l3PtTmp, weight);
  //  SafeHistFill(histMap1D_, getHistName("LepPt", variation.second), l4PtTmp, weight);
  //  SafeHistFill(histMap1D_, getHistName("LepPt1", variation.second), l1PtTmp, weight);
  //  SafeHistFill(histMap1D_, getHistName("LepPt2", variation.second), l2PtTmp, weight);
  //  SafeHistFill(histMap1D_, getHistName("LepPt3", variation.second), l3PtTmp, weight);
  //  SafeHistFill(histMap1D_, getHistName("LepPt4", variation.second), l4PtTmp, weight);
  //  SafeHistFill(histMap1D_, getHistName("e1PtSorted", variation.second), e1PtTmp, weight);
  //  SafeHistFill(histMap1D_, getHistName("e2PtSorted", variation.second), e2PtTmp, weight);
  //  SafeHistFill(histMap1D_, getHistName("LepEta", variation.second), l1Eta, weight);
  //  SafeHistFill(histMap1D_, getHistName("LepEta", variation.second), l2Eta, weight);
  //  SafeHistFill(histMap1D_, getHistName("LepEta", variation.second), l3Eta, weight);
  //  SafeHistFill(histMap1D_, getHistName("LepEta", variation.second), l4Eta, weight);
  SafeHistFill(histMap1D_, getHistName("nJets", variation.second), jetPt->size(), weight);

  int central_nJets = 0;
  for (unsigned int cind = 0; cind < jetPt->size(); cind++)
  {
    if (std::abs(jetEta->at(cind)) < 2.4)
    {
      central_nJets = central_nJets + 1;
    }
  }

  SafeHistFill(histMap1D_, getHistName("nJets_central", variation.second), central_nJets, weight);

  if (jetPt->size() == 0 && jetPt->size() == jetEta->size())
  {
    SafeHistFill(histMap1D_, getHistName("Mass0j", variation.second), Mass, weight);
  }
  else if (jetPt->size() == 1 && jetPt->size() == jetEta->size())
  {
    SafeHistFill(histMap1D_, getHistName("Mass1j", variation.second), Mass, weight);
  }
  else if (jetPt->size() == 2 && jetPt->size() == jetEta->size())
  {
    SafeHistFill(histMap1D_, getHistName("Mass2j", variation.second), Mass, weight);
  }
  else
  {
    if (jetPt->size() == 3 && jetPt->size() == jetEta->size())
    {
      SafeHistFill(histMap1D_, getHistName("Mass3j", variation.second), Mass, weight);
    }

    if (jetPt->size() >= 3 && jetPt->size() == jetEta->size())
    {
      SafeHistFill(histMap1D_, getHistName("Mass34j", variation.second), Mass, weight);
    }

    if (jetPt->size() >= 4 && jetPt->size() == jetEta->size())
    {
      SafeHistFill(histMap1D_, getHistName("Mass4j", variation.second), Mass, weight);
    }
  }

  if (jetPt->size() > 0 && jetPt->size() == jetEta->size())
  {
    SafeHistFill(histMap1D_, getHistName("jetPt[0]", variation.second), jetPt->at(0), weight);
    // if (jetPt->size() ==2 && jetPt->size() == jetEta->size()) {
    SafeHistFill(histMap1D_, getHistName("jetEta[0]", variation.second), jetEta->at(0), weight); //}
    SafeHistFill(histMap1D_, getHistName("absjetEta[0]", variation.second), std::abs(jetEta->at(0)), weight);
    //    SafeHistFill(histMap1D_, getHistName("jetPhi[0]", variation.second), jetPhi->at(0), weight);

    for (unsigned int ind = 0; ind < jetPt->size(); ind++)
    {
      // SafeHistFill(histMap1D_, getHistName("jetEtaAllj", variation.second), jetEta->at(ind), weight);
      // SafeHistFill(histMap1D_, getHistName("absjetEtaAllj", variation.second), std::abs(jetEta->at(ind)), weight);

      if (run < 319077)
      {
        //        SafeHistFill(jetTestMap2D_, getHistName("jetHEM_AB", variation.second), jetPhi->at(ind), jetEta->at(ind), weight);
        //        SafeHistFill(jetTestMap2D_, getHistName("jetHEM2_AB", variation.second), jetPhi->at(ind), jetEta->at(ind), jetPt->at(ind) * weight);
      }
      else
      {
        //        SafeHistFill(jetTestMap2D_, getHistName("jetHEM_CD", variation.second), jetPhi->at(ind), jetEta->at(ind), weight);
        //        SafeHistFill(jetTestMap2D_, getHistName("jetHEM2_CD", variation.second), jetPhi->at(ind), jetEta->at(ind), jetPt->at(ind) * weight);
      }
    }

    if (jetPt->size() == 1 && jetPt->size() == jetEta->size())
    {
      //      SafeHistFill(histMap1D_, getHistName("jetPtN1", variation.second), jetPt->at(0), weight);
      //      SafeHistFill(histMap1D_, getHistName("absjetEtaN1", variation.second), std::abs(jetEta->at(0)), weight);
      //      SafeHistFill(jetTestMap2D_, getHistName("jetPtN1", variation.second), jetPt->at(0), std::abs(jetEta->at(0)), weight);
      if (jetPt->at(0) < 100)
      {
        //        SafeHistFill(histMap1D_, getHistName("absjetEtaN1_100", variation.second), std::abs(jetEta->at(0)), weight);
      }
    }
  }
  if (jetPt->size() > 1 && jetPt->size() == jetEta->size())
  {
    SafeHistFill(histMap1D_, getHistName("jetPt[1]", variation.second), jetPt->at(1), weight);
    // if (jetPt->size() ==2 && jetPt->size() == jetEta->size()) {
    SafeHistFill(histMap1D_, getHistName("jetEta[1]", variation.second), jetEta->at(1), weight); //}
    SafeHistFill(histMap1D_, getHistName("absjetEta[1]", variation.second), std::abs(jetEta->at(1)), weight);
    //    SafeHistFill(histMap1D_, getHistName("jetPhi[1]", variation.second), jetPhi->at(1), weight);
    SafeHistFill(histMap1D_, getHistName("dEtajj", variation.second), dEtajj, weight);
    SafeHistFill(histMap1D_, getHistName("mjj", variation.second), mjj, weight);

    if (jetPt->size() == 2 && jetPt->size() == jetEta->size())
    {
      //      SafeHistFill(histMap1D_, getHistName("jetEta[01]", variation.second), jetEta->at(0), weight);
      //      SafeHistFill(histMap1D_, getHistName("jetEta[01]", variation.second), jetEta->at(1), weight);
      // pt vs eta for 2-jet event lowest pt jet
      //      SafeHistFill(jetTestMap2D_, getHistName("jetPtN2", variation.second), jetPt->at(1), std::abs(jetEta->at(1)), weight);
    }
    //    SafeHistFill(histMap1D_, getHistName("jetPt[01]", variation.second), jetPt->at(0), weight);
    //    SafeHistFill(histMap1D_, getHistName("jetPt[01]", variation.second), jetPt->at(1), weight);
  }

  if (jetPt->size() > 2 && jetPt->size() == jetEta->size())
  {
    //    SafeHistFill(histMap1D_, getHistName("jetPt[2]", variation.second), jetPt->at(2), weight);
    //    SafeHistFill(histMap1D_, getHistName("jetEta[2]", variation.second), jetEta->at(2), weight);
    //    SafeHistFill(histMap1D_, getHistName("jetPhi[2]", variation.second), jetPhi->at(2), weight);

    if (jetPt->size() == 3 && jetPt->size() == jetEta->size())
    {

      // pt vs eta for 3-jet event lowest pt jet
      //      SafeHistFill(jetTestMap2D_, getHistName("jetPtN3", variation.second), jetPt->at(2), std::abs(jetEta->at(2)), weight);
    }
  }

  //    if (!PassesZZjjSelection()){
  //  return;
  //}

  // if (Mass <= 180.0){
  //   return;
  // }

  if (isMC_)
  {
    // std::cout<<run<<":"<<lumi<<":"<<evt<<std::endl;
    // std::cout<<"UpdatedSF:"<<weight<<std::endl;
  }

  // SafeHistFill(histMap1D_, getHistName("Mass", variation.second), Mass,weight);
  // SafeHistFill(histMap1D_, getHistName("dEtajj", variation.second), dEtajj, weight);

  // Summing 12,34 leptons
  // SafeHistFill(histMap1D_, getHistName("Lep12Pt", variation.second), l1Pt, weight);
  // SafeHistFill(histMap1D_, getHistName("Lep12Pt", variation.second), l2Pt, weight);
  // SafeHistFill(histMap1D_, getHistName("Lep34Pt", variation.second), l3Pt, weight);
  // SafeHistFill(histMap1D_, getHistName("Lep34Pt", variation.second), l4Pt, weight);
  // SafeHistFill(histMap1D_, getHistName("Lep12Eta", variation.second), l1Eta, weight);
  // SafeHistFill(histMap1D_, getHistName("Lep12Eta", variation.second), l2Eta, weight);
  // SafeHistFill(histMap1D_, getHistName("Le34Eta", variation.second), l3Eta, weight);
  // SafeHistFill(histMap1D_, getHistName("Lep34Eta", variation.second), l4Eta, weight);
  // SafeHistFill(histMap1D_, getHistName("Z1Mass", variation.second), Z1mass, weight);
  // SafeHistFill(histMap1D_, getHistName("Z2Mass", variation.second), Z2mass, weight);
  // SafeHistFill(histMap1D_, getHistName("Z1Pt", variation.second), Z1pt, weight);
  // SafeHistFill(histMap1D_, getHistName("Z2Pt", variation.second), Z2pt, weight);
  // SafeHistFill(histMap1D_, getHistName("Z1Phi", variation.second), Z1Phi, weight);
  // SafeHistFill(histMap1D_, getHistName("Z2Phi", variation.second), Z2Phi, weight);
  // SafeHistFill(histMap1D_, getHistName("Z1lep1_Pt", variation.second), l1Pt, weight);
  // SafeHistFill(histMap1D_, getHistName("Z1lep1_Eta", variation.second), l1Eta, weight);
  // SafeHistFill(histMap1D_, getHistName("Z1lep1_Phi", variation.second), l1Phi, weight);
  // SafeHistFill(histMap1D_, getHistName("Z1lep1_PdgId", variation.second), l1PdgId, weight);
  // SafeHistFill(histMap1D_, getHistName("Z1lep2_Pt", variation.second), l2Pt, weight);
  // SafeHistFill(histMap1D_, getHistName("Z1lep2_Eta", variation.second), l2Eta, weight);
  // SafeHistFill(histMap1D_, getHistName("Z1lep2_Phi", variation.second), l2Phi, weight);
  // SafeHistFill(histMap1D_, getHistName("Z1lep2_PdgId", variation.second), l2PdgId, weight);
  // SafeHistFill(histMap1D_, getHistName("Z2lep1_Pt", variation.second), l3Pt, weight);
  // SafeHistFill(histMap1D_, getHistName("Z2lep1_Eta", variation.second), l3Eta, weight);
  // SafeHistFill(histMap1D_, getHistName("Z2lep1_Phi", variation.second), l3Phi, weight);
  // SafeHistFill(histMap1D_, getHistName("Z2lep1_PdgId", variation.second), l3PdgId, weight);
  // SafeHistFill(histMap1D_, getHistName("Z2lep2_Pt", variation.second), l4Pt, weight);
  // SafeHistFill(histMap1D_, getHistName("Z2lep2_Eta", variation.second), l4Eta, weight);
  // SafeHistFill(histMap1D_, getHistName("Z2lep2_Phi", variation.second), l4Phi, weight);
  // SafeHistFill(histMap1D_, getHistName("Z2lep2_PdgId", variation.second), l4PdgId, weight);
  // SafeHistFill(hists2D_, getHistName("Z1lep1_Z1lep2_Pt",variation.second),l1Pt,l2Pt,weight);
  // SafeHistFill(hists2D_, getHistName("Z1lep1_Z1lep2_Eta",variation.second),l1Eta,l2Eta,weight);
  // SafeHistFill(hists2D_, getHistName("Z1lep1_Z1lep2_Phi",variation.second),l1Phi,l2Phi,weight);
  // SafeHistFill(hists2D_, getHistName("Z2lep1_Z2lep2_Pt",variation.second),l3Pt,l4Pt,weight);
  // SafeHistFill(hists2D_, getHistName("Z2lep1_Z2lep2_Eta",variation.second),l3Eta,l4Eta,weight);
  // SafeHistFill(hists2D_, getHistName("Z2lep1_Z2lep2_Phi",variation.second),l3Phi,l4Phi,weight);
  ////2D Z1 vs Z2
  // SafeHistFill(hists2D_, getHistName("Z1Mass_Z2Mass",variation.second),Z1mass,Z2mass,weight);

  // if (histMap1D_[getHistName("nvtx", variation.second)] != nullptr) {
  //     b_nvtx->GetEntry(entry);
  //     histMap1D_[getHistName("nvtx", variation.second)]->Fill(nvtx, weight);
  // }
  // if (isMC_)
  //   SafeHistFill(histMap1D_, getHistName("nTruePU", variation.second), nTruePU, weight);
}

void ZZSelector::SetupNewDirectory()
{
  SelectorBase::SetupNewDirectory();
  isaTGC_ = name_.find("atgc") != std::string::npos;
  // std::cout<<"selection in ZZSelector: "<<selection_<<std::endl;
  applyFullSelection_ = (selection_ == ZZselection);
  // std::cout<<applyFullSelection_<<std::endl;
  InitializeHistogramsFromConfig();
  // std::cout<<"Do Histos get initialized"<<std::endl;
}
