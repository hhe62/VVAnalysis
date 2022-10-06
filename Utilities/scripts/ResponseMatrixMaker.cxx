#include "ResponseMatrixMaker.hxx"
#include<cstdlib>
#include<iostream>
#include<cmath> // std::sqrt


namespace
{
  float deltaPhi(float phi1, float phi2)
  {
    const float pi = 3.14159265;
    float out = phi1 - phi2;
    while(out > pi)
      out -= 2. * pi;
    while(out < -1. * pi)
      out += 2. * pi;

    return out;
  }

  float deltaR(float eta1, float phi1, float eta2, float phi2)
  {
    float dPhi = deltaPhi(phi1, phi2);
    float dEta = eta1 - eta2;
    return std::sqrt(dPhi * dPhi + dEta * dEta);
  }

  template<class HType>
  int getBinIndex(HType& h, float x)
  {
    int bin = h.FindBin(x); // root wants a signed int for some reason
    if(h.IsBinUnderflow(bin))
      bin = 1;
    else if(h.IsBinOverflow(bin))
      bin -= 1;

    return bin;
  }

  template<class HType>
  int getBinIndex(HType& h, float x, float y)
  {
    int bin = h.FindBin(x, y);

    if(h.IsBinOverflow(bin))
      {
        int binx, biny, binz;
        h.GetBinXYZ(bin, binx, biny, binz);
        if(binx > h.GetNbinsX())
          binx -= 1;
        if(biny > h.GetNbinsY())
          biny -= 1;

        bin = h.GetBin(binx, biny, binz);
      }
    if(h.IsBinUnderflow(bin))
      {
        int binx, biny, binz;
        h.GetBinXYZ(bin, binx, biny, binz);
        if(!binx)
          binx += 1;
        if(!biny)
          biny += 1;

        bin = h.GetBin(binx, biny, binz);
      }

    return bin;
  }

  template<class HType>
  float getContentFromHist(HType& h, float x)
  {
    int bin = getBinIndex(h, x);
    return h.GetBinContent(bin);
  }

  template<class HType>
  float getContentFromHist(HType& h, float x, float y)
  {
    int bin = getBinIndex(h, x, y);

    return h.GetBinContent(bin);
  }

  template<class HType>
  float getErrorFromHist(HType& h, float x)
  {
    int bin = getBinIndex(h, x);
    return h.GetBinError(bin);
  }

  template<class HType>
  float getErrorFromHist(HType& h, float x, float y)
  {
    int bin = getBinIndex(h, x, y);

    return h.GetBinError(bin);
  }
};


template<typename T>
ResponseMatrixMakerBase<T>::ResponseMatrixMakerBase(const Str& channel,
                                                    const Str& varName,
                                                    const Vec<float>& binning) :
  varName(varName),
  channel(channel),
  binning(binning),
  nPDFVariations(100),
  iAlphaSDn(nPDFVariations), //I switched the up down from the original codes and expect this to be correct
  iAlphaSUp(nPDFVariations+1)
{
  scale = 1.;
  skipSyst = false;

  pdfResponses = TH3D("null", "null", 1, 0., 1., 1, 0., 1., 1, 0., 1.);
}


// template<typename T>
// void ResponseMatrixMakerBase<T>::registerLeptonSFs(const TH2F& h, const Str& lep,
//                                               const Str& upOrDown = "")
// {
//   if(leptonSFHists.find(lep) == UMap<Str, TH2F&>::end)
//     leptonSFHists[lep] = UMap<Str, TH2F&>();
//
//   leptonSFHists[lep][upOrDown] = h;
// }


template<typename T>
const TH2D& ResponseMatrixMakerBase<T>::getResponse(const Str& syst)
{
  if(!responses.size()){
    std::cout<<"When is setup called?"<<std::endl;
    setup();}

  return responses.at(syst);
}


template<typename T>
const TH3D& ResponseMatrixMakerBase<T>::getPDFResponses()
{
  if(!responses.size())
    setup();

  return pdfResponses;
}


template<typename T>
const Vec<TH2D>& ResponseMatrixMakerBase<T>::getScaleResponses()
{
  if(!responses.size())
    setup();

  return scaleResponses;
}


template<typename T>
Vec<Str> ResponseMatrixMakerBase<T>::knownSystematics() const
{
  Vec<Str> out;
  for(auto&& resp : responses)
    {
      if(!resp.first.empty())
        out.push_back(resp.first);
    }

  return out;
}


template<typename T>
void ResponseMatrixMakerBase<T>::setup()
{
  std::cout<<"responsesSize: "<<responses.size()<<std::endl;
  if(responses.size())
    return;

  UPtr<TChain> recoTree(new TChain((getChannel()+"/ntuple").c_str(),
                                   ("recoChain_"+getVar() + "_" + getChannel()).c_str()));
  for(const auto& fn : fileNames)
    {//std::cout<<fn<<std::endl;
      recoTree->Add(fn.c_str());}
  //std::cout<<fileNames<<std::endl;
  //std::cout<<binning<<std::endl;
  // set up lots of things
  Vec<Str> systs = Vec<Str>({"nominal",
        "pu_Up","pu_Down","jes_up","jes_dn","jer_up","jer_dn",
        //"pdf_Up","pdf_Down",
        });
  
  Vec<Str> jsysts = Vec<Str>({"jes_up","jes_dn","jer_up","jer_dn",
        
        });
  //const bool isJetVar = (varName.find("jet") != Str::npos ||
  //                       varName.find("Jet") != Str::npos ||
  //                       varName.find("jj") != Str::npos);
  //if(isJetVar && !skipSyst)
  //  {
  //    systs.push_back("jer_up");
  //    systs.push_back("jer_dn");
  //    systs.push_back("jes_up");
  //    systs.push_back("jes_dn");
  //  }
  bool hasE, hasMu;
  Vec<Str> objects;
  if(channel.find("eeee") != Str::npos)
    {
      objects = Vec<Str>({"e1","e2","e3","e4"});
      hasE = true;
      hasMu = false;
      if(!skipSyst)
        {
          systs.push_back("eEff_Up");
          systs.push_back("eEff_Down");
  //        systs.push_back("eScale_up");
  //        systs.push_back("eScale_dn");
  //        systs.push_back("eRhoRes_up");
  //        systs.push_back("eRhoRes_dn");
  //        systs.push_back("ePhiRes_up");
        }
    }
  else if(channel.find("eemm") != Str::npos)
    {
      objects = Vec<Str>({"e1","e2","m1","m2"});
      hasE = true;
      hasMu = true;
      if(!skipSyst)
        {
          systs.push_back("eEff_Up");
          systs.push_back("eEff_Down");
          systs.push_back("mEff_Up");
          systs.push_back("mEff_Down");
  //        systs.push_back("eScale_up");
  //        systs.push_back("eScale_dn");
  //        systs.push_back("eRhoRes_up");
  //        systs.push_back("eRhoRes_dn");
  //        systs.push_back("ePhiRes_up");
  //        systs.push_back("mClosure_up");
  //        systs.push_back("mClosure_dn");
        }
    }
  else
    {
      objects = Vec<Str>({"m1","m2","m3","m4"});
      hasE = false;
      hasMu = true;
      if(!skipSyst)
        {
          systs.push_back("mEff_Up");
          systs.push_back("mEff_Down");
  //        systs.push_back("mClosure_up");
  //        systs.push_back("mClosure_dn");
        }
    }

  // Scale and PDF systematics only done for samples that have LHE info (e.g.
  // not MCFM)
  // Not doing pdf and scale systematics at the moment. 
  // -> Now include these systematics
  bool hasLHE = bool(recoTree->FindBranch("pdfWeights")); 
  //bool hasLHE = false;

  //see https://twiki.cern.ch/twiki/bin/view/CMS/HowToPDF#How_to_retrieve_LHE_weights_info
  const Vec<size_t> scaleIndicesWeCareAbout = {1,2,3,4,6,8};
  int lastLHEind = 111;
  if (year==2017 || year ==2018){lastLHEind = 111;}
    if (year == 2016) {lastLHEind = 110;}
  if(hasLHE && !skipSyst)
    {//Change the codes to store all scale,PDF+alpha_s variations in a single vector.
       for(size_t i = 0; i <= lastLHEind ; ++i ) //scale 0 to 8, followed by 103 pdf+alpha_s variations
        {
          scaleResponses.push_back(
            TH2D(("scalePDFAlphasVariation"+std::to_string(i)).c_str(), "",
                 binning.size()-1, &binning[0],
                 binning.size()-1, &binning[0]));

    //   for(auto i : scaleIndicesWeCareAbout)
    //     {
    //       scaleResponses.push_back(
    //         TH2D(("scaleVariation"+std::to_string(i)).c_str(), "",
    //              binning.size()-1, &binning[0],
    //              binning.size()-1, &binning[0]));
        }

      //This part remains the same, but the pdfResponses part is no longer used.
      if(!skipSyst)
        {
          Vec<float> iterationBins;
          for(size_t i = 0; i <= nPDFVariations; ++i) //binning points 0 to 100
            iterationBins.push_back(float(i));
          pdfResponses = TH3D("pdfResponses", "",
                              binning.size()-1, &binning[0],
                              binning.size()-1, &binning[0],
                              iterationBins.size()-1, &iterationBins[0]);

          auto scaleWtPtr = &scaleWeights;
          auto pdfWtPtr = &pdfAndAlphaSWeights;
          recoTree->SetBranchAddress("scaleWeights", &scaleWtPtr);
          recoTree->SetBranchAddress("pdfWeights", &pdfWtPtr);

          //systs.push_back("alphaS_up");
         // systs.push_back("alphaS_dn");
        }
    }

  for(auto& s : systs)
    {
      responses[s] = TH2D(s.c_str(), "",
                          binning.size()-1, &binning[0],
                          binning.size()-1, &binning[0]);
    }

  //responses["nominal"] = TH2D("nominal", "",
  //                     binning.size()-1, &binning[0],
  //                     binning.size()-1, &binning[0]);

  // Get gen info
  UPtr<TChain> trueTree = UPtr<TChain>(new TChain((getChannel()+"Gen/ntuple").c_str(),
                                                  ("trueChain_"+getVar() + "_" + getChannel()).c_str()));
  for(const auto& fn : fileNames)
    trueTree->Add(fn.c_str());

  UPtr<UMap<size_t, T> > trueVals = this->getTrueValues(*trueTree.get(),
                                                        objects);

  trueTree.reset(); // gone -- don't use any more

  // Set up common branches
  setCommonBranches(*recoTree, objects);

  this->setRecoBranches(*recoTree.get(), objects);

  bool doPUWt = (puWeightHists.find("") != puWeightHists.end());
  bool doPUWtUp = (puWeightHists.find("Up") != puWeightHists.end());
  bool doPUWtDn = (puWeightHists.find("Down") != puWeightHists.end());
  float upperEdge = binning[binning.size()-1];
  //std::cout<<"upperEdge: "<<upperEdge<<std::endl;
  float lowerEdge = binning[binning.size()-2];
  //std::cout<<"lowerEdge: "<<lowerEdge<<std::endl;
  float binCenter = (lowerEdge)+((upperEdge - lowerEdge)/2.); 
  //std::cout<<"Last bin center: "<<binCenter<<std::endl;
  // Loop through base reco tree, fill most things
  for(size_t row = 0; row < size_t(std::abs(recoTree->GetEntries())); ++row)
    {
      recoTree->GetEntry(row);
      //std::cout<<"genWeight Before: "<<genWeight<<std::endl;

      if (year==2016 ||year==2017){genWeight*=L1prefiringWeight;}

      genWeight *= jetPUSFmulfac;

      //std::cout<<"genWeight After: "<<genWeight<<std::endl;

      // elements needed for event weights
      float puWt = (doPUWt ? ::getContentFromHist(puWeightHists.at(""), truePU) : 1.);
      float puWtUp = (doPUWtUp ? ::getContentFromHist(puWeightHists.at("Up"), truePU) : 1.);
      float puWtDn = (doPUWtDn ? ::getContentFromHist(puWeightHists.at("Down"), truePU) : 1.);

      float lepSF = this->getLepSF(objects);
      float lepSFEUp = 1.;
      float lepSFEDn = 1.;
      float lepSFMUp = 1.;
      float lepSFMDn = 1.;

      if(hasE)
        {
          lepSFEUp = this->getLepSF(objects, 1., 0.);
          lepSFEDn = this->getLepSF(objects, -1., 0.);
        }
      if(hasMu)
        {
          lepSFMUp = this->getLepSF(objects, 0., 1.);
          lepSFMDn = this->getLepSF(objects, 0., -1.);
        }

      // for(size_t i = 0; i < objects.size(); ++i)
      //   {
      //     if(objects.at(i)[0] == "e")
      //       {
      //         lepSF *= ::getContentFromHist(leptonSFHists["e"][""], lPt.at(i), std::abs(lEta.at(i)));
      //         lepSFEUp *= ::getContentFromHist(leptonSFHists["e"]["up"], lPt.at(i), std::abs(lEta.at(i)));
      //         lepSFEDn *= ::getContentFromHist(leptonSFHists["e"]["dn"], lPt.at(i), std::abs(lEta.at(i)));
      //       }
      //     if(objects.at(i)[0] == "m")
      //       {
      //         lepSF *= ::getContentFromHist(leptonSFHists["m"][""], lPt.at(i), lEta.at(i));
      //         lepSFMUp *= ::getContentFromHist(leptonSFHists["m"]["up"], lPt.at(i), lEta.at(i));
      //         lepSFMDn *= ::getContentFromHist(leptonSFHists["m"]["dn"], lPt.at(i), lEta.at(i));
      //       }
      //   }

      auto iTrue = trueVals->find(evt);
      if(iTrue == trueVals->end())
        continue;
      //const T& trueVal = iTrue->second;

      T& trueVal = iTrue->second;
      //if(this->getVar().c_str() == "jetPt0"){trueVal=trueVal.at(0);}

      //Do JES and JER systematics here
      for(auto& js : jsysts){
        if(this->selectEvent(js))
        {  
         
          T val = this->getEventResponse(js);
	
          const float nominalWeight = scale * puWt * lepSF * genWeight;
	 

	        if (val >= upperEdge){val=binCenter;}
	   
          if (trueVal >= upperEdge){trueVal=binCenter;}
	  
         
          this->fillResponse(responses[js], val, trueVal, nominalWeight);

        }
      }

      if(this->selectEvent())
        {  
          // Nominal value
          //const T val = this->getEventResponse();
          //if(this->nJets != this->valuevecfloat->size()) {std::cout<<"==================size unequal========================"<<this->nJets<<" vs "<< this->valuevecfloat->size()<<std::endl;}
          T val = this->getEventResponse();
	  //std::cout<<"========check=================="<<this->getEventResponse()<<"=============="<<val<<std::endl;
	  //std::cout<<"=======upperedge larger?==================="<<upperEdge<<"========"<<(val > upperEdge)<<std::endl;
	  //if(this->getVar().c_str() == "jetPt0"){val=val.at(0);}
          const float nominalWeight = scale * puWt * lepSF * genWeight;
	  //std::cout<<"True Val "<<trueVal<<" RecoVal "<<val<<std::endl;
          //const float nominalWeight = scale * lepSF * genWeight;
          // fill histos that use nominal value but with different weights
         
          //this section is for variables like leppt or all multibranch response matrices
	  
          //for(size_t i = 0; i < val.size(); ++i){
          //  if (val[i] > upperEdge){val[i]=binCenter;}
            //std::cout<<"val: "<<val[i]<<std::endl;
          //}

          //for(size_t i = 0; i < trueVal.size(); ++i){
          //  if (trueVal[i] > upperEdge){trueVal[i]=binCenter;}
            //std::cout<<"trueVal: "<<trueVal[i]<<std::endl;
          //} 
	  
	  
          //this section is for variables like mass or all singlebranch response matrices

	  //std::cout<<"val is of simple float type "<<val<<std::endl;

	  if (val >= upperEdge){val=binCenter;}
	    //std::cout<<"==============value:bincenter"<<val<<" "<<binCenter<<std::endl;}
          if (trueVal >= upperEdge){trueVal=binCenter;}
	  //if (val>=3){
	  //std::cout<<"Test nJets first==============="<<this->getEventResponse()<<std::endl;
	  //std::cout<<"Test nJets val ==============="<<val<<std::endl;

	  //std::cout<<"Test nJets bincenter ==============="<<(val > upperEdge)<<binCenter<<std::endl;
	  //std::cout<<"Test nJets bincenter Truth ==============="<<(trueVal > upperEdge)<<binCenter<<std::endl;}

          //std::cout<<"val: "<<val<<std::endl;
          this->fillResponse(responses["nominal"], val, trueVal, nominalWeight);

          if(!skipSyst)
            {
              this->fillResponse(responses["pu_Up"], val, trueVal, scale * puWtUp * lepSF * genWeight);
              this->fillResponse(responses["pu_Down"], val, trueVal, scale * puWtDn * lepSF * genWeight);
              //Already summed all 100 variations and I have a pdfUp and pdfDown histogram
              //this->fillResponse(responses["pdf_Up"], val, trueVal, scale * puWt * pdfWtUp * lepSF * genWeight);
              //this->fillResponse(responses["pdf_Down"], val, trueVal, scale * puWt * pdfWtDn * lepSF * genWeight);
              if(hasE)
                {
                  this->fillResponse(responses["eEff_Up"], val, trueVal, scale * puWt * lepSFEUp * genWeight);
                  this->fillResponse(responses["eEff_Down"], val, trueVal, scale * puWt * lepSFEDn * genWeight);
                  //this->fillResponse(responses["eEff_Up"], val, trueVal, scale * lepSFEUp * genWeight);
                  //this->fillResponse(responses["eEff_Down"], val, trueVal, scale * lepSFEDn * genWeight);
                }
              if(hasMu)
                {
                  this->fillResponse(responses["mEff_Up"], val, trueVal, scale * puWt * lepSFMUp * genWeight);
                  this->fillResponse(responses["mEff_Down"], val, trueVal, scale * puWt * lepSFMDn * genWeight);
                  //this->fillResponse(responses["mEff_Up"], val, trueVal, scale * lepSFMUp * genWeight);
                  //this->fillResponse(responses["mEff_Down"], val, trueVal, scale * lepSFMDn * genWeight);
                }

	      
	      
              if(hasLHE && pdfAndAlphaSWeights.at(0))
                {
                  // fill once for each scale variation
                  float nominalWeightScaleNorm = nominalWeight / scaleWeights.at(0);
                  for(size_t ind = 0; ind < 9; ++ind)
                    this->fillResponse(scaleResponses.at(ind), val, trueVal,
                                       nominalWeightScaleNorm * scaleWeights.at(ind));
                  /* float nominalWeightScaleNorm = nominalWeight / scaleWeights.at(0);
                  for(size_t ind = 0; ind < scaleIndicesWeCareAbout.size(); ++ind)
                    this->fillResponse(scaleResponses.at(ind), val, trueVal,
                                       nominalWeightScaleNorm * scaleWeights.at(scaleIndicesWeCareAbout.at(ind))); */

	
                  float nominalWeightPDFNorm = nominalWeight / pdfAndAlphaSWeights.at(0); //nominal PDF not included for 2016
                  for(size_t ind = 0; ind <= lastLHEind-9; ++ind)
                    this->fillResponse(scaleResponses.at(ind+9), val, trueVal,
                                        nominalWeightScaleNorm* pdfAndAlphaSWeights.at(ind));
                                       //nominalWeightPDFNorm * pdfAndAlphaSWeights.at(ind));
                                       //use nominalWeightScaleNorm here since the 9th index is not nominal PDF for 2016, but 0th always scale nominal

                  // the last two items in the PDF weight vector are alpha_S variations
                 //In our ntuplizer uwvv codes, scaleWeight corresponds to lhe weights index 0 to 8, and pdf weights 9 to 9999 (used until 111). 

		  // fill the 3-D histogram with one response for each PDF variation
                  //float nominalWeightPDFNorm = nominalWeight / pdfAndAlphaSWeights.at(0);
                  
		  //for(size_t ind = 0; ind < nPDFVariations; ++ind)
		  //{this->fillResponse(pdfResponses, val, trueVal, ind,
                  //                     nominalWeightPDFNorm * pdfAndAlphaSWeights.at(ind));}

                  // the last two items in the PDF weight vector are alpha_S variations, the following fills are not used since we already have them in scaleResponses above
                  //this->fillResponse(responses["alphaS_up"], val, trueVal,
                  //                   nominalWeightPDFNorm * pdfAndAlphaSWeights.at(iAlphaSUp));
                  //this->fillResponse(responses["alphaS_dn"], val, trueVal,
                  //                   nominalWeightPDFNorm * pdfAndAlphaSWeights.at(iAlphaSDn));
                }

              // changes to jet scale/resolution actually change numbers
              //if(isJetVar)
              //  {
              //    if(this->selectEvent("jer_up"))
              //      this->fillResponse(responses["jer_up"], this->getEventResponse("jer_up"),
              //                         trueVal, nominalWeight);
              //    if(this->selectEvent("jer_dn"))
              //      this->fillResponse(responses["jer_dn"], this->getEventResponse("jer_dn"),
              //                         trueVal, nominalWeight);
              //    if(this->selectEvent("jes_up"))
              //      this->fillResponse(responses["jes_up"], this->getEventResponse("jes_up"),
              //                         trueVal, nominalWeight);
              //    if(this->selectEvent("jes_dn"))
              //      this->fillResponse(responses["jes_dn"], this->getEventResponse("jes_dn"),
              //                         trueVal, nominalWeight);
              //  }
            }
        }
    }

  recoTree.reset(); // gone -- don't use any more

  Map<Str,Str> systTreesNeeded;
  //if(hasE && !skipSyst)
  //  {
  //    if(systFileNames.find("eScaleUp") != systFileNames.end())
  //      systTreesNeeded["eScale_up"] = "eScaleUp";
  //    if(systFileNames.find("eScaleDn") != systFileNames.end())
  //      systTreesNeeded["eScale_dn"] = "eScaleDn";
  //    if(systFileNames.find("eRhoResUp") != systFileNames.end())
  //      systTreesNeeded["eRhoRes_up"] = "eRhoResUp";
  //    if(systFileNames.find("eRhoResDn") != systFileNames.end())
  //      systTreesNeeded["eRhoRes_dn"] = "eRhoResDn";
  //    if(systFileNames.find("ePhiResUp") != systFileNames.end())
  //      systTreesNeeded["ePhiRes_up"] = "ePhiResUp";
  //  }
  //if(hasMu && !skipSyst)
  //  {
  //    if(systFileNames.find("mClosureUp") != systFileNames.end())
  //      systTreesNeeded["mClosure_up"] = "mClosureUp";
  //    if(systFileNames.find("mClosureDn") != systFileNames.end())
  //      systTreesNeeded["mClosure_dn"] = "mClosureDn";
  //  }

  // systematics requiring other ntuples
  for(auto& treeInfo : systTreesNeeded)
    {
      const Str& systName = treeInfo.first;
      const Str& treeName = treeInfo.second;

      UPtr<TChain> t = UPtr<TChain>(new TChain((getChannel()+"/ntuple").c_str(),
                                               ("chain_"+getVar() + "_" + getChannel()).c_str()));
      for(const auto& fn : systFileNames[treeName])
        t->Add(fn.c_str());

      setCommonBranches(*t, objects);
      this->setRecoBranches(*t, objects);

      for(size_t row = 0; row < size_t(std::abs(t->GetEntries())); ++row)
        {
          t->GetEntry(row);
          if (year==2016 ||year==2017){genWeight*=L1prefiringWeight;}

          float puWt = (doPUWt ? ::getContentFromHist(puWeightHists.at(""), truePU) : 1.);
          //float puWt = 1.0;
          float lepSF = this->getLepSF(objects);

           // for(size_t i = 0; i < objects.size(); ++i)
           //   {
           //     if(objects.at(i)[0] == "e")
           //       {
           //         lepSF *= ::getContentFromHist(leptonSFHists["e"][""], lPt.at(i), std::abs(lEta.at(i)));
           //       }
           //     if(objects.at(i)[0] == "m")
           //       {
           //         lepSF *= ::getContentFromHist(leptonSFHists["m"][""], lPt.at(i), lEta.at(i));
           //       }
           //   }

          auto iTrue = trueVals->find(evt);

          if(iTrue == trueVals->end())
            continue;

          const T& trueVal = iTrue->second;

          if(this->selectEvent(systName))
            {
              this->fillResponse(responses[systName], this->getEventResponse(),
                                 trueVal, scale * puWt * lepSF * genWeight);
            }
        }
    } // new tree disappears here

  trueVals.reset();
}


template<typename T>
bool ResponseMatrixMakerBase<T>::selectTrueEvent(float mZ1, float mZ2) const
{
  return mZ1 > 60. && mZ1 < 120. && mZ2 > 60. && mZ2 < 120.;
}


template<typename T>
float *
ResponseMatrixMakerBase<T>::getmZPtr(TChain& t,
                                     const Str& obj1, const Str& obj2,
                                     float& maybeUseThis) const
{
  Str bName = obj1+"_"+obj2+"_"+"Mass";

  if(!t.GetBranch(bName.c_str())->GetAddress())
    t.SetBranchAddress(bName.c_str(), &maybeUseThis);

  return (float* const)(t.GetBranch(bName.c_str())->GetAddress());
}


template<typename T>
float ResponseMatrixMakerBase<T>::getLepSF(const Vec<Str>& objects,
                                           float eSyst, float mSyst)
{
  float out = 1.;
  //This functionality gets overwritten by UseSFHist version of this function
  //so hopefully we don't need to access the SFs and errors directly from the trees.
  //for(size_t i = 0; i < objects.size(); ++i)
  //  {
  //    if(objects.at(i).compare(0,1,"e") == 0)
  //      out *= lSF.at(i) + eSyst * lSFErr.at(i);
  //    else if(objects.at(i).compare(0,1,"m") == 0)
  //      out *= lSF.at(i) + mSyst * lSFErr.at(i);
  //  }

  return out;
}


template<typename T>
void ResponseMatrixMakerBase<T>::setCommonBranches(TChain& t, const Vec<Str>& objects)
{
  // if(lPt.size() != objects.size())
  //   lPt.resize(objects.size());
  // if(lEta.size() != objects.size())
  //   lEta.resize(objects.size());

  //Since we provide SF files at the time of running selectors and unfolding script
  //so hopefully we don't need to access the SFs and errors directly from the trees.-TEST
  //if(lSF.size() != objects.size())
  //  lSF.resize(objects.size());
  //if(lSFErr.size() != objects.size())
  //  lSFErr.resize(objects.size());

  //for(size_t i = 0; i < objects.size(); ++i)
  //  {
  //    // t.SetBranchAddress((objects[i]+"Pt").c_str(), &(lPt[i]));
  //    // t.SetBranchAddress((objects[i]+"Eta").c_str(), &(lEta[i]));
  //    t.SetBranchAddress((objects[i]+"EffScaleFactor").c_str(), &(lSF[i]));
  //    t.SetBranchAddress((objects[i]+"EffScaleFactorError").c_str(), &(lSFErr[i]));
  //  }
  // hopefully we don't need to access these SF branches from the trees- TEST
  t.SetBranchAddress("nTruePU", &truePU);
  t.SetBranchAddress("evt", &evt);
  t.SetBranchAddress("genWeight", &genWeight);
  t.SetBranchAddress("jetPUSFmulfac", &jetPUSFmulfac);
  
    if (year ==2016 || year==2017){
    t.SetBranchAddress("L1prefiringWeight", &L1prefiringWeight);}

    if(channel.find("eeee") != Str::npos){    
    t.SetBranchAddress("e1SIP3D", &l1SIP3D);
    t.SetBranchAddress("e2SIP3D", &l2SIP3D);
    t.SetBranchAddress("e3SIP3D", &l3SIP3D);
    t.SetBranchAddress("e4SIP3D", &l4SIP3D);
    
    t.SetBranchAddress("e1ZZTightIDNoVtx", &l1IsTight);
    t.SetBranchAddress("e1ZZIso", &l1IsIso);
    
    t.SetBranchAddress("e2ZZTightIDNoVtx", &l2IsTight);
    t.SetBranchAddress("e2ZZIso", &l2IsIso);
    
    t.SetBranchAddress("e3ZZTightIDNoVtx", &l3IsTight);
    t.SetBranchAddress("e3ZZIso", &l3IsIso);
    
    t.SetBranchAddress("e4ZZTightIDNoVtx", &l4IsTight);
    t.SetBranchAddress("e4ZZIso", &l4IsIso);

    t.SetBranchAddress("e1PdgId", &l1PdgId);
    t.SetBranchAddress("e2PdgId", &l2PdgId);
    t.SetBranchAddress("e3PdgId", &l3PdgId);
    t.SetBranchAddress("e4PdgId", &l4PdgId);
    }
    

    if(channel.find("eemm") != Str::npos){    
    t.SetBranchAddress("e1SIP3D", &l1SIP3D);
    t.SetBranchAddress("e2SIP3D", &l2SIP3D);
    t.SetBranchAddress("m1SIP3D", &l3SIP3D);
    t.SetBranchAddress("m2SIP3D", &l4SIP3D);
    
    t.SetBranchAddress("e1ZZTightIDNoVtx", &l1IsTight);
    t.SetBranchAddress("e1ZZIso", &l1IsIso);
    
    t.SetBranchAddress("e2ZZTightIDNoVtx", &l2IsTight);
    t.SetBranchAddress("e2ZZIso", &l2IsIso);

    t.SetBranchAddress("m1ZZIso", &l3IsIso);
    if (year == 2016 || year == 2017)
    {
      t.SetBranchAddress("m1PASTightIDNoVtx", &l3IsTight);
      t.SetBranchAddress("m2PASTightIDNoVtx", &l4IsTight);
    }
    else
    {
      t.SetBranchAddress("m1ZZTightIDNoVtx", &l3IsTight);
      t.SetBranchAddress("m2ZZTightIDNoVtx", &l4IsTight);
    }
    t.SetBranchAddress("m2ZZIso", &l4IsIso);

    t.SetBranchAddress("e1PdgId", &l1PdgId);
    t.SetBranchAddress("e2PdgId", &l2PdgId);
    t.SetBranchAddress("m1PdgId", &l3PdgId);
    t.SetBranchAddress("m2PdgId", &l4PdgId);
    }
    

    if(channel.find("mmmm") != Str::npos){    
    t.SetBranchAddress("m1SIP3D", &l1SIP3D);
    t.SetBranchAddress("m2SIP3D", &l2SIP3D);
    t.SetBranchAddress("m3SIP3D", &l3SIP3D);
    t.SetBranchAddress("m4SIP3D", &l4SIP3D);
    
     if (year == 2016 || year == 2017)
    {
      t.SetBranchAddress("m1PASTightIDNoVtx", &l1IsTight);
      t.SetBranchAddress("m2PASTightIDNoVtx", &l2IsTight);
      t.SetBranchAddress("m3PASTightIDNoVtx", &l3IsTight);
      t.SetBranchAddress("m4PASTightIDNoVtx", &l4IsTight);
    }
    else
    {
      t.SetBranchAddress("m1ZZTightIDNoVtx", &l1IsTight);
      t.SetBranchAddress("m2ZZTightIDNoVtx", &l2IsTight);
      t.SetBranchAddress("m3ZZTightIDNoVtx", &l3IsTight);
      t.SetBranchAddress("m4ZZTightIDNoVtx", &l4IsTight);
    }
    t.SetBranchAddress("m1ZZIso", &l1IsIso);
    t.SetBranchAddress("m2ZZIso", &l2IsIso);
    t.SetBranchAddress("m3ZZIso", &l3IsIso);
    t.SetBranchAddress("m4ZZIso", &l4IsIso);

    t.SetBranchAddress("m1PdgId", &l1PdgId);
    t.SetBranchAddress("m2PdgId", &l2PdgId);
    t.SetBranchAddress("m3PdgId", &l3PdgId);
    t.SetBranchAddress("m4PdgId", &l4PdgId);
    }
  


}


template<typename T>
BranchValueResponseMatrixMaker<T>::BranchValueResponseMatrixMaker(const Str& channel, const Str& varName,
                                                                  const Vec<float>& binning) :
  SimpleValueResponseMatrixMakerBase<T>(channel, varName, binning)
{;}


template<typename T>
UPtr<UMap<size_t, T> >
BranchValueResponseMatrixMaker<T>::getTrueValues(TChain& t,
                                                 const Vec<Str>& objects,
                                                 const Str& syst) const
{
  UPtr<UMap<size_t, T> > out(new UMap<size_t, T>());

  unsigned long long trueEvt;
  T val;
  t.SetBranchAddress("evt", &trueEvt);
  t.SetBranchAddress(this->getVar().c_str(), &val);

  float _doNotUse1 = 0.;
  float _doNotUse2 = 0.;
  float *const mZ1 = this->getmZPtr(t, objects.at(0), objects.at(1),
                                    _doNotUse1);
  float *const mZ2 = this->getmZPtr(t, objects.at(2), objects.at(3),
                                    _doNotUse2);

  for(size_t row = 0; row < size_t(std::abs(t.GetEntries())); ++row)
    {
      t.GetEntry(row);

      if(!this->selectTrueEvent(*mZ1,*mZ2))
        continue;

      (*out)[trueEvt] = val;
    }

  return std::move(out);
}


template<typename T> void
BranchValueResponseMatrixMaker<T>::setRecoBranches(TChain& t, const Vec<Str>& objects)
{

  if(this->getVar().find("Mass") != Str::npos){
    t.SetBranchAddress("Mass", &value);}
    //std::cout<<"SetRecoBranches called========================================================"<<std::endl;}
  else if (this->getVar().find("nJets") == Str::npos && this->getVar().find("mjj") == Str::npos){t.SetBranchAddress(this->getVar().c_str(), &value);} //nJets and mjj already set elsewhere
}


template<typename T> T
BranchValueResponseMatrixMaker<T>::getEventResponse(const Str& opt) const
{
  return value;
}


template<class R>
AbsValueResponseMatrixMaker<R>::AbsValueResponseMatrixMaker(const Str& channel,
                                                            const Str& varName,
                                                            const Vec<float>& binning) :
  R(channel, varName, binning)
{;}


template<class R>
UPtr<UMap<size_t, typename R::ValType> >
AbsValueResponseMatrixMaker<R>::getTrueValues(TChain& trueTree,
                                              const Vec<Str>& objects,
                                              const Str& syst) const
{
  UPtr<UMap<size_t, typename R::ValType> > out = R::getTrueValues(trueTree,
                                                                  objects,
                                                                  syst);
  for(auto& val : *out)
    this->doAbs(val.second);

  return std::move(out);
}


template<class R>
typename R::ValType
AbsValueResponseMatrixMaker<R>::getEventResponse(const Str& option) const
{
  typename R::ValType out = R::getEventResponse(option);
  this->doAbs(out);
  return out;
}


template<class R>
typename R::ValType&
AbsValueResponseMatrixMaker<R>::doAbs(typename R::ValType& val) const
{
  val = std::abs(val);
  return val;
}


template<typename T>
JetBranchResponseMatrixMakerBase<T>::JetBranchResponseMatrixMakerBase(const Str& channel,
                                                                      const Str& varName,
                                                                      const Vec<float>& binning) :
  BranchValueResponseMatrixMaker<T>(channel, varName, binning)
{;}


template<typename T>
T JetBranchResponseMatrixMakerBase<T>::getEventResponse(const Str& syst) const
{
  Vec<float>* tmp_jvec = valuevecfloat;
  if(syst.find("jes_up") != Str::npos){
    tmp_jvec = value_jesUp;}
  if(syst.find("jes_dn") != Str::npos){
    tmp_jvec = value_jesDn;}
  if(syst.find("jer_up") != Str::npos){
    tmp_jvec = value_jerUp;}
  if(syst.find("jer_dn") != Str::npos){
    tmp_jvec = value_jerDn;}

  if(this->getVar()=="jetPt[0]"){
    //std::cout<<"response value"<<valuevecfloat->at(0)<<std::endl;
    return tmp_jvec->at(0);}
  else if(this->getVar()=="jetPt[1]"){
    //std::cout<<"response value"<<valuevecfloat->at(0)<<std::endl;
    return tmp_jvec->at(1);}
  else if(this->getVar()=="jetEta[0]"){
    //std::cout<<"response value"<<valuevecfloat->at(0)<<std::endl;
    return tmp_jvec->at(0);}
  else if(this->getVar()=="jetEta[1]"){
    //std::cout<<"response value"<<valuevecfloat->at(0)<<std::endl;
    return tmp_jvec->at(1);}
  else if(this->getVar()=="absjetEta[0]"){
    //std::cout<<"response value"<<valuevecfloat->at(0)<<std::endl;
    return std::abs(tmp_jvec->at(0));}
  else if(this->getVar()=="absjetEta[1]"){
    //std::cout<<"response value"<<valuevecfloat->at(0)<<std::endl;
    return std::abs(tmp_jvec->at(1));}
  else if(this->getVar()=="dEtajj"){
    //std::cout<<"response value"<<valuevecfloat->at(0)<<std::endl;
    return std::abs(tmp_jvec->at(1)-tmp_jvec->at(0));}

  else{return BranchValueResponseMatrixMaker<T>::getEventResponse(syst);}
  //return BranchValueResponseMatrixMaker<T>::getEventResponse(syst);
}


template<typename T>
void
JetBranchResponseMatrixMakerBase<T>::setRecoBranches(TChain& t,
                                                     const Vec<Str>& objects)
{

  my_doNotUse1 = 0.;
  my_doNotUse2 = 0.;
  my_mZ1 = this->getmZPtr(t, objects.at(0), objects.at(1),my_doNotUse1);
  my_mZ2 = this->getmZPtr(t, objects.at(2), objects.at(3),my_doNotUse2);

  //std::cout<<this->getVar()<<(this->getVar()=="jetPt0")<<std::endl;
  if(this->getVar().find("jetPt") != Str::npos){

    t.SetBranchAddress("jetPt", &valuevecfloat);
    t.SetBranchAddress("jetPt_jesUp", &value_jesUp);
    t.SetBranchAddress("jetPt_jesDown", &value_jesDn);
    t.SetBranchAddress("jetPt_jerUp", &value_jerUp);
    t.SetBranchAddress("jetPt_jerDown", &value_jerDn);}

  else if(this->getVar().find("jetEta") != Str::npos || this->getVar()=="dEtajj"){

    t.SetBranchAddress("jetEta", &valuevecfloat);
    t.SetBranchAddress("jetEta_jesUp", &value_jesUp);
    t.SetBranchAddress("jetEta_jesDown", &value_jesDn);
    t.SetBranchAddress("jetEta_jerUp", &value_jerUp);
    t.SetBranchAddress("jetEta_jerDown", &value_jerDn);}


  else{BranchValueResponseMatrixMaker<T>::setRecoBranches(t, objects);

    //t.SetBranchAddress((this->getVar()+"_jesUp").c_str(), &value_jesUp);    //value_jesUp etc. are defined as vectors here, but they are not actually used yet so shouldn't matter.
    //t.SetBranchAddress((this->getVar()+"_jesDown").c_str(), &value_jesDn);
    //t.SetBranchAddress((this->getVar()+"_jerUp").c_str(), &value_jerUp);
    //t.SetBranchAddress((this->getVar()+"_jerDown").c_str(), &value_jerDn);
  }
}

template<class T>
DijetBranchResponseMatrixMaker<T>::DijetBranchResponseMatrixMaker(const Str& channel,
                                                               const Str& varName,
                                                               const Vec<float>& binning) :
  JetBranchResponseMatrixMakerBase<T>(channel, varName, binning)
{;}

template<class T>
UPtr<UMap<size_t, T> >
DijetBranchResponseMatrixMaker<T>::getTrueValues(TChain& trueTree,
                                              const Vec<Str>& objects,
                                              const Str& syst) const
{
  UPtr<UMap<size_t, T> > out(new UMap<size_t, T>());

  unsigned long long trueEvt;
  T trueVal;
  Vec<float>* trueValvecfloat=NULL;
  unsigned int trueNJets;
  float trueMjj;
  trueTree.SetBranchAddress("evt", &trueEvt);
  if(this->getVar().find("jetPt") != Str::npos){trueTree.SetBranchAddress("jetPt", &trueValvecfloat);}
  else if(this->getVar().find("jetEta") != Str::npos ||this->getVar()=="dEtajj"){trueTree.SetBranchAddress("jetEta", &trueValvecfloat);}
  else if (this->getVar().find("Mass") != Str::npos)
    {trueTree.SetBranchAddress("Mass", &trueVal);} //All Massnj variables use the Mass address, with nJets selection applied later below
      //std::cout<<"===========================================================================================Mass name corrected for:"<<this->getVar().c_str()<<std::endl;}
  else if (this->getVar().find("nJets") == Str::npos && this->getVar().find("mjj") == Str::npos){trueTree.SetBranchAddress(this->getVar().c_str(), &trueVal);}
  trueTree.SetBranchAddress("nJets", &trueNJets);
  trueTree.SetBranchAddress("mjj", &trueMjj);

  float _doNotUse1 = 0.;
  float _doNotUse2 = 0.;
  float *const mZ1 = this->getmZPtr(trueTree, objects.at(0), objects.at(1),
                                    _doNotUse1);
  float *const mZ2 = this->getmZPtr(trueTree, objects.at(2), objects.at(3),
                                    _doNotUse2);

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      if(!this->selectTrueEvent(*mZ1,*mZ2) && this->getVar().find("Full") == Str::npos)
        {continue;}

      if(this->getVar()=="Mass"){(*out)[trueEvt] = trueVal;}
      else if(this->getVar()=="Mass0j" && trueNJets == 0){(*out)[trueEvt] = trueVal;}
      else if(this->getVar()=="Mass1j" && trueNJets == 1){(*out)[trueEvt] = trueVal;}
	//std::cout<<trueEvt<<" "<<trueVal<<std::endl; }
      else if(this->getVar()=="Mass2j" && trueNJets == 2){(*out)[trueEvt] = trueVal;}
      else if(this->getVar()=="Mass3j" && trueNJets == 3){(*out)[trueEvt] = trueVal;}
      else if(this->getVar()=="Mass34j" && trueNJets >= 3){(*out)[trueEvt] = trueVal;}
      else if(this->getVar()=="Mass4j" && trueNJets >= 4){(*out)[trueEvt] = trueVal;}

       if(this->getVar()=="MassFull"){(*out)[trueEvt] = trueVal;}
      else if(this->getVar()=="Mass0jFull" && trueNJets == 0){(*out)[trueEvt] = trueVal;}
      else if(this->getVar()=="Mass1jFull" && trueNJets == 1){(*out)[trueEvt] = trueVal;}
	//std::cout<<trueEvt<<" "<<trueVal<<std::endl; }
      else if(this->getVar()=="Mass2jFull" && trueNJets == 2){(*out)[trueEvt] = trueVal;}
      else if(this->getVar()=="Mass3jFull" && trueNJets == 3){(*out)[trueEvt] = trueVal;}
      else if(this->getVar()=="Mass34jFull" && trueNJets >= 3){(*out)[trueEvt] = trueVal;}
      else if(this->getVar()=="Mass4jFull" && trueNJets >= 4){(*out)[trueEvt] = trueVal;}

      else if(this->getVar()=="nJets"){(*out)[trueEvt] = trueNJets;}
      else if(this->getVar()=="mjj" && trueNJets >= 2){(*out)[trueEvt] = trueMjj;}
      else if(this->getVar()=="dEtajj" && trueNJets >= 2){(*out)[trueEvt] = std::abs(trueValvecfloat->at(0)-trueValvecfloat->at(1));}
      else if (trueNJets >= 1 && this->getVar()=="jetPt[0]"){(*out)[trueEvt] = trueValvecfloat->at(0);}
      else if (trueNJets >= 2 && this->getVar()=="jetPt[1]"){(*out)[trueEvt] = trueValvecfloat->at(1);}
      else if (trueNJets >= 1 && this->getVar()=="jetEta[0]"){(*out)[trueEvt] = trueValvecfloat->at(0);}
      else if (trueNJets >= 2 && this->getVar()=="jetEta[1]"){(*out)[trueEvt] = trueValvecfloat->at(1);}
      else if (trueNJets >= 1 && this->getVar()=="absjetEta[0]"){(*out)[trueEvt] = std::abs(trueValvecfloat->at(0));}
      else if (trueNJets >= 2 && this->getVar()=="absjetEta[1]"){(*out)[trueEvt] = std::abs(trueValvecfloat->at(1));}
      //else {(*out)[trueEvt] = trueVal;}

      //if(trueNJets >= 1){
	//std::cout<<trueNJets<<":"<<trueValvecfloat->size()<<std::endl;
	//if(this->getVar()=="jetPt0"){(*out)[trueEvt] = trueValvecfloat->at(0);}
	//std::cout<<"selected true values:"<<trueValvecfloat->at(0)<<std::endl;}
        //else{(*out)[trueEvt] = trueVal;}}
    }

  return std::move(out); //move performed on return automatically so this is not actually needed?
}

template<class T>
void
DijetBranchResponseMatrixMaker<T>::setRecoBranches(TChain& t, const Vec<Str>& objects)
{
  JetBranchResponseMatrixMakerBase<T>::setRecoBranches(t, objects);
  
  t.SetBranchAddress("nJets", &nJets);
  t.SetBranchAddress("mjj", &mjj);
  //t.SetBranchAddress("Mass", &Mass);
  t.SetBranchAddress("nJets_jesUp", &nJets_jesUp);
  t.SetBranchAddress("nJets_jesDown", &nJets_jesDn);
  t.SetBranchAddress("nJets_jerUp", &nJets_jerUp);
  t.SetBranchAddress("nJets_jerDown", &nJets_jerDn);
  t.SetBranchAddress("mjj_jesUp", &mjj_jesUp);
  t.SetBranchAddress("mjj_jesDown", &mjj_jesDn);
  t.SetBranchAddress("mjj_jerUp", &mjj_jerUp);
  t.SetBranchAddress("mjj_jerDown", &mjj_jerDn);
}

template<typename T>
T DijetBranchResponseMatrixMaker<T>::getEventResponse(const Str& syst) const
{
  //std::cout<<"getEventResponse overwrite success======================================="<<std::endl;
  
  if(this->getVar()=="nJets")
  {
    unsigned int tmp_val = nJets;
    if(syst.find("jes_up") != Str::npos){
    tmp_val = nJets_jesUp;}
  else if(syst.find("jes_dn") != Str::npos){
    tmp_val= nJets_jesDn;}
  else if(syst.find("jer_up") != Str::npos){
    tmp_val = nJets_jerUp;}
  else if(syst.find("jer_dn") != Str::npos){
    tmp_val= nJets_jerDn;}
  
  return tmp_val;
  }

  else if(this->getVar()=="mjj")
  {
    float tmp_val = mjj;
    if(syst.find("jes_up") != Str::npos){
    tmp_val = mjj_jesUp;}
  else if(syst.find("jes_dn") != Str::npos){
    tmp_val= mjj_jesDn;}
  else if(syst.find("jer_up") != Str::npos){
    tmp_val = mjj_jerUp;}
  else if(syst.find("jer_dn") != Str::npos){
    tmp_val= mjj_jerDn;}
  
  return tmp_val;
  }
   // return mjj;}
  else{return JetBranchResponseMatrixMakerBase<T>::getEventResponse(syst);}
  //return BranchValueResponseMatrixMaker<T>::getEventResponse(syst);
}

//template<typename T>
//void 
//DijetBranchResponseMatrixMaker<T>::fillResponse(TH2D& h, const Vec<float>& val, const Vec<float>& trueVal, float w) const
//{

//  h.Fill(val.at(0), trueVal.at(0), w); //currently only used to handle jetPt case                                                                                                                    


//}

template<class T>
bool
DijetBranchResponseMatrixMaker<T>::selectEvent(const Str& syst) const
{
  float mZ1=*(JetBranchResponseMatrixMakerBase<T>::my_mZ1);
  float mZ2=*(JetBranchResponseMatrixMakerBase<T>::my_mZ2);
  bool  mass_sel = false;
  bool  mass_selFull = false;
  bool  mass_selBase = mZ1 > 60. && mZ1 < 120. && mZ2 > 60. && mZ2 < 120.;
  bool  tight_sel = (this->tightZ1Leptons() && this->tightZ2Leptons());
  
  if (this->year ==2016 || this->year==2017){
    bool sip_sel= this->l1SIP3D < 4.0 && this->l2SIP3D < 4.0 && this->l3SIP3D < 4.0 && this->l4SIP3D < 4.0;
    mass_sel = mass_selBase && tight_sel && sip_sel;
    mass_selFull = tight_sel && sip_sel;
  }
  if (this->year ==2018) {
    mass_sel = mass_selBase && tight_sel;
    mass_selFull = tight_sel;

  }
  //std::cout<<"Confirm mZ1,mZ2:  "<<mZ1<<" "<<mZ2<<std::endl;
  unsigned int tmp_nJets = nJets;
  if(syst.empty()){tmp_nJets = nJets;}
  else if(syst.find("jes_up") != Str::npos)
    {tmp_nJets = nJets_jesUp;}
  else if(syst.find("jes_dn") != Str::npos)
    {tmp_nJets = nJets_jesDn;} 
  else if(syst.find("jer_up") != Str::npos)
    {tmp_nJets = nJets_jerUp ;}
  else if(syst.find("jer_dn") != Str::npos)
    {tmp_nJets = nJets_jerDn ;}

    //std::cout<<"===========================Selection with no syst entered==============="<<std::endl; 
    //std::cout<<this->getVar().c_str()<<" "<<this->getVar()<<"  "<<(this->getVar().c_str()==this->getVar())<<std::endl;
    //std::cout<<(this->getVar().c_str()=="jetEta[0]")<<std::endl;
    //std::cout<<(this->getVar()=="jetEta[0]")<<std::endl;
    //std::cout<<"==================size unequal========================"<<nJets<<" vs "<< this->valuevecfloat->size()<<std::endl;

    if (this->getVar() == "Mass0j") {return tmp_nJets == 0 && mass_sel;}
    if (this->getVar() == "Mass1j") {return tmp_nJets == 1 && mass_sel;}
    if (this->getVar() == "Mass2j") {return tmp_nJets == 2 && mass_sel;}
    if (this->getVar() == "Mass3j") {return tmp_nJets == 3 && mass_sel;}
    if (this->getVar() == "Mass34j") {return tmp_nJets >= 3 && mass_sel;}
    if (this->getVar() == "Mass4j") {return tmp_nJets >= 4 && mass_sel;}
    if (this->getVar() == "nJets") {return mass_sel;}
    if (this->getVar() == "mjj") {return tmp_nJets >= 2 && mass_sel;}
    if (this->getVar() == "dEtajj") {return tmp_nJets >= 2 && mass_sel;}
    if (this->getVar() == "Mass") {return  mass_sel;}// && mjj>100 && Mass>180;
    if (this->getVar() == "jetPt[0]") {return tmp_nJets >= 1 && mass_sel;}
    if (this->getVar() == "jetPt[1]") {return tmp_nJets >= 2 && mass_sel;}
    if (this->getVar() == "jetEta[0]") {return tmp_nJets >= 1 && mass_sel;}
    if (this->getVar() == "jetEta[1]") {return tmp_nJets >= 2 && mass_sel;}
    if (this->getVar() == "absjetEta[0]") {return tmp_nJets >= 1 && mass_sel;}
    if (this->getVar() == "absjetEta[1]") {return tmp_nJets >= 2 && mass_sel;}

    if (this->getVar() == "MassFull") {return  mass_selFull;}
    if (this->getVar() == "Mass0jFull") {return tmp_nJets == 0 && mass_selFull;}
    if (this->getVar() == "Mass1jFull") {return tmp_nJets == 1 && mass_selFull;}
    if (this->getVar() == "Mass2jFull") {return tmp_nJets == 2 && mass_selFull;}
    if (this->getVar() == "Mass3jFull") {return tmp_nJets == 3 && mass_selFull;}
    if (this->getVar() == "Mass34jFull") {return tmp_nJets >= 3 && mass_selFull;}
    if (this->getVar() == "Mass4jFull") {return tmp_nJets >= 4 && mass_selFull;}
  

  

  std::cout<<"======================NO MATCH in selection!====================================="<<std::endl;
  return mass_sel;
}
  //return nJets >= 1 && mass_sel;// && mjj>100 && Mass>180;
  //}

testJets::testJets(const Str& channel,const Str& varName,const Vec<float>& binning):
  DijetBranchResponseMatrixMaker<float>(channel, varName, binning)
{;}


SelectedZResponseMatrixMakerBase::SelectedZResponseMatrixMakerBase(const Str& channel,
                                                                   const Str& varName,
                                                                   const Vec<float>& binning,
                                                                   const Str& compVar) :
  SimpleValueResponseMatrixMakerBase<float>(channel, varName, binning),
  z1CompVarName(constructZVarName(channel, compVar, true)),
  z2CompVarName(constructZVarName(channel, compVar, false)),
  z1RespVarName(constructZVarName(channel, varName, true)),
  z2RespVarName(constructZVarName(channel, varName, false)),
  compIsResp(z1CompVarName == z1RespVarName || z2CompVarName == z2RespVarName)
{;}


UPtr<UMap<size_t, float> >
SelectedZResponseMatrixMakerBase::getTrueValues(TChain& trueTree,
                                                const Vec<Str>& objects,
                                                const Str& syst) const
{
  UPtr<UMap<size_t, float> > out(new UMap<size_t, float>());

  unsigned long long trueEvt;
  float z1CompTrue;
  float z2CompTrue;
  float z1RespTrueValue;
  float z2RespTrueValue;
  float* z1RespTruePtr;
  float* z2RespTruePtr;

  trueTree.SetBranchAddress("evt", &trueEvt);
  trueTree.SetBranchAddress(z1CompVarName.c_str(), &z1CompTrue);
  trueTree.SetBranchAddress(z2CompVarName.c_str(), &z2CompTrue);

  float _doNotUse1 = 0.;
  float _doNotUse2 = 0.;
  float *const mZ1 = this->getmZPtr(trueTree, objects.at(0), objects.at(1),
                                    _doNotUse1);
  float *const mZ2 = this->getmZPtr(trueTree, objects.at(2), objects.at(3),
                                    _doNotUse2);

  if(compIsResp)
    {
      z1RespTruePtr = &z1CompTrue;
      z2RespTruePtr = &z2CompTrue;
    }
  else
    {
      trueTree.SetBranchAddress(z1RespVarName.c_str(), &z1RespTrueValue);
      trueTree.SetBranchAddress(z2RespVarName.c_str(), &z2RespTrueValue);
      z1RespTruePtr = &z1RespTrueValue;
      z2RespTruePtr = &z2RespTrueValue;
    }

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      if(!this->selectTrueEvent(*mZ1,*mZ2))
        continue;

      (*out)[trueEvt] = (this->z1IsBetter(z1CompTrue, z2CompTrue) ?
                         *z1RespTruePtr :
                         *z2RespTruePtr);
    }

  return std::move(out);
}


void
SelectedZResponseMatrixMakerBase::setRecoBranches(TChain& t,
                                                  const Vec<Str>& objects)
{
  t.SetBranchAddress(z1CompVarName.c_str(), &z1CompVar);
  t.SetBranchAddress(z2CompVarName.c_str(), &z2CompVar);
  if(compIsResp)
    {
      z1RespVarPtr = &z1CompVar;
      z2RespVarPtr = &z2CompVar;
    }
  else
    {
      t.SetBranchAddress(z1RespVarName.c_str(), &z1RespVarValue);
      t.SetBranchAddress(z2RespVarName.c_str(), &z2RespVarValue);
      z1RespVarPtr = &z1RespVarValue;
      z2RespVarPtr = &z2RespVarValue;
    }
}


float
SelectedZResponseMatrixMakerBase::getEventResponse(const Str& option) const
{
  if(this->z1IsBetter(z1CompVar, z2CompVar))
    return *z1RespVarPtr;
  return *z2RespVarPtr;
}


Str
SelectedZResponseMatrixMakerBase::constructZVarName(const Str& channel,
                                                    const Str& var,
                                                    bool wantZ1) const
{
  if(channel == "mmmm")
    {
      if(wantZ1)
        return "m1_m2_" + var;
      return "m3_m4_" + var;
    }

  if(wantZ1)
    return "e1_e2_" + var;

  if(channel == "eemm")
    return "m1_m2_" + var;

  return "e3_e4_" + var;
}


Z1ByMassResponseMatrixMaker::Z1ByMassResponseMatrixMaker(const Str& channel,
                                                         const Str& varName,
                                                         const Vec<float>& binning) :
  SelectedZResponseMatrixMakerBase(channel, varName, binning, "Mass"),
  comparisonVarName("Mass")
{;}


bool Z1ByMassResponseMatrixMaker::z1IsBetter(const float z1Comp,
                                             const float z2Comp) const
{
  return std::abs(z1Comp - 91.1876) < std::abs(z2Comp - 91.1876);
}


Z2ByMassResponseMatrixMaker::Z2ByMassResponseMatrixMaker(const Str& channel,
                                                         const Str& varName,
                                                         const Vec<float>& binning) :
  SelectedZResponseMatrixMakerBase(channel, varName, binning, "Mass"),
  comparisonVarName("Mass")
{;}


bool Z2ByMassResponseMatrixMaker::z1IsBetter(const float z1Comp,
                                             const float z2Comp) const
{
  return std::abs(z1Comp - 91.1876) > std::abs(z2Comp - 91.1876);
}


Z1ByPtResponseMatrixMaker::Z1ByPtResponseMatrixMaker(const Str& channel,
                                                     const Str& varName,
                                                     const Vec<float>& binning) :
  SelectedZResponseMatrixMakerBase(channel, varName, binning, "Pt"),
  comparisonVarName("Pt")
{;}


bool Z1ByPtResponseMatrixMaker::z1IsBetter(const float z1Comp,
                                           const float z2Comp) const
{
  return z1Comp > z2Comp;
}


Z2ByPtResponseMatrixMaker::Z2ByPtResponseMatrixMaker(const Str& channel,
                                                     const Str& varName,
                                                     const Vec<float>& binning) :
  SelectedZResponseMatrixMakerBase(channel, varName, binning, "Pt"),
  comparisonVarName("Pt")
{;}


bool Z2ByPtResponseMatrixMaker::z1IsBetter(const float z1Comp,
                                           const float z2Comp) const
{
  return z1Comp < z2Comp;
}


ZZCompositeResponseMatrixMakerBase::ZZCompositeResponseMatrixMakerBase(const Str& channel,
                                                                       const Str& varName,
                                                                       const Vec<float>& binning) :
  SimpleValueResponseMatrixMakerBase<float>(channel, varName, binning),
  z1VarName(constructZVarName(channel, varName, true)),
  z2VarName(constructZVarName(channel, varName, false))
{;}


UPtr<UMap<size_t, float> >
ZZCompositeResponseMatrixMakerBase::getTrueValues(TChain& trueTree,
                                                  const Vec<Str>& objects,
                                                  const Str& syst) const
{
  UPtr<UMap<size_t, float> > out(new UMap<size_t, float>());

  unsigned long long trueEvt;
  float z1VarTrue;
  float z2VarTrue;

  trueTree.SetBranchAddress("evt", &trueEvt);
  trueTree.SetBranchAddress(z1VarName.c_str(), &z1VarTrue);
  trueTree.SetBranchAddress(z2VarName.c_str(), &z2VarTrue);

  float _doNotUse1 = 0.;
  float _doNotUse2 = 0.;
  float *const mZ1 = this->getmZPtr(trueTree, objects.at(0), objects.at(1),
                                    _doNotUse1);
  float *const mZ2 = this->getmZPtr(trueTree, objects.at(2), objects.at(3),
                                    _doNotUse2);

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      if(!this->selectTrueEvent(*mZ1,*mZ2))
        continue;

      (*out)[trueEvt] = this->calculateZZVar(z1VarTrue, z2VarTrue);
    }

  return std::move(out);
}


void
ZZCompositeResponseMatrixMakerBase::setRecoBranches(TChain& t,
                                                    const Vec<Str>& objects)
{
  t.SetBranchAddress(z1VarName.c_str(), &z1InputVar);
  t.SetBranchAddress(z2VarName.c_str(), &z2InputVar);
}


float
ZZCompositeResponseMatrixMakerBase::getEventResponse(const Str& option) const
{
  return this->calculateZZVar(z1InputVar, z2InputVar);
}


Str
ZZCompositeResponseMatrixMakerBase::constructZVarName(const Str& channel,
                                                      const Str& var,
                                                      bool wantZ1) const
{
  if(channel == "mmmm")
    {
      if(wantZ1)
        return "m1_m2_" + var;
      return "m3_m4_" + var;
    }

  if(wantZ1)
    return "e1_e2_" + var;

  if(channel == "eemm")
    return "m1_m2_" + var;

  return "e3_e4_" + var;
}


ZZDeltaPhiResponseMatrixMaker::ZZDeltaPhiResponseMatrixMaker(const Str& channel,
                                                             const Str& varName,
                                                             const Vec<float>& binning) :
  ZZCompositeResponseMatrixMakerBase(channel, "Phi", binning)
{;}


float ZZDeltaPhiResponseMatrixMaker::calculateZZVar(float z1Phi,
                                                    float z2Phi) const
{
  return ::deltaPhi(z1Phi, z2Phi);
}


ZZDeltaRResponseMatrixMaker::ZZDeltaRResponseMatrixMaker(const Str& channel,
                                                         const Str& varName,
                                                         const Vec<float>& binning) :
  SimpleValueResponseMatrixMakerBase<float>(channel, varName, binning),
  z1EtaBranchName(constructZVarName(channel, "Eta", true)),
  z2EtaBranchName(constructZVarName(channel, "Eta", false)),
  z1PhiBranchName(constructZVarName(channel, "Phi", true)),
  z2PhiBranchName(constructZVarName(channel, "Phi", false))
{;}


UPtr<UMap<size_t, float> >
ZZDeltaRResponseMatrixMaker::getTrueValues(TChain& trueTree,
                                           const Vec<Str>& objects,
                                           const Str& syst) const
{
  UPtr<UMap<size_t, float> > out(new UMap<size_t, float>());

  unsigned long long trueEvt;
  float z1EtaTrue;
  float z2EtaTrue;
  float z1PhiTrue;
  float z2PhiTrue;

  trueTree.SetBranchAddress("evt", &trueEvt);
  trueTree.SetBranchAddress(z1EtaBranchName.c_str(), &z1EtaTrue);
  trueTree.SetBranchAddress(z2EtaBranchName.c_str(), &z2EtaTrue);
  trueTree.SetBranchAddress(z1PhiBranchName.c_str(), &z1PhiTrue);
  trueTree.SetBranchAddress(z2PhiBranchName.c_str(), &z2PhiTrue);

  float _doNotUse1 = 0.;
  float _doNotUse2 = 0.;
  float *const mZ1 = this->getmZPtr(trueTree, objects.at(0), objects.at(1),
                                    _doNotUse1);
  float *const mZ2 = this->getmZPtr(trueTree, objects.at(2), objects.at(3),
                                    _doNotUse2);

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      if(!this->selectTrueEvent(*mZ1,*mZ2))
        continue;

      (*out)[trueEvt] = ::deltaR(z1EtaTrue, z1PhiTrue, z2EtaTrue, z2PhiTrue);
    }

  return std::move(out);
}

void
ZZDeltaRResponseMatrixMaker::setRecoBranches(TChain& t, const Vec<Str>& objects)
{
  t.SetBranchAddress(z1EtaBranchName.c_str(), &z1Eta);
  t.SetBranchAddress(z2EtaBranchName.c_str(), &z2Eta);
  t.SetBranchAddress(z1PhiBranchName.c_str(), &z1Phi);
  t.SetBranchAddress(z2PhiBranchName.c_str(), &z2Phi);
}


float
ZZDeltaRResponseMatrixMaker::getEventResponse(const Str& option) const
{
  return ::deltaR(z1Eta, z1Phi, z2Eta, z2Phi);
}


Str
ZZDeltaRResponseMatrixMaker::constructZVarName(const Str& channel,
                                               const Str& var,
                                               bool wantZ1) const
{
  if(channel == "mmmm")
    {
      if(wantZ1)
        return "m1_m2_" + var;
      return "m3_m4_" + var;
    }

  if(wantZ1)
    return "e1_e2_" + var;

  if(channel == "eemm")
    return "m1_m2_" + var;

  return "e3_e4_" + var;
}


template<typename T>
MultiBranchResponseMatrixMakerBase<T>::MultiBranchResponseMatrixMakerBase(const Str& channel,
                                                                          const Str& varName,
                                                                          const Vec<float>& binning) :
  VectorValueResponseMatrixMakerBase<T>(channel, varName, binning)
{;}


template<typename T>
UPtr<UMap<size_t, Vec<T> > >
MultiBranchResponseMatrixMakerBase<T>::getTrueValues(TChain& trueTree,
                                                     const Vec<Str>& objects,
                                                     const Str& syst) const
{
  UPtr<UMap<size_t, Vec<T> > > out(new UMap<size_t, Vec<T> >());

  const Vec<Str> varNames = this->constructVarNames(this->getChannel(),
                                                    this->getVar());

  unsigned long long trueEvt;
  Vec<T> trueValues(varNames.size());

  trueTree.SetBranchAddress("evt", &trueEvt);
  for(size_t i = 0; i < varNames.size(); ++i)
    trueTree.SetBranchAddress(varNames[i].c_str(), &trueValues[i]);

  float _doNotUse1 = 0.;
  float _doNotUse2 = 0.;
  float *const mZ1 = this->getmZPtr(trueTree, objects.at(0), objects.at(1),
                                    _doNotUse1);
  float *const mZ2 = this->getmZPtr(trueTree, objects.at(2), objects.at(3),
                                    _doNotUse2);

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      if(!this->selectTrueEvent(*mZ1,*mZ2))
        continue;

      // copies elements so should be safe
      (*out)[trueEvt] = trueValues;
    }

  return std::move(out);
}


template<typename T>
void
MultiBranchResponseMatrixMakerBase<T>::setRecoBranches(TChain& t,
                                                       const Vec<Str>& objects)
{
  const Vec<Str> varNames = this->constructVarNames(this->getChannel(),
                                                    this->getVar());

  values.resize(varNames.size());

  // this is technically dangerous, but I think it should be ok as long as
  // values never changes size between calls to this function
  for(size_t i = 0; i < varNames.size(); ++i)
    t.SetBranchAddress(varNames[i].c_str(), &values[i]);
}


template<typename T>
Vec<T>
MultiBranchResponseMatrixMakerBase<T>::getEventResponse(const Str& option) const
{
  // copies elements so should be safe
  return values;
}


template<typename T>
Vec<Str>
MultiBranchResponseMatrixMakerBase<T>::constructVarNames(const Str& channel,
                                                         const Str& var) const
{
  Vec<Str> objects = this->constructObjectNames(channel);

  Vec<Str> out;
  for(const auto& obj : objects)
    out.push_back(obj + var);

  return out;
}


AllLeptonBranchResponseMatrixMaker::AllLeptonBranchResponseMatrixMaker(const Str& channel,
                                                                       const Str& varName,
                                                                       const Vec<float>& binning) :
  MultiBranchResponseMatrixMakerBase<float>(channel, varName, binning)
{;}


Vec<Str>
AllLeptonBranchResponseMatrixMaker::constructObjectNames(const Str& channel) const
{
  if(channel == "mmmm")
    return Vec<Str>({"m1","m2","m3","m4"});
  else if(channel == "eemm")
    return Vec<Str>({"m1","m2","e1","e2"});
  else
    return Vec<Str>({"e1","e2","e3","e4"});
}


BothZsBranchResponseMatrixMaker::BothZsBranchResponseMatrixMaker(const Str& channel,
                                                                 const Str& varName,
                                                                 const Vec<float>& binning) :
  MultiBranchResponseMatrixMakerBase<float>(channel, varName, binning)
{;}


Vec<Str>
BothZsBranchResponseMatrixMaker::constructObjectNames(const Str& channel) const
{
  if(channel == "mmmm")
    return Vec<Str>({"m1_m2_","m3_m4_"});
  else if(channel == "eemm")
    return Vec<Str>({"m1_m2_","e1_e2_"});
  else
    return Vec<Str>({"e1_e2_","e3_e4_"});
}


LeptonMaxBranchResponseMatrixMaker::LeptonMaxBranchResponseMatrixMaker(const Str& channel,
                                                                       const Str& varName,
                                                                       const Vec<float>& binning) :
  SimpleValueResponseMatrixMakerBase<float>(channel, varName, binning),
  varNames(constructVarNames(channel, varName))
{;}


UPtr<UMap<size_t, float> >
LeptonMaxBranchResponseMatrixMaker::getTrueValues(TChain& trueTree,
                                                  const Vec<Str>& objects,
                                                  const Str& syst) const
{
  UPtr<UMap<size_t, float> > out(new UMap<size_t, float>());

  unsigned long long trueEvt;
  Vec<float> trueValues(varNames.size());

  trueTree.SetBranchAddress("evt", &trueEvt);
  for(size_t i = 0; i < varNames.size(); ++i)
    trueTree.SetBranchAddress(varNames[i].c_str(), &trueValues[i]);

  float _doNotUse1 = 0.;
  float _doNotUse2 = 0.;
  float *const mZ1 = this->getmZPtr(trueTree, objects.at(0), objects.at(1),
                                    _doNotUse1);
  float *const mZ2 = this->getmZPtr(trueTree, objects.at(2), objects.at(3),
                                    _doNotUse2);

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      if(!this->selectTrueEvent(*mZ1,*mZ2))
        continue;

      float max = -999999999.;
      for(const auto& val : trueValues)
        {
          if(val > max)
            max = val;
        }

      (*out)[trueEvt] = max;
    }

  return std::move(out);
}


void
LeptonMaxBranchResponseMatrixMaker::setRecoBranches(TChain& t,
                                                    const Vec<Str>& objects)
{
  values.resize(varNames.size());

  // this is technically dangerous, but I think it should be ok as long as
  // values never changes size between calls to this function
  for(size_t i = 0; i < varNames.size(); ++i)
    t.SetBranchAddress(varNames[i].c_str(), &values[i]);
}


float
LeptonMaxBranchResponseMatrixMaker::getEventResponse(const Str& option) const
{
  float max = -999999999.;
  for(const auto& val : values)
    {
      if(val > max)
        max = val;
    }

  return max;
}


Vec<Str>
LeptonMaxBranchResponseMatrixMaker::constructVarNames(const Str& channel,
                                                      const Str& var) const
{
  Vec<Str> objects;

  if(channel == "mmmm")
    objects = {"m1","m2","m3","m4"};
  else if(channel == "eemm")
    objects = {"m1","m2","e1","e2"};
  else
    objects = {"e1","e2","e3","e4"};

  Vec<Str> out;
  for(const auto& obj : objects)
    out.push_back(obj + var);

  return out;
}


template<typename T, size_t _N>
NthJetResponseMatrixMaker<T,_N>::NthJetResponseMatrixMaker(const Str& channel, const Str& varName,
                                                           const Vec<float>& binning) :
  SimpleValueResponseMatrixMakerBase<T>(channel, varName, binning)
{
  allJetValues[""] = &allJetValues_object;
  allJetValues["jes_up"] = &allJetValues_jesUp_object;
  allJetValues["jes_dn"] = &allJetValues_jesDn_object;
  allJetValues["jer_up"] = &allJetValues_jerUp_object;
  allJetValues["jer_dn"] = &allJetValues_jerDn_object;

  gROOT->ProcessLine("#include<vector>");
}


template<typename T, size_t _N>
UPtr<UMap<size_t, T> >
NthJetResponseMatrixMaker<T,_N>::getTrueValues(TChain& trueTree,
                                               const Vec<Str>& objects,
                                               const Str& syst) const
{
  UPtr<UMap<size_t, T> > out(new UMap<size_t, T>());

  unsigned long long trueEvt;

  // some pointer bullshit to make ROOT happy
  Vec<T> allJetsTrue_object;
  Vec<T>* allJetsTrue = &allJetsTrue_object;
  trueTree.SetBranchAddress("evt", &trueEvt);
  trueTree.SetBranchAddress(this->getVar().c_str(), &allJetsTrue);

  float _doNotUse1 = 0.;
  float _doNotUse2 = 0.;
  float *const mZ1 = this->getmZPtr(trueTree, objects.at(0), objects.at(1),
                                    _doNotUse1);
  float *const mZ2 = this->getmZPtr(trueTree, objects.at(2), objects.at(3),
                                    _doNotUse2);

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      if(!this->selectTrueEvent(*mZ1,*mZ2))
        continue;

      if(allJetsTrue->size() > _N)
        (*out)[trueEvt] = allJetsTrue->at(_N);
    }

  return std::move(out);
}


template<typename T, size_t _N> void
NthJetResponseMatrixMaker<T,_N>::setRecoBranches(TChain& t, const Vec<Str>& objects)
{
  allJetValues.clear();

  t.SetBranchAddress(this->getVar().c_str(), &allJetValues[""]);
  t.SetBranchAddress((this->getVar()+"_jesUp").c_str(), &allJetValues["jes_up"]);
  t.SetBranchAddress((this->getVar()+"_jesDown").c_str(), &allJetValues["jes_dn"]);
  t.SetBranchAddress((this->getVar()+"_jerUp").c_str(), &allJetValues["jer_up"]);
  t.SetBranchAddress((this->getVar()+"_jerDown").c_str(), &allJetValues["jer_dn"]);
}


template<typename T, size_t _N> T
NthJetResponseMatrixMaker<T,_N>::getEventResponse(const Str& syst) const
{
  // size will already be checked

  auto iVal = allJetValues.find(syst);

  if(iVal == allJetValues.end())
    return allJetValues.at("")->at(_N);

  return iVal->second->at(_N);
}


template<typename T, size_t _N> bool
NthJetResponseMatrixMaker<T,_N>::selectEvent(const Str& syst) const
{
  auto iVal = allJetValues.find(syst);

  if(iVal == allJetValues.end())
    return allJetValues.at("")->size() > _N;

  return iVal->second->size() > _N;
}


template<class R>
UseSFHists<R>::UseSFHists(const Str& channel, const Str& varName,
                          const Vec<float>& binning) :
  R(channel, varName, binning)
{
  // set some defaults so things don't break
  hEleSelSF.reset(new TH2F("eSelSFDefault", "", 1, 0., 1000., 1, 0., 1000.));
  hEleSelSF->SetBinContent(4, 1.);
  hEleSelGapSF.reset(new TH2F("eSelGapSFDefault", "", 1, 0., 1000., 1, 0., 1000.));
  hEleSelGapSF->SetBinContent(4, 1.);
  hEleLowRecoSF.reset(new TH2F("eLowRecoSFDefault", "", 1, 0., 1000., 1, 0., 1000.));
  hEleLowRecoSF->SetBinContent(4, 1.);
  hEleRecoSF.reset(new TH2F("eRecoSFDefault", "", 1, 0., 1000., 1, 0., 1000.));
  hEleRecoSF->SetBinContent(4, 1.);
  hMuSF.reset(new TH2F("mSFDefault", "", 1, 0., 1000., 1, 0., 1000.));
  hMuSF->SetBinContent(4, 1.);
  hMuSFErr.reset(new TH2F("mSFErrDefault", "", 1, 0., 1000., 1, 0., 1000.));
}


template<class R> void
UseSFHists<R>::registerElectronSelectionSFHist(const TH2F& h)
{
  hEleSelSF.reset((TH2F*)h.Clone());
}


template<class R> void
UseSFHists<R>::registerElectronSelectionGapSFHist(const TH2F& h)
{
  hEleSelGapSF.reset((TH2F*)h.Clone());
}


template<class R> void
UseSFHists<R>::registerElectronRecoSFHist(const TH2F& h)
{
  hEleRecoSF.reset((TH2F*)h.Clone());
}

template<class R> void
UseSFHists<R>::registerElectronLowRecoSFHist(const TH2F& h)
{
  hEleLowRecoSF.reset((TH2F*)h.Clone());
}

template<class R> void
UseSFHists<R>::registerMuonSFHist(const TH2F& h)
{
  hMuSF.reset((TH2F*)h.Clone());
}


template<class R> void
UseSFHists<R>::registerMuonSFErrorHist(const TH2F& h)
{
  hMuSFErr.reset((TH2F*)h.Clone());
}


template<class R> void
UseSFHists<R>::setupOneLepton(TChain& t, const Str& obj,
                              float*& ptPtr, float& ptVal,
                              float*& etaPtr, float& etaVal,
                              bool*& isGapPtr, bool& isGapVal)
{
  Str ptName = obj+"Pt";
  if(t.GetBranch(ptName.c_str())->GetAddress())
    ptPtr = (float*)(t.GetBranch(ptName.c_str())->GetAddress());
  else
    {
      ptPtr = &ptVal;
      t.SetBranchAddress(ptName.c_str(), ptPtr);
    }

  Str etaName = "";
  Str isGapName = obj+"IsGap";
  if(obj.compare(0,1,"e") == 0)
    {
      etaName = obj+"SCEta";

      if(t.GetBranch(isGapName.c_str())->GetAddress())
        isGapPtr = (bool*)(t.GetBranch(isGapName.c_str())->GetAddress());
      else
        {
          isGapPtr = &isGapVal;
          t.SetBranchAddress(isGapName.c_str(), isGapPtr);
        }
    }
  else
    {
      isGapPtr = &isGapVal;
      isGapVal = false; // just in case

      etaName = obj+"Eta";
    }

  if(t.GetBranch(etaName.c_str())->GetAddress())
    etaPtr = (float*)(t.GetBranch(etaName.c_str())->GetAddress());
  else
    {
      etaPtr = &etaVal;
      t.SetBranchAddress(etaName.c_str(), etaPtr);
    }
}


template<class R> void
UseSFHists<R>::setRecoBranches(TChain& t, const Vec<Str>& objects)
{
  R::setRecoBranches(t, objects);

  lPtsSF.resize(4);
  lEtasSF.resize(4);
  lIsGapSF.resize(4);

  setupOneLepton(t, objects.at(0),
                 lPtsSF.at(0), l1PtSF_value,
                 lEtasSF.at(0), l1EtaSF_value,
                 lIsGapSF.at(0), l1IsGapSF_value);
  setupOneLepton(t, objects.at(1),
                 lPtsSF.at(1), l2PtSF_value,
                 lEtasSF.at(1), l2EtaSF_value,
                 lIsGapSF.at(1), l2IsGapSF_value);
  setupOneLepton(t, objects.at(2),
                 lPtsSF.at(2), l3PtSF_value,
                 lEtasSF.at(2), l3EtaSF_value,
                 lIsGapSF.at(2), l3IsGapSF_value);
  setupOneLepton(t, objects.at(3),
                 lPtsSF.at(3), l4PtSF_value,
                 lEtasSF.at(3), l4EtaSF_value,
                 lIsGapSF.at(3), l4IsGapSF_value);
}


template<class R> float
UseSFHists<R>::getLepSF(const Vec<Str>& leptons,
                        float eSyst, float mSyst)
{
  float out = 1.;

  for(size_t i = 0; i < leptons.size(); ++i)
    {
      if(leptons.at(i).compare(0,1,"e") == 0)
        {
          float sf = 1.;
          float err = 0.;

          if(*lIsGapSF.at(i))
            {
              sf *= ::getContentFromHist(*hEleSelGapSF, *lEtasSF.at(i),
                                         *lPtsSF.at(i));
              if(eSyst)
                err += ::getErrorFromHist(*hEleSelGapSF, *lEtasSF.at(i),
                                          *lPtsSF.at(i));
            }
          else
            {
              sf *= ::getContentFromHist(*hEleSelSF, *lEtasSF.at(i),
                                         *lPtsSF.at(i));
              if(eSyst)
                err += ::getErrorFromHist(*hEleSelSF, *lEtasSF.at(i),
                                          *lPtsSF.at(i));
            }
          if(*lPtsSF.at(i)<20.)
            { 
            sf *= ::getContentFromHist(*hEleLowRecoSF, *lEtasSF.at(i),
                                     *lPtsSF.at(i));
            if(eSyst)
              {
                float recoErr = ::getErrorFromHist(*hEleLowRecoSF,
                                                   *lEtasSF.at(i),
                                                   *lPtsSF.at(i));
                if(*lPtsSF.at(i) < 10.)
                  recoErr += 0.01;
                err = eSyst * std::sqrt(err*err + recoErr*recoErr);
              }
            }
          else 
          {
            sf *= ::getContentFromHist(*hEleRecoSF, *lEtasSF.at(i),
                                     *lPtsSF.at(i));
            if(eSyst)
              {
                float recoErr = ::getErrorFromHist(*hEleRecoSF,
                                                   *lEtasSF.at(i),
                                                   *lPtsSF.at(i));
                if(*lPtsSF.at(i) > 75.)
                  recoErr += 0.01;
                err = eSyst * std::sqrt(err*err + recoErr*recoErr);
              }
          }
          out *= (sf + err);
        }
      else if(leptons.at(i).compare(0,1,"m") == 0)
        {
          float sf = ::getContentFromHist(*hMuSF, *lEtasSF.at(i),
                                          *lPtsSF.at(i));
          if(mSyst)
              sf += mSyst * ::getContentFromHist(*hMuSFErr,
                                                 *lEtasSF.at(i),
                                                 *lPtsSF.at(i));

          out *= sf;
        }
    }

  return out;
}


template<class R>
RelaxGenZCuts<R>::RelaxGenZCuts(const Str& channel, const Str& varName,
                                const Vec<float>& binning) :
  R(channel, varName, binning)
{;}

typedef Vec<float> Vecfloat;
typedef SimpleValueResponseMatrixMakerBase<float> FloatResponseMatrixMakerBase;
typedef SimpleValueResponseMatrixMakerBase<Vecfloat> test1;
typedef BranchValueResponseMatrixMaker<float> FloatBranchResponseMatrixMaker;
typedef AbsValueResponseMatrixMaker<FloatBranchResponseMatrixMaker> AbsFloatBranchResponseMatrixMaker;
typedef BranchValueResponseMatrixMaker<unsigned int> UIntBranchResponseMatrixMaker;
typedef DijetBranchResponseMatrixMaker<float> dijetmaker;
typedef AbsValueResponseMatrixMaker<dijetmaker> AbsDijetBranchResponseMatrixMaker;
typedef JetBranchResponseMatrixMakerBase<unsigned int> JetUIntBranchResponseMatrixMaker;
typedef JetBranchResponseMatrixMakerBase<float> JetFloatBranchResponseMatrixMaker;
typedef AbsValueResponseMatrixMaker<ZZDeltaPhiResponseMatrixMaker> ZZAbsDeltaPhiResponseMatrixMaker;
typedef NthJetResponseMatrixMaker<float,0> FirstJetFloatResponseMatrixMaker;
typedef NthJetResponseMatrixMaker<float,1> SecondJetFloatResponseMatrixMaker;
typedef AbsValueResponseMatrixMaker<FirstJetFloatResponseMatrixMaker> FirstJetAbsFloatResponseMatrixMaker;
typedef AbsValueResponseMatrixMaker<SecondJetFloatResponseMatrixMaker> SecondJetAbsFloatResponseMatrixMaker;
//typedef DijetBranchResponseMatrixMaker<Vecfloat> testJet; 

typedef UseSFHists<FloatBranchResponseMatrixMaker>       SFHistFloatBranchResponseMatrixMaker;
typedef UseSFHists<AbsFloatBranchResponseMatrixMaker>    SFHistAbsFloatBranchResponseMatrixMaker;
typedef UseSFHists<UIntBranchResponseMatrixMaker>        SFHistUIntBranchResponseMatrixMaker;
typedef UseSFHists<testJets>                             SFHisttestJet;
//typedef UseSFHists<DijetBranchResponseMatrixMaker>       SFHistDijetBranchResponseMatrixMaker;
typedef UseSFHists<AbsDijetBranchResponseMatrixMaker>    SFHistAbsDijetBranchResponseMatrixMaker;
typedef UseSFHists<JetUIntBranchResponseMatrixMaker>     SFHistJetUIntBranchResponseMatrixMaker;
typedef UseSFHists<JetFloatBranchResponseMatrixMaker>    SFHistJetFloatBranchResponseMatrixMaker;
typedef UseSFHists<Z1ByMassResponseMatrixMaker>          SFHistZ1ByMassResponseMatrixMaker;
typedef UseSFHists<Z2ByMassResponseMatrixMaker>          SFHistZ2ByMassResponseMatrixMaker;
typedef UseSFHists<Z1ByPtResponseMatrixMaker>            SFHistZ1ByPtResponseMatrixMaker;
typedef UseSFHists<Z2ByPtResponseMatrixMaker>            SFHistZ2ByPtResponseMatrixMaker;
typedef UseSFHists<ZZDeltaRResponseMatrixMaker>          SFHistZZDeltaRResponseMatrixMaker;
typedef UseSFHists<ZZAbsDeltaPhiResponseMatrixMaker>     SFHistZZAbsDeltaPhiResponseMatrixMaker;
typedef UseSFHists<AllLeptonBranchResponseMatrixMaker>   SFHistAllLeptonBranchResponseMatrixMaker;
typedef UseSFHists<LeptonMaxBranchResponseMatrixMaker>   SFHistLeptonMaxBranchResponseMatrixMaker;
typedef UseSFHists<BothZsBranchResponseMatrixMaker>      SFHistBothZsBranchResponseMatrixMaker;
typedef UseSFHists<FirstJetFloatResponseMatrixMaker>     SFHistFirstJetFloatResponseMatrixMaker;
typedef UseSFHists<SecondJetFloatResponseMatrixMaker>    SFHistSecondJetFloatResponseMatrixMaker;
typedef UseSFHists<FirstJetAbsFloatResponseMatrixMaker>  SFHistFirstJetAbsFloatResponseMatrixMaker;
typedef UseSFHists<SecondJetAbsFloatResponseMatrixMaker> SFHistSecondJetAbsFloatResponseMatrixMaker;

typedef RelaxGenZCuts<FloatBranchResponseMatrixMaker> FullSpectrumFloatResponseMatrixMaker;
typedef UseSFHists<FullSpectrumFloatResponseMatrixMaker> SFHistFullSpectrumFloatResponseMatrixMaker;

#if defined(__ROOTCLING__)
#pragma link C++ class FloatBranchResponseMatrixMaker;
#pragma link C++ class AbsFloatBranchResponseMatrixMaker;
#pragma link C++ class UIntBranchResponseMatrixMaker;
//#pragma link C++ class DijetBranchResponseMatrixMaker;
#pragma link C++ class testJet;
#pragma link C++ class SFHisttestJet;
#pragma link C++ class AbsDijetBranchResponseMatrixMaker;
#pragma link C++ class JetUIntBranchResponseMatrixMaker;
#pragma link C++ class JetFloatBranchResponseMatrixMaker;
#pragma link C++ class Z1ByMassResponseMatrixMaker;
#pragma link C++ class Z2ByMassResponseMatrixMaker;
#pragma link C++ class Z1ByPtResponseMatrixMaker;
#pragma link C++ class Z2ByPtResponseMatrixMaker;
#pragma link C++ class ZZDeltaRResponseMatrixMaker;
#pragma link C++ class ZZAbsDeltaPhiResponseMatrixMaker;
#pragma link C++ class AllLeptonBranchResponseMatrixMaker;
#pragma link C++ class LeptonMaxBranchResponseMatrixMaker;
#pragma link C++ class BothZsBranchResponseMatrixMaker;
#pragma link C++ class FirstJetFloatResponseMatrixMaker;
#pragma link C++ class SecondJetFloatResponseMatrixMaker;
#pragma link C++ class FirstJetAbsFloatResponseMatrixMaker;
#pragma link C++ class SecondJetAbsFloatResponseMatrixMaker;

#pragma link C++ class SFHistFloatBranchResponseMatrixMaker;
#pragma link C++ class SFHistAbsFloatBranchResponseMatrixMaker;
#pragma link C++ class SFHistUIntBranchResponseMatrixMaker;
#pragma link C++ class SFHistDijetBranchResponseMatrixMaker;
#pragma link C++ class SFHistAbsDijetBranchResponseMatrixMaker;
#pragma link C++ class SFHistJetUIntBranchResponseMatrixMaker;
#pragma link C++ class SFHistJetFloatBranchResponseMatrixMaker;
#pragma link C++ class SFHistZ1ByMassResponseMatrixMaker;
#pragma link C++ class SFHistZ2ByMassResponseMatrixMaker;
#pragma link C++ class SFHistZ1ByPtResponseMatrixMaker;
#pragma link C++ class SFHistZ2ByPtResponseMatrixMaker;
#pragma link C++ class SFHistZZDeltaRResponseMatrixMaker;
#pragma link C++ class SFHistZZAbsDeltaPhiResponseMatrixMaker;
#pragma link C++ class SFHistAllLeptonBranchResponseMatrixMaker;
#pragma link C++ class SFHistLeptonMaxBranchResponseMatrixMaker;
#pragma link C++ class SFHistBothZsBranchResponseMatrixMaker;
#pragma link C++ class SFHistFirstJetFloatResponseMatrixMaker;
#pragma link C++ class SFHistSecondJetFloatResponseMatrixMaker;
#pragma link C++ class SFHistFirstJetAbsFloatResponseMatrixMaker;
#pragma link C++ class SFHistSecondJetAbsFloatResponseMatrixMaker;

#pragma link C++ class FullSpectrumFloatResponseMatrixMaker;
#endif
