#include "SimpleAnalysis/AnalysisClass.h"

DefineAnalysis(StrongDilepton2019)

void StrongDilepton2019::Init()
{
  addRegions({
              "SRZ_soft",
              "SRZ_hard",
              "SR_low",
              "SR_medium",
              "SR_high",
//              "SRC",
//              "SRC_MET"
            });

  addRegions({
              "CRZ",
              "CRFS_soft",
              "CRFS_hard",
              "CRT_soft",
              "CRT_hard",
              "CRZ_low",
              "CRZ_medium",
              "CRZ_high",
              "CRFS_low",
              "CRFS_medium",
              "CRFS_high",
//              "CRC",
//              "CRC_MET"
            });

  addRegions({
              "VRZ",
              "VRT_soft",
              "VRT_hard",
              "VRS_soft",
              "VRS_hard",
              "VRFS_soft",
              "VRFS_hard",
              "VR_WZ",
              "VR_ZZ",
              "VR_3L",
              "VR1_low",
              "VR1_medium",
              "VR1_high",
//              "VRC",
//              "VRC_MET"
            });

  addHistogram("SR_low_mll",15,12,312);
  addHistogram("SR_medium_mll",15,12,312);
  addHistogram("SR_high_mll",15,12,312);
  addHistogram("CRFS_low_mll",15,12,312);
  addHistogram("CRFS_medium_mll",15,12,312);
  addHistogram("CRFS_high_mll",15,12,312);
  addHistogram("VR1_low_mll",15,12,312);
  addHistogram("VR1_medium_mll",15,12,312);
  addHistogram("VR1_high_mll",15,12,312);

}

void StrongDilepton2019::ProcessEvent(AnalysisEvent *event)
{
    float gen_met      = event->getGenMET();
    float gen_ht       = event->getGenHT();
    int channel_number = event->getMCNumber();
    bool debug = false;
    if(debug) std::cout << "This is after gen-variables" << std::endl;
    auto baselineElectrons = filterCrack(event->getElectrons(10, 2.47, ELooseBLLH|ED0Sigma5|EZ05mm)); 
    auto baselineMuons     = event->getMuons(10, 2.5, MuMedium|MuZ05mm); 
    auto jets              = event->getJets(20., 2.8);  /// JVT cut not applied
    auto metVec            = event->getMET();
    if(debug) std::cout << "This is after baseine leptons" << std::endl;
    /// Standard SUSY overlap removal      
    auto radiusCalcMuon = [] (const AnalysisObject& muon, const AnalysisObject& ) { return 0.1 + 9.6/muon.Pt(); };
    auto radiusCalcEl   = [] (const AnalysisObject& electron, const AnalysisObject& ) { return 0.1 + 9.6/electron.Pt(); };
    
    jets               = overlapRemoval(jets, baselineElectrons, 0.2, NOT(BTag85MV2c10)); /// not entirely correct
    jets               = overlapRemoval(jets, baselineMuons, 0.4, LessThan3Tracks); 
    baselineElectrons  = overlapRemoval(baselineElectrons, jets, radiusCalcEl); 
    baselineMuons      = overlapRemoval(baselineMuons, jets, radiusCalcMuon); 
    baselineElectrons  = overlapRemoval(baselineElectrons, baselineMuons,0.01);  
    if(debug) std::cout << "This is after overlap removal" << std::endl;
    auto signalElectrons = filterObjects(baselineElectrons, 10, 2.0, EMediumLH|EIsoFCTight);  /// missing ECIDS
    auto signalMuons     = filterObjects(baselineMuons, 10, 2.5, MuD0Sigma3|MuIsoFCTightTrackOnly);

    // require signal jets to be 30 GeV
    auto signalJets      = filterObjects(jets, 30);
  
    // combine into signalLeptons for easy counting
    auto signalLeptons   = signalElectrons + signalMuons;
    sortObjectsByPt( signalLeptons );

    // get b-jets
    auto bjets = filterObjects(signalJets, 30., 2.5, BTag77MV2c10);

    int n_bjets       = bjets.size();
    int n_jets        = signalJets.size();
    int n_leptons     = signalLeptons.size();
    double met        = metVec.Et();

    // require at least 2 signal leptons
    if(n_leptons<2) return;
    
    if(debug) std::cout << "This is after n_leptons cut" << std::endl;
    
    // Leading leptons mll
    TLorentzVector dilepton( signalLeptons[0].Px()+signalLeptons[1].Px() , signalLeptons[0].Py()+signalLeptons[1].Py() , 
                           signalLeptons[0].Pz()+signalLeptons[1].Pz() , signalLeptons[0].E() +signalLeptons[1].E()  );
    
    if(debug) std::cout << "This is after dilepton assignment cut" << std::endl;
    float mll = dilepton.M();
    float pTll = dilepton.Pt();
    // inclusive - HT 
    float HT  = sumObjectsPt(signalJets);
    if(debug) std::cout << "This is after HT assignment" << std::endl;
    // dphimin between leading 4 signal jets and met
    float dphiMin2 = minDphi(metVec, signalJets, 2);
    if(debug) std::cout << "This is after dphiMin2 assignment" << std::endl;
    // Leading two leptons
    float mt2     = calcMT2(signalLeptons[0],signalLeptons[1],metVec);
    if(debug) std::cout << "This is after mT2 assignment" << std::endl;
    // 3rd leading lepton and met
    float mT3 = n_leptons>2 ? calcMT(signalLeptons[2], metVec) : 0.;
    if(debug) std::cout << "This is after mT3 assignment" << std::endl;
    // mm is channel 0; ee is channel 1; em is channel 2; me is channel 3
    int channel = -1;
    if (signalElectrons.size()==0) channel = 0;
    else if (signalMuons.size()==0) channel = 1;
    else if (n_leptons==2) channel = signalElectrons[0].Pt()>signalMuons[0].Pt()?2:3;
    else {
        // 3 or more leptons, need to check momenta
        if      (signalMuons.size()>1     && signalMuons[1].Pt()>signalElectrons[0].Pt()) channel = 0;
        else if (signalElectrons.size()>1 && signalElectrons[1].Pt()>signalMuons[0].Pt()) channel = 1;
        else channel = signalElectrons[0].Pt()>signalMuons[0].Pt()?2:3;
    }
    if(debug) std::cout << "This is after channel assignment" << std::endl;
    // Add OS requirement?
    if (signalLeptons[0].Pt()>50){
        //std::cout << met << " " << HT << " " << n_jets << " " << mll << " " << channel << " " << dphiMin2 << std::endl;
        // "Hard" lepton regions
        // On-shell signal regions
        if (met>400 && HT> 400 && n_jets>1 && 81<mll && mll<101 && channel<2 && dphiMin2>0.4) accept("SRZ_soft");
        if (met>200 && HT>1200 && n_jets>1 && 81<mll && mll<101 && channel<2 && dphiMin2>0.4) accept("SRZ_hard");
        // Edge signal regions
        if (met>250 && HT> 200 && n_jets>1 && mt2>70 && mll>12 && channel<2 && dphiMin2>0.4){ accept("SR_low"); fill("SR_low_mll",mll); }
        if (met>400 && HT> 400 && n_jets>1 && mt2>25 && mll>12 && channel<2 && dphiMin2>0.4){ accept("SR_medium"); fill("SR_medium_mll",mll); }
        if (met>200 && HT>1200 && n_jets>1 &&           mll>12 && channel<2 && dphiMin2>0.4){ accept("SR_high"); fill("SR_high_mll",mll); }
        // On-shell control regions
        if (met<60  && HT> 600 && n_jets>1 && 81<mll && mll<101 && channel<2 && dphiMin2>0.4) accept("CRZ");
        if (met>400 && HT> 400 && n_jets>1 && 61<mll && mll<121 && channel>1 && dphiMin2>0.4) accept("CRFS_soft");
        if (met>200 && HT>1200 && n_jets>1 && 61<mll && mll<121 && channel>1 && dphiMin2>0.4) accept("CRFS_hard");
        if (met>400 && HT> 400 && n_jets>1 && 45<mll && (81>mll || mll>101) && channel<2 && dphiMin2>0.4) accept("CRT_soft");
        if (met>200 && HT>1200 && n_jets>1 && 45<mll && (81>mll || mll>101) && channel<2 && dphiMin2>0.4) accept("CRT_hard");
        // Edge control regions
        if (met<60 && HT> 200 && n_jets>1 && mt2>70 && mll>12 && channel<2 && dphiMin2>0.4) accept("CRZ_low");
        if (met<60 && HT> 400 && n_jets>1 && mt2>25 && mll>12 && channel<2 && dphiMin2>0.4) accept("CRZ_medium");
        if (met<60 && HT>1200 && n_jets>1 &&           mll>12 && channel<2 && dphiMin2>0.4) accept("CRZ_high");
        if (met>250 && HT> 200 && n_jets>1 && mt2>70 && mll>12 && channel>1 && dphiMin2>0.4){ accept("CRFS_low"); fill("CRFS_low_mll",mll); }
        if (met>400 && HT> 400 && n_jets>1 && mt2>25 && mll>12 && channel>1 && dphiMin2>0.4){ accept("CRFS_medium"); fill("CRFS_medium_mll",mll); }
        if (met>200 && HT>1200 && n_jets>1 &&           mll>12 && channel>1 && dphiMin2>0.4){ accept("CRFS_high"); fill("CRFS_high_mll",mll); }
        // On-shell validation regions
        if (met<225 && HT> 600 && n_jets>1 && 81<mll && mll<101 && channel<2 && dphiMin2>0.4) accept("VRZ");
        if (met>100 && met<200 && HT> 400 && n_jets>1 && 45<mll && (81>mll || mll>101) && channel<2 && dphiMin2>0.4) accept("VRT_soft");
        if (met>100 && met<200 && HT>1200 && n_jets>1 && 45<mll && (81>mll || mll>101) && channel<2 && dphiMin2>0.4) accept("VRT_hard");
        if (met>100 && met<200 && HT> 400 && n_jets>1 && 81<mll && mll<101 && channel<2 && dphiMin2>0.4) accept("VRS_soft");
        if (met>100 && met<200 && HT>1200 && n_jets>1 && 81<mll && mll<101 && channel<2 && dphiMin2>0.4) accept("VRS_hard");
        if (met>100 && met<200 && HT> 400 && n_jets>1 && 61<mll && mll<121 && channel>1 && dphiMin2>0.4) accept("VRFS_soft");
        if (met>100 && met<200 && HT>1200 && n_jets>1 && 61<mll && mll<121 && channel>1 && dphiMin2>0.4) accept("VRFS_hard");
        // Diboson validation regions
        if (n_leptons==3 && met>100 && met<200 && mT3<100 && n_bjets==0) accept("VR_WZ");
        if (n_leptons==4 && met<100 && n_bjets==0) accept("VR_ZZ");
        if (n_leptons==3 && met> 60 && met<100 && HT>200 && n_jets>1 && 81<mll && mll<101 && dphiMin2>0.4) accept("VR_3L");
        // Edge validation regions
        if (met>100 && met<200 && HT> 200 && n_jets>1 && mt2>70 && mll>12 && channel<2 && dphiMin2>0.4){ accept("VR1_low"); fill("VR1_low_mll",mll); }
        if (met>100 && met<200 && HT> 400 && n_jets>1 && mt2>25 && mll>12 && channel<2 && dphiMin2>0.4){ accept("VR1_medium"); fill("VR1_medium_mll",mll); }
        if (met>100 && met<200 && HT>1200 && n_jets>1 &&           mll>12 && channel<2 && dphiMin2>0.4){ accept("VR1_high"); fill("VR1_high_mll",mll); }

  }
  if(debug) std::cout << "This is after getting signal cuts" << std::endl;

//  // Soft lepton regions -- need to think about lepton momenta here
//  if (met>250 && pTll<20 && n_jets>1 && mll>12 && channel<2){ accept("SRC"); }
//  if (met>500 && pTll<75 && n_jets>1 && mll>12 && channel<2){ accept("SRC_MET"); }
//  // Soft lepton control regions
//  if (met>250 && pTll<20 && n_jets>1 && mll>12 && channel>1){ accept("CRC"); }
//  if (met>500 && pTll<75 && n_jets>1 && mll>12 && channel>1){ accept("CRC_MET"); }
//  // Soft lepton validation regions
//  if (                                  mll>12 && channel<2){ accept("VRC"); }
//  if (                                  mll>12 && channel<2){ accept("VRC_MET"); }

//   ntupVar("gen_met", gen_met);
//   ntupVar("gen_ht", gen_ht);
//   ntupVar("channel_number", channel_number);
  ntupVar("channel", channel);

  ntupVar("met", met);
  ntupVar("HT", HT);
  ntupVar("mT2", mt2);
  ntupVar("mT3", mT3);
  ntupVar("mll", mll);
  ntupVar("pTll", pTll);
  ntupVar("dphiMin2", dphiMin2);

  ntupVar("n_jets", n_jets);
  ntupVar("n_bjets", n_bjets);
  ntupVar("n_electrons", static_cast<int>(signalElectrons.size()));
  ntupVar("n_muons", static_cast<int>(signalMuons.size()));
  ntupVar("n_leptons", n_leptons);

  return;
}
