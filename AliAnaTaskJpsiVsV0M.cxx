/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// \file   AliAnaTaskJpsiVsV0M.cxx
// \author Sarah Herrmann <sarah.herrmann@cern.ch>
//
// \brief This code makes basic event selection and muon track selection, and
//        fills 3 trees containing V0 multiplicity and muon information
//        for different event selections
//        Intended to work on the grid
// \date 10/10/23

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAnaTaskJpsiVsV0M.h"


#include "AliMuonTrackCuts.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliESDUtils.h"

class AliAnaTaskJpsiVsV0M;    // my analysis class

using namespace std;

ClassImp(AliAnaTaskJpsiVsV0M) // classimp: necessary for root

AliAnaTaskJpsiVsV0M::AliAnaTaskJpsiVsV0M() : AliAnalysisTaskSE(),
    fAOD(0), fOutputList(0),
    fTreeDiMu(0x0),
    fTreeSingleMu(0x0),
    fTreeINT7(0x0),
    fCentMethod("V0M"),
    fMuonTrackCuts(0),
    fTriggerMatchLevelMuon(0),
    fRunN(-1),
    fispileupspd(kFALSE),
    fv0multTot(-1.),
    fv0multcorr(-1.),
    fitsmult(-1.),
    fv0mpercentile(-1.),
    fv0apercentile(-1.),
    fv0cpercentile(-1.),
    fcl1percentile(-1.),
    fspdpercentile(-1.),
    fv0mpercentile7(-1.),
    fv0apercentile7(-1.),
    fv0cpercentile7(-1.),
    fcl1percentile7(-1.),
    fspdpercentile7(-1.),
    fUtils(0x0),
    fzvtx(0.),
    ftracksmu(0x0),
    etaMu(0.),
    phiMu(0.),
    ptMu(0.),
    fv0multTotINT7(0.),
    fv0multCorrINT7(0.),
    hV0MultTot(0),
    hV0MultCorr(0),
    nentries(0)
{
    for(Int_t j=0; j<64; j++)
    {
      fv0mult[j] = 0.0;
      fv0multSMu[j] = 0.0;

      hV0MultPerChannelNoMu[j]= NULL;
    }

}
//_____________________________________________________________________________
AliAnaTaskJpsiVsV0M::AliAnaTaskJpsiVsV0M(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0),
    fTreeDiMu(0x0),
    fTreeSingleMu(0x0),
    fTreeINT7(0x0),
    fCentMethod("V0M"),
    fMuonTrackCuts(new AliMuonTrackCuts("stdMuonCuts","stdMuonCuts")),
    fTriggerMatchLevelMuon(0),
    fRunN(-1),
    fispileupspd(kFALSE),
    fv0multTot(-1.),
    fv0multcorr(-1.),
    fitsmult(-1.),
    fv0mpercentile(-1.),
    fv0apercentile(-1.),
    fv0cpercentile(-1.),
    fcl1percentile(-1.),
    fspdpercentile(-1.),
    fv0mpercentile7(-1.),
    fv0apercentile7(-1.),
    fv0cpercentile7(-1.),
    fcl1percentile7(-1.),
    fspdpercentile7(-1.),
    fUtils(0x0),
    fzvtx(0.),
    ftracksmu(0x0),
    etaMu(0.),
    phiMu(0.),
    ptMu(0.),
    fv0multTotINT7(0.),
    fv0multCorrINT7(0.),
    hV0MultTot(0),
    hV0MultCorr(0),
    nentries(0)
{
    // constructor

    for(Int_t j=0; j<64; j++)
    {
      fv0mult[j] = 0.0;
      fv0multSMu[j] = 0.0;

      hV0MultPerChannelNoMu[j]= NULL;
    }

    fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
    fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

    fTriggerMatchLevelMuon=1;//change that if needs to be configured

    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it,
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the output of the analysis: in this case it's a list of histograms
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
    DefineOutput(2, TTree::Class());

    DefineOutput(3, TTree::Class());

    DefineOutput(4, TTree::Class());
}
//_____________________________________________________________________________
AliAnaTaskJpsiVsV0M::~AliAnaTaskJpsiVsV0M()
{
    // destructor

    delete fMuonTrackCuts;

    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnaTaskJpsiVsV0M::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you create the histograms that you want to use
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)

    //adding the histograms for the kINT7

    hV0MultTot = new TH1F("hV0MultTot", "hV0MultTot", 20000, 0, 10000);
    fOutputList->Add(hV0MultTot);

    hV0MultCorr = new TH1F("hV0MultCorr", "hV0MultCorr", 20000, 0, 10000);
    fOutputList->Add(hV0MultCorr);

    for(Int_t ichannel=0; ichannel<64; ichannel++)
    {
      hV0MultPerChannelNoMu[ichannel] = new TH1F(Form("hV0MultPerChannelNoMu_%d",ichannel),Form("hV0MultPerChannelNoMu_%d",ichannel),500,0,250);
      fOutputList->Add(hV0MultPerChannelNoMu[ichannel]);
    }

    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output



    ftracksmu = new TClonesArray("AliAODTrack",10);
    fUtils = new AliAnalysisUtils();

    fTreeDiMu = new TTree("DiMuonEvents","Dimuon events");
    SetTreeDiMu(fTreeDiMu);
    PostData(2,fTreeDiMu);//the output 2 is defined in the task call of the first constructor

    fTreeSingleMu = new TTree("SingleMuonEvents","SingleMuonEvents");
    SetTreeSingleMu(fTreeSingleMu);
    PostData(3,fTreeSingleMu);

    fTreeINT7 = new TTree("kINT7Events","kINT7Events");
    SetTreeINT7(fTreeINT7);
    PostData(4,fTreeINT7);
}
//_____________________________________________________________________________
void AliAnaTaskJpsiVsV0M::SetTreeDiMu(TTree *tr)
{
  tr->Branch("fRunN",&fRunN);
  //tr->Branch("fispileupspd",&fispileupspd);
  tr->Branch("fv0multTot",&fv0multTot);
  tr->Branch("fv0multcorr",&fv0multcorr);
  tr->Branch("fitsmult",&fitsmult);
  tr->Branch("fv0mpercentile",&fv0mpercentile);
  tr->Branch("fv0apercentile",&fv0apercentile);
  tr->Branch("fv0cpercentile",&fv0cpercentile);
  tr->Branch("fcl1percentile",&fcl1percentile);
  tr->Branch("fspdpercentile",&fspdpercentile);
  tr->Branch("fzvtx",&fzvtx);
  tr->Branch("fv0mult",fv0mult,"fv0mult[64]/F");
  tr->Branch("ftracksmu",&ftracksmu);
}

void AliAnaTaskJpsiVsV0M::SetTreeSingleMu(TTree *tr)
{
  tr->Branch("fv0multSMu",fv0multSMu,"fv0multSMu[64]/F");
  tr->Branch("etaMu",&etaMu);
  tr->Branch("phiMu",&phiMu);
  tr->Branch("ptMu",&ptMu);
}

void AliAnaTaskJpsiVsV0M::SetTreeINT7(TTree *tr)
{
  tr->Branch("fv0multTotINT7",&fv0multTotINT7,"fv0multTotINT7/F");
  tr->Branch("fv0multCorrINT7",&fv0multCorrINT7,"fv0multCorrINT7/F");
  tr->Branch("fv0mpercentile7",&fv0mpercentile7);
  tr->Branch("fv0apercentile7",&fv0apercentile7);
  tr->Branch("fv0cpercentile7",&fv0cpercentile7);
  tr->Branch("fcl1percentile7",&fcl1percentile7);
  tr->Branch("fspdpercentile7",&fspdpercentile7);
}
//==============================================================================
void AliAnaTaskJpsiVsV0M::GetAcceptedTracksMuonArm(AliAODEvent *aodEvent)
{

  // fills the array of muon tracks that pass the cuts

  ftracksmu->Clear();

  Int_t nTracks = aodEvent->GetNumberOfTracks();

  AliAODTrack *track = 0;

  for (Int_t iTrack=0; iTrack<nTracks; iTrack++)
  {
    track = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
    if(!track) AliFatal("Not a standard AOD");

    if (track->IsMuonTrack() && track->GetMatchTrigger()>=fTriggerMatchLevelMuon)
    {


      if (track->Eta() > -4.0 && track->Eta() < -2.5)
      {

        Double_t rabs = track->GetRAtAbsorberEnd();
        if (rabs > 17.6 && rabs < 89.5)
        {

          UInt_t mask = fMuonTrackCuts->GetSelectionMask(track);

          if (mask & AliMuonTrackCuts::kMuPdca)
          {
            new ((*ftracksmu)[ftracksmu->GetEntriesFast()]) AliAODTrack(*track);
          }
        }
      }//if (track->Eta() > -4.0 && track->Eta() < -2.5)
    }//if (track->IsMuonTrack() && track->GetMatchTrigger()>=fTriggerMatchLevelMuon)
  }//for (Int_t iTrack=0; iTrack<nTracks; iTrack++)
}

//==============================================================================
void AliAnaTaskJpsiVsV0M::NotifyRun()
{
  // Notify run, called automatically by the framework
  fMuonTrackCuts->SetRun(fInputHandler);
}

//==============================================================================

void AliAnaTaskJpsiVsV0M::UserExec(Option_t *)
{
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you
    // have access to the current event.
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file

    if(!fAOD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event


    fRunN = fAOD->GetRunNumber();//will get filled into the tree after

    Bool_t isSingleMuonEvent = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMuonSingleLowPt7);
    Bool_t iskINT7Event = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
    Bool_t isDiMuonUSEvent = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt7);//Is dimuon unlike sign event

    if (!(isSingleMuonEvent || iskINT7Event || isDiMuonUSEvent))
    {
      return;
    }


    //Get the centrality
    Double_t percentile;

    AliMultSelection *multSelection = (AliMultSelection *)fAOD->FindListObject("MultSelection");
    if(!multSelection) return;

    if (fUtils->IsPileUpSPD(fAOD))
    {
      //removing pileup
      return;
    }

    //other pileup selection: from Christiane's task

    if(fAOD->IsPileupFromSPDInMultBins())
    {
        //This event is pileUp from AOD
        return;
    }

    if(fUtils->IsPileUpEvent(fAOD))
    {
        //This event is pileUp
        return;
    }

    percentile = multSelection->GetMultiplicityPercentile(fCentMethod.Data());



    // ------- Vertex selection: standard
    const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
    if (spdVtx->GetNContributors()<=0) return;
    Double_t cov[6]={0};
    spdVtx->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (spdVtx->IsFromVertexerZ() && (zRes>0.25)) return;
    if ((spdVtx->GetZ()>10) || (spdVtx->GetZ()<-10)) return;


    //--------- Fill the muon track TClonesArray: next step
    GetAcceptedTracksMuonArm(fAOD);// fills the array of muon tracks ftracksmu
                                           //that pass the cuts
                                           //is necessary for all 3 trigger classes

    if (isDiMuonUSEvent)
    {
      //tracksmu is already filled with the muons I need


      fv0mpercentile = multSelection->GetMultiplicityPercentile("V0M");
      fv0apercentile = multSelection->GetMultiplicityPercentile("V0A");
      fv0cpercentile = multSelection->GetMultiplicityPercentile("V0C");
      fcl1percentile = multSelection->GetMultiplicityPercentile("SPDClusters");
      fspdpercentile = multSelection->GetMultiplicityPercentile("SPDTracklets");

      // pileup detection with SPD
      //fispileupspd   = fUtils->IsPileUpSPD(fAOD);

      //Multiplicity in the ITS
      fitsmult = ((AliVAODHeader*)fAOD->GetHeader())->GetNumberOfITSClusters(0)+
        ((AliVAODHeader*)fAOD->GetHeader())->GetNumberOfITSClusters(1);

      fzvtx = spdVtx->GetZ();//fzvtx contains the z position of the SPD vertex


      //-------------------------- V0 multiplicity -------------------------------
      AliAODVZERO *vzeroAOD = dynamic_cast<AliAODVZERO *>( dynamic_cast<AliAODEvent *>(fAOD)->GetVZEROData());
      Float_t V0AMult = vzeroAOD->GetMTotV0A();
      Float_t V0CMult = vzeroAOD->GetMTotV0C();
      fv0multTot=V0AMult+V0CMult;//total V0 multiplicity in this event

      //------------------------   V0M Correction
      Float_t vzeroMultACorr=V0AMult, vzeroMultCCorr=V0CMult;
      vzeroMultACorr = AliESDUtils::GetCorrV0A(V0AMult,fzvtx);
      vzeroMultCCorr = AliESDUtils::GetCorrV0C(V0CMult,fzvtx);
      fv0multcorr = vzeroMultACorr + vzeroMultCCorr; // corrected V0 multiplicity


      //---- Store V0 signal channel per channel
      for (Int_t iChannel=0; iChannel<64; iChannel++)
      {
        fv0mult[iChannel]= fAOD->GetVZEROData()->GetMultiplicity(iChannel);
      }

      fTreeDiMu->Fill();//fills all the branches with their different TClonesArray or elements

      PostData(2, fTreeDiMu);

    }//if (isDiMuonUSEvent)





    if (isSingleMuonEvent)
    {
      for (Int_t iChannel=0; iChannel<64; iChannel++)
      {
        fv0multSMu[iChannel]= fAOD->GetVZEROData()->GetMultiplicity(iChannel);
      }

      if (!(ftracksmu->GetEntriesFast()==1))
      {
        printf("There are %d muons in this event number %lld, but there should be just 1\n", ftracksmu->GetEntriesFast(), nentries);
        return;
        //note: happens sometimes with 2 and 0 muons
      }

      AliAODTrack *muon = (AliAODTrack*) ftracksmu->UncheckedAt(0);
      etaMu = muon->Eta();
      phiMu = muon->Phi();
      ptMu = muon->Pt();

      fTreeSingleMu->Fill();
      PostData(3, fTreeSingleMu);
    }//if (isSingleMuonEvent)

    if (iskINT7Event)
    {
      fv0mpercentile7 = multSelection->GetMultiplicityPercentile("V0M");
      fv0apercentile7 = multSelection->GetMultiplicityPercentile("V0A");
      fv0cpercentile7 = multSelection->GetMultiplicityPercentile("V0C");
      fcl1percentile7 = multSelection->GetMultiplicityPercentile("SPDClusters");
      fspdpercentile7 = multSelection->GetMultiplicityPercentile("SPDTracklets");

      //-------------------------- V0 multiplicity -------------------------------
      AliAODVZERO *vzeroAOD = dynamic_cast<AliAODVZERO *>( dynamic_cast<AliAODEvent *>(fAOD)->GetVZEROData());
      Float_t V0AMult = vzeroAOD->GetMTotV0A();
      Float_t V0CMult = vzeroAOD->GetMTotV0C();
      fv0multTotINT7=V0AMult+V0CMult;//total V0 multiplicity in this event

      //------------------------   V0M Correction
      Float_t vzeroMultACorr=V0AMult, vzeroMultCCorr=V0CMult;
      vzeroMultACorr = AliESDUtils::GetCorrV0A(V0AMult,fzvtx);
      vzeroMultCCorr = AliESDUtils::GetCorrV0C(V0CMult,fzvtx);
      fv0multCorrINT7 = vzeroMultACorr + vzeroMultCCorr; // corrected V0 multiplicity

      hV0MultTot->Fill(fv0multTotINT7);
      hV0MultCorr->Fill(fv0multCorrINT7);

      if (ftracksmu->GetEntriesFast()==0)
      {
        //fill the 64 histograms with the corresponding mult in V0 channels

        for (Int_t iChannel=0; iChannel<64; iChannel++)
        {
          hV0MultPerChannelNoMu[iChannel]->Fill(fAOD->GetVZEROData()->GetMultiplicity(iChannel));
        }

      }
      //postdata for the histograms
      PostData(1, fOutputList);

      fTreeINT7->Fill();
      PostData(4, fTreeINT7);
    }//if (iskINT7Event)

    nentries++;
}//UserExec


//_____________________________________________________________________________
void AliAnaTaskJpsiVsV0M::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
