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

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnaMCMultEpsilonTask.h"
#include "AliAODMCParticle.h"

#include "AliAODTrack.h"
#include "AliMultSelection.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliEventCuts.h"
#include "AliAODMCHeader.h"
#include <AliHeader.h>
#include "TCanvas.h"
#include "AliESDUtils.h"
#include "AliVertexingHFUtils.h"

class AliAnaMCMultEpsilonTask;

using namespace std;

ClassImp(AliAnaMCMultEpsilonTask) // classimp: necessary for root

AliAnaMCMultEpsilonTask::AliAnaMCMultEpsilonTask() : AliAnalysisTaskSE(),
    fAOD(0), fOutputList(0), fMCEvent(0), hNV0CNchZvtx(0), hNV0ANchZvtx(0), hNV0MNchZvtx(0), hNV0CNchZvtxInMSAcc(0), hNV0CCorrNchZvtx(0), hNV0ACorrNchZvtx(0), hNV0MCorrNchZvtx(0), hNV0CCorrNchZvtxInMSAcc(0), hNV0CNch(0),
    fTreeMC(0x0), fisPSAndkINT7(kTRUE), fisVtxQA(kTRUE),fhasSPDvtx(kTRUE),fhasMultSel(kTRUE),
    fv0multTot(-1.),
    fv0multcorr(-1.),
    fv0cmultTot(-1.),
    fv0cmultcorr(-1.),
    fSPDMult(-1),
    fzvtx(0.),
    fnchInEta1(-1), fnchV0C(-1), fnchV0A(-1), fnchMS(-1)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnaMCMultEpsilonTask::AliAnaMCMultEpsilonTask(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0), fMCEvent(0), hNV0CNchZvtx(0), hNV0ANchZvtx(0), hNV0MNchZvtx(0), hNV0CNchZvtxInMSAcc(0), hNV0CCorrNchZvtx(0), hNV0ACorrNchZvtx(0), hNV0MCorrNchZvtx(0), hNV0CCorrNchZvtxInMSAcc(0), hNV0CNch(0),
    fTreeMC(0x0), fisPSAndkINT7(kTRUE), fisVtxQA(kTRUE),fhasSPDvtx(kTRUE),fhasMultSel(kTRUE),
    fv0multTot(-1.),
    fv0multcorr(-1.),
    fv0cmultTot(-1.),
    fv0cmultcorr(-1.),
    fSPDMult(-1),
    fzvtx(0.),
    fnchInEta1(-1), fnchV0C(-1), fnchV0A(-1), fnchMS(-1)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it,
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!

    DefineOutput(2, TTree::Class());//fTreeMC
}
//_____________________________________________________________________________
AliAnaMCMultEpsilonTask::~AliAnaMCMultEpsilonTask()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnaMCMultEpsilonTask::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)


    //Declaration of histograms

    hNV0CNchZvtx = new TH3F("hNV0CNchZvtx","hNV0CNchZvtx; V0C mult; N_{ch}^{prim}; #it{z}_{vtx} (cm)",1000,0,1000,1000,0,1000,20,-10,10);
    hNV0ANchZvtx = new TH3F("hNV0ANchZvtx","hNV0ANchZvtx; V0A mult; N_{ch}^{prim}; #it{z}_{vtx} (cm)",1000,0,1000,1000,0,1000,20,-10,10);
    hNV0MNchZvtx = new TH3F("hNV0MNchZvtx","hNV0MNchZvtx; V0M mult; N_{ch}^{prim}; #it{z}_{vtx} (cm)",1000,0,1000,1000,0,1000,20,-10,10);
    hNV0CNchZvtxInMSAcc = new TH3F("hNV0CNchZvtxInMSAcc","hNV0CNchZvtxInMSAcc; V0C mult; N_{ch}^{prim}; #it{z}_{vtx} (cm)",1000,0,1000,1000,0,1000,20,-10,10);

    hNV0CCorrNchZvtx = new TH3F("hNV0CCorrNchZvtx","hNV0CCorrNchZvtx; V0C mult corr; N_{ch}^{prim}; #it{z}_{vtx} (cm)",1000,0,1000,1000,0,1000,20,-10,10);
    hNV0ACorrNchZvtx = new TH3F("hNV0ACorrNchZvtx","hNV0ACorrNchZvtx; V0A mult corr; N_{ch}^{prim}; #it{z}_{vtx} (cm)",1000,0,1000,1000,0,1000,20,-10,10);
    hNV0MCorrNchZvtx = new TH3F("hNV0MCorrNchZvtx","hNV0MCorrNchZvtx; V0M mult corr; N_{ch}^{prim}; #it{z}_{vtx} (cm)",1000,0,1000,1000,0,1000,20,-10,10);
    hNV0CCorrNchZvtxInMSAcc = new TH3F("hNV0CCorrNchZvtxInMSAcc","hNV0CCorrNchZvtxInMSAcc; V0C mult corr; N_{ch}^{prim}; #it{z}_{vtx} (cm)",1000,0,1000,1000,0,1000,20,-10,10);

    hNV0CNch = new TH2F("hNV0CNch","hNV0CNch; V0C mult; N_{ch}^{prim}",1000,0,1000,1000,0,1000);

    fOutputList->Add(hNV0CNchZvtx);
    fOutputList->Add(hNV0ANchZvtx);
    fOutputList->Add(hNV0MNchZvtx);
    fOutputList->Add(hNV0CNchZvtxInMSAcc);
    fOutputList->Add(hNV0CCorrNchZvtx);
    fOutputList->Add(hNV0ACorrNchZvtx);
    fOutputList->Add(hNV0MCorrNchZvtx);
    fOutputList->Add(hNV0CCorrNchZvtxInMSAcc);
    fOutputList->Add(hNV0CNch);




    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
    fTreeMC = new TTree("fTreeMC","fTreeMC");
    SetTree(fTreeMC);
    PostData(2,fTreeMC);//the output 2 is defined in the task call of the first constructor
}
//________________________________________________________________________
void AliAnaMCMultEpsilonTask::SetTree(TTree *tr)
{
  tr->Branch("fisPSAndkINT7",&fisPSAndkINT7);
  tr->Branch("fisVtxQA",&fisVtxQA);
  tr->Branch("fhasSPDvtx",&fhasSPDvtx);
  tr->Branch("fhasMultSel",&fhasMultSel);

  tr->Branch("fv0multTot",&fv0multTot);
  tr->Branch("fv0multcorr",&fv0multcorr);
  tr->Branch("fv0cmultTot",&fv0cmultTot);
  tr->Branch("fv0cmultcorr",&fv0cmultcorr);
  tr->Branch("fSPDMult",&fSPDMult);
  tr->Branch("fzvtx",&fzvtx);

  tr->Branch("fnchInEta1",&fnchInEta1);
  tr->Branch("fnchV0C",&fnchV0C);
  tr->Branch("fnchV0A",&fnchV0A);
  tr->Branch("fnchMS",&fnchMS);

}
//_____________________________________________________________________________
void AliAnaMCMultEpsilonTask::UserExec(Option_t *)
{

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
        // check if there actually is an event:
    if(!fAOD) { return; }

    Bool_t iskINT7Event = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
    if (!iskINT7Event)
    {
      //then PS&kINT7 is kFALSE
      fisPSAndkINT7=kFALSE;
    }

    AliMultSelection *multSelection = (AliMultSelection *)fAOD->FindListObject("MultSelection");
    if(!multSelection) fhasMultSel=kFALSE;


    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if(fMCEvent) ProcessMCParticles();
    if(!fMCEvent) { return; }


    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}

void AliAnaMCMultEpsilonTask::ProcessMCParticles()
{
    //MC zvtx
    TArrayF MC_Vtx_true_XYZ(3);
    fMCEvent->GenEventHeader()->PrimaryVertex(MC_Vtx_true_XYZ);
    float VertexZ = MC_Vtx_true_XYZ[2];

    //selections on the SPD vertex
    const AliAODVertex *spdVtx = fAOD->GetPrimaryVertexSPD();
    float spdVtxZ = 0.;
    if(!spdVtx)
    {
      fhasSPDvtx=kFALSE;
      fisVtxQA=kFALSE;
      spdVtxZ=VertexZ;//MC z vtx
    }
    else
    {
      //SPD vertex exists
      spdVtxZ = spdVtx->GetZ();
      if(TMath::Abs(spdVtxZ) > 10) fisVtxQA=kFALSE;
      if (spdVtx->GetNContributors()<=0) fisVtxQA=kFALSE;
      Double_t cov[6]={0};
      spdVtx->GetCovarianceMatrix(cov);
      Double_t zRes = TMath::Sqrt(cov[5]);
      if (spdVtx->IsFromVertexerZ() && (zRes>0.25)) fisVtxQA=kFALSE;
    }


    fzvtx=spdVtxZ;


    TClonesArray *stack = nullptr;
    TList *lst = fAOD->GetList();
    stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
    if (!stack) { return; }

    AliAODMCHeader *mcHeader = 0;
    mcHeader = (AliAODMCHeader*)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEHFTreeCreator::UserExec: MC header branch not found!\n");
      return;
    }

    TClonesArray* AODMCTrackArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
        if (!AODMCTrackArray){
            Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
            this->Dump();
            return;
        }

    //-------------------------- V0 multiplicity -------------------------------
    AliAODVZERO *vzeroAOD = dynamic_cast<AliAODVZERO *>( dynamic_cast<AliAODEvent *>(fAOD)->GetVZEROData());
    Float_t V0AMult = vzeroAOD->GetMTotV0A();
    Float_t V0CMult = vzeroAOD->GetMTotV0C();
    fv0multTot=V0AMult+V0CMult;//total V0 multiplicity in this event
    fv0cmultTot=V0CMult;//total V0 multiplicity in this event

    //------------------------   V0M Correction
    Float_t vzeroMultACorr=V0AMult, vzeroMultCCorr=V0CMult;
    vzeroMultACorr = AliESDUtils::GetCorrV0A(V0AMult,spdVtxZ);
    vzeroMultCCorr = AliESDUtils::GetCorrV0C(V0CMult,spdVtxZ);
    fv0multcorr = vzeroMultACorr + vzeroMultCCorr; // corrected V0 multiplicity
    fv0cmultcorr = vzeroMultCCorr; // corrected V0 multiplicity
    //________________________________________________________________________


    fSPDMult= AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(fAOD, -1, 1);


    Int_t TrueMCParticlesInV0C = 0;
    Int_t TrueMCParticlesInV0A = 0;
    Int_t TrueMCParticlesInV0M = 0;
    Int_t TrueMCParticlesInMSAcc = 0;
    Int_t TrueMCParticlesINELgt0 = 0;

    for (Int_t mcTrack = 0;  mcTrack < (fMCEvent->GetNumberOfTracks()); mcTrack++)
    {

        AliMCParticle *particle = (AliMCParticle*)fMCEvent->GetTrack(mcTrack);
        if (!particle) continue;

        //if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(mcTrack, mcHeader, stack)) { continue; }


        Double_t vz = particle->Zv();

        Bool_t TrackIsPrim = particle->IsPhysicalPrimary();
        if (!TrackIsPrim)continue;

        Double_t pt = particle->Pt();
        //if(pt<.15)continue;
        //if(pt>5)continue;
        //cout<<pt<<endl;
        Double_t trueMCeta = particle->Eta();
        Double_t trueMCcharge = particle->Charge();
        if(trueMCcharge == 0)continue;

        if(trueMCeta >-1 && trueMCeta <1)
        {
            TrueMCParticlesINELgt0++;
        }

        if(trueMCeta >-3.7 && trueMCeta <-1.7)
        {
            TrueMCParticlesInV0C++;
        }

        if(trueMCeta >2.8 && trueMCeta <5.1)
        {
            TrueMCParticlesInV0A++;
        }

        if(trueMCeta >-4 && trueMCeta <-2.5)
        {
            TrueMCParticlesInMSAcc++;
        }



    }//for (Int_t mcTrack = 0;  mcTrack < (fMCEvent->GetNumberOfTracks()); mcTrack++)

    TrueMCParticlesInV0M=TrueMCParticlesInV0C+TrueMCParticlesInV0A;
    fnchInEta1=TrueMCParticlesINELgt0;
    fnchV0C=TrueMCParticlesInV0C;
    fnchV0A=TrueMCParticlesInV0A;
    fnchMS=TrueMCParticlesInMSAcc;

    hNV0CNchZvtx->Fill(fv0cmultTot, TrueMCParticlesInV0C, spdVtxZ);
    hNV0ANchZvtx->Fill(V0AMult, TrueMCParticlesInV0A, spdVtxZ);
    hNV0MNchZvtx->Fill(fv0multTot, TrueMCParticlesInV0M, spdVtxZ);
    hNV0CNchZvtxInMSAcc->Fill(fv0cmultTot, TrueMCParticlesInMSAcc, spdVtxZ);

    hNV0CCorrNchZvtx->Fill(fv0cmultcorr, TrueMCParticlesInV0C, spdVtxZ);
    hNV0ACorrNchZvtx->Fill(vzeroMultACorr, TrueMCParticlesInV0A, spdVtxZ);
    hNV0MCorrNchZvtx->Fill(fv0multcorr, TrueMCParticlesInV0M, spdVtxZ);
    hNV0CCorrNchZvtxInMSAcc->Fill(fv0cmultcorr, TrueMCParticlesInMSAcc, spdVtxZ);

    hNV0CNch->Fill(fv0cmultTot, TrueMCParticlesInV0C);



    fTreeMC->Fill();//fills all the branches with their different elements

    PostData(2, fTreeMC);

}







//_____________________________________________________________________________
void AliAnaMCMultEpsilonTask::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
