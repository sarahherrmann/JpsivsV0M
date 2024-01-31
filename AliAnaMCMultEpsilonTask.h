/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnaMCMultEpsilonTask_H
#define AliAnaMCMultEpsilonTask_H

#include "AliAnalysisTaskSE.h"

#include "AliAnalysisUtils.h"
#include "AliGenEventHeader.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliEventCuts.h"
#include "TH3F.h"
#include "TH2F.h"

class AliAnaMCMultEpsilonTask : public AliAnalysisTaskSE
{
    public:
                                AliAnaMCMultEpsilonTask();
                                AliAnaMCMultEpsilonTask(const char *name);
        virtual                 ~AliAnaMCMultEpsilonTask();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        virtual void            ProcessMCParticles();


        void SetTree(TTree *tr);



    private:

        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list

        TH3F *hNV0CNchZvtx;
        TH3F *hNV0ANchZvtx;
        TH3F *hNV0MNchZvtx;

        TH3F *hNV0CNchZvtxInMSAcc;//in the eta acc of Muon Spectrometer -4,-2.5

        TH3F *hNV0CCorrNchZvtx;
        TH3F *hNV0ACorrNchZvtx;
        TH3F *hNV0MCorrNchZvtx;

        TH3F *hNV0CCorrNchZvtxInMSAcc;//in the eta acc of Muon Spectrometer -4,-2.5

        TH2F *hNV0CNch;


        AliMCEvent* fMCEvent;       //! corresponding MC event


        TTree *fTreeMC;       //! Output tree


        //Tree branches for fTreeMC:
        Bool_t fisPSAndkINT7; //! isPS&kINT7
        Bool_t fisVtxQA; //! is the event selected by vtxQA ?
        Bool_t fhasSPDvtx; //! does it have a SPD vtx ?
        Bool_t fhasMultSel; //! does this event have mult selection

        Float_t fv0multTot; //!
        Float_t fv0multcorr; //!
        Float_t fv0cmultTot; //!
        Float_t fv0cmultcorr; //!
        Double_t fzvtx; //! SPD vertex z coordinate, or MC vtx if fhasSPDvtx=kFALSE
        Int_t fSPDMult; //! N_tracklets in SPD, with eta in [-1,1]

        Int_t fnchInEta1; //! N_ch^prim with eta in [-1,1]
        Int_t fnchV0C; //! N_ch^prim with eta in [-3.7,-1.7]
        Int_t fnchV0A; //! N_ch^prim with eta in [2.8,5.1]
        Int_t fnchMS; //! N_ch^prim with eta in [-4,-2.5]


        AliAnaMCMultEpsilonTask(const AliAnaMCMultEpsilonTask&); // not implemented
        AliAnaMCMultEpsilonTask& operator=(const AliAnaMCMultEpsilonTask&); // not implemented

        ClassDef(AliAnaMCMultEpsilonTask, 1);
};

#endif
