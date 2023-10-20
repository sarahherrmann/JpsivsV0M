/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnaTaskJpsiVsV0M_H
#define AliAnaTaskJpsiVsV0M_H

#include "AliAnalysisTaskSE.h"

class AliMuonTrackCuts;
class AliAnalysisUtils;

class AliAnaTaskJpsiVsV0M : public AliAnalysisTaskSE
{
    public:
                                AliAnaTaskJpsiVsV0M();
                                AliAnaTaskJpsiVsV0M(const char *name);
        virtual                 ~AliAnaTaskJpsiVsV0M();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);


        //----- ANALYSIS
        void NotifyRun();//automatically called
        void SetTreeDiMu(TTree *tr);
        void SetTreeSingleMu(TTree *tr);
        void SetTreeINT7(TTree *tr);
        void GetAcceptedTracksMuonArm(AliAODEvent *aodEvent);


    private:
        AliAODEvent*            fAOD;           //! input event
        TList*                  fOutputList;    //! output list for histograms

        TString fCentMethod;
        AliMuonTrackCuts* fMuonTrackCuts; // Muon Track cuts
        Short_t fTriggerMatchLevelMuon;

        Long64_t nentries;

        TTree *fTreeDiMu;       //! Output tree for dimu unlike sign event trigger
        TTree *fTreeSingleMu;   //! Output tree for single muon low pt ev trigger
        TTree *fTreeINT7;      //! Output tree for kINT7 ev trigger
        AliAnalysisUtils* fUtils;   //! analysis utils to detect pileup

        //Tree branches for fTreeDiMu:
        Int_t   fRunN; //!
        Bool_t fispileupspd; //!
        Float_t fv0multTot; //!
        Float_t fv0multcorr; //!
        Float_t fv0cmultTot; //!
        Float_t fv0cmultcorr; //!
        Double_t fitsmult; //!
        Double_t fv0mpercentile; //!
        Double_t fv0apercentile; //!
        Double_t fv0cpercentile; //!
        Double_t fcl1percentile; //!
        Double_t fspdpercentile; //!
        Double_t fzvtx; //!
        Float_t fv0mult[64]; ///! multiplicity for each channel
        TClonesArray *ftracksmu; //! muon tracks

        //Tree branches for fTreeSingleMu:
        Float_t fv0multSMu[64]; //! multiplicity for each channel
        Float_t etaMu;          //! eta of the single muon
        Float_t phiMu;          //! phi of the single muon
        Float_t ptMu;           //! pt of the single muon


        //Tree branches for fTreeINT7:
        Float_t fv0multTotINT7;  //! total multiplicity in V0
        Float_t fv0multCorrINT7;  //! total corrected multiplicity in V0
        Float_t fv0cmultTotINT7;  //! total multiplicity in V0C
        Float_t fv0cmultCorrINT7;  //! total corrected multiplicity in V0C
        Double_t fv0mpercentile7; //!
        Double_t fv0apercentile7; //!
        Double_t fv0cpercentile7; //!
        Double_t fcl1percentile7; //!
        Double_t fspdpercentile7; //!

        //histograms filled in kINT7 events:
        TH1F* hV0MultTot;        //! total multiplicity in the whole V0
        TH1F* hV0MultCorr;       //! total multiplicity in the whole V0, corrected
        TH1F* hV0CMultTot;        //! total multiplicity in the whole V0C
        TH1F* hV0CMultCorr;       //! total multiplicity in the whole V0C, corrected
        TH1F* hV0MultPerChannelNoMu[64]; //! multiplicity per V0 channel when there is no muon

        AliAnaTaskJpsiVsV0M(const AliAnaTaskJpsiVsV0M&); // not implemented
        AliAnaTaskJpsiVsV0M& operator=(const AliAnaTaskJpsiVsV0M&); // not implemented

        ClassDef(AliAnaTaskJpsiVsV0M, 1);
};

#endif
