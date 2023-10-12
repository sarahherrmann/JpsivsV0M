#pragma once
#include "AliAODTrack.h"

class treeReaderHelper
{
 public:
   treeReaderHelper();
   ~treeReaderHelper() = default;
   //void Clear();


   // Methods definition
   void SetMassCuts(Double_t lowMassCut, Double_t upMassCut) { fLowMassCut = lowMassCut; fUpMassCut = upMassCut; }
   void SetJpsiEtaCuts(Double_t min, Double_t max) { fJpsiEtaMin = min; fJpsiEtaMax = max; }
   void InitHistograms()
   {
     Int_t nBinsMinv = TMath::Nint((fUpMassCut - fLowMassCut)/0.1);
     fJpsiMinv = new TH1D("fJpsiMinv","Minv;M^{inv} (GeV/c^{2})",2*nBinsMinv,fLowMassCut,fUpMassCut);
     fJpsiPt = new TH1D("fJpsiPt",";#it{p}_{T} of J/#psi (GeV/c)",20,0,10);
     fMuonEta = new TH1D("fMuonEta","fMuonEta",50,-4,-2.5);
     fMuonPhi = new TH1D("fMuonPhi","fMuonPhi",100,0,2*TMath::Pi());

     fNV0OverNV0Mean = new TH1F("fNV0OverNV0Mean","fNV0OverNV0Mean",5000,0,50);

     for (int ichannel=0; ichannel<64; ichannel++)//can be reduced to 32, only V0C
     {
       //fV0MultPerChannelNoMu[ichannel] = new TH1F(Form("fV0MultPerChannelNoMu_%d",ichannel),Form("fV0MultPerChannelNoMu_%d",ichannel),500,0,250);
       fV0MultPerChannelOneMu[ichannel] = new TH1F(Form("fV0MultPerChannelOneMu_%d",ichannel),Form("fV0MultPerChannelOneMu_%d",ichannel),500,0,250);
       fV0MultPerChannelTwoMu[ichannel] = new TH1F(Form("fV0MultPerChannelTwoMu_%d",ichannel),Form("fV0MultPerChannelTwoMu_%d",ichannel),500,0,250);
     }

     fEtaPhiV0MeanPerChannelNoMu = new TH2F("fEtaPhiV0MeanPerChannelNoMu","fEtaPhiV0MeanPerChannelNoMu; #eta; #varphi",4,-3.7,-1.7,8,0,2*TMath::Pi());

     for (int ichannel=0; ichannel<32; ichannel++)
     {
       fEtaPhiV0PerChannelOneMu[ichannel] = new TH2F(Form("fEtaPhiV0PerChannelOneMu_%d",ichannel),Form("fEtaPhiV0PerChannelOneMu_%d; #eta; #varphi",ichannel),4,-3.7,-1.7,8,0,2*TMath::Pi());
       fEtaPhiV0MeanPerChannelOneMu[ichannel] = new TProfile2D(Form("fEtaPhiV0MeanPerChannelOneMu_%d",ichannel),Form("fEtaPhiV0MeanPerChannelOneMu_%d; #eta; #varphi",ichannel),4,-3.7,-1.7,8,0,2*TMath::Pi());
       hDiffEtaPhiMeanV0[ichannel] = new TH2F(Form("hDiffEtaPhiMeanV0_%d",ichannel),Form("hDiffEtaPhiMeanV0_%d; #eta; #varphi",ichannel),4,-3.7,-1.7,8,0,2*TMath::Pi());
     }


   }

   void SetTreeAddresses(TTree *trDiMu, TTree *trSingleMu, TTree *trInt7);//set the tree branch adresses

   void GetNV0Histogram(TList *histoList, TTree *trInt7);// get the NV0/<NV0> histogram, to determine binning
   void processV0PerChannelSingleMu(TTree *trSingleMu);//fills the fV0MultPerChannel histograms
   void readEvents(TTree *tr);//reads the events inside the tree and calls methods processDimuons ...

   void DeriveExcessV0MeanPerChannel();//must be called after all the previous methods


   void writeOutput(const char *outputFileName);


   Double_t CalcMinv(AliAODTrack *track1,AliAODTrack *track2) const;
   Double_t CalcY(AliAODTrack *track1,AliAODTrack *track2) const;
   Double_t CalcPt(AliAODTrack *track1,AliAODTrack *track2) const;
   Double_t CalcPhi(AliAODTrack *track1,AliAODTrack *track2) const;
   Double_t CalcCos(AliAODTrack *track1,AliAODTrack *track2) const;


   //tree branches objects for dimuon tree
   Int_t   fRunN; //!
   //Bool_t fispileupspd; //!
   Float_t fv0multTot; //!
   Float_t fv0multcorr; //!
   Double_t fitsmult; //!
   Double_t fv0mpercentile; //!
   Double_t fv0apercentile; //!
   Double_t fv0cpercentile; //!
   Double_t fcl1percentile; //!
   Double_t fspdpercentile; //!
   Double_t fzvtx; //!
   Float_t fv0mult[64]; ///! multiplicity for each channel
   TClonesArray *ftracksmu = new TClonesArray("AliAODTrack",10); //! muon tracks
   //the definition here is needed !

   //Tree branches for fTreeSingleMu:
   Float_t fv0multSMu[64]; //! multiplicity for each channel
   Float_t etaMu;          //! eta of the single muon
   Float_t phiMu;          //! phi of the single muon
   Float_t ptMu;           //! pt of the single muon


   //Tree branches for fTreeINT7:
   Float_t fv0multTotINT7;  //! total multiplicity in V0
   Float_t fv0multCorrINT7;  //! total corrected multiplicity in V0
   Double_t fv0mpercentile7; //!
   Double_t fv0apercentile7; //!
   Double_t fv0cpercentile7; //!
   Double_t fcl1percentile7; //!
   Double_t fspdpercentile7; //!

   //histograms filled in kINT7 events:
   TH1F* hV0MultTot;        //! total multiplicity in the whole V0
   TH1F* hV0MultCorr;       //! total multiplicity in the whole V0, corrected
   TH1F* hV0MultPerChannelNoMu[64]; //! multiplicity per V0 channel when there is no muon

 private:


   // Cut definitions
   Double_t fLowMassCut;
   Double_t fUpMassCut;

   Double_t fJpsiEtaMin;
   Double_t fJpsiEtaMax;

   //for V0 channelMap

   std::map<std::pair<int,int>, int> fChannelMapV0;
   float etaRanges[4][2]={{-3.7,-3.2},{-3.2,-2.7},{-2.7,-2.2},{-2.2,-1.7}};//just for V0C
   float phiRanges[9]={0, TMath::Pi()/4, TMath::Pi()/2, 3*TMath::Pi()/4, TMath::Pi(), 5*TMath::Pi()/4, 6*TMath::Pi()/4, 7*TMath::Pi()/4, 2*TMath::Pi()};
   int GetChannelFromEtaPhi(float eta, float phi);//get the number of the V0 channel in that eta and phi region
   std::vector<int> fChannelsWithMuons;//contains chId of V0 channels with muons

   void processDimuons();//fills some dimuon histograms and applies dimuon cuts





   //Histograms
   TH1D *fJpsiMinv;
   TH1D *fJpsiPt;
   TH1D *fMuonEta;
   TH1D *fMuonPhi;

   TH1F *fNV0OverNV0Mean;

   TH1F *fV0MultPerChannelNoMu[64];//V0 mult in each channel, when no muon hits the channel
   TH1F *fV0MultPerChannelOneMu[64];//V0 mult in each channel, when one muon hits the channel
   TH1F *fV0MultPerChannelTwoMu[64];//V0 mult in each channel, when two muons hit the channel

   TH2F *fEtaPhiV0MeanPerChannelNoMu;//eta phi map of the amount of signal in V0 channels, when no muon in the event
   TH2F *fEtaPhiV0PerChannelOneMu[32];//V0 mult in a channel, when one muon hits the channel, and eta phi distribution of mult
   TProfile2D *fEtaPhiV0MeanPerChannelOneMu[32];//V0 mult in a channel, when one muon hits the channel, and eta phi distribution of mult

   TH2F *hDiffEtaPhiMeanV0[32];

};
#include "treeReaderHelper.cxx" //really needed, else linking problems appear
