#pragma once
#include "AliAODTrack.h"

class treeReaderHelper
{
 public:
   treeReaderHelper();
   ~treeReaderHelper() = default;


   float massBinWidth = 0.05;//bin width for invariant mass

   // Methods definition
   void SetMassCuts(Double_t lowMassCut, Double_t upMassCut) { fLowMassCut = lowMassCut; fUpMassCut = upMassCut; }
   void SetJpsiEtaCuts(Double_t min, Double_t max) { fJpsiEtaMin = min; fJpsiEtaMax = max; }
   void InitHistograms()
   {

     Int_t nBinsMinv = TMath::Nint((fUpMassCut - fLowMassCut)/massBinWidth);
     float binEdgesNV0C[] = {0.,30.,57.5,76.,100.,133.,155.,186.,235.,335.,10000.};
     float binEdgesPt[]={0.,1.,3.,1000};

     float binMinv[]={2.,5.};
     float binEdgesMinv[nBinsMinv+1];
     float valueBinEdge=fLowMassCut;
     int ibin=0;
     while(valueBinEdge<fUpMassCut+massBinWidth)
     {
       binEdgesMinv[ibin]=valueBinEdge;
       valueBinEdge=valueBinEdge+massBinWidth;
       ibin++;
     }

     hMultPtJpsiMinv = new TH3D("hMultPtJpsiMinv",";#it{p}_{T};NV0C;M^{inv} (GeV/c^{2})",(sizeof(binEdgesNV0C) / sizeof(binEdgesNV0C[0]) - 1), binEdgesNV0C,(sizeof(binEdgesPt) / sizeof(binEdgesPt[0]) - 1), binEdgesPt, nBinsMinv,binEdgesMinv);
     hMultPtJpsiMinv->Sumw2();

     hMultPtJpsiMinvCorr = new TH3D("hMultPtJpsiMinvCorr",";#it{p}_{T};NV0C;M^{inv} (GeV/c^{2})",(sizeof(binEdgesNV0C) / sizeof(binEdgesNV0C[0]) - 1), binEdgesNV0C,(sizeof(binEdgesPt) / sizeof(binEdgesPt[0]) - 1), binEdgesPt, nBinsMinv,binEdgesMinv);
     hMultPtJpsiMinvCorr->Sumw2();

     hMultMinusMultCorr = new TH1F("hMultMinusMultCorr","hMultMinusMultCorr",5000,-20,300);


     fJpsiMinv = new TH1D("fJpsiMinv","Minv;M^{inv} (GeV/c^{2})",nBinsMinv,fLowMassCut,fUpMassCut);
     fJpsiPt = new TH1D("fJpsiPt",";#it{p}_{T} of J/#psi (GeV/c)",200,0,20);


     fNV0OverNV0Mean = new TH1F("fNV0OverNV0Mean","fNV0OverNV0Mean",5000,0,50);
     fNV0COverNV0CMean = new TH1F("fNV0COverNV0CMean","fNV0COverNV0CMean",20000,0,200);
     hV0CMultCorrCopy = new TH1F("hV0CMultCorrCopy","hV0CMultCorrCopy",20000,0,10000);


     hV0CMultInSingleMu = new TH1F("hV0CMultInSingleMu","hV0CMultInSingleMu",20000,0,10000);

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

       fDiMuCorrectionWDiMuEv[ichannel] = new TH1F(Form("fDiMuCorrectionWDiMuEv_%d",ichannel),Form("fDiMuCorrectionWDiMuEv_%d",ichannel),5000,-20,300);
       fDiMuCorrectionWSiMuEv[ichannel] = new TH1F(Form("fDiMuCorrectionWSiMuEv_%d",ichannel),Form("fDiMuCorrectionWSiMuEv_%d",ichannel),5000,-20,300);

       fEtaPhiMeanV0DiMuonEv[ichannel] = new TProfile2D(Form("fEtaPhiMeanV0DiMuonEv_%d",ichannel),Form("fEtaPhiMeanV0DiMuonEv_%d; #eta; #varphi",ichannel),4,-3.7,-1.7,8,0,2*TMath::Pi());
     }


   }

   void SetTreeAddresses(TTree *trDiMu, TTree *trSingleMu, TTree *trInt7);//set the tree branch adresses

   void GetNV0Histogram(TList *histoList, TTree *trInt7);// get the NV0/<NV0> histogram, to determine binning
   void processV0PerChannelSingleMu(TTree *trSingleMu);//fills the fV0MultPerChannel histograms
   void readEvents(TTree *tr);//reads the events inside the tree and calls methods processDimuons ...

   void DeriveExcessV0MeanPerChannel();//must be called after all the previous methods
   void fillCorrFactorMap();//must be called after processV0PerChannelSingleMu

   void writeOutput(const char *outputFileName);//write all the histograms to an output root file


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


   TH1F* hV0CMultCorrCopy;

   TH1F* hV0CMultInSingleMu;//total V0C mult in single muon events

 private:


   // Cut definitions
   Double_t fLowMassCut;
   Double_t fUpMassCut;

   Double_t fJpsiEtaMin;
   Double_t fJpsiEtaMax;

   //for V0 channelMap

   std::map<std::pair<int,int>, int> fChannelMapV0;//key:index of eta and phi, value: channel id
   float etaRanges[4][2]={{-3.7,-3.2},{-3.2,-2.7},{-2.7,-2.2},{-2.2,-1.7}};//just for V0C
   float phiRanges[9]={0, TMath::Pi()/4, TMath::Pi()/2, 3*TMath::Pi()/4, TMath::Pi(), 5*TMath::Pi()/4, 6*TMath::Pi()/4, 7*TMath::Pi()/4, 2*TMath::Pi()};
   int GetChannelFromEtaPhi(float eta, float phi);//get the number of the V0 channel in that eta and phi region
   std::vector<int> fChannelsWithMuons;//contains chId of V0 channels with muons

   std::map<int, std::pair<int,int>> fChannelAccrossChannelId;//key:index of channel, value: channel ids of channel at phi+pi/2 and phi-pi/2
   //gives you the indexes of the 2 channels accross from the key channel

   std::map<int, float> fChannelToCorrFactor;//key: index of the channel
   //value: the correction factor deltai for that channel

   void processDimuons(bool IsFirstTime = true);//fills some dimuon histograms and applies dimuon cuts, if IsFirstTime=false, then does the correction for dimuon and fills the corrected histogram
   Float_t correctDimuons(AliAODTrack *track1, AliAODTrack *track2, bool IsCorrWithDiMu = false);//applies the correction for muons on the fv0cmult
   //returns the corrected total fv0c mult, for dimuon contribution and the fzvtx corr from ESDUtils
   //weight is the ratio of the V0C mult in dimuon ev / V0C mult in single muon ev


   void bringTo02Pi(float& phi);//bring any phi angle to [0,2pi]





   //Histograms
   TH1D *fJpsiMinv;
   TH1D *fJpsiPt;

   TH3D *hMultPtJpsiMinv;//3D histogram of invariant mass (z) versus mult in V0C(x) and pt of Jpsi (y)
   TH3D *hMultPtJpsiMinvCorr;//3D histogram of invariant mass (z) versus mult in V0C(x) corrected for dimuons and pt of Jpsi (y)

   TH1F *hMultMinusMultCorr;//difference between uncorrected mult and corrected mult for dimu

   TH1F *fNV0OverNV0Mean;
   TH1F *fNV0COverNV0CMean;

   TH1F *fV0MultPerChannelNoMu[64];//V0 mult in each channel, when no muon hits the channel
   TH1F *fV0MultPerChannelOneMu[64];//V0 mult in each channel, when one muon hits the channel
   TH1F *fV0MultPerChannelTwoMu[64];//V0 mult in each channel, when two muons hit the channel

   TH2F *fEtaPhiV0MeanPerChannelNoMu;//eta phi map of the amount of signal in V0 channels, when no muon in the event
   TH2F *fEtaPhiV0PerChannelOneMu[32];//V0 mult in a channel, when one muon hits the channel, and eta phi distribution of mult
   TProfile2D *fEtaPhiV0MeanPerChannelOneMu[32];//V0 mult in a channel, when one muon hits the channel, and eta phi distribution of mult

   TH2F *hDiffEtaPhiMeanV0[32];

   TH1F *fDiMuCorrectionWSiMuEv[32];//=signal in this channel - pedestal for this channel, for dimuon correction
   TH1F *fDiMuCorrectionWDiMuEv[32];//=signal in this channel - pedestal for this channel, for dimuon correction
   TProfile2D *fEtaPhiMeanV0DiMuonEv[32];//V0 mult in a channel, when one muon hits the channel, in dimuon events and eta phi distribution of mult

};
#include "treeReaderHelper.cxx" //really needed, else linking problems appear
