#pragma once
#include "AliAODTrack.h"

class treeReaderHelper
{
 public:
   treeReaderHelper();
   ~treeReaderHelper() = default;


   float massBinWidth = 0.005;//bin width for invariant mass

   // Methods definition
   void SetRunNbMap(const char *fperiod, int subset=0);
   void SetMassCuts(Double_t lowMassCut, Double_t upMassCut) { fLowMassCut = lowMassCut; fUpMassCut = upMassCut; }
   void SetJpsiEtaCuts(Double_t min, Double_t max) { fJpsiEtaMin = min; fJpsiEtaMax = max; }
   void InitHistograms()
   {

     Int_t nBinsMinv = TMath::Nint((fUpMassCut - fLowMassCut)/massBinWidth);

     //float binEdgesNV0C[] = {0.0, 28.665, 54.715, 72.56, 95.545, 127.145, 148.775, 178.085, 225.30, 324.495, 619.53998};
     float binEdgesNV0C[] = {0.0, 28.665, 54.715, 72.56, 95.545, 127.145, 148.775, 178.085, 225.30, 324.495, 364.435, 403.850, 696.985};
     float binEdgesPt[]={0.,1.,3.,1000};

     float binEdgesMinv[nBinsMinv+1];
     float valueBinEdge=fLowMassCut;
     int ibin=0;
     while(valueBinEdge<fUpMassCut+massBinWidth)
     {
       binEdgesMinv[ibin]=valueBinEdge;
       valueBinEdge=valueBinEdge+massBinWidth;
       ibin++;
     }


     float binEdgesPtFine[] = {0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 5.25, 5.50, 5.75, 6.00, 6.25, 6.50, 6.75, 7.00, 7.25, 7.50, 7.75, 8.00, 8.25, 8.50, 8.75, 9.00, 9.25, 9.50, 9.75, 10.00, 10.50, 11.00, 11.50, 12.00, 12.50, 13.00, 13.50, 14.00, 15.00, 16.00, 17.00, 18.00, 19.00, 20.00};

     hMultPtJpsiMinv = new TH3D("hMultPtJpsiMinv",";NV0C;#it{p}_{T};M^{inv} (GeV/c^{2})",(sizeof(binEdgesNV0C) / sizeof(binEdgesNV0C[0]) - 1), binEdgesNV0C,(sizeof(binEdgesPt) / sizeof(binEdgesPt[0]) - 1), binEdgesPt, nBinsMinv,binEdgesMinv);
     hMultPtJpsiMinv->Sumw2();

     hMultPtJpsiMinvCorr = new TH3D("hMultPtJpsiMinvCorr",";NV0C;#it{p}_{T};M^{inv} (GeV/c^{2})",(sizeof(binEdgesNV0C) / sizeof(binEdgesNV0C[0]) - 1), binEdgesNV0C,(sizeof(binEdgesPt) / sizeof(binEdgesPt[0]) - 1), binEdgesPt, nBinsMinv,binEdgesMinv);
     hMultPtJpsiMinvCorr->Sumw2();

     hMultPtJpsiMinvCorrAxE = new TH3D("hMultPtJpsiMinvCorrAxE",";NV0C;#it{p}_{T};M^{inv} (GeV/c^{2})",(sizeof(binEdgesNV0C) / sizeof(binEdgesNV0C[0]) - 1), binEdgesNV0C,(sizeof(binEdgesPt) / sizeof(binEdgesPt[0]) - 1), binEdgesPt, nBinsMinv,binEdgesMinv);
     hMultPtJpsiMinvCorrAxE->Sumw2();

     hMultPtJpsiMinvCorrAxERenorm = new TH3D("hMultPtJpsiMinvCorrAxERenorm",";NV0C;#it{p}_{T};M^{inv} (GeV/c^{2})",(sizeof(binEdgesNV0C) / sizeof(binEdgesNV0C[0]) - 1), binEdgesNV0C,(sizeof(binEdgesPt) / sizeof(binEdgesPt[0]) - 1), binEdgesPt, nBinsMinv,binEdgesMinv);
     hMultPtJpsiMinvCorrAxERenorm->Sumw2();

     hMultPtJpsiMinvAllCorrY1 = new TH3D("hMultPtJpsiMinvAllCorrY1",";NV0C;#it{p}_{T};M^{inv} (GeV/c^{2})",(sizeof(binEdgesNV0C) / sizeof(binEdgesNV0C[0]) - 1), binEdgesNV0C,(sizeof(binEdgesPt) / sizeof(binEdgesPt[0]) - 1), binEdgesPt, nBinsMinv,binEdgesMinv);
     hMultPtJpsiMinvAllCorrY1->Sumw2();

     hMultPtJpsiMinvAllCorrY2 = new TH3D("hMultPtJpsiMinvAllCorrY2",";NV0C;#it{p}_{T};M^{inv} (GeV/c^{2})",(sizeof(binEdgesNV0C) / sizeof(binEdgesNV0C[0]) - 1), binEdgesNV0C,(sizeof(binEdgesPt) / sizeof(binEdgesPt[0]) - 1), binEdgesPt, nBinsMinv,binEdgesMinv);
     hMultPtJpsiMinvAllCorrY2->Sumw2();

     for (int ibinY=0; ibinY<5; ibinY++)//5 bins in rapidity
     {
       hMultPtJpsiMinvAllCorrY[ibinY] = new TH3D(Form("hMultPtJpsiMinvAllCorrY_%d",ibinY),";NV0C;#it{p}_{T};M^{inv} (GeV/c^{2})",(sizeof(binEdgesNV0C) / sizeof(binEdgesNV0C[0]) - 1), binEdgesNV0C,(sizeof(binEdgesPt) / sizeof(binEdgesPt[0]) - 1), binEdgesPt, nBinsMinv,binEdgesMinv);
       hMultPtJpsiMinvAllCorrY[ibinY]->Sumw2();
     }

     hMultMinusMultCorr = new TH1F("hMultMinusMultCorr","hMultMinusMultCorr",5000,-20,300);

     hMultMinusMultCorr1mu = new TH1F("hMultMinusMultCorr1mu","hMultMinusMultCorr1mu",5000,-20,300);


     fJpsiMinv = new TH1D("fJpsiMinv","Minv;M^{inv} (GeV/c^{2})",nBinsMinv,fLowMassCut,fUpMassCut);
     fJpsiPt = new TH1D("fJpsiPt",";#it{p}_{T} of J/#psi (GeV/c)",200,0,20);
     fJpsiY = new TH1D("fJpsiY",";#it{y} of J/#psi",200,-4.1,-2.4);


     fNV0OverNV0Mean = new TH1F("fNV0OverNV0Mean","fNV0OverNV0Mean",5000,0,50);
     fNV0COverNV0CMean = new TH1F("fNV0COverNV0CMean","fNV0COverNV0CMean",20000,0,200);

     hV0CMultCorrCopy = new TH1F("hV0CMultCorrCopy","hV0CMultCorrCopy",600000,0,3000);
     hV0CAfterNorm = new TH1F("hV0CAfterNorm","hV0CAfterNorm",600000,0,3000);

     hSPDMultVsV0CMult = new TH2F("hSPDMultVsV0CMult","hSPDMultVsV0CMult; V0C mult; SPD mult",1000,0,1000,200,0,200);

     // ------- TEST: Only for 18m

     fNV0COverGlobalMean_292273 = new TH1F("fNV0COverGlobalMean_292273","fNV0COverGlobalMean_292273",20000,0,200);
     fNV0COverGlobalMean_290223 = new TH1F("fNV0COverGlobalMean_290223","fNV0COverGlobalMean_290223",20000,0,200);

     // ---------------------------


     hV0CMultInSingleMu = new TH1F("hV0CMultInSingleMu","hV0CMultInSingleMu",20000,0,10000);

     hMuonExcess = new TH1F("hMuonExcess","hMuonExcess",25,0,25);

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


     //histograms to check the zvtx dependency of V0Cmult
     NevZvtx = new TH1F("NevZvtx","NevZvtx",200,-10,10);
     V0CMultCorrvsZvtx = new TH2F("V0CMultCorrvsZvtx","V0CMultCorrvsZvtx; #it{z}_{vtx} (cm); #it{N}_{V0C}",40,-10,10,700,0,700);
     V0CMultCorrRenormvsZvtx = new TH2F("V0CMultCorrRenormvsZvtx","V0CMultCorrRenormvsZvtx; #it{z}_{vtx} (cm); #it{N}_{V0C}",40,-10,10,700,0,700);
     V0CMultvsZvtx = new TH2F("V0CMultvsZvtx","V0CMultvsZvtx; #it{z}_{vtx} (cm); #it{N}_{V0C}",40,-10,10,700,0,700);

     V0MultCorrvsZvtx = new TH2F("V0MultCorrvsZvtx","V0MultCorrvsZvtx; #it{z}_{vtx} (cm); #it{N}_{V0C}",40,-10,10,1000,0,1000);
     V0MultvsZvtx = new TH2F("V0MultvsZvtx","V0MultvsZvtx; #it{z}_{vtx} (cm); #it{N}_{V0C}",40,-10,10,1000,0,1000);


   }

   void SetListV0Norm(TFile *fileNormV0);
   void SetTreeAddresses(TTree *trDiMu, TTree *trSingleMu, TTree *trInt7);//set the tree branch adresses
   void GetJpsiAxE(TH2F *hAccxEff)//Get the Jpsi AxE from Zaida's file
   {
     hJpsiAxE=(TH2F*)hAccxEff->Clone();
   };

   void GetNV0Histogram(TList *histoList, TTree *trInt7);// get the NV0/<NV0> histogram, to determine binning
   void processV0PerChannelSingleMu(TTree *trSingleMu);//fills the fV0MultPerChannel histograms
   void readEvents(TTree *tr);//reads the events inside the tree and calls methods processDimuons ...

   void DeriveExcessV0MeanPerChannel();//must be called after all the previous methods
   void fillCorrFactorMap();//must be called after processV0PerChannelSingleMu

   Float_t GetV0CMeanCorrection(TList* fV0List, Float_t V0Corr, Float_t refMult_V0C, Int_t runNb, Float_t zvtx);//compute anti-aging correction for V0C

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
   Int_t fSPDMult;//N tracklets in SPD, eta in -1 1

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
   Double_t fzvtxkINT7; //!SPD zvtx in kINT7 events
   Int_t fRunNInt7;//! run number
   Int_t fSPDMult7;//N tracklets in SPD, eta in -1 1

   //histograms filled in kINT7 events:
   TH1F* hV0MultTot;        //! total multiplicity in the whole V0
   TH1F* hV0MultCorr;       //! total multiplicity in the whole V0, corrected
   TH1F* hV0CMultTot;        //! total multiplicity in the whole V0C
   TH1F* hV0CMultCorr;       //! total multiplicity in the whole V0C, corrected
   TH1F* hV0MultPerChannelNoMu[64]; //! multiplicity per V0 channel when there is no muon


   TH1F* hV0CMultCorrCopy;
   TH1F* hV0CAfterNorm;

   TH2F* hSPDMultVsV0CMult;

   TH1F* hV0CMultInSingleMu;//total V0C mult in single muon events
   TH1F* hMuonExcess;//e[i] for each V0C channel

   TH2F* hJpsiAxE;

 private:


   // Cut definitions
   Double_t fLowMassCut;
   Double_t fUpMassCut;

   Double_t fJpsiEtaMin;
   Double_t fJpsiEtaMax;

   //runList definition

   const Int_t *runList;
   Int_t nRuns;

   std::map<int,int> fMapRunList;

   const Float_t *fParamZvtxCorr;//parameters for fZvtxCorr, stored in runList.h

   //For normalisation reference
   Float_t fMeanV0C_corr;//the mean <NV0C> in this period, used as a reference
   Float_t fRefV0C=77.5;

   TList *listV0CNorm;
   TList *listV0ANorm;
   TList *listV0MNorm;

   TF1 *fZvtxCorr;//The function for the zvtx correction

   //for V0 channelMap

   std::map<std::pair<int,int>, int> fChannelMapV0;//key:index of eta and phi, value: channel id
   float etaRanges[4][2]={{-3.7,-3.2},{-3.2,-2.7},{-2.7,-2.2},{-2.2,-1.7}};//just for V0C
   float phiRanges[9]={0, TMath::Pi()/4, TMath::Pi()/2, 3*TMath::Pi()/4, TMath::Pi(), 5*TMath::Pi()/4, 6*TMath::Pi()/4, 7*TMath::Pi()/4, 2*TMath::Pi()};

   std::vector<int> fChannelsWithMuons;//contains chId of V0 channels with muons

   std::map<int, std::pair<int,int>> fChannelAccrossChannelId;//key:index of channel, value: channel ids of channel at phi+pi/2 and phi-pi/2
   //gives you the indexes of the 2 channels accross from the key channel

   std::map<int, float> fChannelToCorrFactor;//key: index of the channel
   //value: the correction factor deltai for that channel

   int GetChannelFromEtaPhi(float eta, float phi);//get the number of the V0 channel in that eta and phi region

   void processDimuons();//fills some dimuon histograms and applies dimuon cuts, if IsFirstTime=false, then does the correction for dimuon and fills the corrected histogram
   Float_t correctDimuons(AliAODTrack *track1, AliAODTrack *track2, int &nMuons);//applies the correction for muons on the fv0cmult
   //returns the corrected total fv0c mult, for dimuon contribution and the fzvtx corr from ESDUtils



   void bringTo02Pi(float& phi);//bring any phi angle to [0,2pi]





   //Histograms
   TH1D *fJpsiMinv;
   TH1D *fJpsiPt;
   TH1D *fJpsiY;

   TH3D *hMultPtJpsiMinv;//3D histogram of invariant mass (z) versus mult in V0C(x) and pt of Jpsi (y)
   TH3D *hMultPtJpsiMinvCorr;//3D histogram of invariant mass (z) versus mult in V0C(x) and pt of Jpsi (y), corrected for dimuons
   TH3D *hMultPtJpsiMinvCorrAxE;//3D histogram of invariant mass (z) versus mult in V0C(x) and pt of Jpsi (y), corrected for dimuons and Jpsi AxE
   TH3D *hMultPtJpsiMinvCorrAxERenorm;//3D histogram of invariant mass (z) versus mult in V0C(x) and pt of Jpsi (y), corrected for dimuons and Jpsi AxE

   TH3D *hMultPtJpsiMinvAllCorrY1;//3D histogram of invariant mass, all corrections and -4.0<y<-3.2
   TH3D *hMultPtJpsiMinvAllCorrY2;//3D histogram of invariant mass, all corrections and -3.2<=y<-2.4

   TH3D *hMultPtJpsiMinvAllCorrY[5];//3D histogram of invariant mass, all corrections and y in 5 different bins


   TH1F *hMultMinusMultCorr;//difference between uncorrected mult and corrected mult for dimu
   TH1F *hMultMinusMultCorr1mu;//difference between uncorrected mult and corrected mult for dimu, only one mu in V0C acc


   TH1F *fNV0OverNV0Mean;
   TH1F *fNV0COverNV0CMean;

   TH1F *fNV0COverGlobalMean_292273;
   TH1F *fNV0COverGlobalMean_290223;

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


   TH1F *NevZvtx;//Zvtx distribution of kINT7 events
   TH2F *V0CMultvsZvtx;//Zvtx distribution of V0C mult, for kINT7 events
   TH2F *V0CMultCorrvsZvtx;//Zvtx distribution of V0C mult, for kINT7 events
   TH2F *V0CMultCorrRenormvsZvtx;//Zvtx distribution of V0C mult, for kINT7 events

   TH2F *V0MultvsZvtx;//Zvtx distribution of V0C mult, for kINT7 events
   TH2F *V0MultCorrvsZvtx;//Zvtx distribution of V0C mult, for kINT7 events

};
#include "treeReaderHelper.cxx" //really needed, else linking problems appear
