#include "treeReaderHelper.h"
#include "AliAODTrack.h"
#include "AliVVZERO.h"
#include "AliESDUtils.h"
#include "runLists.h"

treeReaderHelper::treeReaderHelper()
{
  //fill the map only once, in the initialization
  for (int ichannel=0; ichannel<32; ichannel++)
  {
    float etaMid =  AliVVZERO::GetVZEROEtaMin(ichannel)+0.25;
    float phiMid = AliVVZERO::GetVZEROAvgPhi(ichannel);

    int indexEtaMin = 0;
    int indexPhiMin = 0;
    while ((etaMid>etaRanges[indexEtaMin][1]) && indexEtaMin<3 )
    {
      indexEtaMin++;
    }

    while (phiMid>phiRanges[indexPhiMin])
    {
      indexPhiMin++;
    }
    indexPhiMin--;
    //filling the map
    fChannelMapV0[{indexEtaMin, indexPhiMin}] = ichannel;
  }//for (int ichannel=0; ichannel<32; ichannel++)

  for (int ichannel=0; ichannel<32; ichannel++)
  {
    //-----------------------
    //Other map: channelId to the 2 channels at -+pi/2
    //Needs previous map to be filled !! Hence the different loop

    float etaMid =  AliVVZERO::GetVZEROEtaMin(ichannel)+0.25;
    float phiMid = AliVVZERO::GetVZEROAvgPhi(ichannel);

    float phiPlusPiOver2 = phiMid + TMath::Pi()/2;
    float phiMinusPiOver2 = phiMid - TMath::Pi()/2;

    bringTo02Pi(phiPlusPiOver2);//between 0 and 2pi
    bringTo02Pi(phiMinusPiOver2);//between 0 and 2pi

    int channelPlusPiOver2 = GetChannelFromEtaPhi(etaMid, phiPlusPiOver2);//channel across the channel hit
    int channelMinusPiOver2 = GetChannelFromEtaPhi(etaMid, phiMinusPiOver2);//channel across the channel hit


    //printf("---- channelId = %d, channel+pi/2 = %d, channel-pi/2 = %d\n", ichannel, channelPlusPiOver2, channelMinusPiOver2);

    fChannelAccrossChannelId[ichannel]={channelPlusPiOver2, channelMinusPiOver2};
  }//for (int ichannel=0; ichannel<32; ichannel++)

}

void treeReaderHelper::SetRunNbMap(const char *fperiod, int subset=0)
{

  //Fills the runlist
  //and the nb of runs in this period
  //and set the appropriate parameters for the fZvtx function
  //(that corrects the profile of zvtx of V0C)

  TString periodStr = fperiod;

  fParamZvtxCorr = fParam18m;//default to avoid getting errors the first time I read the file

  if (periodStr.Contains("18b")) { runList = runList18b; nRuns = sizeof(runList18b)/sizeof(runList18b[0]); fParamZvtxCorr = fParam18b; }
  else if (periodStr.Contains("18c")) { runList = runList18c; nRuns = sizeof(runList18c)/sizeof(runList18c[0]); fParamZvtxCorr = fParam18c; }
  else if (periodStr.Contains("18d")) { runList = runList18d; nRuns = sizeof(runList18d)/sizeof(runList18d[0]); fParamZvtxCorr = fParam18d; }
  else if (periodStr.Contains("18e")) { runList = runList18e; nRuns = sizeof(runList18e)/sizeof(runList18e[0]); fParamZvtxCorr = fParam18e; }

  else if (periodStr.Contains("18f") && subset == 1) { runList = runList18f1; nRuns = sizeof(runList18f1)/sizeof(runList18f1[0]); fParamZvtxCorr = fParam18f; }

  else if (periodStr.Contains("18f") && subset <= 0) { runList = runList18f; nRuns = sizeof(runList18f)/sizeof(runList18f[0]); fParamZvtxCorr = fParam18f; }

  else if (periodStr.Contains("18l")) { runList = runList18l; nRuns = sizeof(runList18l)/sizeof(runList18l[0]); fParamZvtxCorr = fParam18l; }
  else if (periodStr.Contains("18m") && subset <= 0) { runList = runList18m; nRuns = sizeof(runList18m)/sizeof(runList18m[0]); fParamZvtxCorr = fParam18m; }
  else if (periodStr.Contains("18m") && subset == 1) { runList = runList18m1; nRuns = sizeof(runList18m1)/sizeof(runList18m1[0]); }
  else if (periodStr.Contains("18m") && subset == 2) { runList = runList18m2; nRuns = sizeof(runList18m2)/sizeof(runList18m2[0]); }
  else if (periodStr.Contains("18o")) { runList = runList18o; nRuns = sizeof(runList18o)/sizeof(runList18o[0]); fParamZvtxCorr = fParam18o; }
  else if (periodStr.Contains("18p")) { runList = runList18p; nRuns = sizeof(runList18p)/sizeof(runList18p[0]); fParamZvtxCorr = fParam18p; }

  else if (periodStr.Contains("17h") && subset <= 0) { runList = runList17h; nRuns = sizeof(runList17h)/sizeof(runList17h[0]); fParamZvtxCorr = fParam17h; }
  else if (periodStr.Contains("17h") && subset == 1) { runList = runList17hManuel; nRuns = sizeof(runList17hManuel)/sizeof(runList17hManuel[0]); }
  else if (periodStr.Contains("17i")) { runList = runList17i; nRuns = sizeof(runList17i)/sizeof(runList17i[0]); fParamZvtxCorr = fParam17i; }
  else if (periodStr.Contains("17k")) { runList = runList17k; nRuns = sizeof(runList17k)/sizeof(runList17k[0]); fParamZvtxCorr = fParam17k; }
  else if (periodStr.Contains("17l")) { runList = runList17l; nRuns = sizeof(runList17l)/sizeof(runList17l[0]); fParamZvtxCorr = fParam17l; }
  else if (periodStr.Contains("17m")) { runList = runList17m; nRuns = sizeof(runList17m)/sizeof(runList17m[0]); fParamZvtxCorr = fParam17m; }
  else if (periodStr.Contains("17o1")) { runList = runList17o1; nRuns = sizeof(runList17o1)/sizeof(runList17o1[0]); fParamZvtxCorr = fParam17o1; }
  else if (periodStr.Contains("17o2")) { runList = runList17o2; nRuns = sizeof(runList17o2)/sizeof(runList17o2[0]); fParamZvtxCorr = fParam17o2; }
  else if (periodStr.Contains("17o3")) { runList = runList17o3; nRuns = sizeof(runList17o3)/sizeof(runList17o3[0]); fParamZvtxCorr = fParam17o3; }
  else if (periodStr.Contains("17r")) { runList = runList17r; nRuns = sizeof(runList17r)/sizeof(runList17r[0]); fParamZvtxCorr = fParam17r; }

  else if (periodStr.Contains("16h")) { runList = runList16h; nRuns = sizeof(runList16h)/sizeof(runList16h[0]); fParamZvtxCorr = fParam16h; }
  else if (periodStr.Contains("16j")) { runList = runList16j; nRuns = sizeof(runList16j)/sizeof(runList16j[0]); fParamZvtxCorr = fParam16j; }
  else if (periodStr.Contains("16k1")) { runList = runList16k1; nRuns = sizeof(runList16k1)/sizeof(runList16k1[0]); fParamZvtxCorr = fParam16k1; }
  else if (periodStr.Contains("16k2")) { runList = runList16k2; nRuns = sizeof(runList16k2)/sizeof(runList16k2[0]); fParamZvtxCorr = fParam16k2; }
  else if (periodStr.Contains("16o")) { runList = runList16o; nRuns = sizeof(runList16o)/sizeof(runList16o[0]); fParamZvtxCorr = fParam16o; }
  else if (periodStr.Contains("16p")) { runList = runList16p; nRuns = sizeof(runList16p)/sizeof(runList16p[0]); fParamZvtxCorr = fParam16p; }
  else {
    printf("unknown period: %s",fperiod);
    return;
  }


}

void treeReaderHelper::SetListV0Norm(TFile *fileNormV0)
{
  listV0CNorm = (TList *)fileNormV0->Get("V0C_mean_zvtx_List");

  //the function for the zvtx correction of V0C
  fZvtxCorr = new TF1("f1","(1.+[0]*x+[1]*x*x+[2]*x*x*x)",-10,10);
  //set its parameters, derived from a fit done on the TProfile mean V0C versus zvtx
  //fit done on each period
  fZvtxCorr->SetParameter(0,fParamZvtxCorr[0]);
  fZvtxCorr->SetParameter(1,fParamZvtxCorr[1]);
  fZvtxCorr->SetParameter(2, fParamZvtxCorr[2]);

}

void treeReaderHelper::SetTreeAddresses(TTree *trDiMu, TTree *trSingleMu, TTree *trInt7)
{
  trDiMu->SetBranchAddress("fRunN",&fRunN);
  //tr->SetBranchAddress("fispileupspd",&fispileupspd);
  trDiMu->SetBranchAddress("fv0multTot",&fv0multTot);
  trDiMu->SetBranchAddress("fv0multcorr",&fv0multcorr);
  trDiMu->SetBranchAddress("fv0cmultTot",&fv0cmultTot);
  trDiMu->SetBranchAddress("fv0cmultcorr",&fv0cmultcorr);
  trDiMu->SetBranchAddress("fitsmult",&fitsmult);
  trDiMu->SetBranchAddress("fv0mpercentile",&fv0mpercentile);
  //trDiMu->SetBranchAddress("fv0apercentile",&fv0apercentile);
  //trDiMu->SetBranchAddress("fv0cpercentile",&fv0cpercentile);
  trDiMu->SetBranchAddress("fcl1percentile",&fcl1percentile);
  trDiMu->SetBranchAddress("fspdpercentile",&fspdpercentile);
  trDiMu->SetBranchAddress("fzvtx",&fzvtx);
  trDiMu->SetBranchAddress("ftracksmu",&ftracksmu);
  trDiMu->SetBranchAddress("fv0mult",fv0mult);
  trDiMu->SetBranchAddress("fSPDMult",&fSPDMult);

  trSingleMu->SetBranchAddress("fv0multSMu",fv0multSMu);
  trSingleMu->SetBranchAddress("etaMu",&etaMu);
  trSingleMu->SetBranchAddress("phiMu",&phiMu);
  trSingleMu->SetBranchAddress("ptMu",&ptMu);


  trInt7->SetBranchAddress("fv0multTotINT7",&fv0multTotINT7);
  trInt7->SetBranchAddress("fv0multCorrINT7",&fv0multCorrINT7);
  trInt7->SetBranchAddress("fv0cmultTotINT7",&fv0cmultTotINT7);
  trInt7->SetBranchAddress("fv0cmultCorrINT7",&fv0cmultCorrINT7);
  trInt7->SetBranchAddress("fv0mpercentile7",&fv0mpercentile7);
  trInt7->SetBranchAddress("fv0apercentile7",&fv0apercentile7);
  trInt7->SetBranchAddress("fv0cpercentile7",&fv0cpercentile7);
  trInt7->SetBranchAddress("fcl1percentile7",&fcl1percentile7);
  trInt7->SetBranchAddress("fspdpercentile7",&fspdpercentile7);
  trInt7->SetBranchAddress("fzvtxkINT7",&fzvtxkINT7);
  trInt7->SetBranchAddress("fRunNInt7",&fRunNInt7);
  trInt7->SetBranchAddress("fSPDMult7",&fSPDMult7);
}

void treeReaderHelper::GetNV0Histogram(TList *histoList, TTree *trInt7)
{
  //gets the NV0/<NV0> histogram
  //obtains the V0 signal per channel when no muon in the event ,
  //and fill fV0MultPerChannelNoMu[ichannel] and fEtaPhiV0MeanPerChannelNoMu

  hV0MultTot=(TH1F *)histoList->FindObject("hV0MultTot");//get the histograms from the histo list
  hV0MultCorr=(TH1F *)histoList->FindObject("hV0MultCorr");

  hV0CMultTot=(TH1F *)histoList->FindObject("hV0CMultTot");//get the histograms from the histo list
  hV0CMultCorr=(TH1F *)histoList->FindObject("hV0CMultCorr");

  float meanNV0Tot = hV0MultTot->GetMean();//=<NV0>
  float meanNV0Corr = hV0MultCorr->GetMean();//=<NV0>

  fMeanV0C_corr=hV0CMultCorr->GetMean();//=<NV0C>


  float NV0OverNV0Mean, NV0COverNV0CMean;

  Long64_t nentries = trInt7->GetEntries();
   for (Long64_t i=0;i<nentries;i++)
   {
     trInt7->GetEntry(i);
     NV0OverNV0Mean=(fv0multCorrINT7*1.0)/meanNV0Corr;
     fNV0OverNV0Mean->Fill(NV0OverNV0Mean);

     NV0COverNV0CMean=(fv0cmultCorrINT7*1.0)/fMeanV0C_corr;
     fNV0COverNV0CMean->Fill(NV0COverNV0CMean);

     hV0CMultCorrCopy->Fill(fv0cmultCorrINT7);

     NevZvtx->Fill(fzvtxkINT7);
     V0CMultvsZvtx->Fill(fzvtxkINT7,fv0cmultTotINT7);
     V0CMultCorrvsZvtx->Fill(fzvtxkINT7,fv0cmultCorrINT7);

     V0MultvsZvtx->Fill(fzvtxkINT7,fv0multTotINT7);
     V0MultCorrvsZvtx->Fill(fzvtxkINT7,fv0multCorrINT7);


     Float_t corrNormV0C = GetV0CMeanCorrection(listV0CNorm, fv0cmultCorrINT7, fRefV0C, fRunNInt7, fzvtxkINT7);
     hV0CAfterNorm->Fill(corrNormV0C);

     V0CMultCorrRenormvsZvtx->Fill(fzvtxkINT7,corrNormV0C);


     hSPDMultVsV0CMult->Fill(corrNormV0C, fSPDMult7);


     //TEST: only for 18m
     if (fRunNInt7==290223)
     {
       fNV0COverGlobalMean_290223->Fill((corrNormV0C*1.0)/fMeanV0C_corr);
     }

     if (fRunNInt7==292273)
     {
       fNV0COverGlobalMean_292273->Fill((corrNormV0C*1.0)/fMeanV0C_corr);
     }

     //End of TEST
   }

   //obtains the V0 signal per channel when no muon in the event ,
   //and fill fV0MultPerChannelNoMu[ichannel] and fEtaPhiV0MeanPerChannelNoMu

   for (int ichannel=0; ichannel<32; ichannel++)//saving only the histograms with V0C channels
   {
     //V0 mult in each channel when there is no muon in the event
     fV0MultPerChannelNoMu[ichannel]=(TH1F *)histoList->FindObject(Form("hV0MultPerChannelNoMu_%d",ichannel));

     //no muon in event
     float etaMi = AliVVZERO::GetVZEROEtaMin(ichannel) + 0.25;
     float phiMi = AliVVZERO::GetVZEROAvgPhi(ichannel);

     float meanV0MultPerChannel = fV0MultPerChannelNoMu[ichannel]->GetMean();

     fEtaPhiV0MeanPerChannelNoMu->Fill(etaMi, phiMi, meanV0MultPerChannel);//the mean multiplicity in each channel
   }


}


void treeReaderHelper::DeriveExcessV0MeanPerChannel()
{
  //computes the difference fEtaPhiV0MeanPerChannelOneMu[k]-fEtaPhiV0MeanPerChannelNoMu
  //and fills it in the hDiffEtaPhiMeanV0[k] histograms

  float excess =0.0;

  for (int ichannel=0; ichannel<32; ichannel++)
  {
    for (int etabin=1; etabin<fEtaPhiV0MeanPerChannelNoMu->GetNbinsX()+1; etabin++)
    {
      for (int phibin=1; phibin<fEtaPhiV0MeanPerChannelNoMu->GetNbinsY()+1; phibin++)
      {
        excess = fEtaPhiV0MeanPerChannelOneMu[ichannel]->GetBinContent(etabin, phibin) - fEtaPhiV0MeanPerChannelNoMu->GetBinContent(etabin, phibin);

        hDiffEtaPhiMeanV0[ichannel]->SetBinContent(etabin, phibin, excess);
      }
    }
  }


}

void treeReaderHelper::readEvents(TTree *tr)
{
  fillCorrFactorMap();
  //read events one by one and calls the different processFunctions for each events
  Long64_t nentries = tr->GetEntries();
   for (Long64_t i=0;i<nentries;i++)
   {
     tr->GetEntry(i);

     processDimuons();//fills some dimuon histograms and applies dimuon cuts

   }

}//readTree

void treeReaderHelper::writeOutput(const char *outputFileName)
{
  TFile* outputfile = new TFile(outputFileName, "RECREATE");
  fJpsiMinv->Write();
  fJpsiPt->Write();
  fJpsiY->Write();

  hMultPtJpsiMinv->Write();
  hMultPtJpsiMinvCorr->Write();
  hMultPtJpsiMinvCorrAxE->Write();
  hMultPtJpsiMinvCorrAxERenorm->Write();

  hMultPtJpsiMinvAllCorrY1->Write();
  hMultPtJpsiMinvAllCorrY2->Write();

  for (int ibinY=0; ibinY<5; ibinY++)//5 bins in rapidity
  {
    hMultPtJpsiMinvAllCorrY[ibinY]->Write();
  }

  hMultMinusMultCorr->Write();
  hMultMinusMultCorr1mu->Write();

  fNV0OverNV0Mean->Write();
  fNV0COverNV0CMean->Write();

  hV0CMultCorrCopy->Write();
  hV0CAfterNorm->Write();

  hSPDMultVsV0CMult->Write();

  fNV0COverGlobalMean_292273->Write();
  fNV0COverGlobalMean_290223->Write();

  hV0CMultInSingleMu->Write();

  hMuonExcess->Write();

  V0CMultvsZvtx->Write();
  V0CMultCorrvsZvtx->Write();
  V0CMultCorrRenormvsZvtx->Write();

  V0MultvsZvtx->Write();
  V0MultCorrvsZvtx->Write();

  NevZvtx->Write();

  gDirectory->mkdir("V0MultPerChannel");
  outputfile->cd("V0MultPerChannel");

  for (int ichannel=0; ichannel<32; ichannel++)//saving only the histograms with V0C channels
  {
    fV0MultPerChannelNoMu[ichannel]->Write();//same for one muon and two muons
    fV0MultPerChannelOneMu[ichannel]->Write();
    //fV0MultPerChannelTwoMu[ichannel]->Write();
    //fEtaPhiV0PerChannelOneMu[ichannel]->Write();//useless now, since we have the mean below
    fEtaPhiV0MeanPerChannelOneMu[ichannel]->Write();

    //hDiffEtaPhiMeanV0[ichannel]->Write();
    fDiMuCorrectionWDiMuEv[ichannel]->Write();
    fDiMuCorrectionWSiMuEv[ichannel]->Write();

    fEtaPhiMeanV0DiMuonEv[ichannel]->Write();
  }
  fEtaPhiV0MeanPerChannelNoMu->Write();

  outputfile->Close();
}

void treeReaderHelper::processV0PerChannelSingleMu(TTree *trSingleMu)
{
  //obtains the V0 channels which have muon signal, and fills the histograms
  //fV0MultPerChannelOneMu and fEtaPhiV0PerChannelOneMu

  Long64_t nentries = trSingleMu->GetEntries();
   for (Long64_t i=0;i<nentries;i++)
   {
     trSingleMu->GetEntry(i);

     int channelMu = GetChannelFromEtaPhi(etaMu,phiMu);//V0 channel where the single muon has hit
     if (channelMu<0)
     {
       //the GetChannelFromEtaPhi returned -1, eta or phi out of range
       continue;
     }

     fV0MultPerChannelOneMu[channelMu]->Fill(fv0multSMu[channelMu]);

     //temporary: fill the raw tot mult in V0C in single muon events (with one muon in V0C acc)
     float totMultSingleMu=0;

     //now we fill the histogram fEtaPhiV0PerChannelOneMu[idChannelMu] with eta, phi and v0mult for all channels
     float etaMid, phiMid; // eta and phi of the middle of each channel
     for (int ichannel=0; ichannel<32; ichannel++)//if ichannel>31, the signal is in V0A, no muons there
     {//could be made faster with a map filled in the beginning
       etaMid = AliVVZERO::GetVZEROEtaMin(ichannel) + 0.25;
       phiMid = AliVVZERO::GetVZEROAvgPhi(ichannel);
       fEtaPhiV0PerChannelOneMu[channelMu]->Fill(etaMid, phiMid, fv0multSMu[ichannel]);//now useless

       fEtaPhiV0MeanPerChannelOneMu[channelMu]->Fill(etaMid, phiMid, fv0multSMu[ichannel]);//since it is a TProfile2D, it will compute the mean V0 mult in each x,y bin
       //by doing Z(x,y)=sum(z(x,y))/nb of entries in (x,y bin) I guess ??

       totMultSingleMu+=fv0multSMu[ichannel];
     }

     hV0CMultInSingleMu->Fill(totMultSingleMu);

   }




  //after that we have to take the mean in each eta,phi bin to fill the 32 fEtaPhiV0MeanPerChannelOneMu
  //histograms--> automatically done with the TProfile2D

}

void treeReaderHelper::fillCorrFactorMap()
{
  //fills the map of correction factor (derived with single muon histos)
  //for each channel Id
	float corrFactor=0.;
  std::pair<int, int> channelsAcross;
  float meanPedestal;
  float meanSignal;
	for (int ichannel=0; ichannel<32; ichannel++)
  {
    float etaMi = AliVVZERO::GetVZEROEtaMin(ichannel) + 0.25;
    float phiMi = AliVVZERO::GetVZEROAvgPhi(ichannel);

    channelsAcross=fChannelAccrossChannelId[ichannel];


    meanPedestal = fEtaPhiV0MeanPerChannelOneMu[channelsAcross.first]->GetBinContent(fEtaPhiV0MeanPerChannelOneMu[channelsAcross.first]->FindBin(etaMi,phiMi));
    meanPedestal+= fEtaPhiV0MeanPerChannelOneMu[channelsAcross.second]->GetBinContent(fEtaPhiV0MeanPerChannelOneMu[channelsAcross.second]->FindBin(etaMi,phiMi));

    meanPedestal=meanPedestal/2.0;

    meanSignal=fEtaPhiV0MeanPerChannelOneMu[ichannel]->GetBinContent(fEtaPhiV0MeanPerChannelOneMu[ichannel]->FindBin(etaMi,phiMi));

    corrFactor=meanSignal-meanPedestal;
    //computed the corrfactor and now store it in fChannelToCorrFactor[ichannel]
	 fChannelToCorrFactor[ichannel]=corrFactor;
   hMuonExcess->Fill(ichannel, corrFactor);
  }

}



void treeReaderHelper::processDimuons()
{
  //================================= Dimuons =================================

  for (Int_t iTrMA1=0; iTrMA1<(ftracksmu->GetEntriesFast()-1); iTrMA1++)
  { // first muon
    AliAODTrack *track1 = (AliAODTrack*) ftracksmu->UncheckedAt(iTrMA1);
    for (Int_t iTrMA2=iTrMA1+1; iTrMA2<ftracksmu->GetEntriesFast(); iTrMA2++)
    { // second muon
      AliAODTrack *track2 = (AliAODTrack*) ftracksmu->UncheckedAt(iTrMA2);

      if ((track1->Charge()*track2->Charge()) >= 0) continue;//same charge muons, not interesting

      Double_t minv = CalcMinv(track1,track2);//invariant mass
      if (minv<fLowMassCut || minv>fUpMassCut) continue;//cut on invariant mass



      Double_t pT = CalcPt(track1,track2);//pT of dimuon

      Double_t y = CalcY(track1,track2);//rapidity of dimuon
      if (y < fJpsiEtaMin || y > fJpsiEtaMax) continue;//cut on rapidity


      //al the cuts have been applied, now fill histograms

      fJpsiPt->Fill(pT);
      fJpsiMinv->Fill(minv);
      fJpsiY->Fill(y);

      hMultPtJpsiMinv->Fill(fv0cmultcorr,pT,minv);

      int nMuonHitV0C = -1; //number of muon from the Jpsi decay hitting the V0C


      Float_t correctedfv0cmultcorr = correctDimuons(track1,track2, nMuonHitV0C);

      Float_t correctedrenorm_fv0cmultcorr=GetV0CMeanCorrection(listV0CNorm, correctedfv0cmultcorr, fRefV0C, fRunN, fzvtx);

      hMultPtJpsiMinvCorr->Fill(correctedfv0cmultcorr,pT,minv);

      Float_t jpsiAxE = hJpsiAxE->GetBinContent(hJpsiAxE->FindBin(pT,y));
      if(jpsiAxE<=0)
      {
        jpsiAxE=1e-6;
      }

      hMultPtJpsiMinvCorrAxE->Fill(correctedfv0cmultcorr,pT,minv, 1./jpsiAxE);

      hMultMinusMultCorr->Fill(fv0cmultcorr-correctedfv0cmultcorr);

      hMultPtJpsiMinvCorrAxERenorm->Fill(correctedrenorm_fv0cmultcorr,pT,minv, 1./jpsiAxE);

      if(y<-3.2)
      {
        hMultPtJpsiMinvAllCorrY1->Fill(correctedrenorm_fv0cmultcorr,pT,minv, 1./jpsiAxE);
      }
      else
      {
        hMultPtJpsiMinvAllCorrY2->Fill(correctedrenorm_fv0cmultcorr,pT,minv, 1./jpsiAxE);
      }

      //filling the inv-mass 3D histograms with different y ranges for J/psi

      float yMinBin = -4.0;

      for (int ibinY=0; ibinY<5; ibinY++)//5 bins in rapidity
      {
        if(y>yMinBin && y<=(yMinBin+0.3))
        {
          hMultPtJpsiMinvAllCorrY[ibinY]->Fill(correctedrenorm_fv0cmultcorr,pT,minv, 1./jpsiAxE);
        }
        yMinBin=yMinBin+0.3;

      }



    }// second muon
  }// first muon

}

Float_t treeReaderHelper::correctDimuons(AliAODTrack *track1, AliAODTrack *track2, int &nMuons)
{

  Float_t correctedformuons_fv0mult;

  //---------------------------------- muon 1 ----------------------------------
  float eta1 = track1->Eta();
  float phi1 = track1->Phi();

  int channel1 = GetChannelFromEtaPhi(eta1, phi1);//channel hit by muon 1
  //if not -1

  //Note: the histograms of the mean value of mult for each fv0 channels are already filled
  //the histograms named fEtaPhiV0MeanPerChannelOneMu[ichannel], these are used in
  //fChannelToCorrFactor

  float corrFactor1=0.;
  if (channel1 != -1)
  {
    corrFactor1 = fChannelToCorrFactor[channel1];//delta_i
  }

  //---------------------------------- muon 2 ----------------------------------
  float eta2 = track2->Eta();
  float phi2 = track2->Phi();

  int channel2 = GetChannelFromEtaPhi(eta2, phi2);//channel hit by muon 2
  //if not -1

  float corrFactor2=0.;
  if (channel2 != -1)
  {
    corrFactor2 = fChannelToCorrFactor[channel2];//delta_i
  }

  //----------------------------------------------------------------------------

  int type=0;//1 is muon in channel 1 only, 2 is muon in channel 2 only, 3 is muon in channel 1 & 2

  if (channel1 == -1 && channel2 == -1)//no muon hit the V0C
  {
    nMuons=0;
    return fv0cmultcorr;//return the mult without additionnal correction
  }

  float V0CMult=fv0cmultTot;
  //recorrect with ESD utils the V0C mult, with fzvtx

  //since fv0cmultTot is just the sum of the signal in each channel, I can just substract the corrections


  if (channel1 != -1 && channel2 == -1)
  {
    //muon 1 only hit channel
    fDiMuCorrectionWSiMuEv[channel1]->Fill(fv0mult[channel1]-corrFactor1);
    V0CMult=V0CMult-corrFactor1;
    hMultMinusMultCorr1mu->Fill(corrFactor1);

    nMuons=1;
  }

  if (channel1 == -1 && channel2 != -1)
  {
    //muon 2 only hit channel
    fDiMuCorrectionWSiMuEv[channel2]->Fill(fv0mult[channel2]-corrFactor2);
    V0CMult=V0CMult-corrFactor2;
    hMultMinusMultCorr1mu->Fill(corrFactor2);

    nMuons=1;
  }

  if (channel1 != -1 && channel2 != -1)//both muons hit
  {
    fDiMuCorrectionWSiMuEv[channel1]->Fill(fv0mult[channel1]-corrFactor1);
    fDiMuCorrectionWSiMuEv[channel2]->Fill(fv0mult[channel2]-corrFactor2);

    V0CMult=V0CMult-corrFactor1-corrFactor2;

    nMuons=2;
  }

  correctedformuons_fv0mult = AliESDUtils::GetCorrV0C(V0CMult,fzvtx);


  return correctedformuons_fv0mult;
}

Float_t treeReaderHelper::GetV0CMeanCorrection(TList* fV0List, Float_t V0Corr, Float_t refMult_V0C, Int_t runNb, Float_t zvtx)
{
  //fV0List is the list containing the histos per run
  //V0CCorr is th V0C mult corrected by ESDUtils
  TProfile *profileV0=(TProfile *)fV0List->FindObject(Form("%d_V0C_mean_zvtx",runNb));



  if (fV0List->FindObject(Form("%d_V0C_mean_zvtx",runNb))==0)
  {
    //the profile is not available for this run
    throw std::exception();
  }

  Float_t meanV0InRun = profileV0->GetMean(2);


  Float_t V0CorrNorm = V0Corr*(refMult_V0C/meanV0InRun);
  return V0CorrNorm/(fZvtxCorr->Eval(zvtx));
}


//===============================================================================

Double_t treeReaderHelper::CalcMinv(AliAODTrack *track1, AliAODTrack *track2) const
{
  return TMath::Sqrt((track1->E()+track2->E())*(track1->E()+track2->E())-
		     ((track1->Px()+track2->Px())*(track1->Px()+track2->Px())+
		      (track1->Py()+track2->Py())*(track1->Py()+track2->Py())+
		      (track1->Pz()+track2->Pz())*(track1->Pz()+track2->Pz())));
}

Double_t treeReaderHelper::CalcY(AliAODTrack *track1, AliAODTrack *track2) const
{
  Double_t e = track1->E()+track2->E();
  Double_t pz = track1->Pz()+track2->Pz();
  Double_t y;
  if(TMath::Abs(e-pz)>1e-7)
    y = 0.5*TMath::Log((e+pz)/(e-pz));
  else
    y = -1e6;

  return y;
}


Double_t treeReaderHelper::CalcPhi(AliAODTrack *track1, AliAODTrack *track2) const
{
  Double_t px = track1->Px() + track2->Px();
  Double_t py = track1->Py() + track2->Py();
  Double_t phi = ((px*px+py*py) != 0) ? TMath::Pi()+TMath::ATan2(-py, -px) : -1e6; // phi

  return phi;
}

Double_t treeReaderHelper::CalcPt(AliAODTrack *track1, AliAODTrack *track2) const
{
  return TMath::Sqrt((track1->Px()+track2->Px())*
		     (track1->Px()+track2->Px())+
		     (track1->Py()+track2->Py())*
		     (track1->Py()+track2->Py()));
}

Double_t treeReaderHelper::CalcCos(AliAODTrack *track1, AliAODTrack *track2) const
{
  Double_t pt1 = TMath::Sqrt(track1->Px()*track1->Px()+track1->Py()*track1->Py());
  Double_t pt2 = TMath::Sqrt(track2->Px()*track2->Px()+track2->Py()*track2->Py());
  if (pt1 > 0. && pt2 > 0.)
    return (track1->Px()*track2->Px()+track1->Py()*track2->Py())/pt1/pt2;
  else
    return -100.;
}

int treeReaderHelper::GetChannelFromEtaPhi(float eta, float phi)
{
  if (eta>-1.7 || eta<-3.7)
  {
    //eta out of range
    return -1;
  }

  if (phi>2*TMath::Pi() || phi<0)
  {
    //phi out of range
    printf(">>>>>>> WARNING: phi = %f out of range\n", phi);
    return -1;
  }

  int i = 0;
  float realEtaMin = -3.7;
  while (eta>etaRanges[i][1])
  {
    i++;
  }

  realEtaMin=etaRanges[i][0];

  int j = 0;
  float realPhiMin = 0;
  while (phi>phiRanges[j])
  {
    j++;
  }
  realPhiMin = phiRanges[j-1];

  if (j-1 < 0)
  {
    return -1;
  }
  return fChannelMapV0[{i,j-1}];
}

//=============================================================================

void treeReaderHelper::bringTo02Pi(float& phi)
{
  // ensure angle in [0:2pi] for the any input angle
  while (phi < 0.f) {
    phi += 2*TMath::Pi();
  }
  while (phi > 2*TMath::Pi()) {
    phi -= 2*TMath::Pi();
  }
}
