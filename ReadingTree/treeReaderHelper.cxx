#include "treeReaderHelper.h"
#include "AliAODTrack.h"
#include "AliVVZERO.h"

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

  float meanNV0CCorr = hV0CMultCorr->GetMean();//=<NV0>


  float NV0OverNV0Mean, NV0COverNV0CMean;

  Long64_t nentries = trInt7->GetEntries();
   for (Long64_t i=0;i<nentries;i++)
   {
     trInt7->GetEntry(i);
     NV0OverNV0Mean=(fv0multCorrINT7*1.0)/meanNV0Corr;
     fNV0OverNV0Mean->Fill(NV0OverNV0Mean);

     NV0COverNV0CMean=(fv0cmultCorrINT7*1.0)/meanNV0CCorr;
     fNV0COverNV0CMean->Fill(NV0COverNV0CMean);

     hV0CMultTotCopy->Fill(fv0cmultTotINT7);
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
  //read events one by one and calls the different processFunctions for each events
  Long64_t nentries = tr->GetEntries();
   for (Long64_t i=0;i<nentries;i++)
   {
     tr->GetEntry(i);

     //printf(">>>>>>> Event %lld fv0cmultcorr = %f\n", i, fv0cmultcorr);

     processDimuons();//fills some dimuon histograms and applies dimuon cuts

     //processV0PerChannel();

   }


}//readTree

void treeReaderHelper::writeOutput(const char *outputFileName)
{
  TFile* outputfile = new TFile(outputFileName, "RECREATE");
  fJpsiMinv->Write();
  fJpsiPt->Write();

  hMultPtJpsiMinv->Write();

  fNV0OverNV0Mean->Write();
  fNV0COverNV0CMean->Write();

  

  gDirectory->mkdir("V0MultPerChannel");
  outputfile->cd("V0MultPerChannel");

  for (int ichannel=0; ichannel<32; ichannel++)//saving only the histograms with V0C channels
  {
    fV0MultPerChannelNoMu[ichannel]->Write();//same for one muon and two muons
    fV0MultPerChannelOneMu[ichannel]->Write();
    //fV0MultPerChannelTwoMu[ichannel]->Write();
    //fEtaPhiV0PerChannelOneMu[ichannel]->Write();//useless now, since we have the mean below
    fEtaPhiV0MeanPerChannelOneMu[ichannel]->Write();

    hDiffEtaPhiMeanV0[ichannel]->Write();
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

     //now we fill the histogram fEtaPhiV0PerChannelOneMu[idChannelMu] with eta, phi and v0mult for all channels
     float etaMid, phiMid; // eta and phi of the middle of each channel
     for (int ichannel=0; ichannel<32; ichannel++)//if ichannel>31, the signal is in V0A, no muons there
     {//could be made faster with a map filled in the beginning
       etaMid = AliVVZERO::GetVZEROEtaMin(ichannel) + 0.25;
       phiMid = AliVVZERO::GetVZEROAvgPhi(ichannel);
       fEtaPhiV0PerChannelOneMu[channelMu]->Fill(etaMid, phiMid, fv0multSMu[ichannel]);//now useless

       fEtaPhiV0MeanPerChannelOneMu[channelMu]->Fill(etaMid, phiMid, fv0multSMu[ichannel]);//since it is a TProfile2D, it will compute the mean V0 mult in each x,y bin
       //by doing Z(x,y)=sum(z(x,y))/nb of entries in (x,y bin) I guess ??
     }

   }




  //after that we have to take the mean in each eta,phi bin to fill the 32 fEtaPhiV0MeanPerChannelOneMu
  //histograms--> automatically done with the TProfile2D

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
      fJpsiPt->Fill(pT);
      // Int_t ptBinJpsi = fJpsiPtAxis -> FindBin(pT);
      // if (ptBinJpsi<1 || ptBinJpsi>fNbinsJpsiPt) continue;//underflow or overflow , ignore
      Double_t y = CalcY(track1,track2);//rapidity of dimuon
      if (y < fJpsiEtaMin || y > fJpsiEtaMax) continue;//cut on rapidity


      //al the cuts have been applied, now fill histograms
      fJpsiMinv->Fill(minv);

      //printf("---------------fv0cmultcorr = %f, pT = %f, minv = %f\n", fv0cmultcorr, pT, minv);

      hMultPtJpsiMinv->Fill(fv0cmultcorr,pT,minv);



    }// second muon
  }// first muon

}

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
  if (eta>-1.7 && eta<-3.7)
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
