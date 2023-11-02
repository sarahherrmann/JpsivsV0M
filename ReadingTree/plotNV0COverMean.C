void plotNV0COverMean()
{
  const char* fInName[4]={"JpsiRead_18d.root", "JpsiRead_18e.root", "JpsiRead_18f.root", "JpsiRead_18l.root"};
  const char* periodName[4]={"18d","18e","18f","18l"};
  int colors[4]={kMagenta, kBlue, kGreen, kRed};
  TFile *fIn[4];

  TH1F *fNV0COverNV0CMean[4];


  for (int i=0; i<4; i++)
  {
    fIn[i] = TFile::Open(fInName[i]);

    fNV0COverNV0CMean[i]= (TH1F*)fIn[i]->Get("fNV0COverNV0CMean");
    fNV0COverNV0CMean[i]->SetName(Form("fNV0COverNV0CMean%s",periodName[i]));
    fNV0COverNV0CMean[i]->SetLineColor(colors[i]);
  }

  TCanvas *cPercentile = new TCanvas("cPercentile","cPercentile",700,500);

  TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
  legend->SetBorderSize(0);

  gStyle->SetOptStat("");

  for (int i=0; i<4; i++)
  {
    if (i==0)
    {
      fNV0COverNV0CMean[i]->DrawNormalized("le");
    }
    else{
      fNV0COverNV0CMean[i]->DrawNormalized("sameLE");
    }

    legend->AddEntry(fNV0COverNV0CMean[i],Form("LHC%s",periodName[i]),"le");
  }





  legend->Draw();



}
