int colors[3]={kRed, kBlue, kMagenta};

void deriveResults(const char *inputFile = "../fitJpsi/FitResults_18d.root")
{

  gStyle->SetOptStat("");
  TFile* inputf = TFile::Open(inputFile,"READ");

  TH1F* hV0CMultCorr= (TH1F*)inputf->Get("hV0CMultCorr");
  TH2D* hIntJpsi= (TH2D*)inputf->Get("hIntJpsi");
  TH1D* hIntJpsiMean= (TH1D*)inputf->Get("hIntJpsiMean");

  const int nBinsMult = hIntJpsi->GetNbinsX();
  const int nBinsPt = hIntJpsi->GetNbinsY();

  std::vector<double>* vrelatV0C[nBinsPt];//vrelatV0C=NV0/<NV0> averaged in each ranges
  std::vector<double>* vrelatNJpsi[nBinsPt];//vrelatNJpsi
  std::vector<double>* verrors[nBinsPt];//errors on Jpsi
  std::vector<double>* vwidths[nBinsPt];//widths on the x axis

  //Fill the values of mean NV0C/<NV0C> for each mult class, and the corresponding Jpsi yield




  TGraphErrors* result[nBinsPt];//TGraph containing the result for each pt bin

  double njpsi;//number of jpsi in one mult bin
  double Njpsi;//number of jpsi integrated in mult
  double nevkInt7;//number of kINT7 events in one mult bin
  double NevkInt7;//number of kINT7 events integrated in mult

  double lowerEdgeNV0C;//lower edge of one mult bin
  double upperEdgeNV0C;//upper edge of one mult bin

  for (int j=0;j<nBinsPt;j++)
  {
    //j+1: index for pt jpsi bin
    vrelatV0C[j]=new std::vector<double>;
    vrelatNJpsi[j]=new std::vector<double>;
    verrors[j]=new std::vector<double>;
    vwidths[j]=new std::vector<double>;

    Njpsi=hIntJpsiMean->GetBinContent(j+1);
    NevkInt7=hV0CMultCorr->Integral();

    printf("--------------- >>> PT BIN %d\n", j+1);


    for (int i=0;i<nBinsMult;i++)
    {
      //i+1: index for mult bin


      njpsi=hIntJpsi->GetBinContent(i+1,j+1);

      int globalBin=hIntJpsi->GetBin(i+1,j+1);

      lowerEdgeNV0C=hIntJpsi->GetXaxis()->GetBinLowEdge(i+1);//mult low edge limit for NV0C
      upperEdgeNV0C=hIntJpsi->GetXaxis()->GetBinLowEdge(i+1)+ hIntJpsi->GetXaxis()->GetBinWidth(i+1);

      nevkInt7=hV0CMultCorr->Integral(hV0CMultCorr->FindBin(lowerEdgeNV0C), hV0CMultCorr->FindBin(upperEdgeNV0C));

      double sumNV0CInRange=0;

      for(int k=hV0CMultCorr->FindBin(lowerEdgeNV0C);k<=hV0CMultCorr->FindBin(upperEdgeNV0C);k++)
      {
        Float_t cont = hV0CMultCorr->GetBinContent(k);
        Float_t mult = hV0CMultCorr->GetXaxis()->GetBinCenter(k);
        sumNV0CInRange += (cont*mult);
      }

      double binwidthNV0COverMeanNV0C=(upperEdgeNV0C-lowerEdgeNV0C)/hV0CMultCorr->GetMean();



      float errnJpsi=hIntJpsi->GetBinError(i+1,j+1);
      float errNJpsi=hIntJpsiMean->GetBinError(j+1);

      double nV0Crel=(sumNV0CInRange/nevkInt7)/(hV0CMultCorr->GetMean());
      double relatJpsiYield=(njpsi/nevkInt7)/(Njpsi/NevkInt7);


      //compute the statistical error on relative Jpsi yield
      float err=relatJpsiYield*sqrt(pow(errnJpsi/njpsi,2)+pow(errNJpsi/Njpsi,2));

      if (binwidthNV0COverMeanNV0C>100)//for the last point
      {
        binwidthNV0COverMeanNV0C=(nV0Crel-lowerEdgeNV0C/hV0CMultCorr->GetMean())*2;
      }

      vrelatV0C[j]->emplace_back((sumNV0CInRange/nevkInt7)/(hV0CMultCorr->GetMean()));
      vrelatNJpsi[j]->emplace_back((njpsi/nevkInt7)/(Njpsi/NevkInt7));
      verrors[j]->emplace_back(err);
      vwidths[j]->emplace_back( binwidthNV0COverMeanNV0C/ 2.);

      printf("--- multBin = %d , nV0Crel = %f, relatJpsiYield = %f, error =%f\n", i+1, nV0Crel, relatJpsiYield, err);


    }


    result[j]= new TGraphErrors(vrelatV0C[j]->size(), vrelatV0C[j]->data(), vrelatNJpsi[j]->data(), vwidths[j]->data(), verrors[j]->data());
    result[j]->SetName(Form("JpsivsV0M_ptBin_%d",j));
    result[j]->SetTitle("dN_{J/#psi}");
    result[j]->SetMarkerStyle(8);
    result[j]->SetMarkerSize(1.);
    result[j]->SetMarkerColor(colors[j]);
  }




  auto* axes = new TH1F("plot","",36, 0, 7);
  axes->GetYaxis()->SetRangeUser(0,9);
  axes->GetXaxis()->SetTitle("#frac{#it{N_{V0C}}}{<#it{N_{V0C}}>}");
  axes->GetYaxis()->SetTitle("#frac{dN_{J/#psi}/d#it{N_{V0C}}}{<dN_{J/#psi}/d#it{N_{V0C}}>}");
  axes->GetXaxis()->SetTitleOffset(1.2);


  TCanvas* canvasResult = new TCanvas("canvasResult", "canvasResult",0,53,800,600);
  canvasResult->Range(-0.8708919,-1.408696,7.870892,9.841304);
  canvasResult->SetGrid();
  canvasResult->SetTopMargin(0.07478261);
  canvasResult->SetBottomMargin(0.1252174);

  axes->DrawCopy();

  TLegend *legend=new TLegend(0.1453634,0.6747826,0.4260652,0.8747826,NULL,"brNDC");
  legend->SetBorderSize(0);

  for (int iresult=0; iresult<nBinsPt; iresult++)
  {
    result[iresult]->Draw("P same");

    float lowerEdgePtBin=hIntJpsi->GetYaxis()->GetBinLowEdge(iresult+1);
    float upperEdgePtBin=hIntJpsi->GetYaxis()->GetBinLowEdge(iresult+1)+ hIntJpsi->GetYaxis()->GetBinWidth(iresult+1);

    if (iresult+1 == nBinsPt)
    {
      //last pt bin, pt > 3 GeV/c
      legend->AddEntry(result[iresult],Form("#it{p_{T}} > %.1f GeV/#it{c}",lowerEdgePtBin),"lpe");
    }
    else
    {
      legend->AddEntry(result[iresult],Form("%.1f #leq #it{p_{T}} #leq %.1f GeV/#it{c}",lowerEdgePtBin,upperEdgePtBin),"lpe");
    }
  }

  legend->Draw();

}
