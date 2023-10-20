void testProjection(const char *inputFile = "../ReadingTree/JpsiRead_18c.root")
{


  TFile* inputf = TFile::Open(inputFile,"READ");
  //TH1D* fJpsiMinv = (TH1D*)inputf->Get((ddir + "/EventSelection").c_str());//later: for each bin of NV0/<NV0>

  TH1D* hJpsiMinv = (TH1D*)inputf->Get("fJpsiMinv");

  TH3D* hJpsi3D = (TH3D*)inputf->Get("hMultPtJpsiMinv");
  TH1D* hJpsi3D_pz = hJpsi3D->ProjectionZ("_pz",1,1,1,1);

  TH1D* hJpsi3D_py = (TH1D*) hJpsi3D->Project3D("z");


  printf("NbinsZ = %d\n", hJpsi3D->GetNbinsZ() );

  TCanvas *cProj = new TCanvas("cProj","cProj",10,10,700,500);
  //hJpsi3D_pz->Draw();//doesn't draw anything
  hJpsi3D_py->Draw();


  float massBinWidth=0.05;
  float binEdgesMinv[61];
  float valueBinEdge=2;
  int ibin=0;
  while(valueBinEdge<5.05)//fUpMassCut+massBinWidth
  {
    binEdgesMinv[ibin]=valueBinEdge;
    valueBinEdge=valueBinEdge+massBinWidth;
    ibin++;
  }

  printf("binEdgesMinv %s", "=");
  for (int i=0;i<61;i++){
    printf("%f, ", binEdgesMinv[i]);
  }
}
