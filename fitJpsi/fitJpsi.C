const int npx = 5000;//number of points for TF1
enum paramCBE {kNorm, kMean, kSigma, kAlpha, kN, kAlphap, kNp, kNparamsCBE};
int nParamsPol2Exp = 4;
bool fRejectFitPoints = true;
double fFitRejectRangeLow = 2.7;
double fFitRejectRangeHigh =  3.3;

// From 13 TeV pp data fit
Double_t CB2alpha = 9.79765e-01;
Double_t CB2n = 6.96799e+00;
Double_t CB2alpha2 = 1.86328e+00;
Double_t CB2n2 = 1.49963e+01;

Double_t functionSignalCrystalBallExtended(Double_t *x,Double_t *par);
Double_t functionBackgroundPol2Exp(Double_t *x, Double_t *par);
Double_t functionSignalCBEBackgroundPol2Exp(Double_t *x, Double_t *par);
Double_t functionSignal2CBEBackgroundPol2Exp(Double_t *x, Double_t *par);

void fitJpsi(const char *inputFile = "../ReadingTree/JpsiRead_18c.root")
{


  TFile* inputf = TFile::Open(inputFile,"READ");

  TH1D* hJpsiMinv = (TH1D*)inputf->Get("fJpsiMinv");

  TH3D* hJpsi3D = (TH3D*)inputf->Get("hMultPtJpsiMinv");
  TH1D* hJpsi3D_pz = hJpsi3D->ProjectionZ("_pz",1,1,1,1);



  //--------------------- Fit background function ----------------------------

  //---- step 1: get the parameters for the left and right fit

  TCanvas *cBKG = new TCanvas("cBKG","Fitting background func",10,10,700,500);

  TF1 *fitBkgLeft = new TF1("fitBkgLeft",functionBackgroundPol2Exp,2,fFitRejectRangeLow,4);//xmin, xmax, nparameters
  fitBkgLeft->SetNpx(npx);
  fitBkgLeft->SetLineWidth(4);
  fitBkgLeft->SetLineColor(kBlack);

  hJpsiMinv->Fit("fitBkgLeft","0","ep", 2, fFitRejectRangeLow);

  double parBkgLeft[4];
  fitBkgLeft->GetParameters(parBkgLeft);

  TF1 *fitBkgRight = new TF1("fitBkgRight",functionBackgroundPol2Exp,fFitRejectRangeHigh,5,4);//xmin, xmax, nparameters
  fitBkgRight->SetNpx(npx);
  fitBkgRight->SetLineWidth(4);
  fitBkgRight->SetLineColor(kBlack);

  hJpsiMinv->Fit("fitBkgLeft","0","ep", fFitRejectRangeHigh, 5);

  double parBkgRight[4];
  fitBkgRight->GetParameters(parBkgRight);

  //---- step 2: combined fit on left and right (if done without previous fit parameters, may not converge)

  fRejectFitPoints=true;

  TF1 *fitBkgPrelim = new TF1("fitBkgPrelim",functionBackgroundPol2Exp,2,5,4);//xmin, xmax, nparameters
  fitBkgPrelim->SetNpx(npx);
  fitBkgPrelim->SetLineWidth(4);
  fitBkgPrelim->SetLineColor(kCyan);

  //add the pol x Exp parameters

  fitBkgPrelim->SetParameter(0,(parBkgRight[0]+parBkgLeft[0])/2.0);
  fitBkgPrelim->SetParameter(1,(parBkgRight[1]+parBkgLeft[1])/2.0);
  fitBkgPrelim->SetParameter(2,(parBkgRight[2]+parBkgLeft[2])/2.0);
  fitBkgPrelim->SetParameter(3,(parBkgRight[3]+parBkgLeft[3])/2.0);

  hJpsiMinv->Fit("fitBkgPrelim","V+","ep", 2, 5);

  fRejectFitPoints=false;

  double parBkgPrelim[4];
  fitBkgPrelim->GetParameters(parBkgPrelim);

  //--------------------- With added CBE for psi(2S) -------------------------

  TF1 *fitSignalAndBkg2 = new TF1("fitSignalAndBkg2",functionSignal2CBEBackgroundPol2Exp,2,5,8);//xmin, xmax, nparameters
  fitSignalAndBkg2->SetNpx(npx);
  fitSignalAndBkg2->SetLineWidth(4);
  fitSignalAndBkg2->SetLineColor(kBlack);

  //For Jpsi
  fitSignalAndBkg2->SetParameter(nParamsPol2Exp+kNorm,hJpsiMinv->GetBinContent(hJpsiMinv->FindBin(3.096)));//from previous fit fitSignalAndBkg
  fitSignalAndBkg2->SetParameter(nParamsPol2Exp+kMean,3.096); // mean jpsi
  fitSignalAndBkg2->SetParLimits(nParamsPol2Exp+kMean,2.9,3.3);// limits of mean jpsi
  fitSignalAndBkg2->SetParameter(nParamsPol2Exp+kSigma,0.070);   // sigma Jpsi
  fitSignalAndBkg2->SetParLimits(nParamsPol2Exp+kSigma,0.035,0.14);// limits of sigma jpsi

  //add the pol x Exp parameters

  fitSignalAndBkg2->SetParameter(0,parBkgPrelim[0]);
  fitSignalAndBkg2->SetParameter(1,parBkgPrelim[1]);
  fitSignalAndBkg2->SetParameter(2,parBkgPrelim[2]);
  fitSignalAndBkg2->SetParameter(3,parBkgPrelim[3]);

  fitSignalAndBkg2->SetParameter(7,0.5*hJpsiMinv->GetBinContent(hJpsiMinv->FindBin(3.69))); //normalisation of Psi2s


  TCanvas *c3 = new TCanvas("c3","Fitting CBE for Jpsi and psi2S and bkg",10,10,700,500);
  c3->SetGrid();

  hJpsiMinv->Fit("fitSignalAndBkg2","V+","ep");

  //--------------------------get the different background and signal functions

  //Define the 3 functions, bkg + jpsi signal + psi(2S) signal

  TF1 *fitBkg = new TF1("fitBkg",functionBackgroundPol2Exp,2,5,4);//xmin, xmax, nparameters
  fitBkg->SetNpx(npx);
  fitBkg->SetLineWidth(4);
  fitBkg->SetLineColor(kBlue);

  TF1 *fitJpsiSignal = new TF1("fitJpsiSignal",functionSignalCrystalBallExtended,2,5,7);//7 parameters
  fitJpsiSignal->SetLineColor(kRed);
  fitJpsiSignal->SetLineWidth(4);
  fitJpsiSignal->SetNpx(npx);

  TF1 *fitPsi2SSignal = new TF1("fitPsi2SSignal",functionSignalCrystalBallExtended,2,5,7);//7 parameters
  fitPsi2SSignal->SetLineColor(kGreen);
  fitPsi2SSignal->SetLineWidth(4);
  fitPsi2SSignal->SetNpx(npx);




  double parAll[8];//all parameters of the total fit
  fitSignalAndBkg2->GetParameters(parAll);

  Double_t parJpsi[7]={parAll[4], parAll[5], parAll[6], CB2alpha, CB2n, CB2alpha2, CB2n2};
  Double_t parPsi2S[7]={parAll[7], parAll[5]+0.589, parAll[6]*1.05, CB2alpha, CB2n, CB2alpha2, CB2n2};


  fitBkg->SetParameters(parAll);//will take only the first 4 parameters
  fitBkg->Draw("same");

  fitJpsiSignal->SetParameters(parJpsi);
  fitJpsiSignal->Draw("same");

  fitPsi2SSignal->SetParameters(parPsi2S);
  fitPsi2SSignal->Draw("same");


  // calculate the Jpsi integral and error
  Double_t intJpsi = fitJpsiSignal->Integral(2,5)/hJpsiMinv->GetBinWidth(1);

  TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->AddEntry(hJpsiMinv,"Data","lpe");
  legend->AddEntry(fitBkg,"Background fit","l");
  legend->AddEntry(fitJpsiSignal,"J/#psi fit","l");
  legend->AddEntry(fitPsi2SSignal,"#psi(2S) fit","l");
  legend->AddEntry(fitSignalAndBkg2,"Global fit","l");
  legend->Draw();






}

//I guess it is CB2
Double_t functionSignalCrystalBallExtended(Double_t *x,Double_t *par)
{
  // Extended Crystal Ball : 7 parameters
  //
  // par[0] = Normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = alpha
  // par[4] = n
  // par[5] = alpha'
  // par[6] = n'

  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
    return par[0]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail
  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;
    return par[0]*(a/TMath::Power(b - t, par[4]));
  }

  if (t >= absAlpha2) //right tail
  {

    Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
    Double_t d = par[6]/absAlpha2 - absAlpha2;
    return par[0]*(c/TMath::Power(d + t, par[6]));
  }

  return 0. ;
}
// taken from https://github.com/alisw/AliPhysics/blob/master/PWG/muondep/AliAnalysisMuMuJpsiResult.cxx

Double_t functionBackgroundPol2Exp(Double_t *x, Double_t *par)
{
  // pol2 x exp : 4 params

  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }
  return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0])*TMath::Exp(par[3]*x[0]);
}
// taken from https://github.com/alisw/AliPhysics/blob/master/PWG/muondep/AliAnalysisMuMuJpsiResult.cxx without point rejection

Double_t functionSignalCBEBackgroundPol2Exp(Double_t *x, Double_t *par)
{
  //4 params for the pol2exp and 7 for CBE
  return functionBackgroundPol2Exp(x,par) + functionSignalCrystalBallExtended(x,&par[4]);
}


Double_t functionSignal2CBEBackgroundPol2Exp(Double_t *x, Double_t *par)
{
  //4 params for the pol2exp 3 for the first CBE (Jpsi) and 1 for the normalisation of psi(2S)
  //8 parameters in total
  //p4 = norm Jpsi
  //p5 = mean Jpsi
  //p6= sigma Jpsi
  Double_t parJpsi[7]={par[4], par[5], par[6], CB2alpha, CB2n, CB2alpha2, CB2n2};
  Double_t parPsi2S[7]={par[7], par[5]+0.589, par[6]*1.05, CB2alpha, CB2n, CB2alpha2, CB2n2};
  //Mean psi(2S) = Mean(Jpsi)+ 0.589 GeV/c^2
  //Sigma psi(2S) = Sigma(Jpsi)*1.05
  return functionBackgroundPol2Exp(x,par) + functionSignalCrystalBallExtended(x,parJpsi) + functionSignalCrystalBallExtended(x,parPsi2S);
  //bckg + jpsi + psi(2S)
}
