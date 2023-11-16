const int npx = 5000;//number of points for TF1
enum paramCBE {kNorm, kMean, kSigma, kAlpha, kN, kAlphap, kNp, kNparamsCBE};
int nParamsPol2Exp = 4;//same for fVWG
bool fRejectFitPoints = true;//reject point for background functions
double fFitRejectRangeLow = 2.7;
double fFitRejectRangeHigh =  3.3;

double fFitRangeLow;
double fFitRangeHigh;

enum bkgFunction {kPol2Exp, kVWG};

// From 13 TeV pp data fit
// Double_t CB2alpha = 9.79765e-01;
// Double_t CB2n = 6.96799e+00;
// Double_t CB2alpha2 = 1.86328e+00;
// Double_t CB2n2 = 1.49963e+01;

//From the 13 TeV pp fit of Manuel Guittiere's AN
Double_t CB2alpha;
Double_t CB2n;
Double_t CB2alpha2;
Double_t CB2n2;

Double_t functionSignalCrystalBallExtended(Double_t *x,Double_t *par);//signal
Double_t functionBackgroundPol2Exp(Double_t *x, Double_t *par);
Double_t functionSignal2CBEBackgroundPol2Exp(Double_t *x, Double_t *par);//signal+bkground
Double_t fVWG(Double_t *x, Double_t *par);//Variable width gaussain, bkground funtion
Double_t functionSignalCB2BackgroundVWG(Double_t *x, Double_t *par);

void fitJpsi1histogram(int i, int j,TH1D* hJpsiMinv, unsigned int bkgFunction);


TCanvas *cJpsiResults[12][5];
TCanvas *cJpsiResultsMean[5];//mean values of Jpsi (integrated in mult)

double intJpsi[12][5];//contains the number of jpsi for each mult and pt bins
double intJpsiMean[5];//contains the number of jpsi for each pt bins (integrated in mult)
TH2D *hIntJpsi;//the weight is the jpsi integral, for the corresponding (mult,pt) bin
TH1D *hIntJpsiMean;//the weight is the jpsi integral, for the corresponding pt bin

void initFitParams(unsigned int tailParam, unsigned int fitRange)//bkgFunction not used here
{
  //for the systematic errors, different tail params, bkg functions and fit ranges can be chosen
  if (!tailParam)//tailParam==0, tail parameters from data
  {
    //From the 13 TeV pp data fit of Manuel Guittiere's AN
    CB2alpha = 0.883;
    CB2n = 9.940;
    CB2alpha2 = 1.832;
    CB2n2 = 15.323;
  }
  else
  {
    //tailParam==1, tail parameters from MC
    //From the 13 TeV pp MC fit of Manuel Guittiere's AN
    CB2alpha = 0.993;
    CB2n = 2.9075;
    CB2alpha2 = 2.182;
    CB2n2 = 3.122;
  }

  if (!fitRange)//fitRange==0, first fit range
  {
    fFitRangeLow=2.0;
    fFitRangeHigh=4.9;
  }
  else
  {
    //fitRange==1, second fit range
    fFitRangeLow=2.1;
    fFitRangeHigh=5.0;
  }


}

void fitJpsi(unsigned int tailParam, unsigned int bkgFunction, unsigned int fitRange, const char *inputFile = "../ReadingTree/JpsiRead_18c.root", const char *outputFile = "FitResults_18c.root", const char *histname="hMultPtJpsiMinv")
{


  TFile* inputf = TFile::Open(inputFile,"READ");

  TH1F* hV0CMultCorr= (TH1F*)inputf->Get("hV0CMultCorrCopy");

  TH3D* hJpsi3D = (TH3D*)inputf->Get(histname);

  //use hJpsi3D binning to create the binning for hIntJpsi
  TAxis *xaxis = hJpsi3D->GetXaxis();//mult bins from hJpsi3D
  TAxis *yaxis = hJpsi3D->GetYaxis();//pT bins from hJpsi3D


  hIntJpsi = new TH2D("hIntJpsi", "hIntJpsi", hJpsi3D->GetNbinsX(), xaxis->GetXbins()->GetArray(), hJpsi3D->GetNbinsY(), yaxis->GetXbins()->GetArray());
  //number of jpsi versus mult and pt

  hIntJpsiMean = new TH1D("hIntJpsiMean", "hIntJpsiMean", hJpsi3D->GetNbinsY(), yaxis->GetXbins()->GetArray());//number of Jpsi vs pt

  const int nBinsMult = hJpsi3D->GetNbinsX();
  const int nBinsPt = hJpsi3D->GetNbinsY();

  TH1D* hJpsi3D_pz[nBinsMult][nBinsPt];

  TH1D* hJpsi3D_pz_mean[nBinsPt];//used to calculate the mean value of Jpsi integrated over all mult classes

  initFitParams(tailParam,fitRange);

  TFile* outputfile = new TFile(outputFile, "RECREATE");


  for (int i=0;i<nBinsMult;i++)
  {
    //i: index for mult bin
    for (int j=0;j<nBinsPt;j++)
    {
      //j: index for pt jpsi bin
      hJpsi3D_pz[i][j]=hJpsi3D->ProjectionZ(Form("_pz_%d_%d",i+1,j+1),i+1,i+1,j+1,j+1);//i and j are bin numbers, they begin at 1
      hJpsi3D_pz[i][j]->SetTitle("");
      printf("Fitting minv histogram for multbin = %d, ptbin = %d\n", i+1,j+1);
      fitJpsi1histogram(i,j,hJpsi3D_pz[i][j], bkgFunction);






      cJpsiResults[i][j]->Write();
    }
  }



  hIntJpsi->Write();


  //computing the mean values, integrated over mult
  int i = -1;

  for (int j=0;j<nBinsPt;j++)
  {
    //j: index for pt jpsi bin
    hJpsi3D_pz_mean[j]=hJpsi3D->ProjectionZ(Form("_pz_mean_%d",j+1),1,nBinsMult,j+1,j+1);//i and j are bin numbers, they begin at 1
    hJpsi3D_pz_mean[j]->SetTitle("");
    printf("Fitting minv histogram for integrated multbins, ptbin = %d\n",j+1);
    fitJpsi1histogram(-1,j,hJpsi3D_pz_mean[j], bkgFunction);

    cJpsiResultsMean[j]->Write();
  }



  hIntJpsiMean->Write();

  hV0CMultCorr->Write("hV0CMultCorr");

  outputfile->Close();

}//void fitJpsi

void fitJpsi1histogram(int i, int j,TH1D* hJpsiMinv, unsigned int bkgFunction)//fit the histogram corresponding to bin i in mult and j in pt
{
  //the bkground function can change between 0 (Pol2Exp) and 1 (VWG)

  //--------------------- Fit background function ----------------------------

  TF1 *fitBkgLeft, *fitBkgRight, *fitBkgPrelim;
  if(!bkgFunction)
  {
    //bkgFunction==0, will be a Pol2Exp
    fitBkgLeft = new TF1("fitBkgLeft",functionBackgroundPol2Exp,fFitRangeLow,fFitRejectRangeLow,4);//xmin, xmax, nparameters
    fitBkgRight = new TF1("fitBkgRight",functionBackgroundPol2Exp,fFitRejectRangeHigh,fFitRangeHigh,4);//xmin, xmax, nparameters
    fitBkgPrelim = new TF1("fitBkgPrelim",functionBackgroundPol2Exp,fFitRangeLow,fFitRangeHigh,4);//xmin, xmax, nparameters
  }
  else
  {
    //bkgFunction==1, will be a VWG
    fitBkgLeft = new TF1("fitBkgLeft",fVWG,2,fFitRejectRangeLow,4);//xmin, xmax, nparameters
    fitBkgRight = new TF1("fitBkgRight",fVWG,fFitRejectRangeHigh,fFitRangeHigh,4);//xmin, xmax, nparameters
    fitBkgPrelim = new TF1("fitBkgPrelim",fVWG,2,5,4);//xmin, xmax, nparameters

    fitBkgLeft->SetParameter(0,2*hJpsiMinv->GetBinContent(hJpsiMinv->FindBin(2.0)));//normalisation
    fitBkgLeft->SetParameter(1,-1.0);//mean
    fitBkgLeft->SetParLimits(1,-10.0,2.0);//mean limits
    fitBkgLeft->SetParameter(2,0.5);//sigma A
    fitBkgLeft->SetParLimits(0,0.5*hJpsiMinv->GetBinContent(hJpsiMinv->FindBin(2.0)),200*hJpsiMinv->GetBinContent(hJpsiMinv->FindBin(2.0)));// limits of N

    fitBkgRight->SetParameter(0,2*hJpsiMinv->GetBinContent(hJpsiMinv->FindBin(2.0)));//normalisation
    fitBkgRight->SetParameter(1,-1.0);//mean
    fitBkgRight->SetParLimits(1,-10.0,2.0);//mean limits
    fitBkgRight->SetParameter(2,0.5);//sigma A
    fitBkgRight->SetParLimits(0,0.5*hJpsiMinv->GetBinContent(hJpsiMinv->FindBin(2.0)),200*hJpsiMinv->GetBinContent(hJpsiMinv->FindBin(2.0)));// limits of N

  }

  //---- step 1: get the parameters for the left and right fit



  fitBkgLeft->SetNpx(npx);
  fitBkgLeft->SetLineWidth(4);
  fitBkgLeft->SetLineColor(kBlack);

  hJpsiMinv->Fit("fitBkgLeft","0q","ep", 2, fFitRejectRangeLow);

  double parBkgLeft[4];
  fitBkgLeft->GetParameters(parBkgLeft);


  fitBkgRight->SetNpx(npx);
  fitBkgRight->SetLineWidth(4);
  fitBkgRight->SetLineColor(kBlack);

  hJpsiMinv->Fit("fitBkgRight","0q","ep", fFitRejectRangeHigh, 5);

  double parBkgRight[4];
  fitBkgRight->GetParameters(parBkgRight);

  //---- step 2: combined fit on left and right (if done without previous fit parameters, may not converge)

  fRejectFitPoints=true;


  fitBkgPrelim->SetNpx(npx);
  fitBkgPrelim->SetLineWidth(4);
  fitBkgPrelim->SetLineColor(kCyan);

  //add the pol x Exp parameters

  fitBkgPrelim->SetParameter(0,(parBkgRight[0]+parBkgLeft[0])/2.0);
  fitBkgPrelim->SetParameter(1,(parBkgRight[1]+parBkgLeft[1])/2.0);
  fitBkgPrelim->SetParameter(2,(parBkgRight[2]+parBkgLeft[2])/2.0);
  fitBkgPrelim->SetParameter(3,(parBkgRight[3]+parBkgLeft[3])/2.0);

  hJpsiMinv->Fit("fitBkgPrelim","0q","ep", 2, 5);//hJpsiMinv->Fit("fitBkgPrelim","V+","ep", 2, 5);

  fRejectFitPoints=false;

  double parBkgPrelim[4];
  fitBkgPrelim->GetParameters(parBkgPrelim);

  //--------------------- With added CBE for psi(2S) -------------------------

  TF1 *fitSignalAndBkg2;
  if(!bkgFunction)
  {
    //bkgFunction==0, will be a Pol2Exp
    fitSignalAndBkg2 = new TF1("fitSignalAndBkg2",functionSignal2CBEBackgroundPol2Exp,fFitRangeLow,fFitRangeHigh,8);//xmin, xmax, nparameters
  }
  else
  {
    //bkgFunction==1, will be a VWG
    fitSignalAndBkg2 = new TF1("fitSignalAndBkg2",functionSignalCB2BackgroundVWG,fFitRangeLow,fFitRangeHigh,8);//xmin, xmax, nparameters
  }

  //TF1 *fitSignalAndBkg2 = new TF1("fitSignalAndBkg2",functionSignal2CBEBackgroundPol2Exp,fFitRangeLow,fFitRangeHigh,8);//xmin, xmax, nparameters
  fitSignalAndBkg2->SetNpx(npx);
  fitSignalAndBkg2->SetLineWidth(4);
  fitSignalAndBkg2->SetLineColor(kBlack);

  //For Jpsi
  fitSignalAndBkg2->SetParameter(nParamsPol2Exp+kNorm,hJpsiMinv->GetBinContent(hJpsiMinv->FindBin(3.096)));//from previous fit fitSignalAndBkg
  fitSignalAndBkg2->SetParameter(nParamsPol2Exp+kMean,3.096); // mean jpsi
  fitSignalAndBkg2->SetParLimits(nParamsPol2Exp+kMean,2.9,3.3);// limits of mean jpsi
  fitSignalAndBkg2->SetParameter(nParamsPol2Exp+kSigma,0.070);   // sigma Jpsi
  fitSignalAndBkg2->SetParLimits(nParamsPol2Exp+kSigma,0.035,0.14);// limits of sigma jpsi

  //add the background function parameters

  // fitSignalAndBkg2->SetParameter(0,parBkgPrelim[0]);
  // fitSignalAndBkg2->SetParameter(1,parBkgPrelim[1]);
  // fitSignalAndBkg2->SetParameter(2,parBkgPrelim[2]);
  // fitSignalAndBkg2->SetParameter(3,parBkgPrelim[3]);

  fitSignalAndBkg2->FixParameter(0,parBkgPrelim[0]);
  fitSignalAndBkg2->FixParameter(1,parBkgPrelim[1]);
  fitSignalAndBkg2->FixParameter(2,parBkgPrelim[2]);
  fitSignalAndBkg2->FixParameter(3,parBkgPrelim[3]);

  fitSignalAndBkg2->SetParameter(7,0.5*hJpsiMinv->GetBinContent(hJpsiMinv->FindBin(3.69))); //normalisation of Psi2s

  if (i==-1)
  {
    cJpsiResultsMean[j]= new TCanvas(Form("cJpsiResultsMean_%d",j+1),"Fitting CBE for Jpsi and psi2S and bkg",10,10,700,500);
    cJpsiResultsMean[j]->SetGrid();
  }
  else
  {
    cJpsiResults[i][j] = new TCanvas(Form("cJpsiResults_%d_%d",i+1,j+1),"Fitting CBE for Jpsi and psi2S and bkg",10,10,700,500);
    cJpsiResults[i][j]->SetGrid();
  }


  //hJpsiMinv->Fit("fitSignalAndBkg2","V+","ep");
  TFitResultPtr resultFit = hJpsiMinv->Fit("fitSignalAndBkg2","qS","ep",fFitRangeLow,fFitRangeHigh);//S to store the fitResults in the pointer
  int statusfit = resultFit->Status();

  if(statusfit==4)
  {
    resultFit = hJpsiMinv->Fit("fitSignalAndBkg2","S","ep",fFitRangeLow,fFitRangeHigh);//S to store the fitResults in the pointer
  }

  printf("------------------------------------ FIT RESULT STATUS %d, mul bin %d, pt bin %d\n", resultFit->Status(), i+1, j+1);
  //--------------------------get the different background and signal functions

  //Define the 3 functions, bkg + jpsi signal + psi(2S) signal
  TF1 *fitBkg;
  if(!bkgFunction)
  {
    //bkgFunction==0, will be a Pol2Exp
    fitBkg = new TF1("fitBkg",functionBackgroundPol2Exp,fFitRangeLow,fFitRangeHigh,4);//xmin, xmax, nparameters
  }
  else
  {
    //bkgFunction==1, will be a VWG
    fitBkg = new TF1("fitBkg",fVWG,fFitRangeLow,fFitRangeHigh,4);//xmin, xmax, nparameters
  }

  //TF1 *fitBkg = new TF1("fitBkg",functionBackgroundPol2Exp,fFitRangeLow,fFitRangeHigh,4);//xmin, xmax, nparameters
  fitBkg->SetNpx(npx);
  fitBkg->SetLineWidth(4);
  fitBkg->SetLineColor(kBlue);

  TF1 *fitJpsiSignal = new TF1("fitJpsiSignal",functionSignalCrystalBallExtended,fFitRangeLow,fFitRangeHigh,7);//7 parameters
  fitJpsiSignal->SetLineColor(kRed);
  fitJpsiSignal->SetLineWidth(4);
  fitJpsiSignal->SetNpx(npx);

  TF1 *fitPsi2SSignal = new TF1("fitPsi2SSignal",functionSignalCrystalBallExtended,fFitRangeLow,fFitRangeHigh,7);//7 parameters
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
  float errJpsi;
  float intJpsiTemp = fitJpsiSignal->Integral(2,4)/hJpsiMinv->GetBinWidth(1);
  //used for the error in case the fit doesn't converge

  //covariant matrix for the error: always symmetric
  TMatrixDSym cv = resultFit->GetCovarianceMatrix();
  Double_t covmat[7][7];//size = n params for the jpsi fit
  for(Int_t i = 0; i < 7; ++i)
  {
    for(Int_t j = 0; j < 7; ++j)
    {
      if(i<=2 && j<=2)
      {
        covmat[i][j] = cv(4+i,4+j);//why the 4 ?? -- only the parameters we are interested about
        //ie parameters for the Jpsi CBE fit
      }
      else
      {
        //the alpha and n parameters, which are fixed, so their covariances is 0
        covmat[i][j] = 0;
      }
      //printf("%.3e  ",covmat[i][j]);


    }
    //printf("\n");
  }

  if (covmat[0][0] == 0) {
     // abnormal fit termination
    printf("Abnormal fit termination.\n");
    errJpsi=TMath::Sqrt(intJpsiTemp);
  }
  else
    errJpsi = fitJpsiSignal->IntegralError(fFitRangeLow,4,&parJpsi[7],&covmat[0][0])/hJpsiMinv->GetBinWidth(1);
    if (isnan(errJpsi))//happens when the integrand is not behaving as expected
    {
      errJpsi = fitJpsiSignal->IntegralError(2.0,3.4,&parJpsi[7],&covmat[0][0])/hJpsiMinv->GetBinWidth(1);
      printf("---- Error was nan, range changed\n");
      if (isnan(errJpsi))//happens when the integrand is not behaving as expected
      {
        errJpsi=TMath::Sqrt(intJpsiTemp);
        printf("---- Error still nan, give up\n");
      }

    }


  if (i ==-1)//we are integrating over the mult
  {
    intJpsiMean[j] = fitJpsiSignal->Integral(2,4)/hJpsiMinv->GetBinWidth(1);
    hIntJpsiMean->SetBinContent(j+1,fitJpsiSignal->Integral(2,4)/hJpsiMinv->GetBinWidth(1));
    hIntJpsiMean->SetBinError(j+1,errJpsi);
  }
  else//diff in mult and pt
  {
    intJpsi[i][j] = fitJpsiSignal->Integral(2,4)/hJpsiMinv->GetBinWidth(1);
    hIntJpsi->SetBinContent(i+1,j+1,fitJpsiSignal->Integral(2,4)/hJpsiMinv->GetBinWidth(1));
    hIntJpsi->SetBinError(i+1,j+1,errJpsi);
  }


  TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  legend->AddEntry(hJpsiMinv,"Data","lpe");
  legend->AddEntry(fitBkg,"Background fit","l");
  legend->AddEntry(fitJpsiSignal,"J/#psi fit","l");
  legend->AddEntry(fitPsi2SSignal,"#psi(2S) fit","l");
  legend->AddEntry(fitSignalAndBkg2,"Global fit","l");
  legend->Draw();

  TLatex* tl;

  if(i==-1)
  {
    tl = new TLatex(0.6022945,0.2536998,Form("#it{N_{J/#psi}} = %f", intJpsiMean[j]));
  }
  else{
    tl = new TLatex(0.6022945,0.2536998,Form("#it{N_{J/#psi}} = %f", intJpsi[i][j]));
  }



  tl->SetNDC();
  tl->SetTextFont(42);
  tl->SetTextSize(0.0422833);
  tl->SetLineWidth(2);
  tl->Draw();



}

//CB2
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

Double_t fVWG(Double_t *x, Double_t *par)
{
  //Variable Width Gaussian with pol1 sigma
  //4 parameters
  Double_t m = x[0];
  Double_t N = par[0];
  Double_t mean = par[1];
  Double_t A = par[2];
  Double_t B = par[3];

  Double_t t = (m-mean)/mean;
  Double_t sigma = A + B*t;

  if (fRejectFitPoints &&  x[0] > fFitRejectRangeLow && x[0] < fFitRejectRangeHigh )
  {
    TF1::RejectPoint();
    return 0.;
  }

  return N*TMath::Exp(-(m-mean)*(m-mean)/(2.*sigma*sigma));
}

Double_t functionSignalCB2BackgroundVWG(Double_t *x, Double_t *par)
{
  //4 params for the VWG, 3 for the first CBE (Jpsi) and 1 for the normalisation of psi(2S)
  //8 parameters in total
  //p4 = norm Jpsi
  //p5 = mean Jpsi
  //p6= sigma Jpsi
  Double_t parJpsi[7]={par[4], par[5], par[6], CB2alpha, CB2n, CB2alpha2, CB2n2};
  Double_t parPsi2S[7]={par[7], par[5]+0.589, par[6]*1.05, CB2alpha, CB2n, CB2alpha2, CB2n2};
  //Mean psi(2S) = Mean(Jpsi)+ 0.589 GeV/c^2
  //Sigma psi(2S) = Sigma(Jpsi)*1.05
  return fVWG(x,par) + functionSignalCrystalBallExtended(x,parJpsi) + functionSignalCrystalBallExtended(x,parPsi2S);
  //bckg + jpsi + psi(2S)
}
