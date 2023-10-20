#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TDatabasePDG.h>
#include <TH2F.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TString.h>
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TFractionFitter.h>
#include "Fit/FitConfig.h"
#include "Fit/Fitter.h"
using namespace std;
#include <fstream>


TH1D *hStatic[3][4]; // static histos for MC components used in the global fit, 3 sets, one for each fit - gen, mc and data
TH1D *hPtSingleMu = 0;
Int_t ReadDataOrMC(const char *fpath, Bool_t isMC,
		   Double_t mumuPtmin, Double_t mumuPtmax, Double_t muPtmin,
		   Int_t n1, Int_t n2,
		   TH2D *hThetaRec, TH2D *hPhiRec, TH2D *hTildaRec, TH2D *hPtRec, TH2D *hYRec, TH1D *hMRec,
		   Double_t thetaMax,
		   TH2D *hThetaRec2 = 0x0, TH2D *hPhiRec2 = 0x0, TH2D *hTildaRec2 = 0x0,
		   TH1D *hTheta1True = 0x0, TH1D *hPhi1True = 0x0, TH1D *hTilda1True = 0x0, TH1D *hPt1True = 0x0,
		   TH1D *hTheta2True = 0x0, TH1D *hPhi2True = 0x0, TH1D *hTilda2True = 0x0, TH1D *hPt2True = 0x0,
		   TH2D *hThetaRes = 0x0, TH2D *hPhiRes = 0x0, TH2D *hTildaRes = 0x0, TH2D *hPtRes = 0x0, TH2D *hYRes = 0x0,
		   RooUnfoldResponse *responseTheta = 0x0, RooUnfoldResponse *responsePhi = 0x0, RooUnfoldResponse *responseTilda = 0x0,
		   Int_t iteration = 0, TF1 **fPt = 0x0, TF1 **fY = 0x0,
		   TH2D *hDecompThetaRec2 = 0x0, TH2D *hDecompPhiRec2 = 0x0, TH2D *hDecompTildaRec2 = 0x0,
		   TF1 *fUndo = 0x0,
		   Bool_t weight2015MC = kFALSE,
		 	Bool_t reweightTrig = kFALSE,
			Bool_t weightUndo = kTRUE);
Int_t ReadGen(const char *fpath,Double_t mumuPtmin, Double_t mumuPtmax,
	      Int_t n1, Int_t n2,
	      TH1D *hThetaGen1, TH1D *hPhiGen1, TH1D *hTildaGen1, TH1D *hPtGen1, TH1D *hYGen1,
	      TH1D *hThetaGen2, TH1D *hPhiGen2, TH1D *hTildaGen2, TH1D *hPtGen2, TH1D *hYGen2,
	      Int_t iteration, TF1 **fPt, TF1 **fY,
	      TH2D *hSL,
	      TH2D *hDecompThetaGen2, TH2D *hDecompPhiGen2, TH2D *hDecompTildaGen2,
	      TF1 *fUndo,
	      Bool_t weight2015MC, Bool_t weightUndo);
Double_t CalcCosThetaPhiHE(TLorentzVector &p1Mom, TLorentzVector &p2Mom, Double_t &phiHE, Double_t &tildaHE);
TH1D *ExtractJpsi(TH2D *hInputData, TH2D *hInputMC, Int_t extractOpt, Bool_t isMC, TList *loutput);
TF1 *FitJpsiMC(TH1D *histoToFit);
Double_t FitJpsi(TH1D *histoToFit, Int_t extractOpt, Double_t &errJpsi, TF1 *fJpsi, TH1D *hFitResults, const  char* nameV);
Double_t DiCrystalBall(double* ax, double* ap);
Double_t BackGround(Double_t *ax, Double_t *ap);
Double_t Pol1Pol2(Double_t *ax, Double_t *ap);
Double_t Pol2Exp(Double_t *ax, Double_t *ap);
Double_t TotalData(Double_t *x, Double_t *ap);
Double_t TotalData1(Double_t *x, Double_t *ap);
Double_t TotalData2(Double_t *x, Double_t *ap);
Double_t triggerError(Double_t *ax, Double_t *ap);
TH1D *GetCorrected(TH1D *hRaw, RooUnfoldResponse *response, TH1D *hAcc, TList *loutput);
TH1D *GetCorrectedWoUnf(TH1D *hRaw, TH1D *hAcc, TList *loutput);
Double_t GetWeight(Int_t iteration, Double_t pT, Double_t Y, TF1 **fPt, TF1 **fY);
void FitPolarizationParams(TH1D *hTheta, TH1D *hPhi, TH1D *hTilda, Double_t thetaMax);
TF1 *UndoPolarisation();
void Decomp(TH2D *hTheta, TH2D *hPhi, TH2D *hTilda,
	    Double_t theta, Double_t phi, Double_t tilda,
	    Double_t thetaTrue, Double_t phiTrue, Double_t tildaTrue,
	    Double_t w, Double_t wUndo);
void FractionFitting(const char *hname,
		     TH1D *hThetaData, TH1D *hPhiData, TH1D *hTildaData,
		     TH2D *hThetaMC, TH2D *hPhiMC, TH2D *hTildaMC,
		     TList *loutput,
		     Double_t thetaMax);
Double_t fracFCN(Double_t *x, Double_t *par);
void MyFractionFitting(const char *hname,
		       TH1D *hThetaData, TH1D *hPhiData, TH1D *hTildaData,
		       TH2D *hThetaMC, TH2D *hPhiMC, TH2D *hTildaMC,
		       TList *loutput,
		       Double_t thetaMax,
		       Int_t fitIndex);
TF1 *fTrig = 0;
// Explanation of input parameters (all params are used in the ame of output list with histos
// fdatapath, fmcpath - paths to files with data and MC trees
// mumuPtmin, mumuPtmax - dimuon pT range
// muPtmin - min muon pT cut
// extractOpt - option to be used in the extraction of Jpsi in data - fit functions, range, alterntive methods (to be defined as needed)
// iteration - the MC pT,y spectra differs from data, especially pT. So iteratively weight the MC in order to reproduce the spectra in data. iteration=0 is just the raw MC as it comes.
// n1, n2 - to be used in MC closure tests. n1 out of n2 MC tree entries will be used as 'data' and the rest as MC (for acceptance, response etc). In case n2 < 0, all MC entries are used as MC (this is the case when one studies data and MC closure).
void analyze(const char *fdatapath = "All_Data", const char *fmcpath = "New_petit",
	     Double_t mumuPtmin = 0.35, Double_t mumuPtmax = 2.00, Double_t muPtmin = 0.0,
	     Int_t extractOpt = 0,
	     Int_t iteration = 0,
	     Int_t n1 = 0, Int_t n2 = -1, Bool_t reweightTrig = kFALSE, const char *typeMCfit="cc")
{
  gSystem->Load("libRooUnfold");

  std::cout<<"ExtractOpt Choice"<<std::endl;
  if(extractOpt==0 || extractOpt==10 ||extractOpt==20 )    std::cout<<"Option Fit VWG"<<std::endl;
  else if(extractOpt==1|| extractOpt==11 ||extractOpt==21 )   std::cout<<"Option Fit Pol1Pol2"<<std::endl;
	else if(extractOpt==2 || extractOpt==12 ||extractOpt==22 )   std::cout<<"Option Fit Pol2Exp"<<std::endl;
	else{
    std::cout<<"No Background function known"<<std::endl;
    std::cout<<"Options possible :"<<std::endl;
    std::cout<< "0 VWG 2.2-6 "<<std::endl;
    std::cout<< "1 Pol1Pol2 2.2-6 "<<std::endl;
    std::cout<< "2 Pol2Exp 2.2-6 "<<std::endl;
    std::cout<< "10 VWG 2-5 "<<std::endl;
    std::cout<< "11 Pol1Pol2 2-5 "<<std::endl;
    std::cout<< "12 Pol2Exp 2-5 "<<std::endl;
		std::cout<< "20 VWG 2.5-4.5 "<<std::endl;
		std::cout<< "21 Pol1Pol2 2.5-4.5 "<<std::endl;
		std::cout<< "22 Pol2Exp 2.5-4.5 "<<std::endl;
    return;
  }

  Int_t nBinsTheta = 25, nBinsPhi = 24;
  Int_t nBinsPt =TMath::Nint((mumuPtmax-mumuPtmin)/0.1);
  Double_t mMin = 2.0, mMax = 6.0; // mass range
  Int_t nBinsM = TMath::Nint((mMax-mMin)/0.05);
  Double_t yMin = -4.15, yMax = -2.35;
  Int_t nBinsY = 12;
  Double_t thetaMax = 0.68;

  // Check if it is incoherent MC in order to apply a weight of 2 for 2015
  Bool_t weight2015MC = kFALSE;
  TString mcstr(fmcpath);
  if (mcstr.Contains("Incoh",TString::kIgnoreCase)) {
    printf("Incoherent MC detected, will apply a weight of 2 for 2015 run period.\n");
    weight2015MC = kTRUE;
    nBinsPhi = 24; // 16 or 20 by default? ...
    thetaMax = 0.68; // check this ...
  }
	Bool_t weightUndo = kTRUE;
	if (mcstr.Contains("New",TString::kIgnoreCase)) {
		printf("New MC detected.\n");
		weightUndo = kFALSE;
	}

  // read Pt,Y MC weights
  TF1 **fPt = 0x0;
  TF1 **fY = 0x0;
  if (iteration > 0) {
		TFile *fweights = 0;
		if(string(typeMCfit)=="cc"){
			fweights=TFile::Open("mcweights.root");
			cout<<"cc"<<endl;
		}else if(string(typeMCfit)=="uu"){
			fweights=TFile::Open("mcweights_uu.root");
			cout<<"uu"<<endl;
		}else if (string(typeMCfit)=="ll") {
			fweights=TFile::Open("mcweights_ll.root");
			cout<<"ll"<<endl;
		}else if (string(typeMCfit)=="lu") {
			fweights=TFile::Open("mcweights_lu.root");
			cout<<"lu"<<endl;
		}else if (string(typeMCfit)=="ul") {
			fweights=TFile::Open("mcweights_ul.root");
			cout<<"ul"<<endl;
		}else if (string(typeMCfit)=="cl") {
			fweights=TFile::Open("mcweights_cl.root");
			cout<<"cl"<<endl;
		}else if (string(typeMCfit)=="cu") {
			fweights=TFile::Open("mcweights_cu.root");
			cout<<"cu"<<endl;
		}else if (string(typeMCfit)=="lc") {
			fweights=TFile::Open("mcweights_lc.root");
			cout<<"lc"<<endl;
		}else if (string(typeMCfit)=="uc") {
			fweights=TFile::Open("mcweights_uc.root");
			cout<<"uc"<<endl;
		}

		fPt = new TF1*[iteration];
    fY = new TF1*[iteration];
    for(Int_t ii = 0; ii < iteration; ++ii) {
      fPt[ii] = (TF1 *)fweights->Get(Form("fPt_%s_%s_%.2f_%.2f_%.2f_%d_%d_%d_%d_%d_%s",
					  fdatapath,fmcpath,
					  mumuPtmin,mumuPtmax,muPtmin,
					  extractOpt,ii,n1,n2,reweightTrig,typeMCfit));
      fY[ii] = (TF1 *)fweights->Get(Form("fY_%s_%s_%.2f_%.2f_%.2f_%d_%d_%d_%d_%d_%s",
					 fdatapath,fmcpath,
					 mumuPtmin,mumuPtmax,muPtmin,
					 extractOpt,ii,n1,n2,reweightTrig,typeMCfit));
      if (fPt[ii] == 0x0 || fY[ii] == 0x0) {
	printf("MC weights for iteration %d not found. Exit.\n",ii);
	exit(0);
      }
    }
  }

  // Create function with which one undoes the polarisation in Starlight
  TF1 *fUndo = UndoPolarisation();

  TH1D *hThetaGen1 = new TH1D("hThetaGen1","",nBinsTheta,-1,1); hThetaGen1->Sumw2();
  TH1D *hPhiGen1 = new TH1D("hPhiGen1","",nBinsPhi,-TMath::Pi(),TMath::Pi()); hPhiGen1->Sumw2();
  TH1D *hTildaGen1 = new TH1D("hTildaGen1","",nBinsPhi,0,TMath::TwoPi()); hTildaGen1->Sumw2();
  TH1D *hPtGen1 = new TH1D("hPtGen1","",nBinsPt,mumuPtmin,mumuPtmax); hPtGen1->Sumw2();
  TH1D *hYGen1 = new TH1D("hYGen1","",nBinsY,yMin,yMax); hYGen1->Sumw2();
  TH1D *hThetaGen2 = new TH1D("hThetaGen2","",nBinsTheta,-1,1); hThetaGen2->Sumw2();
  TH1D *hPhiGen2 = new TH1D("hPhiGen2","",nBinsPhi,-TMath::Pi(),TMath::Pi()); hPhiGen2->Sumw2();
  TH1D *hTildaGen2 = new TH1D("hTildaGen2","",nBinsPhi,0,TMath::TwoPi()); hTildaGen2->Sumw2();
  TH1D *hPtGen2 = new TH1D("hPtGen2","",nBinsPt,mumuPtmin,mumuPtmax); hPtGen2->Sumw2();
  TH1D *hYGen2 = new TH1D("hYGen2","",nBinsY,yMin,yMax); hYGen2->Sumw2();
  TH2D *hSL = new TH2D("hSL","generated phi vs cos(theta) in Starlight",100,-1,1,100,-TMath::Pi(),TMath::Pi()); hSL->Sumw2();
  TH2D *hDecompThetaGen2 = new TH2D("hDecompThetaGen2","Polarisation terms vs cos(theta)",nBinsTheta,-1,1,7,0.5,7.5); hDecompThetaGen2->Sumw2();
  TH2D *hDecompPhiGen2 = new TH2D("hDecompPhiGen2","Polarisation terms vs phi",nBinsPhi,-TMath::Pi(),TMath::Pi(),7,0.5,7.5); hDecompPhiGen2->Sumw2();
  TH2D *hDecompTildaGen2 = new TH2D("hDecompTildaGen2","Polarisation terms vs phi^~",nBinsPhi,0,TMath::TwoPi(),7,0.5,7.5); hDecompTildaGen2->Sumw2();
	hPtSingleMu = new TH1D("hPtSingleMu","hPtSingleMu",100,0,10); hPtSingleMu->Sumw2();

  Int_t nmumugen = ReadGen(fmcpath,mumuPtmin,mumuPtmax,n1,n2,
			   hThetaGen1,hPhiGen1,hTildaGen1,hPtGen1,hYGen1,
			   hThetaGen2,hPhiGen2,hTildaGen2,hPtGen2,hYGen2,
			   iteration,fPt,fY,
			   hSL,
			   hDecompThetaGen2,hDecompPhiGen2,hDecompTildaGen2,
			   fUndo,
			   weight2015MC,weightUndo);

  TH2D *hThetaData = new TH2D("hThetaData","",nBinsTheta,-1,1,nBinsM,mMin,mMax); hThetaData->Sumw2();
  TH2D *hPhiData = new TH2D("hPhiData","",nBinsPhi,-TMath::Pi(),TMath::Pi(),nBinsM,mMin,mMax); hPhiData->Sumw2();
  TH2D *hTildaData = new TH2D("hTildaData","",nBinsPhi,0,TMath::TwoPi(),nBinsM,mMin,mMax); hTildaData->Sumw2();
  TH2D *hPtData = new TH2D("hPtData","",nBinsPt,mumuPtmin,mumuPtmax,nBinsM,mMin,mMax); hPtData->Sumw2();
  TH2D *hYData = new TH2D("hYData","",nBinsY,yMin,yMax,nBinsM,mMin,mMax); hYData->Sumw2();
  TH1D *hMData = new TH1D("hMData","",nBinsM,mMin,mMax); hMData->Sumw2();
  TH2D *hThetaMC = new TH2D("hThetaMC","",nBinsTheta,-1,1,nBinsM,mMin,mMax); hThetaMC->Sumw2();
  TH2D *hPhiMC = new TH2D("hPhiMC","",nBinsPhi,-TMath::Pi(),TMath::Pi(),nBinsM,mMin,mMax); hPhiMC->Sumw2();
  TH2D *hTildaMC = new TH2D("hTildaMC","",nBinsPhi,0,TMath::TwoPi(),nBinsM,mMin,mMax); hTildaMC->Sumw2();
  TH2D *hPtMC = new TH2D("hPtMC","",nBinsPt,mumuPtmin,mumuPtmax,nBinsM,mMin,mMax); hPtMC->Sumw2();
  TH2D *hYMC = new TH2D("hYMC","",nBinsY,yMin,yMax,nBinsM,mMin,mMax); hYMC->Sumw2();
  TH1D *hMMC = new TH1D("hMMC","",nBinsM,mMin,mMax); hMMC->Sumw2();
  TH2D *hThetaMC2 = new TH2D("hThetaMC2","",nBinsTheta,-1,1,nBinsM,mMin,mMax); hThetaMC2->Sumw2();
  TH2D *hPhiMC2 = new TH2D("hPhiMC2","",nBinsPhi,-TMath::Pi(),TMath::Pi(),nBinsM,mMin,mMax); hPhiMC2->Sumw2();
  TH2D *hTildaMC2 = new TH2D("hTildaMC2","",nBinsPhi,0,TMath::TwoPi(),nBinsM,mMin,mMax); hTildaMC2->Sumw2();
  TH2D *hDecompThetaMC2 = new TH2D("hDecompThetaMC2","Polarisation terms vs cos(theta)",nBinsTheta,-1,1,7,0.5,7.5); hDecompThetaMC2->Sumw2();
  TH2D *hDecompPhiMC2 = new TH2D("hDecompPhiMC2","Polarisation terms vs phi",nBinsPhi,-TMath::Pi(),TMath::Pi(),7,0.5,7.5); hDecompPhiMC2->Sumw2();
  TH2D *hDecompTildaMC2 = new TH2D("hDecompTildaMC2","Polarisation terms vs phi^~",nBinsPhi,0,TMath::TwoPi(),7,0.5,7.5); hDecompTildaMC2->Sumw2();

  TH1D *hThetaMCtrue1 = new TH1D("hThetaMCtrue1","",nBinsTheta,-1,1); hThetaMCtrue1->Sumw2();
  TH1D *hPhiMCtrue1 = new TH1D("hPhiMCtrue1","",nBinsPhi,-TMath::Pi(),TMath::Pi()); hPhiMCtrue1->Sumw2();
  TH1D *hTildaMCtrue1 = new TH1D("hTildaMCtrue1","",nBinsPhi,0,TMath::TwoPi()); hTildaMCtrue1->Sumw2();
  TH1D *hPtMCtrue1 = new TH1D("hPtMCtrue1","",nBinsPt,mumuPtmin,mumuPtmax); hPtMCtrue1->Sumw2();
  TH1D *hThetaMCtrue2 = new TH1D("hThetaMCtrue2","",nBinsTheta,-1,1); hThetaMCtrue2->Sumw2();
  TH1D *hPhiMCtrue2 = new TH1D("hPhiMCtrue2","",nBinsPhi,-TMath::Pi(),TMath::Pi()); hPhiMCtrue2->Sumw2();
  TH1D *hTildaMCtrue2 = new TH1D("hTildaMCtrue2","",nBinsPhi,0,TMath::TwoPi()); hTildaMCtrue2->Sumw2();
  TH1D *hPtMCtrue2 = new TH1D("hPtMCtrue2","",nBinsPt,mumuPtmin,mumuPtmax); hPtMCtrue2->Sumw2();
  TH2D *hThetaRes = new TH2D("hThetaRes","True vs reconstructed cos(Theta)",nBinsTheta,-1,1,nBinsTheta,-1,1); hThetaRes->Sumw2();
  TH2D *hPhiRes = new TH2D("hPhiRes","True vs reconstructed phi",nBinsPhi,-TMath::Pi(),TMath::Pi(),nBinsPhi,-TMath::Pi(),TMath::Pi()); hPhiRes->Sumw2();
  TH2D *hTildaRes = new TH2D("hTildaRes","True vs reconstructed phi^tilda",nBinsPhi,0,TMath::TwoPi(),nBinsPhi,0,TMath::TwoPi()); hTildaRes->Sumw2();
  TH2D *hPtRes = new TH2D("hPtRes","True vs reconstructed pT",nBinsPt,mumuPtmin,mumuPtmax,nBinsPt,mumuPtmin,mumuPtmax); hPtRes->Sumw2();
  TH2D *hYRes = new TH2D("hYRes","True vs reconstructed Y",nBinsY,yMin,yMax,nBinsY,yMin,yMax); hYRes->Sumw2();

  // setup response matrices
  RooUnfoldResponse *responseTheta = new RooUnfoldResponse;
  responseTheta->Setup(hThetaMCtrue1,hThetaMCtrue2);
  RooUnfoldResponse *responsePhi = new RooUnfoldResponse;
  responsePhi->Setup(hPhiMCtrue1,hPhiMCtrue2);
  RooUnfoldResponse *responseTilda = new RooUnfoldResponse;
  responseTilda->Setup(hTildaMCtrue1,hTildaMCtrue2);

  Int_t nmumudata = ReadDataOrMC(fdatapath,kFALSE,mumuPtmin,mumuPtmax,muPtmin,0,-1,
				 hThetaData,hPhiData,hTildaData,hPtData,hYData,hMData,
				 thetaMax);

  Int_t nmumumc = ReadDataOrMC(fmcpath,kTRUE,mumuPtmin,mumuPtmax,muPtmin,n1,n2,
			       hThetaMC,hPhiMC,hTildaMC,hPtMC,hYMC,hMMC,
			       thetaMax,
			       hThetaMC2,hPhiMC2,hTildaMC2,
			       hThetaMCtrue1,hPhiMCtrue1,hTildaMCtrue1,hPtMCtrue1,
			       hThetaMCtrue2,hPhiMCtrue2,hTildaMCtrue2,hPtMCtrue2,
			       hThetaRes,hPhiRes,hTildaRes,hPtRes,hYRes,
			       responseTheta,responsePhi,responseTilda,
			       iteration,fPt,fY,
			       hDecompThetaMC2,hDecompPhiMC2,hDecompTildaMC2,
			       fUndo,
			       weight2015MC,
					 	reweightTrig,weightUndo);


  // save all histos in a list with a unique name that contains all the options and values used
  TFile foutput("results.root","UPDATE");
  TList *loutput = new TList;
  loutput->SetName(Form("%s_%s_%.2f_%.2f_%.2f_%d_%d_%d_%d_%d_%s",
			fdatapath,fmcpath,
			mumuPtmin,mumuPtmax,muPtmin,
			extractOpt,iteration,n1,n2,reweightTrig,typeMCfit));


  // First fit integrated inv mass distributions in data and MC, for xcheck of fitting functions
  TF1 *fJpsiMC = FitJpsiMC(hMMC);
  TH1D *hFitResults = new TH1D(Form("%s_res",hMData->GetName()),"",10,0,10);
  Double_t errData;
  Double_t intData = FitJpsi(hMData,extractOpt,errData,fJpsiMC,hFitResults,hMData->GetName());
  loutput->Add(hMMC);
  loutput->Add(hMData);
  loutput->Add(hFitResults);
  //return;
  // fit Jpsi
  TH1D *hThetaDataJpsi = ExtractJpsi(hThetaData,hThetaMC,extractOpt,kFALSE,loutput);
  //return;
  TH1D *hPhiDataJpsi = ExtractJpsi(hPhiData,hPhiMC,extractOpt,kFALSE,loutput);
  TH1D *hTildaDataJpsi = ExtractJpsi(hTildaData,hTildaMC,extractOpt,kFALSE,loutput);
	//return;
  TH1D *hPtDataJpsi = ExtractJpsi(hPtData,hPtMC,extractOpt,kFALSE,loutput);
  TH1D *hYDataJpsi = ExtractJpsi(hYData,hYMC,extractOpt,kFALSE,loutput);
  //return;
  // also for MC, just by projecting to X axis
  TH1D *hThetaMCJpsi = ExtractJpsi(hThetaMC,0x0,-1,kTRUE,loutput);
  TH1D *hPhiMCJpsi = ExtractJpsi(hPhiMC,0x0,-1,kTRUE,loutput);
  TH1D *hTildaMCJpsi = ExtractJpsi(hTildaMC,0x0,-1,kTRUE,loutput);
  TH1D *hPtMCJpsi = ExtractJpsi(hPtMC,0x0,-1,kTRUE,loutput);
  TH1D *hYMCJpsi = ExtractJpsi(hYMC,0x0,-1,kTRUE,loutput);
  TH1D *hThetaMCJpsi2 = ExtractJpsi(hThetaMC2,0x0,-1,kTRUE,loutput);
  TH1D *hPhiMCJpsi2 = ExtractJpsi(hPhiMC2,0x0,-1,kTRUE,loutput);
  TH1D *hTildaMCJpsi2 = ExtractJpsi(hTildaMC2,0x0,-1,kTRUE,loutput);

  //  hMData->

  // Calculate acceptance
  // acc/gen
  TH1D *hAccTheta = (TH1D*)hThetaMCtrue2->Clone("hAccTheta");
  hAccTheta->Reset();
  hAccTheta->Divide(hThetaMCtrue2,hThetaGen2,1,1,"b");

  TH1D *hAccPhi = (TH1D*)hPhiMCtrue2->Clone("hAccPhi");
  hAccPhi->Reset();
  hAccPhi->Divide(hPhiMCtrue2,hPhiGen2,1,1,"b");

  TH1D *hAccTilda = (TH1D*)hTildaMCtrue2->Clone("hAccTilda");
  hAccTilda->Reset();
  hAccTilda->Divide(hTildaMCtrue2,hTildaGen2,1,1,"b");

  // rec/gen
  TH1D *hAccTheta2 = (TH1D*)hThetaMCJpsi2->Clone("hAccTheta2");
  hAccTheta2->Reset();
  hAccTheta2->Divide(hThetaMCJpsi2,hThetaGen2,1,1,"b");

  TH1D *hAccPhi2 = (TH1D*)hPhiMCJpsi2->Clone("hAccPhi2");
  hAccPhi2->Reset();
  hAccPhi2->Divide(hPhiMCJpsi2,hPhiGen2,1,1,"b");

  TH1D *hAccTilda2 = (TH1D*)hTildaMCJpsi2->Clone("hAccTilda2");
  hAccTilda2->Reset();
  hAccTilda2->Divide(hTildaMCJpsi2,hTildaGen2,1,1,"b");

  // Get the fully corrected cos(theta) and phi distributions
  // Unfolding + acceptance corrected
  TH1D *hThetaUnf = GetCorrected(hThetaDataJpsi,responseTheta,hAccTheta,loutput);
  TH1D *hThetaWoUnf = GetCorrectedWoUnf(hThetaDataJpsi,hAccTheta2,loutput);
  TH1D *hPhiUnf = GetCorrected(hPhiDataJpsi,responsePhi,hAccPhi,loutput);
  TH1D *hPhiWoUnf = GetCorrectedWoUnf(hPhiDataJpsi,hAccPhi2,loutput);
  TH1D *hTildaUnf = GetCorrected(hTildaDataJpsi,responseTilda,hAccTilda,loutput);
  TH1D *hTildaWoUnf = GetCorrectedWoUnf(hTildaDataJpsi,hAccTilda2,loutput);

  TH1D *hThetaUnfMC = GetCorrected(hThetaMCJpsi,responseTheta,hAccTheta,loutput);
  TH1D *hThetaWoUnfMC = GetCorrectedWoUnf(hThetaMCJpsi,hAccTheta2,loutput);
  TH1D *hPhiUnfMC = GetCorrected(hPhiMCJpsi,responsePhi,hAccPhi,loutput);
  TH1D *hPhiWoUnfMC = GetCorrectedWoUnf(hPhiMCJpsi,hAccPhi2,loutput);
  TH1D *hTildaUnfMC = GetCorrected(hTildaMCJpsi,responseTilda,hAccTilda,loutput);
  TH1D *hTildaWoUnfMC = GetCorrectedWoUnf(hTildaMCJpsi,hAccTilda2,loutput);

	// save all histos in a list with a unique name that contains all the options and values used

  // // fit fully corrected distributions to extract polarization params
  FitPolarizationParams(hThetaUnf,hPhiUnf,hTildaUnf,thetaMax);
  FitPolarizationParams(hThetaWoUnf,hPhiWoUnf,hTildaWoUnf,thetaMax);
  FitPolarizationParams(hThetaUnfMC,hPhiUnfMC,hTildaUnfMC,thetaMax);
  FitPolarizationParams(hThetaWoUnfMC,hPhiWoUnfMC,hTildaWoUnfMC,thetaMax);
  FitPolarizationParams(hThetaGen1,hPhiGen1,hTildaGen1,thetaMax);
	//
  // FractionFitting("Gen",
	// 	  hThetaGen1,hPhiGen1,hTildaGen1,
	// 	  hDecompThetaGen2,hDecompPhiGen2,hDecompTildaGen2,
	// 	  loutput,
	// 	  thetaMax);
  // FractionFitting("MC",
	// 	  hThetaMCJpsi,hPhiMCJpsi,hTildaMCJpsi,
	// 	  hDecompThetaMC2,hDecompPhiMC2,hDecompTildaMC2,
	// 	  loutput,
	// 	  thetaMax);
  // FractionFitting("Data",
	// 	  hThetaDataJpsi,hPhiDataJpsi,hTildaDataJpsi,
	// 	  hDecompThetaMC2,hDecompPhiMC2,hDecompTildaMC2,
	// 	  loutput,
	// 	  thetaMax);
	//
  // // Fit with chi^2 using MC templates
  MyFractionFitting("MyGen",
  		    hThetaGen1,hPhiGen1,hTildaGen1,
  		    hDecompThetaGen2,hDecompPhiGen2,hDecompTildaGen2,
  		    loutput,
  		    thetaMax,0);
  MyFractionFitting("MyMC",
  		    hThetaMCJpsi,hPhiMCJpsi,hTildaMCJpsi,
  		    hDecompThetaMC2,hDecompPhiMC2,hDecompTildaMC2,
  		    loutput,
  		    thetaMax,1);
  MyFractionFitting("MyData",
  		    hThetaDataJpsi,hPhiDataJpsi,hTildaDataJpsi,
  		    hDecompThetaMC2,hDecompPhiMC2,hDecompTildaMC2,
  		    loutput,
  		    thetaMax,2);


  loutput->Add(hThetaGen1);
  loutput->Add(hPhiGen1);
  loutput->Add(hTildaGen1);
  loutput->Add(hPtGen1);
  loutput->Add(hYGen1);
  loutput->Add(hThetaGen2);
  loutput->Add(hPhiGen2);
  loutput->Add(hTildaGen2);
  loutput->Add(hPtGen2);
  loutput->Add(hYGen2);
  loutput->Add(hThetaData);
  loutput->Add(hPhiData);
  loutput->Add(hTildaData);
  loutput->Add(hPtData);
  loutput->Add(hYData);
  loutput->Add(hThetaMC);
  loutput->Add(hPhiMC);
  loutput->Add(hTildaMC);
  loutput->Add(hPtMC);
  loutput->Add(hYMC);
  loutput->Add(hThetaMC2);
  loutput->Add(hPhiMC2);
  loutput->Add(hTildaMC2);
  loutput->Add(hThetaMCtrue1);
  loutput->Add(hPhiMCtrue1);
  loutput->Add(hTildaMCtrue1);
  loutput->Add(hPtMCtrue1);
  loutput->Add(hThetaMCtrue2);
  loutput->Add(hPhiMCtrue2);
  loutput->Add(hTildaMCtrue2);
  loutput->Add(hPtMCtrue2);
  loutput->Add(hThetaRes);
  loutput->Add(hPhiRes);
  loutput->Add(hTildaRes);
  loutput->Add(hPtRes);
  loutput->Add(hYRes);
  loutput->Add(hThetaDataJpsi);
  loutput->Add(hPhiDataJpsi);
  loutput->Add(hTildaDataJpsi);
  loutput->Add(hPtDataJpsi);
  loutput->Add(hYDataJpsi);
  loutput->Add(hThetaMCJpsi);
  loutput->Add(hPhiMCJpsi);
  loutput->Add(hTildaMCJpsi);
  loutput->Add(hPtMCJpsi);
  loutput->Add(hYMCJpsi);
  loutput->Add(hThetaMCJpsi2);
  loutput->Add(hPhiMCJpsi2);
  loutput->Add(hTildaMCJpsi2);
  loutput->Add(hAccTheta);
  loutput->Add(hAccPhi);
  loutput->Add(hAccTilda);
  loutput->Add(hAccTheta2);
  loutput->Add(hAccPhi2);
  loutput->Add(hAccTilda2);
  loutput->Add(hSL);
  loutput->Add(hDecompThetaGen2);
  loutput->Add(hDecompPhiGen2);
  loutput->Add(hDecompTildaGen2);
  loutput->Add(hDecompThetaMC2);
  loutput->Add(hDecompPhiMC2);
  loutput->Add(hDecompTildaMC2);
  loutput->Add(hPtSingleMu);
  loutput->Write(loutput->GetName(),TObject::kSingleKey | TObject::kOverwrite);
  foutput.Close();

}

Int_t ReadDataOrMC(const char *fpath, Bool_t isMC,
		   Double_t mumuPtmin, Double_t mumuPtmax, Double_t muPtmin,
		   Int_t n1, Int_t n2,
		   TH2D *hThetaRec, TH2D *hPhiRec, TH2D *hTildaRec, TH2D *hPtRec, TH2D *hYRec, TH1D *hMRec,
		   Double_t thetaMax,
		   TH2D *hThetaRec2, TH2D *hPhiRec2, TH2D *hTildaRec2,
		   TH1D *hThetaTrue1, TH1D *hPhiTrue1, TH1D *hTildaTrue1, TH1D *hPtTrue1,
		   TH1D *hThetaTrue2, TH1D *hPhiTrue2, TH1D *hTildaTrue2, TH1D *hPtTrue2,
		   TH2D *hThetaRes, TH2D *hPhiRes, TH2D *hTildaRes, TH2D *hPtRes, TH2D *hYRes,
		   RooUnfoldResponse *responseTheta, RooUnfoldResponse *responsePhi, RooUnfoldResponse *responseTilda,
		   Int_t iteration, TF1 **fPt, TF1 **fY,
		   TH2D *hDecompThetaRec2, TH2D *hDecompPhiRec2, TH2D *hDecompTildaRec2,
		   TF1 *fUndo,
		   Bool_t weight2015MC,
		 	Bool_t reweightTrig,
			Bool_t weightUndo)
{
  if (reweightTrig){
    TFile *fTrigger = TFile::Open("ratio_trigger_function.root");
    fTrig = (TF1 *)fTrigger->Get("fRatioTFR");
    printf("Ratio TFR \n");
  }
  Double_t massmu = TDatabasePDG::Instance()->GetParticle(13)->Mass();

  TFile *file = new TFile(Form("%s/AnalysisResults.root",fpath));
  TTree *treeRec = (TTree*)file->Get("MyTask/fTree");
  auto nentriesRec = treeRec->GetEntries();
  Int_t runNum;
  Double_t fMuPt1rec,fMuPt2rec,fMuEta1rec,fMuEta2rec,fMuPhi1rec,fMuPhi2rec;
  Int_t fMuCharge1rec,fMuCharge2rec,fMuChargerec;
  treeRec->SetBranchAddress("fRunNum",&runNum);
  treeRec->SetBranchAddress("fMuPt1", &fMuPt1rec);
  treeRec->SetBranchAddress("fMuPt2", &fMuPt2rec);
  treeRec->SetBranchAddress("fMuEta1", &fMuEta1rec);
  treeRec->SetBranchAddress("fMuEta2", &fMuEta2rec);
  treeRec->SetBranchAddress("fMuPhi1", &fMuPhi1rec);
  treeRec->SetBranchAddress("fMuPhi2", &fMuPhi2rec);
  treeRec->SetBranchAddress("fMuCharge1", &fMuCharge1rec);
  treeRec->SetBranchAddress("fMuCharge2", &fMuCharge2rec);
  Double_t fMuPt1true=0,fMuPt2true=0,fMuEta1true=0,fMuEta2true=0,fMuPhi1true=0,fMuPhi2true=0;
  Int_t fMuCharge1true=0,fMuCharge2true=0,fMuChargetrue=0;
  if (isMC) {
    treeRec->SetBranchAddress("fMCMuPt1", &fMuPt1true);
    treeRec->SetBranchAddress("fMCMuPt2", &fMuPt2true);
    treeRec->SetBranchAddress("fMCMuEta1", &fMuEta1true);
    treeRec->SetBranchAddress("fMCMuEta2", &fMuEta2true);
    treeRec->SetBranchAddress("fMCMuPhi1", &fMuPhi1true);
    treeRec->SetBranchAddress("fMCMuPhi2", &fMuPhi2true);
    treeRec->SetBranchAddress("fMCMuCharge1", &fMuCharge1true);
    treeRec->SetBranchAddress("fMCMuCharge2", &fMuCharge2true);
  }
  Int_t nmumu = 0, nrec = 0, nacc = 0, nres = 0;
  Double_t wmumu = 0, wrec = 0, wacc = 0, wres = 0;
  Double_t wmumu2 = 0, wrec2 = 0, wacc2 = 0, wres2 = 0;
  TLorentzVector fMuVec1rec,fMuVec2rec;
  TLorentzVector fMuMuRec;
  TLorentzVector fMuVec1true,fMuVec2true;
  TLorentzVector fMuMuTrue;
  for(Int_t ien = 0; ien < nentriesRec; ++ien) {
    // Initialize all vectors (to be safe they can not go uninitialized in the following steps)
    fMuVec1rec.SetPtEtaPhiM(0,0,0,0);
    fMuVec2rec.SetPtEtaPhiM(0,0,0,0);
    fMuMuRec.SetPtEtaPhiM(0,0,0,0);
    fMuVec1true.SetPtEtaPhiM(0,0,0,0);
    fMuVec2true.SetPtEtaPhiM(0,0,0,0);
    fMuMuTrue.SetPtEtaPhiM(0,0,0,0);

    treeRec->GetEntry(ien);
    if (runNum < 240000) {
      printf("Bad run number: %d\n",runNum);
      exit(0);
    }
    Double_t w2015 = 1.0;
    if (runNum < 270000) { // 2015 run period
      if (weight2015MC) w2015 = 2.0;
    }
    if (fMuCharge1rec == 0 || fMuCharge2rec == 0) continue; // invalid entries from the filter task

    if (fMuPt1rec < muPtmin || fMuPt2rec < muPtmin) continue; // tracks below Pt threshold
    if (fMuEta1rec < -4.0 || fMuEta1rec > -2.5) continue; // not applied in the filter task
    if (fMuEta2rec < -4.0 || fMuEta2rec > -2.5) continue; // not applied in the filter task
    if (fMuCharge1rec>0) {
      fMuVec1rec.SetPtEtaPhiM(fMuPt1rec,fMuEta1rec,fMuPhi1rec,massmu);
      fMuVec2rec.SetPtEtaPhiM(fMuPt2rec,fMuEta2rec,fMuPhi2rec,massmu);
    }
    else {
      fMuVec2rec.SetPtEtaPhiM(fMuPt1rec,fMuEta1rec,fMuPhi1rec,massmu);
      fMuVec1rec.SetPtEtaPhiM(fMuPt2rec,fMuEta2rec,fMuPhi2rec,massmu);
    }
    fMuMuRec = fMuVec1rec +fMuVec2rec;
    fMuChargerec = fMuCharge1rec +fMuCharge2rec;
    if (fMuChargerec != 0)			continue; //only mu+mu-
    Double_t mumuPtrec = fMuMuRec.Pt();
    if (mumuPtrec < mumuPtmin || mumuPtrec > mumuPtmax) continue; // dimuon pT cuts

    Double_t mumuMrec = fMuMuRec.M();
    Double_t mumuYrec = fMuMuRec.Rapidity();
    if (mumuYrec < -4.0 || mumuYrec > -2.5) continue; // not sure why this is not applied in filter task already?


    Double_t thetaHErec,phiHErec,tildaHErec;
    thetaHErec = CalcCosThetaPhiHE(fMuVec1rec,fMuVec2rec,phiHErec,tildaHErec);
    if (thetaHErec < -thetaMax || thetaHErec > thetaMax) continue; // cut on cos(theta), to be studied further ...

    Double_t wTrigger1 = 1.0, wTrigger2 = 1.0;
    Double_t wTrigger = wTrigger1*wTrigger2;
    if (isMC) {
      if (reweightTrig){
	wTrigger1 = fTrig->Eval(fMuPt1rec);
	wTrigger2 = fTrig->Eval(fMuPt2rec); //weight from Chun-L
	//      if(fMuPt1rec<0.5 || fMuPt2rec<0.5) cout<<wTrigger1<<"  "<<wTrigger2<<endl;

      }


      hPtSingleMu->Fill(fMuPt1rec,wTrigger1);
      hPtSingleMu->Fill(fMuPt2rec,wTrigger2);
      //if(wTrigger1>2.||wTrigger2>2.)cout<<wTrigger1<<"  "<<wTrigger2<<endl;
    }
    // Get weights in case of MC
    Double_t w = 1.0;
    if (isMC) {
      if (fMuCharge1true>0) {
	fMuVec1true.SetPtEtaPhiM(fMuPt1true,fMuEta1true,fMuPhi1true,massmu);
	fMuVec2true.SetPtEtaPhiM(fMuPt2true,fMuEta2true,fMuPhi2true,massmu);
      }
      else {
	fMuVec2true.SetPtEtaPhiM(fMuPt1true,fMuEta1true,fMuPhi1true,massmu);
	fMuVec1true.SetPtEtaPhiM(fMuPt2true,fMuEta2true,fMuPhi2true,massmu);
      }
      fMuMuTrue = fMuVec1true +fMuVec2true;
      // use MC weight
      Double_t mumuPttrue = fMuMuTrue.Pt();
      Double_t mumuYtrue = fMuMuTrue.Rapidity();
      w = GetWeight(iteration,mumuPttrue,mumuYtrue,fPt,fY);
      w *= w2015;
      w *= wTrigger;
    }

    // Fill data
    Bool_t accept = kFALSE;
    if (!isMC || (n2 < 0 || (ien%n2 < n1))) { // in case of MC, take n1 out of n2 entries as reco data and rest for response martix
      accept = kTRUE;
      hThetaRec->Fill(thetaHErec,mumuMrec,w);
      hPhiRec->Fill(phiHErec,mumuMrec,w);
      hTildaRec->Fill(tildaHErec,mumuMrec,w);
      hPtRec->Fill(mumuPtrec,mumuMrec,w);
      hYRec->Fill(mumuYrec,mumuMrec,w);
      hMRec->Fill(mumuMrec,w);
      nrec++;
      wrec += w2015;
      wrec2 += w;
    }
    // Fill MC
    if (isMC) {
      Double_t mumuPttrue = fMuMuTrue.Pt();
      Double_t mumuYtrue = fMuMuTrue.Rapidity();
      Double_t thetaHEtrue,phiHEtrue,tildaHEtrue;
      thetaHEtrue = CalcCosThetaPhiHE(fMuVec1true,fMuVec2true,phiHEtrue,tildaHEtrue);
      if (accept || n2 < 0) {
	hThetaTrue1->Fill(thetaHEtrue,w);
	hPhiTrue1->Fill(phiHEtrue,w);
	hTildaTrue1->Fill(tildaHEtrue,w);
	hPtTrue1->Fill(mumuPttrue,w);
	nacc++;
	wacc += w2015;
	wacc2 += w;
      }
      if (!accept || n2 < 0) {
	hThetaRec2->Fill(thetaHErec,mumuMrec,w);
	hPhiRec2->Fill(phiHErec,mumuMrec,w);
	hTildaRec2->Fill(tildaHErec,mumuMrec,w);
	hThetaTrue2->Fill(thetaHEtrue,w);
	hPhiTrue2->Fill(phiHEtrue,w);
	hTildaTrue2->Fill(tildaHEtrue,w);
	hPtTrue2->Fill(mumuPttrue,w);
	hThetaRes->Fill(thetaHErec,thetaHEtrue,w);
	hPhiRes->Fill(phiHErec,phiHEtrue,w);
	hTildaRes->Fill(tildaHErec,tildaHEtrue,w);
	hPtRes->Fill(mumuPtrec,mumuPttrue,w);
	hYRes->Fill(mumuYrec,mumuYtrue,w);
	responseTheta->Fill(thetaHErec,thetaHEtrue,w);
	responsePhi->Fill(phiHErec,phiHEtrue,w);
	responseTilda->Fill(tildaHErec,tildaHEtrue,w);
	Double_t wUndo = fUndo->Eval(thetaHEtrue);
	if (weightUndo == kTRUE) {
	  Decomp(hDecompThetaRec2,hDecompPhiRec2,hDecompTildaRec2,
		 thetaHErec,phiHErec,tildaHErec,
		 thetaHEtrue,phiHEtrue,tildaHEtrue,
		 w,wUndo);
	}else{
	  Decomp(hDecompThetaRec2,hDecompPhiRec2,hDecompTildaRec2,
		 thetaHErec,phiHErec,tildaHErec,
		 thetaHEtrue,phiHEtrue,tildaHEtrue,
		 w,1);
	}
	nres++;
	wres += w2015;
	wres2 += w;
      }
    }
    nmumu++;
    wmumu += w2015;
    wmumu2 += w;
  }
  printf("ReadDataOrMC: total:%d reco:%d mc1:%d mc2:%d\n",nmumu,nrec,nacc,nres);
  printf("ReadDataOrMC(weighted2015): total:%.0f reco:%.0f mc1:%.0f mc2:%.0f\n",wmumu,wrec,wacc,wres);
  printf("ReadDataOrMC(weighted(pT,y)): total:%.0f reco:%.0f mc1:%.0f mc2:%.0f\n",wmumu2,wrec2,wacc2,wres2);
  return nmumu;
}

Int_t ReadGen(const char *fpath,Double_t mumuPtmin, Double_t mumuPtmax,
	      Int_t n1, Int_t n2,
	      TH1D *hThetaGen1, TH1D *hPhiGen1, TH1D *hTildaGen1, TH1D *hPtGen1, TH1D *hYGen1,
	      TH1D *hThetaGen2, TH1D *hPhiGen2, TH1D *hTildaGen2, TH1D *hPtGen2, TH1D *hYGen2,
	      Int_t iteration, TF1 **fPt, TF1 **fY,
	      TH2D *hSL,
	      TH2D *hDecompThetaGen2, TH2D *hDecompPhiGen2, TH2D *hDecompTildaGen2,
	      TF1 *fUndo,
	      Bool_t weight2015MC,	Bool_t weightUndo)
{
  Double_t massmu = TDatabasePDG::Instance()->GetParticle(13)->Mass();

  TFile *file = new TFile(Form("%s/AnalysisResults.root",fpath));
  TTree *treeGen = (TTree*)file->Get("MyTask/fTreeGen");
  auto nentriesGen = treeGen->GetEntries();
  Int_t runNum;
  Double_t fMuPt1gen,fMuPt2gen,fMuEta1gen,fMuEta2gen,fMuPhi1gen,fMuPhi2gen;
  Int_t fMuCharge1gen,fMuCharge2gen;
  treeGen->SetBranchAddress("fMCRunNumGen",&runNum);
  treeGen->SetBranchAddress("fMCMuPt1Gen", &fMuPt1gen);
  treeGen->SetBranchAddress("fMCMuPt2Gen", &fMuPt2gen);
  treeGen->SetBranchAddress("fMCMuEta1Gen", &fMuEta1gen);
  treeGen->SetBranchAddress("fMCMuEta2Gen", &fMuEta2gen);
  treeGen->SetBranchAddress("fMCMuPhi1Gen", &fMuPhi1gen);
  treeGen->SetBranchAddress("fMCMuPhi2Gen", &fMuPhi2gen);
  treeGen->SetBranchAddress("fMCMuCharge1Gen", &fMuCharge1gen);
  treeGen->SetBranchAddress("fMCMuCharge2Gen", &fMuCharge2gen);
  TLorentzVector fMuVec1gen,fMuVec2gen;
  TLorentzVector fMuMuGen;
  Int_t nmumu = 0;
  Double_t wmumu = 0, wmumu2 = 0;
  for(Int_t ien = 0; ien < nentriesGen; ++ien) {
    treeGen->GetEntry(ien);
    if (runNum < 240000) {
      printf("Bad run number: %d\n",runNum);
      exit(0);
    }
    Double_t w2015 = 1.0;
    if (runNum < 270000) { // 2015 run period
      if (weight2015MC) w2015 = 2.0;
    }
    if (fMuCharge1gen>0) {
      fMuVec1gen.SetPtEtaPhiM(fMuPt1gen,fMuEta1gen,fMuPhi1gen,massmu);
      fMuVec2gen.SetPtEtaPhiM(fMuPt2gen,fMuEta2gen,fMuPhi2gen,massmu);
    }
    else {
      fMuVec2gen.SetPtEtaPhiM(fMuPt1gen,fMuEta1gen,fMuPhi1gen,massmu);
      fMuVec1gen.SetPtEtaPhiM(fMuPt2gen,fMuEta2gen,fMuPhi2gen,massmu);
    }
    fMuMuGen = fMuVec1gen +fMuVec2gen;
    Double_t mumuPtgen = fMuMuGen.Pt();
    if (mumuPtgen < mumuPtmin || mumuPtgen > mumuPtmax) continue; // dimuon pT cuts
    Double_t thetaHEgen,phiHEgen,tildaHEgen;
    // cut on Y
    Double_t mumuYgen = fMuMuGen.Rapidity();
    //    if (mumuYgen < -4.0 || mumuYgen > -2.5) continue; // for the moment w/o this cut ...
    thetaHEgen = CalcCosThetaPhiHE(fMuVec1gen,fMuVec2gen,phiHEgen,tildaHEgen);
    hSL->Fill(thetaHEgen,phiHEgen);
    // use MC weight
    Double_t w = GetWeight(iteration,mumuPtgen,mumuYgen,fPt,fY);
    w *= w2015;
    if (n2 < 0 || (ien%n2 < n1)) { // take n1 out of n2 entries
      hThetaGen1->Fill(thetaHEgen,w);
      hPhiGen1->Fill(phiHEgen,w);
      hTildaGen1->Fill(tildaHEgen,w);
      hPtGen1->Fill(mumuPtgen,w);
      hYGen1->Fill(mumuYgen,w);
    }
    if (n2 < 0 || (ien%n2 >= n1)) {
      hThetaGen2->Fill(thetaHEgen,w);
      hPhiGen2->Fill(phiHEgen,w);
      hTildaGen2->Fill(tildaHEgen,w);
      hPtGen2->Fill(mumuPtgen,w);
      hYGen2->Fill(mumuYgen,w);
			Double_t wUndo = fUndo->Eval(thetaHEgen);
			if (weightUndo == kTRUE) {
				Decomp(hDecompThetaGen2,hDecompPhiGen2,hDecompTildaGen2,
				 thetaHEgen,phiHEgen,tildaHEgen,
				 thetaHEgen,phiHEgen,tildaHEgen,
				 w,wUndo);
			}else{
				Decomp(hDecompThetaGen2,hDecompPhiGen2,hDecompTildaGen2,
				 thetaHEgen,phiHEgen,tildaHEgen,
				 thetaHEgen,phiHEgen,tildaHEgen,
				 w,1);
			}
    }
    nmumu++;
    wmumu += w2015;
    wmumu2 += w;
  }
  printf("ReadGen: %d  weighted(2015):%.0f  weighted(pT,y):%.0f\n",nmumu,wmumu,wmumu2);
  return nmumu;
}

Double_t CalcCosThetaPhiHE(TLorentzVector &p1Mom, TLorentzVector &p2Mom, Double_t &phiHE, Double_t &tildaHE)
{
  //
  // Calculate phi, theta and phi^tilda in helicity coordinate frame
  //
  // first & second daughter 4-mom as input

  // J/Psi 4-momentum vector
  TLorentzVector motherMom=p1Mom+p2Mom;

  // boost 4-mom vectors to the mother rest frame
  TVector3 beta = (-1.0/motherMom.E())*motherMom.Vect();
  p1Mom.Boost(beta);
  p2Mom.Boost(beta);

  // x,y,z axes
  TVector3 zAxisHE = (motherMom.Vect()).Unit();
  TVector3 beam(0,0,1);
  //  TVector3 yAxisHE = (zAxisHE.Cross(beam)).Unit(); // invert Y axis to be compatible with convention used by Simone
  TVector3 yAxisHE = (beam.Cross(zAxisHE)).Unit();
  TVector3 xAxisHE = (yAxisHE.Cross(zAxisHE)).Unit();

  // calc theta and phi
  Double_t thetaHE = zAxisHE.Dot((p1Mom.Vect()).Unit());

  phiHE = TMath::ATan2((p1Mom.Vect()).Dot(yAxisHE), (p1Mom.Vect()).Dot(xAxisHE));

  if (thetaHE>0)
    tildaHE = phiHE - 1./4.*TMath::Pi();
  else
    tildaHE = phiHE - 3./4.*TMath::Pi();

  if (tildaHE < 0) tildaHE += TMath::TwoPi();

  return thetaHE;
}

TH1D *ExtractJpsi(TH2D *hInputData, TH2D *hInputMC, Int_t extractOpt, Bool_t isMC, TList *loutput)
{
  // hInput's are defined as 2D histo with the invariant mass on Y and the variable of interest (cos(theta),phi,pT,..) on X
  TH1D *hOut = hInputData->ProjectionX(Form("%s_Jpsi",hInputData->GetName()),1,hInputData->GetNbinsY());
  if (!isMC) {
    hOut->Reset();
    for(Int_t ibin = 1; ibin <= hOut->GetNbinsX(); ++ibin) {
      TH1D *hMass = hInputData->ProjectionY(Form("%s_M_%d",hInputData->GetName(),ibin),ibin,ibin);
      TH1D *hMassMC = hInputMC->ProjectionY(Form("%s_M_%d",hInputMC->GetName(),ibin),ibin,ibin);
      // create histogram to store some fit results
      TH1D *hFitResults = new TH1D(Form("%s_res",hMass->GetName()),"",10,0,10);
      hFitResults->SetBinContent(1,hInputData->GetXaxis()->GetBinLowEdge(ibin));
      hFitResults->SetBinContent(2,hInputData->GetXaxis()->GetBinUpEdge(ibin));
      loutput->Add(hMass);
      loutput->Add(hMassMC);
      loutput->Add(hFitResults);
      if (hMass->GetEntries() < 50) continue; // do not fit if too few entries
      Double_t integral,err;
      // fit jpsi
      TF1 *fJpsi = FitJpsiMC(hMassMC); // first fit MC to extract crystal-ball params
      integral = FitJpsi(hMass,extractOpt,err,fJpsi,hFitResults,hMass->GetName());
      hOut->SetBinContent(ibin,integral);
      hOut->SetBinError(ibin,err);

    }
  }
  return hOut;
}

TF1 *FitJpsiMC(TH1D *histoToFit)
{
  Double_t fitMin = 2.2,fitMax = 6.0;
  TF1 *fitDoubleCrystal = new TF1(Form("%s_fitJPsi",histoToFit->GetName()),DiCrystalBall,fitMin,fitMax,7);
  fitDoubleCrystal -> SetParameter(0,histoToFit->GetBinContent(histoToFit->FindBin(3.09))); //normalisation of J/psi
  fitDoubleCrystal -> SetParameter(1,3.09); //mean of J/psi
  fitDoubleCrystal -> SetParameter(2,0.1); //sigma
  fitDoubleCrystal -> SetParameter(3,1.1); //alphaLow of J/psi
  fitDoubleCrystal -> SetParameter(4,2.7); //npowL of J/psi
  fitDoubleCrystal -> SetParameter(5,2.0); //alphaHigh of J/psi
  fitDoubleCrystal -> SetParameter(6,2.7); //npowH of J/psi
  fitDoubleCrystal -> SetNpx(2000);
  fitDoubleCrystal ->SetParNames("Normalisation","Mean","#sigma","#alpha low","n low","#alpha high","n high");

  for (int j = 0; j < 7; ++j) fitDoubleCrystal->FixParameter(j,fitDoubleCrystal->GetParameter(j));
  fitDoubleCrystal->ReleaseParameter(0);
  histoToFit->Fit(fitDoubleCrystal,"QWLI","",fitMin,fitMax);
  fitDoubleCrystal->ReleaseParameter(1);
  fitDoubleCrystal -> SetParLimits(1,2.9,3.3);
  histoToFit->Fit(fitDoubleCrystal,"QWLI","",fitMin,fitMax);

  for(int i=2;i<7;i++) fitDoubleCrystal->FixParameter(i,fitDoubleCrystal->GetParameter(i));
  fitDoubleCrystal->ReleaseParameter(2);
  fitDoubleCrystal -> SetParLimits(2,0.03,0.2);
  histoToFit->Fit(fitDoubleCrystal,"QWLI","",fitMin,fitMax);

  for(int i=3;i<7;i++) fitDoubleCrystal->FixParameter(i,fitDoubleCrystal->GetParameter(i));
  fitDoubleCrystal->ReleaseParameter(3);
  fitDoubleCrystal -> SetParLimits(3,0.5,3.0);
  fitDoubleCrystal->ReleaseParameter(5);
  fitDoubleCrystal -> SetParLimits(5,1.0,4.0);
  histoToFit->Fit(fitDoubleCrystal,"QWLI","",fitMin,fitMax);

  fitDoubleCrystal->ReleaseParameter(4);
  fitDoubleCrystal -> SetParLimits(4,1.1,12.0);
  fitDoubleCrystal->ReleaseParameter(6);
  fitDoubleCrystal -> SetParLimits(6,1.1,10.0);
  histoToFit->Fit(fitDoubleCrystal,"QWLI","",fitMin,fitMax);

  return fitDoubleCrystal;
}

Double_t FitJpsi(TH1D *histoToFit, Int_t extractOpt, Double_t &errJpsi, TF1 *fJpsi, TH1D *hFitResults,const  char* nameV="")
{
  // Use extractOpt to vary fit functions, range, extraction method (for example side-bands)
  Double_t fitMin, fitMax;
  if (extractOpt<10) {
    fitMin = 2.2,fitMax = 6.0;
  }else if (extractOpt>=10 && extractOpt<20){
    fitMin = 2.0,fitMax = 5.0;
  }else{
    fitMin = 2.2, fitMax = 4.8;
  }
  TString varStr(nameV);
  Double_t lowM = histoToFit->GetXaxis()->GetXmin();
  Double_t highM = histoToFit->GetXaxis()->GetXmax();
  TF1 *FitTotal = 0;
  if(extractOpt==0 || extractOpt==10 || extractOpt==20){
    FitTotal=new TF1("TotalFit",TotalData,fitMin,fitMax,12);
    FitTotal->SetNpx(2000);

    FitTotal -> SetParameter(0,histoToFit->GetBinContent(histoToFit->FindBin(2.3))); //normalisation background (taken at M=2.3 GeV/c^2)
  //  FitTotal -> SetParLimits(1,2.0,3.0);
    FitTotal -> SetParameter(1,2.3); //mean of background VWG
    FitTotal -> SetParameter(2,0.7); //A parameter of background VWG
    FitTotal -> SetParameter(3,0.7); //B parameter of background VWG
    FitTotal -> SetParameter(4,histoToFit->GetBinContent(histoToFit->FindBin(3.09))); //normalisation of J/psi
    for(int i = 0;i<6;i++) FitTotal -> SetParameter((5+i),fJpsi->GetParameter(i+1)); //parameters fitJPsi taken from fit to MC
    FitTotal -> SetParameter(11,0.5*histoToFit->GetBinContent(histoToFit->FindBin(3.69))); //normalisation of Psi2s

  }else if(extractOpt == 1 || extractOpt == 11 || extractOpt==21){

    FitTotal=new TF1("TotalFit",TotalData1,fitMin,fitMax,12);
    FitTotal->SetNpx(2000);
    FitTotal -> SetParameter(0,histoToFit->GetBinContent(histoToFit->FindBin(2.3))); //normalisation background (taken at M=2.2 GeV/c^2)
    //    FitTotal -> SetParLimits(1,2.0,3.0);
    //  if(extractOpt == 1 || extractOpt == 11 ){
    // FitTotal -> SetParameter(1,0.005); //par pol1
    // FitTotal -> SetParameter(2,-0.005); //1st par pol2
    // FitTotal -> SetParameter(3,-.001); //2nd par pol2
    //	}else{

    // if (varStr.Contains("Tilda",TString::kIgnoreCase) || varStr.Contains("Phi",TString::kIgnoreCase) ) {
    // 	FitTotal -> SetParameter(1,-0.15); //par pol1
    // 	FitTotal -> SetParameter(2,10); //1st par pol2
    // 	FitTotal -> SetParameter(3,300); //2nd par pol2
    // }else{
    FitTotal -> SetParameter(1,0.08); //par pol1
    FitTotal -> SetParameter(2,0.7); //1st par pol2
    FitTotal -> SetParameter(3,0.1); //2nd par pol2
    //	}
    FitTotal -> SetParameter(4,histoToFit->GetBinContent(histoToFit->FindBin(3.09))); //normalisation of J/psi
    for(int i = 0;i<6;i++) FitTotal -> SetParameter((5+i),fJpsi->GetParameter(i+1)); //parameters fitJPsi taken from fit to MC
    FitTotal -> SetParameter(11,0.5*histoToFit->GetBinContent(histoToFit->FindBin(3.69))); //normalisation of Psi2s
  }else if (extractOpt ==2 || extractOpt==12 || extractOpt == 22) {
    FitTotal=new TF1("TotalFit",TotalData2,fitMin,fitMax,12);
    FitTotal->SetNpx(2000);
    FitTotal -> SetParameter(0,histoToFit->GetBinContent(histoToFit->FindBin(2.3))); //normalisation background (taken at M=2.3 GeV/c^2)
    //FitTotal -> SetParLimits(1,2.0,3.0);
    FitTotal -> SetParameter(1,0.7); //par exp
    FitTotal -> SetParameter(2,0.4); //1st par pol2
    FitTotal -> SetParameter(3,0.3); //2nd par pol2
    FitTotal -> SetParameter(4,histoToFit->GetBinContent(histoToFit->FindBin(3.09))); //normalisation of J/psi
    for(int i = 0;i<6;i++) FitTotal -> SetParameter((5+i),fJpsi->GetParameter(i+1)); //parameters fitJPsi taken from fit to MC
    FitTotal -> SetParameter(11,0.5*histoToFit->GetBinContent(histoToFit->FindBin(3.69))); //normalisation of Psi2s
  }
  FitTotal -> SetParName(0,"Normalisation BG");
  FitTotal -> SetParName(1,"Mean BG");
  FitTotal -> SetParName(2,"Parameter A BG");
  FitTotal -> SetParName(3,"Parameter B BG");
  FitTotal -> SetParName(4,"Normalisation J/#psi");
  FitTotal -> SetParName(5,"Mean J/#psi");
  FitTotal -> SetParName(6,"#sigma J/#psi");
  FitTotal -> SetParName(11,"Normalisation #psi2s");
  for (int j = 0; j < 12; ++j)  FitTotal->FixParameter(j,FitTotal->GetParameter(j));
  histoToFit->Fit(FitTotal,"QIL","",fitMin,fitMax);
  FitTotal->ReleaseParameter(0);
  FitTotal->ReleaseParameter(4);
  //  FitTotal->ReleaseParameter(11);
  //  FitTotal->SetParLimits(11,0.,10.*histoToFit->GetBinContent(histoToFit->FindBin(3.69)));
  histoToFit->Fit(FitTotal,"QIL","",fitMin,fitMax);
  //	FitTotal->SetParLimits(11,0,1e+2);
  //if (extractOpt == 1 || extractOpt == 11 || extractOpt==21) {
  // if(histoToFit->GetBinContent(histoToFit->FindBin(3.69))<8){
  //   cout<<endl;
  //   cout<<endl;
  //   cout<<"pic Psi2s bloque "<<nameV<<endl;
  //   cout<<endl;
  //   cout<<endl;
  // }else{
  FitTotal->ReleaseParameter(11);
  histoToFit->Fit(FitTotal,"QIL","",fitMin,fitMax);
  //    }
  //}else{
  //   FitTotal->ReleaseParameter(11);
    // cout<<endl;
    // cout<<endl;
    // cout<<"2 "<<nameV<<endl;
    // cout<<endl;
    // cout<<endl;
    // histoToFit->Fit(FitTotal,"QLI","",fitMin,fitMax);
//}
  FitTotal->ReleaseParameter(5);
   histoToFit->Fit(FitTotal,"QIL","",fitMin,fitMax);

   FitTotal->ReleaseParameter(6);
   histoToFit->Fit(FitTotal,"QIL","",fitMin,fitMax);
   //if(varStr.Contains("hMData",TString::kIgnoreCase)){
   Int_t cutOpt=0;
   if (extractOpt == 1 || extractOpt == 11 || extractOpt==21) cutOpt=12;
   if (extractOpt == 2 || extractOpt == 12 || extractOpt==22) cutOpt=13;
	 if (extractOpt == 0 || extractOpt == 10 || extractOpt==20) cutOpt=4;
   if(histoToFit->GetBinContent(histoToFit->FindBin(2.2))<cutOpt){
    cout<<"pic BG bloque "<<nameV<<endl;
   }else{
    FitTotal->ReleaseParameter(1);
    histoToFit->Fit(FitTotal,"QIL","",fitMin,fitMax);
    if (extractOpt==0 || extractOpt==10 || extractOpt==20) {
      FitTotal->ReleaseParameter(2);
      FitTotal->ReleaseParameter(3);
      histoToFit->Fit(FitTotal,"QIL","",fitMin,fitMax);
    }else{
      FitTotal->ReleaseParameter(2);
      histoToFit->Fit(FitTotal,"QIL","",fitMin,fitMax);

      FitTotal->ReleaseParameter(3);
      histoToFit->Fit(FitTotal,"QIL","",fitMin,fitMax);
    }
   }

  TFitResultPtr r = histoToFit->Fit(FitTotal,"QILS","",fitMin,fitMax);


  // construct functions for j/psi, psi' and background, and store them in the histogram
  TF1 *fitCrystalJPsi = new TF1("fitJPsiData",DiCrystalBall,fitMin,fitMax,7);
  TF1 *fitCrystalPsi2s = new TF1("fitPsi2sData",DiCrystalBall,fitMin,fitMax,7);
  TF1 *backGround = 0;
  if(extractOpt==0 || extractOpt==10 || extractOpt == 20) backGround =  new TF1("fitBackgroundData",BackGround,fitMin,fitMax,4);
  else if(extractOpt==1 || extractOpt==11|| extractOpt == 21 ) backGround = new TF1("fitBackgroundData",Pol1Pol2,fitMin,fitMax,4);
  else if(extractOpt==2 || extractOpt==12 || extractOpt == 22) backGround = new TF1("fitBackgroundData",Pol2Exp,fitMin,fitMax,4);
  fitCrystalJPsi->SetNpx(2000);
  fitCrystalPsi2s->SetNpx(2000);
  backGround->SetNpx(2000);
  FitTotal->SetLineColor(kBlack);
  backGround->SetLineColor(kBlue);
  fitCrystalJPsi->SetLineColor(kRed);
  fitCrystalPsi2s->SetLineColor(kGreen);

  const double *params = r->GetParams();
  backGround->SetParameters(&params[0]);
  fitCrystalJPsi->SetParameters(&params[4]);
  Double_t params2[7];
  params2[0] = params[11];
  params2[1] = params[5]+0.589; // TDatabasePDG::Instance()->GetParticle(100443)->Mass()-TDatabasePDG::Instance()->GetParticle(443)->Mass();
  params2[2] = 1.09*params[6]; // from AN of Simone Ragoni
  params2[3] = params[7];
  params2[4] = params[8];
  params2[5] = params[9];
  params2[6] = params[10];
  fitCrystalPsi2s->SetParameters(&params2[0]);

  histoToFit->GetListOfFunctions()->Add(backGround);
  histoToFit->GetListOfFunctions()->Add(fitCrystalJPsi);
  histoToFit->GetListOfFunctions()->Add(fitCrystalPsi2s);

  // calculate the Jpsi integral and error
  Double_t intJpsi = fitCrystalJPsi->Integral(lowM,highM)/histoToFit->GetBinWidth(1);
  TMatrixDSym cv = r->GetCovarianceMatrix();
  Double_t covmat[7][7];
  for(Int_t i = 0; i < 7; ++i) {
    for(Int_t j = 0; j < 7; ++j) {
      covmat[i][j] = cv(4+i,4+j);
      //      printf("%f ",covmat[i][j]);
    }
  }
  //  printf("\n");
  if (covmat[0][0] == 0) { // abnormal fit termination
    printf("Abnormal fit termination.\n");
    //errJpsi = 0.; //
    errJpsi=TMath::Sqrt(intJpsi);
  }
  else
    errJpsi = fitCrystalJPsi->IntegralError(lowM,highM,&params[4],&covmat[0][0])/histoToFit->GetBinWidth(1);
  // calculate background under Jpsi in +- 3 sigma
  Double_t intBkg = backGround->Integral(fitCrystalJPsi->GetParameter(1)-3.*fitCrystalJPsi->GetParameter(2),
					 fitCrystalJPsi->GetParameter(1)+3.*fitCrystalJPsi->GetParameter(2))
    /histoToFit->GetBinWidth(1);
  Double_t covmat2[4][4];
  for(Int_t i = 0; i < 4; ++i) {
    for(Int_t j = 0; j < 4; ++j) {
      covmat2[i][j] = cv(i,j);
    }
  }
  Double_t errBkg;
  if (covmat2[0][0] == 0) { // abnormal fit termination
    printf("Abnormal fit termination.\n");
    //    errBkg = 0.;
    errBkg=TMath::Sqrt(intBkg);
  }
  else
    errBkg = backGround->IntegralError(fitCrystalJPsi->GetParameter(1)-3.*fitCrystalJPsi->GetParameter(2),
				       fitCrystalJPsi->GetParameter(1)+3.*fitCrystalJPsi->GetParameter(2),
				       &params[0],&covmat2[0][0])/histoToFit->GetBinWidth(1);

  // store some fit results in form of histogram
  hFitResults->SetBinContent(3,intJpsi);
  hFitResults->SetBinError(3,errJpsi);
  hFitResults->SetBinContent(4,FitTotal->GetChisquare());
  hFitResults->SetBinContent(5,FitTotal->GetNDF());
  hFitResults->SetBinContent(6,intBkg);
  hFitResults->SetBinError(6,errBkg);

  return intJpsi;
}

Double_t DiCrystalBall(double* ax, double* ap){
  double m      = ax[0];
  double norm   = ap[0];  // norm Number of events in the peak.
  double m0     = ap[1];  // Nominal position
  double sigma  = ap[2];  // sigma
  double alphaL = ap[3];  // Low-side alpha.
  double nL     = ap[4];  // Low-side npow
  double alphaR = ap[5];  // High-side alpha.
  double nR     = ap[6];  // High-side npow

  Double_t t = (m-m0)/sigma;
  Double_t absAlphaL = fabs((Double_t)alphaL);
  Double_t absAlphaR = fabs((Double_t)alphaR);

  if (t >= -absAlphaL && t <= absAlphaR) {
    return norm*exp(-0.5*t*t);
  }
  else if (t < -absAlphaL) {
    Double_t a = TMath::Power(nL/absAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    Double_t b = nL/absAlphaL - absAlphaL;
    return norm*a/TMath::Power(b - t, nL);
  }
  else {
    Double_t a = TMath::Power(nR/absAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    Double_t b = nR/absAlphaR - absAlphaR;
    return norm*a/TMath::Power(b + t, nR);
  }
}

Double_t BackGround(Double_t *ax, Double_t *ap){
  double m = ax[0];
  double N = ap[0];
  double m0 = ap[1];
  double A = ap[2];
  double B = ap[3];

  Double_t t = (m-m0);
  Double_t sigma = A + B*t/m0;

  return N*TMath::Gaus(m,m0,sigma);
}

Double_t Pol1Pol2(Double_t *ax, Double_t *ap){
  double x = ax[0];
  return ap[0]*(1+ap[1]*x)/(1+ap[2]*x+ap[3]*x*x);

}
Double_t Pol2Exp(Double_t *ax, Double_t *ap){
  double x = ax[0];
  double xch = 4.;
  double result = 0.;
  if(x<=xch)    result= ap[0]*exp(ap[1]*x)*(1+ap[2]*(x-xch)+ap[3]*(x-xch)*(x-xch));
  else if(x>xch)     result = ap[0]*exp(ap[1]*x);
  return result;
}

Double_t TotalData(Double_t *x, Double_t *ap){
  // Combined background + signal
  Double_t ap2[7];
  ap2[0] = ap[11];
  ap2[1] = ap[5]+0.589; // TDatabasePDG::Instance()->GetParticle(100443)->Mass()-TDatabasePDG::Instance()->GetParticle(443)->Mass();
  ap2[2] = 1.09*ap[6]; // from AN of Simone Ragoni
  ap2[3] = ap[7];
  ap2[4] = ap[8];
  ap2[5] = ap[9];
  ap2[6] = ap[10];
  return BackGround(x,&ap[0]) + DiCrystalBall(x,&ap[4])+DiCrystalBall(x,&ap2[0]);
}

Double_t TotalData1(Double_t *x, Double_t *ap){
  // Combined background + signal
  Double_t ap2[7];
  ap2[0] = ap[11];
  ap2[1] = ap[5]+0.589; // TDatabasePDG::Instance()->GetParticle(100443)->Mass()-TDatabasePDG::Instance()->GetParticle(443)->Mass();
  ap2[2] = 1.09*ap[6]; // from AN of Simone Ragoni
  ap2[3] = ap[7];
  ap2[4] = ap[8];
  ap2[5] = ap[9];
  ap2[6] = ap[10];
  return Pol1Pol2(x,&ap[0]) + DiCrystalBall(x,&ap[4])+DiCrystalBall(x,&ap2[0]);
}

Double_t TotalData2(Double_t *x, Double_t *ap){
  // Combined background + signal
  Double_t ap2[7];
  ap2[0] = ap[11];
  ap2[1] = ap[5]+0.589; // TDatabasePDG::Instance()->GetParticle(100443)->Mass()-TDatabasePDG::Instance()->GetParticle(443)->Mass();
  ap2[2] = 1.09*ap[6]; // from AN of Simone Ragoni
  ap2[3] = ap[7];
  ap2[4] = ap[8];
  ap2[5] = ap[9];
  ap2[6] = ap[10];
  return Pol2Exp(x,&ap[0]) + DiCrystalBall(x,&ap[4])+DiCrystalBall(x,&ap2[0]);
}

Double_t triggerError(Double_t *ax, Double_t *ap)
{
  //From Chun-Lu thesis
  Double_t pt = ax[0];
  Double_t p0 = ap[0];
  Double_t p1 = ap[1];
  Double_t p2 = ap[2];
  Double_t p3 = ap[3];
  Double_t p4 = ap[4];
  Double_t p5 = ap[5];
  Double_t p6 = ap[6];
  Double_t p7 = ap[7];
  Double_t p8 = ap[8];
  Double_t p9 = ap[9];
  Double_t p10 = ap[10];
  Double_t p11 = ap[11];

  Double_t function = 0;
  Double_t RF =  p7 + p0*(TMath::Erf((TMath::Max(pt,p6)-p1)/TMath::Sqrt(2.)/p2)-1.);
  if (pt>p6) {
    function = RF;
  }else if (pt<p6) {
    function = (RF+p3*(TMath::Erf((-TMath::Min(pt,p6)-p4)/TMath::Sqrt(2.)/p5)-TMath::Erf((-p6-p4)/TMath::Sqrt(2.)/p5)));
  }else if (pt>2.0){
    function = p8+(p9/(1+TMath::Exp(-p10*(pt-p11))));
  }
  return function;
}
TH1D *GetCorrected(TH1D *hRaw, RooUnfoldResponse *response, TH1D *hAcc, TList *loutput)
{
  // unfolding
  RooUnfoldBayes unfold(response,hRaw,4,kFALSE,Form("%s_unf",hRaw->GetName()));
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
  TH1D* hUnf= (TH1D*)unfold.Hreco(errorTreatment);
  TH1D* hFold= (TH1D*)response->ApplyToTruth(hUnf,Form("%s_fold",hRaw->GetName()));
  loutput->Add(hUnf);
  loutput->Add(hFold);
  // acceptance correction
  TH1D *hCorr = (TH1D *)hUnf->Clone(Form("%s_corr",hRaw->GetName()));
  hCorr->Reset();
  hCorr->Divide(hUnf,hAcc);
  loutput->Add(hCorr);

  return hCorr;
}

TH1D *GetCorrectedWoUnf(TH1D *hRaw, TH1D *hAcc, TList *loutput)
{
  // only acceptance correction (no unfolding)
  TH1D *hCorrWoUnf = (TH1D *)hRaw->Clone(Form("%s_corrWoUnf",hRaw->GetName()));
  hCorrWoUnf->Reset();
  hCorrWoUnf->Divide(hRaw,hAcc);
  loutput->Add(hCorrWoUnf);

  return hCorrWoUnf;
}

Double_t GetWeight(Int_t iteration, Double_t pT, Double_t Y, TF1 **fPt, TF1 **fY)
{
  // get MC weights
  Double_t w = 1.0;
  for(Int_t ii = 0; ii < iteration; ++ii) {
    w *= fPt[ii]->Eval(pT);
    w *= fY[ii]->Eval(Y);
  }
  return w;
}

void FitPolarizationParams(TH1D *hTheta, TH1D *hPhi, TH1D *hTilda, Double_t thetaMax)
{
  // Fit with 1D functions, see AN of Simone Ragoni
  TF1 *fTheta = new TF1("fTheta","[0]*(1.0+[1]*x*x)",-thetaMax,thetaMax);
  fTheta->SetParameter(0,hTheta->GetBinContent(hTheta->FindBin(0.0)));
  fTheta->SetParameter(1,1.0);
  hTheta->Fit(fTheta,"R");
  TF1 *fPhi = new TF1("fPhi","[0]*(1.0+2.0*[1]*TMath::Cos(2.0*x))");
  fPhi->SetParameter(0,hPhi->GetBinContent(hPhi->FindBin(TMath::Pi()/4.)));
  fPhi->SetParameter(1,0.0);
  hPhi->Fit(fPhi);
  TF1 *fTilda = new TF1("fTilda","[0]*(1.0+TMath::Sqrt(2.0)*[1]*TMath::Cos(x))");
  fTilda->SetParameter(0,hTilda->GetBinContent(hTilda->FindBin(TMath::Pi()/2.)));
  fTilda->SetParameter(1,0.0);
  hTilda->Fit(fTilda);
  // here we have to divide by (3+Lambda_theta) and propage the unc
  // ...
}

TF1 *UndoPolarisation()
{
  // Build function to undo transverse polarisation used by default SL MC
  TF1 *fUndo = new TF1("fUndo","(TMath::Abs(x)<[0])?[1]*(1.+[2]*x*x+[3]*x*x*x*x):[1]*[4]*(1.+[5]*x*x)",-1,1);
  fUndo->SetParameter(0,7.74067e-01);
  fUndo->SetParameter(1,1.0);
  fUndo->SetParameter(2,3.50956e-01);
  fUndo->SetParameter(3,9.58266e-01);
  fUndo->SetParameter(4,9.92048e-01);
  fUndo->SetParameter(5,1.00696e+00);
  return fUndo;
}

void Decomp(TH2D *hTheta, TH2D *hPhi, TH2D *hTilda,
	    Double_t theta, Double_t phi, Double_t tilda,
	    Double_t thetaTrue, Double_t phiTrue, Double_t tildaTrue,
	    Double_t w, Double_t wUndo)
{
  hTheta->Fill(theta,1,w/wUndo*1.0); // const term
  hPhi->Fill(phi,1,w/wUndo*1.0); // const term
  hTilda->Fill(tilda,1,w/wUndo*1.0); // const term

  Double_t w2 = (1.0+thetaTrue*thetaTrue); // 1+cos(theta)^2 term
  hTheta->Fill(theta,2,w/wUndo*w2);
  hPhi->Fill(phi,2,w/wUndo*w2);
  hTilda->Fill(tilda,2,w/wUndo*w2);

  Double_t w3 = (1.0+(1.0-thetaTrue*thetaTrue)*TMath::Cos(2.0*phiTrue)); // 1+sin(theta)^2*cos(2phi)
  hTheta->Fill(theta,3,w/wUndo*w3);
  hPhi->Fill(phi,3,w/wUndo*w3);
  hTilda->Fill(tilda,3,w/wUndo*w3);

  Double_t w4 = (1.0+2.0*TMath::Sqrt(1.0-thetaTrue*thetaTrue)*thetaTrue*TMath::Cos(phiTrue)); // 1+sin(2theta)*cos(phi)
  hTheta->Fill(theta,4,w/wUndo*w4);
  hPhi->Fill(phi,4,w/wUndo*w4);
  hTilda->Fill(tilda,4,w/wUndo*w4);

  Double_t w5 = thetaTrue*thetaTrue; // cos(theta)^2 term
  hTheta->Fill(theta,5,w/wUndo*w5);
  hPhi->Fill(phi,5,w/wUndo*w5);
  hTilda->Fill(tilda,5,w/wUndo*w5);

  Double_t w6 = (1.0-thetaTrue*thetaTrue)*TMath::Cos(2.0*phiTrue); // sin(theta)^2*cos(2phi)
  hTheta->Fill(theta,6,w/wUndo*w6);
  hPhi->Fill(phi,6,w/wUndo*w6);
  hTilda->Fill(tilda,6,w/wUndo*w6);

  Double_t w7 = 2.0*TMath::Sqrt(1.0-thetaTrue*thetaTrue)*thetaTrue*TMath::Cos(phiTrue); // sin(2theta)*cos(phi)
  hTheta->Fill(theta,7,w/wUndo*w7);
  hPhi->Fill(phi,7,w/wUndo*w7);
  hTilda->Fill(tilda,7,w/wUndo*w7);
}

void FractionFitting(const char *hname,
		     TH1D *hThetaData, TH1D *hPhiData, TH1D *hTildaData,
		     TH2D *hThetaMC, TH2D *hPhiMC, TH2D *hTildaMC,
		     TList *loutput,
		     Double_t thetaMax)
{
  //  Double_t thetaMax = ;
  Int_t bin1 = hThetaData->FindBin(-thetaMax+1e-6);
  Int_t bin2 = hThetaData->FindBin(thetaMax-1e-6);
  Int_t nBins = bin2-bin1+1 + hPhiData->GetNbinsX() + hTildaData->GetNbinsX();
  TH1D *hData = new TH1D(Form("%s_FracData",hname),"",nBins,0.5,0.5+Double_t(nBins));
  TH1D *hMC[4];
  for(Int_t j = 0; j < 4; ++j)
    hMC[j] = new TH1D(Form("%s_FracMC%d",hname,j),"",nBins,0.5,0.5+Double_t(nBins));

  Int_t ibin = 0;
  Double_t integralData = 0.0;
  Double_t integralMC[4];
  for(Int_t j = 0; j < 4; ++j) integralMC[j] = 0.0;

  for(Int_t i = bin1; i <= bin2; ++i) {
    ibin++;
    hData->SetBinContent(ibin,hThetaData->GetBinContent(i));
    hData->SetBinError(ibin,hThetaData->GetBinError(i));
    integralData += hThetaData->GetBinContent(i);
    for(Int_t j = 0; j < 4; ++j) {
      hMC[j]->SetBinContent(ibin,hThetaMC->GetBinContent(i,j+1));
      hMC[j]->SetBinError(ibin,hThetaMC->GetBinError(i,j+1));
      integralMC[j] += hThetaMC->GetBinContent(i,j+1);
    }
  }
  for(Int_t i = 1; i <= hPhiData->GetNbinsX(); ++i) {
    ibin++;
    hData->SetBinContent(ibin,hPhiData->GetBinContent(i));
    hData->SetBinError(ibin,hPhiData->GetBinError(i));
    integralData += hPhiData->GetBinContent(i);
    for(Int_t j = 0; j < 4; ++j) {
      hMC[j]->SetBinContent(ibin,hPhiMC->GetBinContent(i,j+1));
      hMC[j]->SetBinError(ibin,hPhiMC->GetBinError(i,j+1));
      integralMC[j] += hPhiMC->GetBinContent(i,j+1);
    }
  }
  for(Int_t i = 1; i <= hTildaData->GetNbinsX(); ++i) {
    ibin++;
    hData->SetBinContent(ibin,hTildaData->GetBinContent(i));
    hData->SetBinError(ibin,hTildaData->GetBinError(i));
    integralData += hTildaData->GetBinContent(i);
    for(Int_t j = 0; j < 4; ++j) {
      hMC[j]->SetBinContent(ibin,hTildaMC->GetBinContent(i,j+1));
      hMC[j]->SetBinError(ibin,hTildaMC->GetBinError(i,j+1));
      integralMC[j] += hTildaMC->GetBinContent(i,j+1);
    }
  }

  TObjArray *mc = new TObjArray(4);        // MC histograms are put in this array
  for(Int_t j = 0; j < 4; ++j) mc->Add(hMC[j]);

  TFractionFitter* fit = new TFractionFitter(hData, mc, "V"); // initialise
  // ROOT::Fit::Fitter* fitter = fit->GetFitter();
  // fitter->Config().MinimizerOptions().SetStrategy(2); // tests ...
  // fitter->Config().SetMinosErrors(kTRUE); // tests ...

  //   fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
  //  fit->SetRangeX(1,bin2-bin1+1);                    // use only the first N bins in the fit
  Int_t status = fit->Fit();               // perform the fit
  std::cout << "fit status: " << status << std::endl;
  // We have to work more on this: what to do in case fit fails and we want to redo via diff algo
  // For the moment, there is no solid solution ... comment out these lines
  // if (status == 4) {
  //   ROOT::Fit::Fitter* fitter = fit->GetFitter();
  //   fitter->Config().SetMinimizer("Minuit2","Simplex");
  //   // Double_t tempval[4], temperr[4];
  //   // for(Int_t j = 0; j < 4; ++j) fit->GetResult(j,tempval[j],temperr[j]);
  //   // for(Int_t j = 0; j < 4; ++j) fitter->Config().ParSettings(j).Set(Form("frac%d",j),tempval[j],0.01);
  //   status = fit->Fit();
  //   std::cout << "fit status: " << status << std::endl;
  // }
  if (status == 0 || status == 4) {                       // check on fit status ...
    TH1* result = (TH1*) fit->GetPlot();
    // hData->Draw("Ep");
    // result->Draw("same");
    loutput->Add(hData);
    loutput->Add(result);
    TH1D *hParams = new TH1D(Form("%s_FracParams",hname),"",7,0,7);
    for(Int_t j = 0; j < 4; ++j) {
      Double_t value, error;
      fit->GetResult(j,value,error);
      hParams->SetBinContent(j+1,value);
      hParams->SetBinError(j+1,error);
      TH1 *hFrac = (TH1*)(fit->GetMCPrediction(j)->Clone(Form("%s_FracResult%d",hname,j)));
      hFrac->Scale(value*integralData/integralMC[j]);
      loutput->Add(hMC[j]);
      loutput->Add(hFrac);
    }
    hParams->SetBinContent(5,fit->GetChisquare());
    hParams->SetBinContent(6,fit->GetNDF());
    hParams->SetBinContent(7,status);
    loutput->Add(hParams);
  }

}

void MyFractionFitting(const char *hname,
		       TH1D *hThetaData, TH1D *hPhiData, TH1D *hTildaData,
		       TH2D *hThetaMC, TH2D *hPhiMC, TH2D *hTildaMC,
		       TList *loutput,
		       Double_t thetaMax,
		       Int_t fitIndex)
{
  Int_t bin1 = hThetaData->FindBin(-thetaMax+1e-6);
  Int_t bin2 = hThetaData->FindBin(thetaMax-1e-6);
  Int_t nBins = bin2-bin1+1 + hPhiData->GetNbinsX() + hTildaData->GetNbinsX();
  TH1D *hData = new TH1D(Form("%s_FracData",hname),"",nBins,0.5,0.5+Double_t(nBins));
  TH1D *hMC[4];
  for(Int_t j = 0; j < 4; ++j)
    hMC[j] = new TH1D(Form("%s_FracMC%d",hname,j),"",nBins,0.5,0.5+Double_t(nBins));

  Int_t ibin = 0;
  Double_t integralData = 0.0;
  Double_t integralMC[4];
  for(Int_t j = 0; j < 4; ++j) integralMC[j] = 0.0;

  Int_t index[4] = {1,5,6,7}; // 1-const,5-cos^2theta,6-...
  for(Int_t i = bin1; i <= bin2; ++i) {
    ibin++;
    hData->SetBinContent(ibin,hThetaData->GetBinContent(i));
    hData->SetBinError(ibin,hThetaData->GetBinError(i));
    integralData += hThetaData->GetBinContent(i);
    for(Int_t j = 0; j < 4; ++j) {
      hMC[j]->SetBinContent(ibin,hThetaMC->GetBinContent(i,index[j]));
      hMC[j]->SetBinError(ibin,hThetaMC->GetBinError(i,index[j]));
      integralMC[j] += hThetaMC->GetBinContent(i,index[j]);
    }
  }
  for(Int_t i = 1; i <= hPhiData->GetNbinsX(); ++i) {
    ibin++;
    hData->SetBinContent(ibin,hPhiData->GetBinContent(i));
    hData->SetBinError(ibin,hPhiData->GetBinError(i));
    integralData += hPhiData->GetBinContent(i);
    for(Int_t j = 0; j < 4; ++j) {
      hMC[j]->SetBinContent(ibin,hPhiMC->GetBinContent(i,index[j]));
      hMC[j]->SetBinError(ibin,hPhiMC->GetBinError(i,index[j]));
      integralMC[j] += hPhiMC->GetBinContent(i,index[j]);
    }
  }
  for(Int_t i = 1; i <= hTildaData->GetNbinsX(); ++i) {
    ibin++;
    hData->SetBinContent(ibin,hTildaData->GetBinContent(i));
    hData->SetBinError(ibin,hTildaData->GetBinError(i));
    integralData += hTildaData->GetBinContent(i);
    for(Int_t j = 0; j < 4; ++j) {
      hMC[j]->SetBinContent(ibin,hTildaMC->GetBinContent(i,index[j]));
      hMC[j]->SetBinError(ibin,hTildaMC->GetBinError(i,index[j]));
      integralMC[j] += hTildaMC->GetBinContent(i,index[j]);
    }
  }
  {
    hData->SetStats(0);
    hData->SetLabelColor(0);
    hData->SetLineWidth(2);
    TF1 *fFracFunc = new TF1(Form("%s_fFracFunc",hname),fracFCN,
			     hData->GetXaxis()->GetXmin(),
			     hData->GetXaxis()->GetXmax(),5);
    fFracFunc->FixParameter(4,Double_t(fitIndex)); // transmit the histo index via parameter 4
    fFracFunc->SetNpx(5000);
    fFracFunc->SetLineWidth(2);
    fFracFunc->SetLineColor(kYellow-2);
    for(Int_t j = 0; j < 4; ++j)      hStatic[fitIndex][j] = hMC[j];
    TFitResultPtr r = hData->Fit(fFracFunc,"S");
    TMatrixD cov = r->GetCorrelationMatrix();
    TMatrixD cor = r->GetCovarianceMatrix();
    // cov.Print();
    // cor.Print();

    loutput->Add(hData);
    for(Int_t j = 0; j < 4; ++j) loutput->Add(hMC[j]);

	  TGraph *gr[3][2];
		TGraph *gr2[3][2];
	  for(Int_t i = 0; i < 3; ++i){
	    for(Int_t j = 0; j < 2; ++j){
				gr[i][j] = new TGraph(80);
				gr2[i][j] = new TGraph(80);
			}
			gr2[i][0]->SetName(Form("%s_N_%d_1sigma",fFracFunc->GetName(),i+1));
	  r->Contour(0,i+1,gr2[i][0],0.683);
		gr2[i][1]->SetName(Form("%s_N_%d_2sigma",fFracFunc->GetName(),i+1));
		r->Contour(0,i+1,gr2[i][0],0.683);
		}


	  gr[0][0]->SetName(Form("%s_%d_%d_1sigma",fFracFunc->GetName(),1,2));
	  r->Contour(1,2,gr[0][0],0.683);
    gr[0][1]->SetName(Form("%s_%d_%d_2sigma",fFracFunc->GetName(),1,2));
    r->Contour(1,2,gr[0][1],0.954);
    gr[1][0]->SetName(Form("%s_%d_%d_1sigma",fFracFunc->GetName(),1,3));
    r->Contour(1,3,gr[1][0],0.683);
    gr[1][1]->SetName(Form("%s_%d_%d_2sigma",fFracFunc->GetName(),1,3));
    r->Contour(1,3,gr[1][1],0.954);
    gr[2][0]->SetName(Form("%s_%d_%d_1sigma",fFracFunc->GetName(),2,3));
    r->Contour(2,3,gr[2][0],0.683);
    gr[2][1]->SetName(Form("%s_%d_%d_2sigma",fFracFunc->GetName(),2,3));
    r->Contour(2,3,gr[2][1],0.954);
    for(Int_t i = 0; i < 3; ++i){
      for(Int_t j = 0; j < 2; ++j){
	loutput->Add(gr[i][j]);
	loutput->Add(gr2[i][j]);

}
}
  }
}

Double_t fracFCN(Double_t *x, Double_t *par)
{
  Int_t fitInd = TMath::Nint(par[4]);
  Int_t iBin = hStatic[fitInd][0]->FindBin(x[0]);
  Double_t integral = (hStatic[fitInd][0]->Integral(1,hStatic[fitInd][0]->GetNbinsX())+
		       par[1]*hStatic[fitInd][1]->Integral(1,hStatic[fitInd][1]->GetNbinsX())+
		       par[2]*hStatic[fitInd][2]->Integral(1,hStatic[fitInd][2]->GetNbinsX())+
		       par[3]*hStatic[fitInd][3]->Integral(1,hStatic[fitInd][3]->GetNbinsX()));
  // std::cout<<hStatic[fitInd][0]->GetBinContent(iBin)<<std::endl;

  return (par[0]/integral*(hStatic[fitInd][0]->GetBinContent(iBin)+
			   par[1]*hStatic[fitInd][1]->GetBinContent(iBin)+
			   par[2]*hStatic[fitInd][2]->GetBinContent(iBin)+
			   par[3]*hStatic[fitInd][3]->GetBinContent(iBin)));
}
