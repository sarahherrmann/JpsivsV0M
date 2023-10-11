#include "treeReaderHelper.h"
#include "AliAODTrack.h"

void readTree(const char *inputFile = "../local_AnalysisResults.root", const char *fperiod = "18b")
{
  TFile *fIn = TFile::Open(inputFile);

  TTree *trDiMu = (TTree *)fIn->Get("JpsiVsV0M/DiMuonEvents");
  TTree *trSingleMu = (TTree *)fIn->Get("JpsiVsV0M/SingleMuonEvents");
  TTree *trInt7 = (TTree *)fIn->Get("JpsiVsV0M/kINT7Events");
  TList *histoList = (TList *)fIn->Get("JpsiVsV0M/kINT7Histograms");

  treeReaderHelper treeReader;

  treeReader.SetMassCuts(2.0,5.0);
  treeReader.SetJpsiEtaCuts(-4.0,-2.5);
  treeReader.InitHistograms();
  treeReader.SetTreeAddresses(trDiMu, trSingleMu, trInt7);
  treeReader.GetNV0Histogram(histoList, trInt7);
  treeReader.processV0PerChannelSingleMu(trSingleMu);
  treeReader.DeriveExcessV0MeanPerChannel();
  treeReader.readEvents(trDiMu);
  treeReader.writeOutput(Form("JpsiRead_%s.root",fperiod));


}//readTree
