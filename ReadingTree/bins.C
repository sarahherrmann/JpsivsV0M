void bins(const char *fn = "JpsiRead_18c.root")
{
  TFile *f = TFile::Open(fn);
  TH1F *h = (TH1F*)f->Get("fNV0OverNV0Mean");
  Float_t totalEntries = h->GetEntries();

  // Example works for 10 mult intervals

  //Float_t targetMult[] = {0.38,1.01,1.54,2.03,2.56,3.23,3.92,4.68,5.67,7.28}; // last one is not really needed, in practice it is 'inf'

  //Float_t targetMult[] = {0.4,0.76,1.41,2.26,3.03,3.92,4.33,4.96,5.67,7.28}; // from https://arxiv.org/pdf/2005.11123v2.pdf + 2 last points from previous

  Float_t targetMult[] = {0.4,0.76,1.41,2.26,3.03,3.92,4.33,4.96,5.85,6.55}; // from https://arxiv.org/pdf/2005.11123v2.pdf + 2 last points custom

  Float_t sum = 0, entries = 0;
  Int_t currentInt = 0;
  Double_t V0values[10];//Double_t V0values[10];
  Double_t V0limits[9];
  Double_t V0percentiles[10];
  for(Int_t ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
    Float_t cont = h->GetBinContent(ibin);
    Float_t mult = h->GetXaxis()->GetBinCenter(ibin);
    entries += cont;
    sum += (cont*mult);
    if ((currentInt < 9) && (sum/entries) > targetMult[currentInt]) {
      V0values[currentInt] = sum/entries;
      V0limits[currentInt] = h->GetXaxis()->GetBinUpEdge(ibin);
      V0percentiles[currentInt] = entries/totalEntries;
      sum = entries = 0;
      currentInt++;
    }
  }
  // fill the values for the last interval
  V0values[9] = sum/entries;
  V0percentiles[9] = entries/totalEntries;

  for(Int_t i = 0; i < 10; ++i) {
    printf("Interval %d:  [%f -> %f]   NV0/<NV0> = %f  percentile = %f\n",
	   i,(i==0)?0:V0limits[i-1],(i==9)?-1:V0limits[i],
	   V0values[i],V0percentiles[i]);
  }
}
