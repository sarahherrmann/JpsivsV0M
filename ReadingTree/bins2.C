void bins2(const char *fn = "JpsiRead_18c.root")
{
  TFile *f = TFile::Open(fn);
  //TH1F *h = (TH1F*)f->Get("fNV0COverNV0CMean");

  TH1F *h = (TH1F*)f->Get("hV0CMultTotCopy");
  Float_t totalEntries = h->GetEntries();

  // Example works for 10 mult intervals
  Float_t targetPercentiles[] = {0.298,0.198,0.098,0.098,0.098,0.048,0.048,0.048,0.038,0.01};
  //corresponding percentiles ranges 100-70%, 70-50, 50-40, 40-30, 30-20, 20-15, 15-10, 10-5, 5-1, 1-0

  Float_t sum = 0, entries = 0;
  Int_t currentInt = 0;
  Double_t V0values[10];
  Double_t V0limits[9];
  Double_t V0percentiles[10];
  for(Int_t ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
    Float_t cont = h->GetBinContent(ibin);
    Float_t mult = h->GetXaxis()->GetBinCenter(ibin);
    entries += cont;
    sum += (cont*mult);
    if ((currentInt < 9) && (entries/totalEntries) > targetPercentiles[currentInt]) {
      V0values[currentInt] = sum/entries;//average in that bin
      V0limits[currentInt] = h->GetXaxis()->GetBinUpEdge(ibin);
      V0percentiles[currentInt] = entries/totalEntries;
      sum = entries = 0;
      currentInt++;
    }
    //    printf("
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
