void binsV0C(const char *fn = "JpsiRead_18c.root", const char *hname = "hV0CMultTotCopy", const char *listname = NULL)
{
  TFile *f = TFile::Open(fn);
  //TH1F *h = (TH1F*)f->Get("fNV0COverNV0CMean");

  TH1F *h;
  if (*listname)//the pointer is not NULL
  {
    //then the histogram is actually in a histolist, we have to retrive the list first
    TList *histoList = (TList *)f->Get(listname);
    h = (TH1F *)histoList->FindObject(hname);//and then retrieve the histo
  }
  else
  {
    h = (TH1F*)f->Get(hname);
  }

  Float_t totalEntries = h->GetEntries();

  // Example works for n mult intervals
  Float_t targetPercentiles[] = {0.298,0.199,0.099,0.099,0.099,0.0498,0.0498,0.0498,0.0398,0.099};
  //corresponding percentiles ranges 100-70%, 70-50, 50-40, 40-30, 30-20, 20-15, 15-10, 10-5, 5-1, 1-0

  //Float_t targetPercentiles[] = {0.298,0.199,0.099,0.099,0.099,0.0498,0.0498,0.0498,0.0398,0.0049,0.00249,0.00249};
  //corresponding percentiles ranges 100-70%, 70-50, 50-40, 40-30, 30-20, 20-15, 15-10, 10-5, 5-1, 1-0.5, 0.5-0.25, 0.25-0

  const int n = sizeof(targetPercentiles)/sizeof(targetPercentiles[0]);
  printf("Number of mult intervals %d\n", n);

  Float_t NV0mean = h->GetMean();
  Float_t binWidth = h->GetBinWidth(1);

  Float_t sum = 0, entries = 0;
  Int_t currentInt = 0;
  Double_t V0values[n];
  Double_t V0limits[n-1];
  Double_t V0percentiles[n];
  for(Int_t ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
    Float_t cont = h->GetBinContent(ibin);
    Float_t mult = h->GetXaxis()->GetBinCenter(ibin);
    entries += cont;
    sum += (cont*mult);
    if ((currentInt < n-1) && (entries/totalEntries) > targetPercentiles[currentInt]) {
      V0values[currentInt] = sum/entries;//average in that bin
      V0limits[currentInt] = h->GetXaxis()->GetBinUpEdge(ibin);
      V0percentiles[currentInt] = entries/totalEntries;
      sum = entries = 0;
      currentInt++;
    }
    //    printf("
  }
  // fill the values for the last interval
  V0values[n-1] = sum/entries;
  V0percentiles[n-1] = entries/totalEntries;

  printf("<NV0> = %f, binwidth = %f \n", NV0mean, binWidth);

  for(Int_t i = 0; i < n; ++i) {
    //printf("Interval %d:  [%f -> %f]   NV0/<NV0> = %f  percentile = %f\n",
    printf("Interval %d:  [%f -> %f]   NV0 = %f  percentile = %f\n",
	   i,(i==0)?0:V0limits[i-1],(i==n-1)?-1:V0limits[i],
	   V0values[i],V0percentiles[i]);



     //NV0= NV0 average in corresponding interval
  }
  for(Int_t i = 0; i < n; ++i) {
  printf("Corresponding Interval %d:  [%f -> %f]   NV0/<NV0> = %f  percentile = %f\n",
  i,(i==0)?0:(V0limits[i-1]/NV0mean),(i==n-1)?-1:(V0limits[i]/NV0mean),
  (V0values[i]/NV0mean),V0percentiles[i]);
  }
}
