fEtaPhiV0MeanPerChannelOneMu[channelMu] contains, in each (eta, phi) bin (that corresponds to 1 V0 channel i),
the mean V0 multiplicity (signal) in channel i

Each TProfile2D fEtaPhiV0MeanPerChannelOneMu[j] contains what is above, when the muon falls into the channel j
allows us to see how far from the muon the V0 channel mult are increased


To derive the excess one should do
for each channel k

fEtaPhiV0MeanPerChannelOneMu[k]-fEtaPhiV0MeanPerChannelNoMu (but you can't substract one TProfile to a TH2F)
GetBinContent (Int_t binx, Int_t biny)//works for both
GetNbinsX
GetNbinsY

hDiffEtaPhiMeanV0[k]->SetBinContent(etabin, phibin, diff);//TH2F containing the excess mean mult

this will be in DeriveExcessV0MeanPerChannel()
