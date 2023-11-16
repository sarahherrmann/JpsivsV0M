Inside analyzeOfLucrezia.C, have a look at the function FitJpsi()


.L fitJpsi_v2.C
fitJpsi("../ReadingTree/JpsiRead_18d.root", "FitResults_18d_CorrDimu_DimuEve.root", "hMultPtJpsiMinvCorr")


fitJpsi("../ReadingTree/JpsiRead_18d.root", "FitResults_18d.root")


02/11/2023 : correct single muon correction (all computed in signle muon events)


fitJpsi("../ReadingTree/JpsiRead_18d.root", "FitResults_18d_CorrDimu_SingleMu.root", "hMultPtJpsiMinvCorr")


---
.L fitJpsi_sys.C
fitJpsi(0,1,0,"../ReadingTree/JpsiRead_18d_vNov2.root","testSys_18d_vNov2.root")
