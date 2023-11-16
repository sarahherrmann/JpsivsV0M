18d:
Uncorrected for dimuons:

.L deriveResults.C
deriveResults()

Corrected for dimuons:
deriveResults("../fitJpsi/FitResults_18d_CorrDimu.root")

deriveResults("../fitJpsi/FitResults_18d_CorrDimu_DimuEve.root")


02/11/23:
Corrected correction, all with single muon events:
deriveResults("../fitJpsi/FitResults_18d_CorrDimu_SingleMu.root")


---
16/11/2023
With working physics selection + test on the systematics


deriveResults("../fitJpsi/testSys_18d_vNov2.root")

What I sent to Cvetan that was strange was from .L fitJpsi_sys.C
root [1] fitJpsi(0,1,0,"../ReadingTree/JpsiRead_18d_vNov2.root","testSys_18d_vNov2.root")
