Local tests:

aliroot "runAnalysisONGRID_full.C(\"18c\", 0, kFALSE, kTRUE, kTRUE)"


Submitted on the grid:

aliroot "runAnalysisONGRID_full.C(\"18c\", 0, kFALSE)"
merge:
aliroot "runAnalysisONGRID_full.C(\"18c\", 0, kTRUE)"
retrieve AnalysisResults.root
aliroot "runAnalysisONGRID_full.C(\"18c\", 0, kTRUE, kFALSE)"


---
18b

aliroot "runAnalysisONGRID_full.C(\"18b\", 0, kFALSE)"
aliroot "runAnalysisONGRID_full.C(\"18b\", 0, kTRUE)"
aliroot "runAnalysisONGRID_full.C(\"18b\", 0, kTRUE, kFALSE)"
## For git
git push -u origin main


Note: 12/10/23 11:15 modification of the task code, added percentile for kINT7 events
+ added more pileup selections, new folder name ("LHC%sJpsiV0M_v2",fperiod)

Should redo 18b and 18c
---
LHC18d

aliroot "runAnalysisONGRID_full.C(\"18d\", 0, kFALSE)"
