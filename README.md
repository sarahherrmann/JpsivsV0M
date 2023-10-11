Local tests:

aliroot "runAnalysisONGRID_full.C(\"18c\", 0, kFALSE, kTRUE, kTRUE)"


Submitted on the grid:

aliroot "runAnalysisONGRID_full.C(\"18c\", 0, kFALSE)"
merge:
aliroot "runAnalysisONGRID_full.C(\"18c\", 0, kTRUE)"
retrieve AnalysisResults.root
aliroot "runAnalysisONGRID_full.C(\"18c\", 0, kTRUE, kFALSE)"

## For git
git push -u origin main
