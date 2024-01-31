// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly
#include "AliAnaMCMultEpsilonTask.h"
#include "/Users/herrmann/aliceAli/AliPhysics/OADB/macros/AddTaskPhysicsSelection.C"
#include "/Users/herrmann/aliceAli/AliPhysics/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"


const Int_t runList17d18[] = { 264076, 264078, 264082, 264085, 264086, 264109, 264110, 264129, 264137, 264138, 264139, 264164, 264168, 264188, 264190, 264194, 264197, 264198, 264232, 264233, 264235, 264238, 264259, 264260, 264261, 264262, 264264, 264265, 264266, 264267, 264273, 264277, 264279, 264281, 264305, 264306, 264312, 264336, 264341, 264345, 264346, 264347};
const Int_t runList19h11a[] = { 286380, 286426, 286427, 286428, 286454, 286455, 286482, 286502, 286508, 286509, 286511, 286566, 286567 ,286568, 286569, 286591, 286592, 286653, 286661, 286695, 286731, 286799, 286801, 286805, 286809, 286846, 286850, 286852, 286874, 286876, 286877, 286907, 286910, 286911, 286930, 286931, 286932, 286933, 286936, 286937, 287000, 287021, 287063, 287064, 287066, 287071, 287072, 287077, 287137, 287155, 287185, 287201, 287202, 287203, 287204, 287208, 287209, 287248, 287249, 287250, 287251, 287254, 287283, 287323, 287324, 287325, 287343, 287344, 287346, 287347, 287349, 287353, 287355, 287356, 287360, 287380, 287381, 287385, 287387, 287388, 287389, 287413, 287451, 287480, 287481, 287484, 287486, 287513, 287516, 287517, 287518, 287521, 287524, 287575, 287578, 287654, 287656, 287657, 287658};

void runAnalysisMCMultEpsilon(const char *fperiod = "17d18", Bool_t merge = kFALSE, Bool_t mergeJDL = kTRUE, Bool_t local = kFALSE)
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    //Bool_t local = kFALSE;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t gridTest = kFALSE;

    //________________Define the periods, the run lists and the path of the data
    TString periodStr = fperiod;
    Int_t nRuns = 0;
    const Int_t *runList;
    if (periodStr.Contains("17d18")) { runList = runList17d18; nRuns = sizeof(runList17d18)/sizeof(runList17d18[0]); }
    else if (periodStr.Contains("19h11a_extra")) { runList = runList19h11a; nRuns = sizeof(runList19h11a)/sizeof(runList19h11a[0]); }
    else {
      printf("unknown period: %s",fperiod);
      return;
    }

    printf("%d runs will be analyzed in %s\n",nRuns,fperiod);

    TString dataDir;
    if (periodStr.Contains("17")) dataDir = Form("/alice/sim/2017/LHC%s",fperiod);
    else if (periodStr.Contains("16")) dataDir = Form("/alice/sim/2016/LHC%s",fperiod);
    else if (periodStr.Contains("19")) dataDir = Form("/alice/sim/2019/LHC%s",fperiod);
    else {
      printf("unknown year: %s",fperiod);
      return;
    }

    TString dataPattern = "/AOD235/*/AliAOD.root";

    //if (periodStr.Contains("17h")) { dataPattern = "/muon_calo_pass2/AOD/*/AliAOD.root"; }

    TString workingDir = Form("LHC%s_MCMultEpsilon",fperiod);

    //___________EOF: Define the periods, the run lists and the path of the data


    // since we will compile a class, tell root where to look for headers
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    printf("Loading physics selection task...\n");
    //  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kTRUE,kFALSE,0,kFALSE);//2nd kFALSE= no pileup cuts

    printf("Loading multiplicity selection task...\n");
    //  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask * multTask = AddTaskMultSelection(kFALSE); // user mode:

    // compile the class and load the add task macro
    // here we have to differentiate between using the just-in-time compiler
    // from root6, or the interpreter of root5
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AliAnaMCMultEpsilonTask.cxx++g");
    AliAnaMCMultEpsilonTask *task = reinterpret_cast<AliAnaMCMultEpsilonTask*>(gInterpreter->ExecuteMacro("AddMCMultEpsilonTask.C"));
#else
    gROOT->LoadMacro("AliAnaMCMultEpsilonTask.cxx++g");
    gROOT->LoadMacro("AddMCMultEpsilonTask.C");
    AliAnaMCMultEpsilonTask *task = AddMCMultEpsilonTask();
#endif


    if(!mgr->InitAnalysis()) return;
    //mgr->SetDebugLevel(2);

    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);
    mgr->SetDebugLevel(5);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        // add a few files to the chain (change this so that your local files are added)
        chain->Add("../MCTask/AliAOD.root");
        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);
    } else {

        if (!TGrid::Connect("alien://")) return;

        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        TString includes_str = "-Wno-deprecated -I$. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include";
        alienHandler->AddIncludePath(includes_str.Data()); // for grid running
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnaMCMultEpsilonTask.cxx AliAnaMCMultEpsilonTask.h");
        alienHandler->SetAnalysisSource("AliAnaMCMultEpsilonTask.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!

        alienHandler->SetMergeViaJDL(mergeJDL);
        if (!merge) {
          alienHandler->SetRunMode("full");
        }
        else {
          alienHandler->SetRunMode("terminate");
        }

        alienHandler->SetAliPhysicsVersion("vAN-20231002_O2-1");//VO_ALICE@AliPhysics::vAN-20231002_O2-1
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        // select the input data
        alienHandler->SetGridDataDir(dataDir.Data());
        alienHandler->SetDataPattern(dataPattern.Data());

        //TString dataPattern = "/muon_calo_pass1/AOD/*/AliAOD.root";
        // MC has no prefix, data has prefix 000
        alienHandler->SetRunPrefix("");
        // runnumber
        for (Int_t i=0;i<nRuns;i++)  alienHandler->AddRunNumber(runList[i]);
        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(75);
        alienHandler->SetExecutable("myTask.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(64800);
        alienHandler->SetJDLName("myTask.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);

        alienHandler->SetNrunsPerMaster(1);

        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate")
        // (see below) mode, set SetMergeViaJDL(kFALSE)
        // to collect final results
        alienHandler->SetMaxMergeStages(2);//alienHandler->SetMaxMergeStages(1);
        //alienHandler->SetMergeViaJDL(kTRUE); //do not retrieve final merged output
        //alienHandler->SetMergeViaJDL(kFALSE);  //retrieve final merged output

        // define the output folders
        alienHandler->SetGridWorkingDir(workingDir.Data());
        alienHandler->SetGridOutputDir("Output");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);


        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(21);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            //alienHandler->SetRunMode("full"); //not merging -step1
            //alienHandler->SetRunMode("terminate"); //merging
            mgr->StartAnalysis("grid");
        }
    }
}
