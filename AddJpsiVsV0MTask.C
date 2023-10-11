AliAnaTaskJpsiVsV0M* AddJpsiVsV0MTask(TString name = "name")
{
    // get the manager via the static access member. since it's static, you don't need
    // to create an instance of the class here to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    // get the input event handler, again via a static method.
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":JpsiVsV0M";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnaTaskJpsiVsV0M* task = new AliAnaTaskJpsiVsV0M(name.Data());
    if(!task) return 0x0;
    task->SelectCollisionCandidates(AliVEvent::kAny);//all triggers accepted
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("kINT7Histograms", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    //added by me for the tree
    mgr->ConnectOutput(task,2,mgr->CreateContainer("DiMuonEvents", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));//same name as the tree

    mgr->ConnectOutput(task,3,mgr->CreateContainer("SingleMuonEvents", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));//same name as the tree

    mgr->ConnectOutput(task,4,mgr->CreateContainer("kINT7Events", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));//same name as the tree
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
