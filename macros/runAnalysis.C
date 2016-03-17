/**************************************************************************
 * Copyright(c) 2013-2014, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifdef __ECLIPSE_IDE
//  few includes and external declarations just for the IDE
#include "Riostream.h"
#include "TSystem.h"
#include "TChain.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliAnalysisTaskPIDqa.h"
#include "AliPhysicsSelectionTask.h"
#include "AliCentralitySelectionTask.h"
#include "AliTaskCDBconnect.h"
#include "AliAnalysisTaskPIDCombined.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMultSelectionTask.h"
extern AliAnalysisGrid* CreateAlienHandler(const char *, Bool_t);
extern AliCentralitySelectionTask *AddTaskCentrality(Bool_t fillHistos=kTRUE, Bool_t aod=kFALSE);
extern AliAnalysisTask *AddTaskPIDResponse(Bool_t, Bool_t,Bool_t,Int_t);
//extern AliAnalysisTask *AddTaskPIDqa(const char *);
extern TChain* CreateESDChain(
  const char* aDataDir = "ESDfiles.txt",
  Int_t aRuns          = 20,
  Int_t offset         = 0,
  Bool_t addFileName   = kFALSE,
  Bool_t addFriend     = kFALSE,
  const char* check    = 0);
extern TChain* CreateAODChain(
    const char* aDataDir = "AODfiles.txt",
    Int_t aRuns          = 20,
    Int_t offset         = 0,
    Bool_t addFileName   = kFALSE,
    const char* friends  = "",
    const char* check    = 0);
AliPhysicsSelectionTask* AddTaskPhysicsSelection(
    Bool_t mCAnalysisFlag = kFALSE,
    Bool_t deprecatedFlag = kTRUE,
    UInt_t computeBG = 0,
    Bool_t useSpecialOutput=kFALSE);
extern AliAnalysisTask *AddTaskTender(Bool_t,Bool_t,Bool_t,Bool_t,Bool_t,Bool_t,Bool_t,Bool_t,Bool_t);
extern AliTaskCDBconnect* AddTaskCDBconnect(const char *path="raw://", Int_t run=0);
extern AliAnalysisTaskPhiCorrelations *AddTaskPhiCorrelations(Int_t analysisMode = 0,
    Bool_t ppRun = kFALSE,
    const char* outputFileName = 0,
    Bool_t eventMixing = kTRUE,
    Int_t zVtxAxis = 0,
    const char* containerName = "histosPhiCorrelations",
    const char* folderName = "PWG4_PhiCorrelations");
AliAnalysisDataContainer* AddTaskFlowQnVectorCorrections(const char *inputHistogramFileName);
#endif // ifdef __ECLIPSE_IDE declaration and includes for the ECLIPSE IDE

#include "runAnalysisCriteria.H"
#include "runAnalysis.H"

void runAnalysis(const char *inputHistogramFileName = "", const char *sRunMode = "full", Bool_t gridMerge = kTRUE) {
  /* WARNING!!! HANDLE WITH CARE!!! precedes GRID and MC */
  Bool_t bASCIIoutput           = kFALSE;

  Bool_t bUseESD                     = kTRUE;
  Bool_t bUseAOD                     = kFALSE;

  Bool_t useTender              = kTRUE;
  Bool_t useCDB                 = kFALSE;
  Bool_t isPhysicsSelection     = kTRUE;

  /* PIDresponse is not optional any longer so,
   * do NOT change this flag at all */
  Bool_t usePIDResponse         = kTRUE;

  Bool_t bRunPIDCombinedTask = kFALSE;

  Bool_t isHeavyIon     = kTRUE;

  // Centrality
  Bool_t UseCentrality  = kTRUE;
  Bool_t UseMultiplicity = kTRUE;

  gROOT->LoadMacro("loadRunOptions.C");
  loadRunOptions();
  gROOT->LoadMacro("runs.C");

  /* SANITY CHECKS! DON'T CHANGE THIS! */
  if (bASCIIoutput) {
    bGRIDPlugin = kFALSE;
  }

  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");

  AliAnalysisGrid      *alienHandler   =   NULL;
  AliAnalysisManager *mgr = new AliAnalysisManager("Flow Qn vector corrections");

  // Enable debug printouts
  mgr->SetDebugLevel(AliLog::kError);
  AliLog::SetGlobalLogLevel(AliLog::kFatal);


  if (bUseESD) {
    AliESDInputHandler* esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
  //  esdH->SetReadFriends(kFALSE);
  }
  else if (bUseAOD) {
    AliAODInputHandler* aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
  }
  else {
    cout << "Neither AOD nor ESD data specified. ABORTING!!!" << endl;
    return;
  }

  if(bGRIDPlugin) {
    gROOT->LoadMacro("CreateAlienHandler.C");
    alienHandler = CreateAlienHandler(sRunMode,gridMerge);
    if (!alienHandler) return;

    mgr->SetGridHandler(alienHandler);
  }

  if(bMC ){
      AliMCEventHandler* mcH = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mcH);
  }

  if( usePIDResponse ) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(bMC,kTRUE,kTRUE,2);
  }

  if( useTender ) {
    gROOT->LoadMacro("$ALICE_PHYSICS/TENDER/TenderSupplies/AddTaskTender.C");
    AliAnalysisTaskSE *tender = AddTaskTender(kFALSE,  /* V0 */
                                              kTRUE,   /* TPC */
                                              kTRUE,   /* TOF */
                                              kFALSE,  /* TRD */
                                              kTRUE,   /* PID */
                                              kFALSE,  /* VTX */
                                              kTRUE,   /* T0 */
                                              kFALSE,   /* EMCAL, no EMCAL for the time being */
                                              kFALSE); /* TRACKFIX */
  }

  if(bUseESD) {
    if(isPhysicsSelection) {
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(bMC);
    }

    if( UseCentrality  ) {
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
      AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
      if (bMC) {
        taskCentrality->SetMCInput();
      }
    }

    if (UseMultiplicity) {
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
      AliMultSelectionTask * task = AddTaskMultSelection(kFALSE); // user mode:
    }
  }


  if (useCDB) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    AliTaskCDBconnect* cdbTask = AddTaskCDBconnect();
  }

  AliAnalysisDataContainer *cinput1 = NULL;
  cinput1 = mgr->GetCommonInputContainer();

  if (bRunPIDCombinedTask == kTRUE){
    AliAnalysisTaskPIDCombined *taskPIDCombined = new AliAnalysisTaskPIDCombined("PID Combined");
    mgr->AddTask(taskPIDCombined);

    AliAnalysisDataContainer *coutput = mgr->CreateContainer("PID Combined", TList::Class(),  AliAnalysisManager::kOutputContainer, "Viscosity.root");
    mgr->ConnectInput( taskPIDCombined, 0,  cinput1);
    mgr->ConnectOutput(taskPIDCombined, 1,  coutput);
    cout << "Info: created PID combined task" << endl;
  }


  /* load the runs list we will analyse */
  loadRuns();
  loadReasons();

  TString debugString="+g";

  gROOT->LoadMacro("QnCorrectionsLog.cxx"+debugString);
  gROOT->LoadMacro("QnCorrectionsEventClasses.cxx"+debugString);
  gROOT->LoadMacro("QnCorrectionsCuts.cxx"+debugString);
  gROOT->LoadMacro("QnCorrectionsHistograms.cxx"+debugString);
  gROOT->LoadMacro("QnCorrectionsDataVector.cxx"+debugString);
  gROOT->LoadMacro("QnCorrectionsQnVector.cxx"+debugString);
  gROOT->LoadMacro("QnCorrectionsCorrectionSteps.cxx"+debugString);
  gROOT->LoadMacro("QnCorrectionsDetector.cxx"+debugString);
  gROOT->LoadMacro("QnCorrectionsManager.cxx"+debugString);
  gROOT->LoadMacro("QnCorrectionsInputGainEqualization.cxx"+debugString);


  gROOT->LoadMacro("QnCorrectionsVarManagerTask.cxx"+debugString);
  gROOT->LoadMacro("QnCorrectionsFillEventTask.cxx"+debugString);
  gROOT->LoadMacro("AnalysisTaskFlowVectorCorrections.cxx"+debugString);


  gROOT->LoadMacro("AddTaskFlowQnVectorCorrections.C"+debugString);
  AliAnalysisDataContainer *corrTask = AddTaskFlowQnVectorCorrections(inputHistogramFileName);

  TChain* chain = 0;

  if( ! bGRIDPlugin ) {
    Int_t numFiles = 0;
    ifstream dataInStream;
    dataInStream.open(szLocalFileList.Data());
    if ( !dataInStream ) {
      cout<<"Data list file does not exist: "<<szLocalFileList.Data()<<endl;
      return;
    }

    string line;

    while ( !dataInStream.eof() ) {
      getline(dataInStream, line);
      if(line.compare("") != 0) {//checks if there is an empty line in the data list
        numFiles++;
      }
    }
    // No need to create a chain - this is handled by the plugin
    if (bUseESD) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
      chain = CreateESDChain(szLocalFileList,numFiles);
    }
    else if (bUseAOD) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
      chain = CreateAODChain(szLocalFileList,numFiles);
    }
  }

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();
  // Start analysis in grid.
  if (bGRIDPlugin){
    mgr->StartAnalysis("grid");
  }
  else{
    mgr->StartAnalysis("local",chain);
  }
}

Int_t setupPar(const char* pararchivename) {
  ///////////////////
  // Setup PAR File//
  ///////////////////
  if (pararchivename) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchivename);

    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");

      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
        Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
        return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
    gSystem->ChangeDirectory("../");
  }
  return 1;
}



