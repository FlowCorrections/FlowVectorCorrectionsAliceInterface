/***************************************************************************
 * Package:       FlowVectorCorrections ALICE glue                         *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2014-2016                                                *
 ***************************************************************************/
//////////////////////////////////////////////////////////////
//
//    Flow Qn Vector correction options:
//
//
//    Calibration method: 0: Calibration of Q/sqrt(M)
//                        1: Calibration of Q/M
//                        2: Calibration of Q/|Q|
//                        3: Calibration of Q
//
//    Calibration step  : 0: Raw
//                        1: Equalization
//                        2: Recentering
//                        3: Twist
//                        4: Scaling
//
//    Equalization method : 0: M/<M>
//                          1: 1+(M-<M>)/sigma(M)
//
//    Channel list      : Array of channel numbers that are included in Q-vector calculation
//
//    Twist and Scaling method: to be implemented
//                              0: Double harmonic track wise (advisable for TPC)
//                              1: Double harmonic Q wise
//                              2: Correlations
//
///////////////////////////////////////////////////////////////


#include <TTree.h>
#include <TSystem.h>
#include <TMath.h>
#include <TROOT.h>
#include <Riostream.h>
#include <TObject.h>
#include <TFile.h>
#include <AliLog.h>
#include "AliAODHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisManager.h"
#include "AliQnCorrectionsHistos.h"

#include "QnCorrectionsEventClasses.h"
#include "QnCorrectionsCuts.h"
#include "QnCorrectionsHistograms.h"
#include "QnCorrectionsDataVector.h"
#include "QnCorrectionsQnVector.h"
#include "QnCorrectionsDetector.h"
#include "QnCorrectionsManager.h"
#include "QnCorrectionsInputGainEqualization.h"
#include "QnCorrectionsQnVectorRecentering.h"
#include "AnalysisTaskFlowVectorCorrections.h"

#include "runAnalisis.H"

#define VAR QnCorrectionsVarManagerTask

void DefineHistograms(QnCorrectionsManager* QnManager, AliQnCorrectionsHistos* histos, TString histClass);

void AddVZERO(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager);
void AddTPC(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager);
void AddTZERO(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager);
void AddFMD(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager);
void AddRawFMD(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager);
void AddZDC(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager);
void AddSPD(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager);

Int_t varForEventMultiplicity;

AliAnalysisDataContainer* AddTaskFlowQnVectorCorrections(const char *inputCalibrationFilename) {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFlowQnVectorCorrections", "No analysis manager found.");
    return 0;
  }

  QnCorrectionsManager *QnManager = new QnCorrectionsManager();
  AnalysisTaskFlowVectorCorrections *taskQnCorrections = new AnalysisTaskFlowVectorCorrections("FlowQnVectorCorrections");

  /* let's establish the event cuts for event selection */
  QnCorrectionsCutsSet *eventCuts = new QnCorrectionsCutsSet();
  eventCuts->Add(new QnCorrectionsCutWithin(VAR::kVtxZ,zvertexMin,zvertexMin));
  if (bUseMultiplicity) {
    varForEventMultiplicity = VAR::kVZEROMultPercentile;
  }
  else {
    varForEventMultiplicity = VAR::kCentVZERO;
  }
  eventCuts->Add(new QnCorrectionsCutWithin(varForEventMultiplicity,centralityMin,centralityMax));
  taskQnCorrections->SetEventCuts(eventCuts);

  /* and the physics selection also */
  if (!b2015DataSet) {
    taskQnCorrections->SelectCollisionCandidates(AliVEvent::kMB);  // Events passing trigger and physics selection for analysis
  }
  else
    taskQnCorrections->SelectCollisionCandidates(AliVEvent::kMB|AliVEvent::kINT7);  // Events passing trigger and physics selection for analysis

  taskQnCorrections->SetUseTPCStandaloneTracks(kFALSE);  // Use of TPC standalone tracks or Global tracks (only for ESD analysis)

  TString histClass = "";
  histClass += "Event_NoCuts;";
  histClass += "Event_Analysis;";
  histClass+= "TrackQA_NoCuts;";

  if (b2015DataSet) {
    AddTPC(taskQnCorrections, QnManager);
    histClass+= "TrackQA_TPC;";
    AddSPD(taskQnCorrections, QnManager);
    histClass+= "TrackletQA_SPD;";
    AddVZERO(taskQnCorrections, QnManager);
    AddTZERO(taskQnCorrections, QnManager);
    AddRawFMD(taskQnCorrections, QnManager);
    AddZDC(taskQnCorrections, QnManager);
  }
  else {
    AddTPC(taskQnCorrections, QnManager);
    histClass+= "TrackQA_TPC;";
//    AddSPD(taskQnCorrections, QnManager);
    histClass+= "TrackletQA_SPD;";
    AddVZERO(taskQnCorrections, QnManager);
    AddTZERO(taskQnCorrections, QnManager);
    AddFMD(taskQnCorrections, ,QnManager);
//    AddRawFMD(taskQnCorrections, QnManager);
    AddZDC(taskQnCorrections, QnManager);
  }


  QnManager->SetShouldFillQnVectorTree(kFALSE);
  QnManager->SetShouldFillQAHistograms(kTRUE);
  QnManager->SetShouldFillOutputHistograms(kTRUE);

  taskQnCorrections->SetFillExchangeContainerWithQvectors(kTRUE);
  taskQnCorrections->SetFillEventQA(kTRUE);

  taskQnCorrections->SetQnCorrectionsManager(QnManager);
  taskQnCorrections->DefineInOutput();
  taskQnCorrections->SetRunsLabels(&listOfRuns);

  /* let's get the calibration file */
  TString inputfilename = inputCalibrationFilename;
  if (inputfilename.Length() != 0) {
    TFile *inputfile = new TFile(inputCalibrationFilename);
    if (inputfile != NULL) {
      taskQnCorrections->SetCalibrationHistograms(inputfile);
    }
  }

  AliQnCorrectionsHistos* hists = taskQnCorrections->GetEventHistograms();
  DefineHistograms(QnManager, hists, histClass);


  mgr->AddTask(taskQnCorrections);

  mgr->ConnectInput(taskQnCorrections,  0, mgr->GetCommonInputContainer());

  //create output containers
  if (QnManager->GetShouldFillOutputHistograms()) {
    AliAnalysisDataContainer *cOutputHist =
      mgr->CreateContainer(QnManager->GetCalibrationHistogramsContainerName(),
          TList::Class(),
          AliAnalysisManager::kOutputContainer,
          "CalibrationHistograms.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotHistQn(), cOutputHist );
  }

  if (QnManager->GetShouldFillQnVectorTree()) {
    AliAnalysisDataContainer *cOutputQvec =
      mgr->CreateContainer("CalibratedQvector",
          TTree::Class(),
          AliAnalysisManager::kOutputContainer,
          "QvectorsTree.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotTree(), cOutputQvec );
  }

  if (QnManager->GetShouldFillQAHistograms()) {
    AliAnalysisDataContainer *cOutputHistQA =
      mgr->CreateContainer(QnManager->GetCalibrationQAHistogramsContainerName(),
          TList::Class(),
          AliAnalysisManager::kOutputContainer,
          "CalibrationQA.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotHistQA(), cOutputHistQA );
  }

  if (taskQnCorrections->GetFillEventQA()) {
    AliAnalysisDataContainer *cOutputQnEventQA =
      mgr->CreateContainer("QnEventQA",
          TList::Class(),
          AliAnalysisManager::kOutputContainer,
          "QnEventQA.root");
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotEventQA(), cOutputQnEventQA );
  }

  AliAnalysisDataContainer *cOutputQvecList =
    mgr->CreateContainer("CalibratedQvectorList",
        TList::Class(),
        AliAnalysisManager::kExchangeContainer,
        "QvectorsList.root");

  if (taskQnCorrections->GetFillExchangeContainerWithQvectors())
    mgr->ConnectOutput(taskQnCorrections, taskQnCorrections->OutputSlotGetListQnVectors(), cOutputQvecList );

  return cOutputQvecList;
}




void AddVZERO(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager){


  Bool_t VZEROchannels[4][64];
  for(Int_t iv0=0; iv0<4; iv0++) for(Int_t ich=0; ich<64; ich++) VZEROchannels[iv0][ich] = kFALSE;

  for(Int_t ich=32; ich<64; ich++) VZEROchannels[0][ich] = kTRUE;  // channel list: kTRUE if channel should be used
  for(Int_t ich=0; ich<32; ich++) VZEROchannels[1][ich] = kTRUE;


  Int_t channelGroups[64];
  for(Int_t ich=0; ich<64; ich++) channelGroups[ich] = Int_t(ich / 8);

  //-----------------------------------------------------------
  // Our event classes for V0
  //
  const Int_t nVZEROdim = 2;
  QnCorrectionsEventClassVariablesSet *CorrEventClasses = new QnCorrectionsEventClassVariablesSet(nVZEROdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(VAR::kVtxZ,
      task->VarName(VAR::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the VZERO detector */
  QnCorrectionsDetector *VZERO = new QnCorrectionsDetector("VZERO", VAR::kVZERO);

  /* the VZEROA detector configuration */
  QnCorrectionsDetectorConfigurationChannels *VZEROAconf =
      new QnCorrectionsDetectorConfigurationChannels(
          "VZEROA",
          CorrEventClasses,
          64, /* number of channels */
          4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROAconf->SetChannelsScheme(VZEROchannels[0], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROAconf->SetQVectorNormalizationMethod(QVNORM_QoverM);
  /* lets configure the equalization of input data */
  QnCorrectionsInputGainEqualization *eqA = new QnCorrectionsInputGainEqualization();
  eqA->SetEqualizationMethod(GEQUAL_averageEqualization);
  eqA->SetShift(1.0);
  eqA->SetScale(0.1);
  eqA->SetUseChannelGroupsWeights(kTRUE);
  VZEROAconf->AddCorrectionOnInputData(eqA);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROAconf->AddCorrectionOnQnVector(new QnCorrectionsQnVectorRecentering());
  /* lets configrure the QA histograms */
  VZEROAconf->SetQACentralityVar(varForEventMultiplicity);
  VZEROAconf->SetQAMultiplicityAxis(100, 0.0, 500.0);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROAconf);

  /* the VZEROC detector configuration */
  QnCorrectionsDetectorConfigurationChannels *VZEROCconf =
      new QnCorrectionsDetectorConfigurationChannels(
          "VZEROC",
          CorrEventClasses,
          64, /* number of channels */
          4); /* number of harmonics: 1, 2, 3 and 4 */
  VZEROCconf->SetChannelsScheme(VZEROchannels[1], channelGroups);
  /* let's configure the Q vector calibration */
  VZEROCconf->SetQVectorNormalizationMethod(QVNORM_QoverM);
  /* lets configure the equalization of input data */
  QnCorrectionsInputGainEqualization *eqC = new QnCorrectionsInputGainEqualization();
  eqC->SetEqualizationMethod(GEQUAL_averageEqualization);
  eqC->SetShift(1.0);
  eqC->SetScale(0.1);
  eqC->SetUseChannelGroupsWeights(kTRUE);
  VZEROCconf->AddCorrectionOnInputData(eqC);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  VZEROCconf->AddCorrectionOnQnVector(new QnCorrectionsQnVectorRecentering());
  /* lets configrure the QA histograms */
  VZEROCconf->SetQACentralityVar(varForEventMultiplicity);
  VZEROCconf->SetQAMultiplicityAxis(100, 0.0, 500.0);

  /* add the configuration to the detector */
  VZERO->AddDetectorConfiguration(VZEROCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(VZERO);
}

void AddTPC(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager){

  /////////////// Add TPC subdetectors ///////////////////

  //-----------------------------------------------------------
  // Our event classes for TPC
  //
  const Int_t nTPCdim = 2;
  QnCorrectionsEventClassVariablesSet *CorrEventClasses = new QnCorrectionsEventClassVariablesSet(nTPCdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(VAR::kVtxZ,
      task->VarName(VAR::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the TPC  detector */
  QnCorrectionsDetector *TPC = new QnCorrectionsDetector("TPC", VAR::kTPC);

  /* the TPC detector configuration */
  QnCorrectionsDetectorConfigurationTracks *TPCconf =
      new QnCorrectionsDetectorConfigurationTracks(
          "TPC",
          CorrEventClasses,
          4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  TPCconf->SetQVectorNormalizationMethod(QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TPCconf->AddCorrectionOnQnVector(new QnCorrectionsQnVectorRecentering());

  /* define the cuts to apply */
  QnCorrectionsCutsSet *cutsTPC = new QnCorrectionsCutsSet();
  cutsTPC->Add(new QnCorrectionsCutWithin(VAR::kDcaXY,-0.3,0.3));
  cutsTPC->Add(new QnCorrectionsCutWithin(VAR::kDcaZ,-0.3,0.3));
  cutsTPC->Add(new QnCorrectionsCutWithin(VAR::kEta,-0.8,0.8));
  cutsTPC->Add(new QnCorrectionsCutWithin(VAR::kPt,0.2,5.));
  cutsTPC->Add(new QnCorrectionsCutWithin(VAR::kTPCncls,70.0,161.0));
  cutsTPC->Add(new QnCorrectionsCutWithin(VAR::kTPCchi2,0.2,4.0));
  TPCconf->SetCuts(cutsTPC);

  /* add the configuration to the detector */
  TPC->AddDetectorConfiguration(TPCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(TPC);
}


void AddSPD(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager){

  /////////////// Add SPD subdetectors ///////////////////

  //-----------------------------------------------------------
  // Our event classes for SPD
  //
  const Int_t nSPDdim = 2;
  QnCorrectionsEventClassVariablesSet *CorrEventClasses = new QnCorrectionsEventClassVariablesSet(nSPDdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(VAR::kVtxZ,
      task->VarName(VAR::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the SPD detector */
  QnCorrectionsDetector *SPD = new QnCorrectionsDetector("SPD", VAR::kSPD);

  /* the SPD detector configuration */
  QnCorrectionsDetectorConfigurationTracks *SPDconf =
      new QnCorrectionsDetectorConfigurationTracks(
          "SPD",
          CorrEventClasses,
          4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  SPDconf->SetQVectorNormalizationMethod(QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  SPDconf->AddCorrectionOnQnVector(new QnCorrectionsQnVectorRecentering());

  /* add the configuration to the detector */
  SPD->AddDetectorConfiguration(SPDconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(SPD);
}

void AddTZERO(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager){

  /////////////// Add TZERO subdetectors ///////////////////

  Bool_t TZEROchannels[2][24];
  for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<23; ich++) TZEROchannels[iv0][ich] = kFALSE;

  for(Int_t ich=12; ich<24; ich++) TZEROchannels[0][ich] = kTRUE;  // channel list: value 1 if channel should be used
  for(Int_t ich=0; ich<12; ich++) TZEROchannels[1][ich] = kTRUE;

  Int_t channelGroups[24];
  for(Int_t ich=0; ich<24; ich++) channelGroups[ich] = Int_t(ich / 12);

  //-----------------------------------------------------------
  // Our event classes for TZERO
  //
  const Int_t nTZEROdim = 2;
  QnCorrectionsEventClassVariablesSet *CorrEventClasses = new QnCorrectionsEventClassVariablesSet(nTZEROdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(VAR::kVtxZ,
      task->VarName(VAR::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the TZERO detector */
  QnCorrectionsDetector *TZERO = new QnCorrectionsDetector("TZERO", VAR::kTZERO);

  /* the TZEROA detector configuration */
  QnCorrectionsDetectorConfigurationChannels *TZEROAconf =
      new QnCorrectionsDetectorConfigurationChannels(
          "TZEROA",
          CorrEventClasses,
          24, /* number of channels */
          4); /* number of harmonics: 1, 2, 3 and 4 */
  TZEROAconf->SetChannelsScheme(TZEROchannels[0], channelGroups);
  /* let's configure the Q vector calibration */
  TZEROAconf->SetQVectorNormalizationMethod(QVNORM_QoverM);
  /* lets configure the equalization of input data */
  QnCorrectionsInputGainEqualization *eqA = new QnCorrectionsInputGainEqualization();
  eqA->SetEqualizationMethod(GEQUAL_averageEqualization);
  eqA->SetShift(1.0);
  eqA->SetScale(0.1);
  eqA->SetUseChannelGroupsWeights(kFALSE);
  TZEROAconf->AddCorrectionOnInputData(eqA);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TZEROAconf->AddCorrectionOnQnVector(new QnCorrectionsQnVectorRecentering());
  /* let's configure the QA histograms */
  TZEROAconf->SetQACentralityVar(varForEventMultiplicity);
  TZEROAconf->SetQAMultiplicityAxis(100, 0.0, 150.0);

  /* add the configuration to the detector */
  TZERO->AddDetectorConfiguration(TZEROAconf);

  /* the TZEROC detector configuration */
  QnCorrectionsDetectorConfigurationChannels *TZEROCconf =
      new QnCorrectionsDetectorConfigurationChannels(
          "TZEROC",
          CorrEventClasses,
          24, /* number of channels */
          4); /* number of harmonics: 1, 2, 3 and 4 */
  TZEROCconf->SetChannelsScheme(TZEROchannels[1], channelGroups);
  /* let's configure the Q vector calibration */
  TZEROCconf->SetQVectorNormalizationMethod(QVNORM_QoverM);
  /* lets configure the equalization of input data */
  QnCorrectionsInputGainEqualization *eqC = new QnCorrectionsInputGainEqualization();
  eqC->SetEqualizationMethod(GEQUAL_averageEqualization);
  eqC->SetShift(1.0);
  eqC->SetScale(0.1);
  eqC->SetUseChannelGroupsWeights(kFALSE);
  TZEROCconf->AddCorrectionOnInputData(eqC);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  TZEROCconf->AddCorrectionOnQnVector(new QnCorrectionsQnVectorRecentering());
  /* let's configure the QA histograms */
  TZEROCconf->SetQACentralityVar(varForEventMultiplicity);
  TZEROCconf->SetQAMultiplicityAxis(100, 0.0, 150.0);

  /* add the configuration to the detector */
  TZERO->AddDetectorConfiguration(TZEROCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(TZERO);
}


void AddZDC(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager){
  /////////////// Add ZDC subdetectors ///////////////////

  Bool_t ZDCchannels[2][10];
  for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<10; ich++) ZDCchannels[iv0][ich] = kFALSE;

  for(Int_t ich=6; ich<10; ich++) ZDCchannels[0][ich] = kTRUE;
  for(Int_t ich=1; ich<5; ich++)  ZDCchannels[1][ich] = kTRUE;

  //-----------------------------------------------------------
  // Our event classes for ZDC
  //
  const Int_t nZDCdim = 3;
  QnCorrectionsEventClassVariablesSet *CorrEventClasses = new QnCorrectionsEventClassVariablesSet(nZDCdim);
  Double_t VtxXbinning[][2] = {{ -0.3, 2}, {0.3, 10 }};
  Double_t VtxYbinning[][2] = {{ -0.3, 2}, {0.3, 10 }};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(VAR::kVtxX,
      task->VarName(VAR::kVtxX), VtxXbinning));
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(VAR::kVtxY,
      task->VarName(VAR::kVtxY), VtxYbinning));
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the ZDC detector */
  QnCorrectionsDetector *ZDC = new QnCorrectionsDetector("ZDC", VAR::kZDC);

  /* the ZDCA detector configuration */
  QnCorrectionsDetectorConfigurationChannels *ZDCAconf =
      new QnCorrectionsDetectorConfigurationChannels(
          "ZDCA",
          CorrEventClasses,
          10, /* number of channels */
          3); /* number of harmonics: 1, 2 and 3 */
  ZDCAconf->SetChannelsScheme(ZDCchannels[0], NULL /* no groups */);
  /* let's configure the Q vector calibration */
  ZDCAconf->SetQVectorNormalizationMethod(QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  ZDCAconf->AddCorrectionOnQnVector(new QnCorrectionsQnVectorRecentering());

  /* add the configuration to the detector */
  ZDC->AddDetectorConfiguration(ZDCAconf);

  /* the ZDCC detector configuration */
  QnCorrectionsDetectorConfigurationChannels *ZDCCconf =
      new QnCorrectionsDetectorConfigurationChannels(
          "ZDCC",
          CorrEventClasses,
          10, /* number of channels */
          3); /* number of harmonics: 1, 2 and 3 */
  ZDCCconf->SetChannelsScheme(ZDCchannels[1], NULL /* no groups */);
  /* let's configure the Q vector calibration */
  ZDCCconf->SetQVectorNormalizationMethod(QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  ZDCCconf->AddCorrectionOnQnVector(new QnCorrectionsQnVectorRecentering());

  /* add the configuration to the detector */
  ZDC->AddDetectorConfiguration(ZDCCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(ZDC);
}


void AddFMD(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  gSystem->Load("libPWGLFforward2");  // for FMD

  // Create the FMD task and add it to the manager
  //===========================================================================
  //--- AOD output handler -----------------------------------------
  AliAODHandler* ret = new AliAODHandler();

  ret->SetOutputFileName("AliAOD.pass2.root");
  mgr->SetOutputEventHandler(ret);

  gSystem->Load("libESDfilter.so");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C");
  AliAnalysisTaskESDfilter *esdfilter = AddTaskESDFilter(kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (esdfilter,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (esdfilter,  0, mgr->GetCommonOutputContainer());

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/AddTaskForwardMult.C");

  ULong_t run = 0; // 0: get from data???
  UShort_t sys = 0; // 0: get from data, 1: pp, 2: AA
  UShort_t sNN = 0; // 0: get from data, otherwise center of mass energy (per nucleon pair)
  Short_t  fld = 0; // 0: get from data, otherwise L3 field in kG

  const Char_t* config = "$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/ForwardAODConfig2.C";
  AliAnalysisTask *taskFmd  = AddTaskForwardMult(bMC, run, sys, sNN, fld, config);

  // --- Make the output container and connect it --------------------
  AliAnalysisDataContainer* histOut =
    mgr->CreateContainer("Forward", TList::Class(),
        AliAnalysisManager::kExchangeContainer,
        "Forward");

  AliAnalysisDataContainer *output =
    mgr->CreateContainer("ForwardResultsP", TList::Class(),
        AliAnalysisManager::kParamContainer,
        "ForwardResultsP");

  mgr->ConnectInput(taskFmd, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskFmd, 1, histOut);

  Bool_t FMDchannels[2][4000];
  for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<4000; ich++) FMDchannels[iv0][ich] = kFALSE;

  for(Int_t ich=2000; ich<4000; ich++) FMDchannels[0][ich] = kTRUE;  // channel list: value 1 if channel should be used
  for(Int_t ich=0; ich<2000; ich++) FMDchannels[1][ich] = kTRUE;

  //-----------------------------------------------------------
  // Our event classes for FMD
  //
  const Int_t nFMDdim = 2;
  QnCorrectionsEventClassVariablesSet *CorrEventClasses = new QnCorrectionsEventClassVariablesSet(nFMDdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(VAR::kVtxZ,
      task->VarName(VAR::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  /* the FMD detector */
  QnCorrectionsDetector *FMD = new QnCorrectionsDetector("FMD", VAR::kFMD);

  /* the FMDA detector configuration */
  QnCorrectionsDetectorConfigurationChannels *FMDAconf =
      new QnCorrectionsDetectorConfigurationChannels(
          "FMDA",
          CorrEventClasses,
          4000, /* number of channels */
          4); /* number of harmonics: 1, 2 and 3 */
  FMDAconf->SetChannelsScheme(FMDchannels[0], NULL /* no groups */);
  /* let's configure the Q vector calibration */
  FMDAconf->SetQVectorNormalizationMethod(QVNORM_QoverM);

  /* add the configuration to the detector */
  FMD->AddDetectorConfiguration(FMDAconf);

  /* the FMDC detector configuration */
  QnCorrectionsDetectorConfigurationChannels *FMDCconf =
      new QnCorrectionsDetectorConfigurationChannels(
          "FMDC",
          CorrEventClasses,
          4000, /* number of channels */
          4); /* number of harmonics: 1, 2, 3 and 4 */
  FMDCconf->SetChannelsScheme(FMDchannels[1], NULL /* no groups */);
  /* let's configure the Q vector calibration */
  FMDCconf->SetQVectorNormalizationMethod(QVNORM_QoverM);

  /* add the configuration to the detector */
  FMD->AddDetectorConfiguration(FMDCconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(FMD);
}


void AddRawFMD(AnalysisTaskFlowVectorCorrections *task, QnCorrectionsManager* QnManager){

  //-----------------------------------------------------------
  // Our event classes for FMD
  //
  const Int_t nFMDdim = 2;
  QnCorrectionsEventClassVariablesSet *CorrEventClasses = new QnCorrectionsEventClassVariablesSet(nFMDdim);
  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }};
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(VAR::kVtxZ,
      task->VarName(VAR::kVtxZ), VtxZbinning));
  CorrEventClasses->Add(new QnCorrectionsEventClassVariable(varForEventMultiplicity,
      Form("Centrality (%s)", task->VarName(varForEventMultiplicity)), Ctbinning));
  ////////// end of binning

  QnCorrectionsCutsSet *cutFMDA = new QnCorrectionsCutsSet();
  cutFMDA->Add(new QnCorrectionsCutWithin(VAR::kFMDEta,0.0,6.0));

  QnCorrectionsCutsSet *cutFMDC = new QnCorrectionsCutsSet();
  cutFMDC->Add(new QnCorrectionsCutWithin(VAR::kFMDEta,-6.0,0.0));

  /* the FMD detector */
  QnCorrectionsDetector *FMDraw = new QnCorrectionsDetector("FMDraw", VAR::kFMDraw);

  /* the FMDAraw detector configuration */
  QnCorrectionsDetectorConfigurationTracks *FMDArawconf =
      new QnCorrectionsDetectorConfigurationTracks(
          "FMDAraw",
          CorrEventClasses,
          4); /* number of harmonics: 1, 2 and 3 */
  /* let's configure the Q vector calibration */
  FMDArawconf->SetQVectorNormalizationMethod(QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  FMDArawconf->AddCorrectionOnQnVector(new QnCorrectionsQnVectorRecentering());
  /* and add the cuts */
  FMDArawconf->SetCuts(cutFMDA);

  /* add the configuration to the detector */
  FMDraw->AddDetectorConfiguration(FMDArawconf);

  /* the FMDCraw detector configuration */
  QnCorrectionsDetectorConfigurationTracks *FMDCrawconf =
      new QnCorrectionsDetectorConfigurationTracks(
          "FMDCraw",
          CorrEventClasses,
          4); /* number of harmonics: 1, 2, 3 and 4 */
  /* let's configure the Q vector calibration */
  FMDCrawconf->SetQVectorNormalizationMethod(QVNORM_QoverM);
  /* let's add the Q vector recentering correction step */
  /* we don't configure it, so we create it anonymous */
  FMDCrawconf->AddCorrectionOnQnVector(new QnCorrectionsQnVectorRecentering());
  /* and add the cuts */
  FMDCrawconf->SetCuts(cutFMDC);

  /* add the configuration to the detector */
  FMDraw->AddDetectorConfiguration(FMDCrawconf);

  /* finally add the detector to the framework manager */
  QnManager->AddDetector(FMDraw);
}

//__________________________________________________________________
void DefineHistograms(QnCorrectionsManager* QnManager, AliQnCorrectionsHistos* histos, TString histClass) {
  //
  // define the histograms
  //
  const Char_t* histClasses = histClass.Data();

  cout << "Defining histograms ..." << endl;
  cout << "histogram classes: " << histClass<< endl;

  //fHistosFile=new TFile(output,"RECREATE");

  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");

  const Int_t kNRunBins = 3000;
  Double_t runHistRange[2] = {137000.,140000.};

  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    cout << "hist class: " << classStr.Data() << endl;

    // Event wise histograms
    if(classStr.Contains("Event")) {
      histos->AddHistClass(classStr.Data());
      histos->AddHistogram(classStr.Data(),"RunNo","Run numbers;Run", kFALSE, kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo);
      histos->AddHistogram(classStr.Data(),"BC","Bunch crossing;BC", kFALSE,3000,0.,3000.,VAR::kBC);
      histos->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag;;", kFALSE,
          2,-0.5,1.5,VAR::kIsPhysicsSelection, 0,0.0,0.0,VAR::kNothing, 0,0.0,0.0,VAR::kNothing, "off;on");

      histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,-30.0,30.0,VAR::kVtxZ);
      //histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,-15.,15.,VAR::kVtxZ);
      histos->AddHistogram(classStr.Data(),"VtxX","Vtx X;vtx X (cm)", kFALSE,300,-1.,1.,VAR::kVtxX);
      histos->AddHistogram(classStr.Data(),"VtxY","Vtx Y;vtx Y (cm)", kFALSE,300,-1.,1.,VAR::kVtxY);


      histos->AddHistogram(classStr.Data(),"CentVZEROvsMultPVZERO","Multiplicity percentile (VZERO);multiplicity VZERO (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kVZEROMultPercentile, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentVZEROvsCentSPD","Centrality(VZERO);centrality VZERO (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentVZERO, 100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentSPD","Centrality(TPC);centrality TPC (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentTPC, 100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentVZERO","Centrality(TPC);centrality TPC (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentTPC, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentZDC","Centrality(TPC);centrality TPC (percents);centrality ZDC (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentTPC, 100, 0.0, 100.0, VAR::kCentZDC);
      histos->AddHistogram(classStr.Data(),"CentZDCvsCentVZERO","Centrality(ZDC);centrality ZDC (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentZDC, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentZDCvsCentSPD","Centrality(ZDC);centrality ZDC (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentZDC, 100, 0.0, 100.0, VAR::kCentSPD);

      histos->AddHistogram(classStr.Data(),"MultVZEROvsCentVZERO","Multiplicity;multiplicity VZERO;VZERO centrality", kFALSE,
          100, 0.0, 32000.0, VAR::kVZEROTotalMult, 100,0.,100., VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"MultSPDvsCentPSD","Multiplicity;SPD tracklets;SPD centrality", kFALSE,
          100, 0.0, 3000.0, VAR::kSPDntracklets, 100,0.,100., VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"MultTPCvsCentSPD","Multiplicity;TPC selected tracks;TPC centrality", kFALSE,
          100, 0.0, 3500.0, VAR::kNtracksSelected, 100,0.,100., VAR::kCentTPC);
      histos->AddHistogram(classStr.Data(),"MultZDCvsCentZDC","Multiplicity;multiplicity ZDC;ZDC centrality", kFALSE,
          100, 0.0, 300000.0, VAR::kZDCTotalEnergy, 100,0.,100., VAR::kCentZDC);


      histos->AddHistogram(classStr.Data(),"MultTPCvsMultVZERO","Multiplicity;tracks TPC;multiplicity VZERO", kFALSE,
          100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 32000.0, VAR::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(),"MultTPCvsMultSPD","Multiplicity;tracklets SPD;tracks TPC", kFALSE,
          100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 3000.0, VAR::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"MultSPDvsMultVZERO","Multiplicity;tracklets SPD;multiplicity VZERO", kFALSE,
          100, 0.0, 32000.0, VAR::kVZEROTotalMult, 100, 0.0, 3000.0, VAR::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"MultTPCvsMultZDC","Multiplicity;tracks TPC;energy ZDC", kFALSE,
          100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultVZEROvsMultZDC","Multiplicity;multiplicity VZERO;energy ZDC", kFALSE,
          100, 0.0, 32000.0, VAR::kVZEROTotalMult, 100, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultSPDvsMultZDC","Multiplicity;tracklets SPD;energy ZDC", kFALSE,
          100, 0.0, 3000.0, VAR::kSPDntracklets, 100, 0.0, 300000.0, VAR::kZDCTotalEnergy);



      histos->AddHistogram(classStr.Data(),"MultVZERO","Multiplicity;multiplicity VZERO", kFALSE,
          320, 0.0, 25000.0, VAR::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(),"MultVZEROA","Multiplicity;multiplicity VZEROA", kFALSE,
          250, 0.0, 9500.0, VAR::kVZEROATotalMult);//10000.0
      histos->AddHistogram(classStr.Data(),"MultVZEROC","Multiplicity;multiplicity VZEROC", kFALSE,
          250, 0.0, 16000.0, VAR::kVZEROCTotalMult);//15000.0
      histos->AddHistogram(classStr.Data(),"MultZDC","Multiplicity;multiplicity ZDC", kFALSE,
          200, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultZDCA","Multiplicity;multiplicity ZDCA", kFALSE,
          200, 0.0, 150000.0, VAR::kZDCATotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultZDCC","Multiplicity;multiplicity ZDCC", kFALSE,
          200, 0.0, 150000.0, VAR::kZDCCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultFMD1","Multiplicity;multiplicity FMD1", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD1TotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD2I","Multiplicity;multiplicity FMD2I", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD2ITotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD2O","Multiplicity;multiplicity FMD2O", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD2OTotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD3I","Multiplicity;multiplicity FMD3I", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD3ITotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD3O","Multiplicity;multiplicity FMD3O", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD3OTotalMult);
      histos->AddHistogram(classStr.Data(),"MultTZEROA","Multiplicity;multiplicity TZEROA", kFALSE,
          300, 0.0, 3000.0, VAR::kTZEROATotalMult);
      histos->AddHistogram(classStr.Data(),"MultTZEROC","Multiplicity;multiplicity TZEROC", kFALSE,
          300, 0.0, 3000.0, VAR::kTZEROCTotalMult);




      histos->AddHistogram(classStr.Data(),"MultPercentVZERO","Multiplicity percentile (VZERO);multiplicity VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kVZEROMultPercentile);
      histos->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC);centrality TPC (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentTPC);
      histos->AddHistogram(classStr.Data(),"CentZDC","Centrality(ZDC);centrality ZDC (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentZDC);

      histos->AddHistogram(classStr.Data(),"CentQuality","Centrality quality;centrality quality", kFALSE,
          100, -50.5, 49.5, VAR::kCentQuality);
      histos->AddHistogram(classStr.Data(),"CentVZERO_Run_prof","<Centrality(VZERO)> vs run;Run; centrality VZERO (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentSPD_Run_prof","<Centrality(SPD)> vs run;Run; centrality SPD (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPC_Run_prof","<Centrality(TPC)> vs run;Run; centrality TPC (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentTPC);
      histos->AddHistogram(classStr.Data(),"CentZDC_Run_prof","<Centrality(ZDC)> vs run;Run; centrality ZDC (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentZDC);


      histos->AddHistogram(classStr.Data(),"NV0sTotal","Number of V0 candidates per event;# pairs", kFALSE,
          1000,0.,30000.,VAR::kNV0total);
      histos->AddHistogram(classStr.Data(),"NV0sSelected","Number of selected V0 candidates per event;# pairs", kFALSE,
          1000,0.,10000.,VAR::kNV0selected);
      histos->AddHistogram(classStr.Data(),"NPairs","Number of candidates per event;# pairs", kFALSE,
          5000,0.,5000.,VAR::kNdielectrons);
      histos->AddHistogram(classStr.Data(),"NPairsSelected", "Number of selected pairs per event; #pairs", kFALSE,
          5000,0.,5000.,VAR::kNpairsSelected);
      histos->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event;# tracks", kFALSE,
          1000,0.,30000.,VAR::kNtracksTotal);
      histos->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event;# tracks", kFALSE,
          1000,0.,30000.,VAR::kNtracksSelected);
      histos->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets; tracklets", kFALSE,
          3000, -0.5, 2999.5, VAR::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"SPDnSingleClusters", "SPD #single clusters; tracklets", kFALSE,
          3000, -0.5, 2999.5, VAR::kSPDnSingleClusters);

      histos->AddHistogram(classStr.Data(),"NV0total_Run_prof", "<Number of total V0s> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNV0total);
      histos->AddHistogram(classStr.Data(),"NV0selected_Run_prof", "<Number of selected V0s> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNV0selected);
      histos->AddHistogram(classStr.Data(),"Ndielectrons_Run_prof", "<Number of dielectrons> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNdielectrons);
      histos->AddHistogram(classStr.Data(),"NpairsSelected_Run_prof", "<Number of selected pairs> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNpairsSelected);
      histos->AddHistogram(classStr.Data(),"NTracksTotal_Run_prof", "<Number of tracks> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNtracksTotal);
      histos->AddHistogram(classStr.Data(),"NTracksSelected_Run_prof", "<Number of selected tracks> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNtracksSelected);
      histos->AddHistogram(classStr.Data(),"SPDntracklets_Run_prof", "<SPD ntracklets> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kSPDntracklets);

      histos->AddHistogram(classStr.Data(),"VtxZ_CentVZERO","Centrality(VZERO) vs vtx. Z;vtx Z (cm); centrality VZERO (%)", kFALSE,
          300,-15.,15.,VAR::kVtxZ, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"VtxZ_CentSPD","Centrality(SPD) vs vtx. Z;vtx Z (cm); centrality SPD (%)", kFALSE,
          300,-15.,15.,VAR::kVtxZ, 100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"VtxZ_CentTPC","Centrality(TPC) vs vtx. Z;vtx Z (cm); centrality TPC (%)", kFALSE,
          300,-15.,15.,VAR::kVtxZ, 100, 0.0, 100.0, VAR::kCentTPC);
      continue;
    }  // end if className contains "Event"


    // Offline trigger histograms
    if(classStr.Contains("OfflineTriggers")) {
      histos->AddHistClass(classStr.Data());

      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += Form("%s",VAR::fOfflineTriggerNames[i]); triggerNames+=";";}

      histos->AddHistogram(classStr.Data(), "Triggers", "Offline triggers fired; ; ;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTrigger, 2, -0.5, 1.5, VAR::kOfflineTriggerFired, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data(), "off;on");
      histos->AddHistogram(classStr.Data(), "Triggers2", "Offline triggers fired; ; ;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 0, 0.0, 0.0, VAR::kNothing, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "Offline triggers fired vs centrality VZERO; ; centrality VZERO;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentVZERO, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentTPC_Triggers2", "Offline triggers fired vs centrality TPC; ; centrality TPC;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentTPC, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentSPD_Triggers2", "Offline triggers fired vs centrality SPD; ; centrality SPD;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentSPD, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentZDC_Triggers2", "Offline triggers fired vs centrality ZDC; ; centrality ZDC;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentZDC, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "VtxZ_Triggers2", "Offline triggers fired vs vtxZ; ; vtx Z (cm.);", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 200, -20.0, 20.0, VAR::kVtxZ, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      continue;
    }

    // Track histograms
    if(classStr.Contains("Tracks")) {
      histos->AddHistClass(classStr.Data());
      for(Int_t ih=0; ih<6; ++ih) {
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
            20, 0.0, 100.0, VAR::kCentVZERO, 500, -1., 1., VAR::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
            20, 0.0, 100.0, VAR::kCentVZERO, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
      }
    }


    // Track histograms
    if(classStr.Contains("TrackQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution; p_{T} (GeV/c^{2});", kFALSE,
          1000, 0.0, 50.0, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
          1000, -1.5, 1.5, VAR::kEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
          1000, 0.0, 6.3, VAR::kPhi);
      histos->AddHistogram(classStr.Data(), "DCAxy", "DCAxy; DCAxy (cm.)", kFALSE,
          1000, -10.0, 10.0, VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz", "DCAz; DCAz (cm.)", kFALSE,
          1000, -10.0, 10.0, VAR::kDcaZ);
      histos->AddHistogram(classStr.Data(), "TPCncls", "TPCncls; TPCncls", kFALSE,
          160, 0.0, 160.0, VAR::kTPCncls);

      // run dependence
      histos->AddHistogram(classStr.Data(), "Pt_Run", "<p_{T}> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, 0.0, 50.0, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Run", "<#eta> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, -1.5, 1.5, VAR::kEta);
      histos->AddHistogram(classStr.Data(), "Phi_Run", "<#varphi> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, 0.0, 6.3, VAR::kPhi);
      histos->AddHistogram(classStr.Data(), "DCAxy_Run", "<DCAxy> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, -10.0, 10.0, VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Run", "<DCAz> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, -10.0, 10.0, VAR::kDcaZ);

      // correlations between parameters
      histos->AddHistogram(classStr.Data(), "Eta_Pt_prof", "<p_{T}> vs #eta; #eta; p_{T} (GeV/c);", kTRUE,
          300, -1.5, +1.5, VAR::kEta, 100, 0.0, 10.0, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt", "p_{T} vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kFALSE,
          300, -0.01, 6.3, VAR::kPhi, 100, 0.0, 2.2, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt_prof", "<p_{T}> vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kTRUE,
          300, 0.0, 6.3, VAR::kPhi, 100, 0.0, 10.0, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
          200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi);
      histos->AddHistogram(classStr.Data(), "TPCncls_Eta_Phi_prof", "<TPC ncls> vs #varphi vs #eta; #eta; #varphi (rad.);TPC ncls", kTRUE,
          200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi, 10, 0.0, 200., VAR::kTPCncls);
      histos->AddHistogram(classStr.Data(), "DCAxy_Eta_Phi_prof", "<DCAxy> vs #varphi vs #eta; #eta; #varphi (rad.);DCAxy (cm)", kTRUE,
          200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi, 10, 0.0, 200., VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Eta_Phi_prof", "<DCAz> vs #varphi vs #eta; #eta; #varphi (rad.);DCAz (cm)", kTRUE,
          200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi, 10, 0.0, 200., VAR::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Pt_DCAxy", "DCAxy vs p_{T}; p_{T} (GeV/c); DCA_{xy} (cm)", kFALSE,
          100, 0.0, 10.0, VAR::kPt, 500, -2.0, 2.0, VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Pt_DCAz", "DCAz vs p_{T}; p_{T} (GeV/c); DCA_{z} (cm)", kFALSE,
          100, 0.0, 10.0, VAR::kPt, 500, -2.0, 2.0, VAR::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Eta_DCAxy", "DCAxy vs #eta; #eta; DCA_{xy} (cm)", kFALSE,
          100, -1.0, 1.0, VAR::kEta, 500, -2.0, 2.0, VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Eta_DCAz", "DCAz vs #eta; #eta; DCA_{z} (cm)", kFALSE,
          100, -1.0, 1.0, VAR::kEta, 500, -2.0, 2.0, VAR::kDcaZ);

      for(Int_t ih=0; ih<6; ++ih) {
        //histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, VAR::kCentVZERO, 500, -1., 1., VAR::kCosNPhi+ih);
        //histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, VAR::kCentVZERO, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_Pt_Eta",ih+1), Form("<cos%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <cos%d #varphi>", ih+1, ih+1), kTRUE,
            20, -1.0, 1.0, VAR::kEta, 30, 0.0, 3.0, VAR::kPt, 500, -1.0, 1.0, VAR::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_Pt_Eta",ih+1), Form("<sin%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <sin%d #varphi>", ih+1, ih+1), kTRUE,
            20, -1.0, 1.0, VAR::kEta, 30, 0.0, 3.0, VAR::kPt, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO_VtxZ",ih+1), Form("<cos%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
            30, -15.0, 15.0, VAR::kVtxZ, 20, 0.0, 100.0, VAR::kCentVZERO, 500, -1., 1., VAR::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO_VtxZ",ih+1), Form("<sin%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
            30, -15.0, 15.0, VAR::kVtxZ, 20, 0.0, 100.0, VAR::kCentVZERO, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
      }
    }

    // Tracklet histograms
    if(classStr.Contains("TrackletQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
          1000, -3.0, 3.0, VAR::kSPDtrackletEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
          300, -0.01, 6.3, VAR::kSPDtrackletPhi);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
          200, -3.0, +3.0, VAR::kSPDtrackletEta, 100, 0.0, 6.3, VAR::kSPDtrackletPhi);
    }

  }

  cout << " done" << endl;
}

