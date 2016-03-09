/*
 ***********************************************************
 Manager for event plane corrections framework
Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
Instructions in AddTask_EPcorrectionsExample.C
2014/12/10
 *********************************************************
 */

#include <Riostream.h>

#include <TROOT.h>
#include <TTimeStamp.h>
#include <TStopwatch.h>
#include <TChain.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliESDEvent.h>
#include "QnCorrectionsCuts.h"
#include "QnCorrectionsManager.h"
#include "AliQnCorrectionsHistos.h"
#include "AliLog.h"

#include "AnalysisTaskFlowVectorCorrections.h"

ClassImp(AnalysisTaskFlowVectorCorrections)


AnalysisTaskFlowVectorCorrections::AnalysisTaskFlowVectorCorrections() :
QnCorrectionsFillEventTask(),
fCalibrateByRun(kTRUE),
fTriggerMask(0),
fListInputHistogramsQnCorrections(0x0),
fEventQAList(0x0),
fEventCuts(NULL),
fLabel(""),
fQAhistograms(""),
fFillEventQA(kTRUE),
fProvideQnVectorsList(kTRUE),
fOutputSlotEventQA(-1),
fOutputSlotHistQA(-1),
fOutputSlotHistQn(-1),
fOutputSlotQnVectorsList(-1),
fOutputSlotTree(-1)
{
  //
  // Default constructor
  //
}

//_________________________________________________________________________________
AnalysisTaskFlowVectorCorrections::AnalysisTaskFlowVectorCorrections(const char* name) :
    QnCorrectionsFillEventTask(name),
fCalibrateByRun(kTRUE),
fTriggerMask(0),
fListInputHistogramsQnCorrections(0x0),
fEventQAList(0x0),
fEventCuts(NULL),
fLabel(""),
fQAhistograms(""),
fFillEventQA(kTRUE),
fProvideQnVectorsList(kTRUE),
fOutputSlotEventQA(-1),
fOutputSlotHistQA(-1),
fOutputSlotHistQn(-1),
fOutputSlotQnVectorsList(-1),
fOutputSlotTree(-1)
{
  //
  // Constructor
  //

  fEventQAList = new TList();
  fEventQAList->SetName("EventQA");
  fEventQAList->SetOwner(kTRUE);

  fEventHistos = new AliQnCorrectionsHistos();
}

//_________________________________________________________________________________
void AnalysisTaskFlowVectorCorrections::DefineInOutput(){

  if(!fQnCorrectionsManager) {
    AliFatal("First configure QnCorrecionsManager!!\n");
    return;
  }

  DefineInput(0,TChain::Class());
  Int_t outputSlot = 1;
  // Calibration histograms
  if (fQnCorrectionsManager->GetShouldFillOutputHistograms()) {
    DefineOutput(outputSlot, TList::Class());
    fOutputSlotHistQn = outputSlot++;
  }
  // Calibrated qvector tree
  if (fQnCorrectionsManager->GetShouldFillQnVectorTree()) {
    DefineOutput(outputSlot, TTree::Class());
    fOutputSlotTree = outputSlot++;
  }
  // Qvector QA histograms
  if (fQnCorrectionsManager->GetShouldFillQAHistograms()) {
    DefineOutput(outputSlot, TList::Class());
    fOutputSlotHistQA = outputSlot++;
  }
  // Calibrated qvector list
  if (fProvideQnVectorsList) {
    DefineOutput(outputSlot, TList::Class());
    fOutputSlotQnVectorsList=outputSlot++;
  }
  // Event QA histograms
  if (fFillEventQA) {
    DefineOutput(outputSlot, TList::Class());
    fOutputSlotEventQA=outputSlot++;
  }
}


//_________________________________________________________________________________
void AnalysisTaskFlowVectorCorrections::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //
  fQnCorrectionsManager->InitializeQnCorrectionsFramework();

  if (fQnCorrectionsManager->GetShouldFillOutputHistograms())
    PostData(fOutputSlotHistQn, fQnCorrectionsManager->GetOutputHistogramsList());
  if (fQnCorrectionsManager->GetShouldFillQnVectorTree())
    PostData(fOutputSlotTree, fQnCorrectionsManager->GetQnVectorTree());
  if (fQnCorrectionsManager->GetShouldFillQAHistograms())
    PostData(fOutputSlotHistQA, fQnCorrectionsManager->GetQAHistogramsList());
  if (fFillEventQA)
    PostData(fOutputSlotEventQA, fEventQAList);
}

/// The current run has changed. Usually it is only sent before
/// the first event is handled.
/// Notify the framework manager that the current label has changed.
void AnalysisTaskFlowVectorCorrections::NotifyRun() {

  if (fCalibrateByRun) fQnCorrectionsManager->SetCurrentProcessListName(Form("%d", this->fCurrentRunNumber));
}

void AnalysisTaskFlowVectorCorrections::UserExec(Option_t *){
  //
  // Main loop. Called for every event
  //

  fEvent = InputEvent();
  fQnCorrectionsManager->ClearEvent();

  fDataBank = fQnCorrectionsManager->GetDataContainer();

  FillEventData();

  fEventHistos->FillHistClass("Event_NoCuts", fDataBank);

  if (IsEventSelected(fDataBank)) {
    fEventHistos->FillHistClass("Event_Analysis", fDataBank);

    fQnCorrectionsManager->ProcessEvent();
  }  // end if event selection

  if(fProvideQnVectorsList)
    PostData(fOutputSlotQnVectorsList, fQnCorrectionsManager->GetQnVectorList());
}  // end loop over events


void AnalysisTaskFlowVectorCorrections::FinishTaskOutput()
{
  //
  // Finish Task
  //
  fQnCorrectionsManager->FinalizeQnCorrectionsFramework();

  THashList* hList = (THashList*) fEventHistos->HistList();
  for(Int_t i=0; i<hList->GetEntries(); ++i) {
    THashList* list = (THashList*)hList->At(i);
    fEventQAList->Add(list);
  }
}

Bool_t AnalysisTaskFlowVectorCorrections::IsEventSelected(Float_t* values) {

  if(!fEventCuts) return kTRUE;
  return fEventCuts->IsSelected(values);
}


