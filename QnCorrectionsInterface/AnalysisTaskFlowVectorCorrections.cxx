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
#include "QnCorrectionsFillEvent.h"
#include "QnCorrectionsCuts.h"
#include "QnCorrectionsVarManager.h"
#include "QnCorrectionsManager.h"
#include "AliQnCorrectionsHistos.h"
#include "AliLog.h"

#include "AnalysisTaskFlowVectorCorrections.h"

ClassImp(AnalysisTaskFlowVectorCorrections)


AnalysisTaskFlowVectorCorrections::AnalysisTaskFlowVectorCorrections() :
AliAnalysisTaskSE(),
fCalibrateByRun(kTRUE),
fTriggerMask(0),
fListInputHistogramsQnCorrections(0x0),
fEventQAList(0x0),
fQnCorrectionsManager(NULL),
fEventCuts(NULL),
fEventFiller(NULL),
fEventHistos(NULL),
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
            AliAnalysisTaskSE(),
            fCalibrateByRun(kTRUE),
            fTriggerMask(0),
            fListInputHistogramsQnCorrections(0x0),
            fEventQAList(0x0),
            fQnCorrectionsManager(NULL),
            fEventCuts(NULL),
            fEventFiller(NULL),
            fEventHistos(NULL),
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

  fEventFiller = new QnCorrectionsFillEvent();

  fEventQAList = new TList();
  fEventQAList->SetName("EventQA");
  fEventQAList->SetOwner(kTRUE);

  fQnCorrectionsQAHistos = new AliQnCorrectionsHistos();
  fEventFiller->SetQnCorrectionsQAHistos(fQnCorrectionsQAHistos);
}

//_________________________________________________________________________________
void AnalysisTaskFlowVectorCorrections::DefineInOutput(){

  if(!fQnCorrectionsManager) {
    AliFatal("First configure EventPlaneManager!!\n");
    return;
  }

  DefineInput(0,TChain::Class());
  Int_t outputSlot = 1;
  // Calibration histograms
  if (fQnCorrectionsManager->ShouldFillHistogramsQnCorrections()) {
    DefineOutput(outputSlot, TList::Class());
    fOutputSlotHistQn = outputSlot++;
  }
  // Calibrated qvector tree
  if (fQnCorrectionsManager->ShouldFillTreeQnVectors()) {
    DefineOutput(outputSlot, TTree::Class());
    fOutputSlotTree = outputSlot++;
  }
  // Qvector QA histograms
  if (fQnCorrectionsManager->ShouldFillHistogramsQA()) {
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

  fEventFiller->SetQnCorrectionsManager(fQnCorrectionsManager);

  fQnCorrectionsManager->Initialize();

  if (fQnCorrectionsManager->ShouldFillHistogramsQnCorrections())
    PostData(fOutputSlotHistQn, fQnCorrectionsManager->GetListOutputHistogramsQnCorrections());
  if (fQnCorrectionsManager->ShouldFillTreeQnVectors())
    PostData(fOutputSlotTree, fQnCorrectionsManager->GetTreeQnVectors());
  if (fQnCorrectionsManager->ShouldFillHistogramsQA())
    PostData(fOutputSlotHistQA, fQnCorrectionsManager->GetListHistogramsQA());
  if (fFillEventQA)
    PostData(fOutputSlotEventQA, fEventQAList);
}

void AnalysisTaskFlowVectorCorrections::UserExec(Option_t *){
  //
  // Main loop. Called for every event
  //

  AliVEvent* event = InputEvent();
  fQnCorrectionsManager->ClearEvent();

  Float_t* values = fQnCorrectionsManager->GetDataContainer();

  fEventFiller->Process((AliAnalysisTaskSE*) this, event, values);

  if (fCalibrateByRun) fQnCorrectionsManager->SetCurrentProcessListName(Form("%d",(Int_t) (values[QnCorrectionsVarManager::kRunNo])));

  fEventHistos->FillHistClass("Event_NoCuts", values);

  if (IsEventSelected(values)) {
    fEventHistos->FillHistClass("Event_Analysis", values);

    fQnCorrectionsManager->ProcessEvent();
  }  // end if event selection

  if(fProvideQnVectorsList)
    PostData(fOutputSlotQnVectorsList, fQnCorrectionsManager->GetListQnVectors());
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


