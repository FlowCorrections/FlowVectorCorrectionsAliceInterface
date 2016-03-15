#ifndef ANALYSISTASKFLOWVECTORCORRECTION_H
#define ANALYSISTASKFLOWVECTORCORRECTION_H

/***************************************************************************
 * Package:       FlowVectorCorrections ALICE glue                         *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2014-2016                                                *
 ***************************************************************************/

#include <TObject.h>
#include "Rtypes.h"

#include "TFile.h"
#include "TTree.h"

#include "AliAnalysisTaskSE.h"
#include "QnCorrectionsFillEventTask.h"

class AliAnalysis;
class QnCorrectionsManager;
class QnCorrectionsCutsSet;
class AliQnCorrectionsHistos;

class AnalysisTaskFlowVectorCorrections : public QnCorrectionsFillEventTask {

public:
  AnalysisTaskFlowVectorCorrections();
  AnalysisTaskFlowVectorCorrections(const char *name);
  virtual ~AnalysisTaskFlowVectorCorrections(){}


  virtual void UserExec(Option_t *);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  virtual void NotifyRun();


  void SetRunByRunCalibration(Bool_t enable) { fCalibrateByRun = enable; }
  void SetQnCorrectionsManager(QnCorrectionsManager* QnManager)  {fQnCorrectionsManager = QnManager;}
  void SetEventCuts(QnCorrectionsCutsSet *cuts)  {fEventCuts = cuts;}
  void SetFillExchangeContainerWithQvectors(Bool_t enable = kTRUE) { fProvideQnVectorsList = enable; }
  void SetFillEventQA(Bool_t enable = kTRUE) { fFillEventQA = enable; }
  void SetTrigger(UInt_t triggerbit) {fTriggerMask=triggerbit;}
  void AddHistogramClass(TString hist) {fQAhistograms+=hist+";";}
  void SetCalibrationHistograms(TList* input) {fListInputHistogramsQnCorrections = input;}
  void DefineInOutput();

  QnCorrectionsManager *GetQnCorrectionsManager() {return fQnCorrectionsManager;}
  AliQnCorrectionsHistos* GetEventHistograms() {return fEventHistos;}
  QnCorrectionsCutsSet* GetEventCuts()  const {return fEventCuts;}
  Int_t OutputSlotEventQA()        const {return fOutputSlotEventQA;}
  Int_t OutputSlotHistQA()        const {return fOutputSlotHistQA;}
  Int_t OutputSlotHistQn()        const {return fOutputSlotHistQn;}
  Int_t OutputSlotGetListQnVectors() const {return fOutputSlotQnVectorsList;}
  Int_t OutputSlotTree()          const {return fOutputSlotTree;}
  Bool_t IsEventSelected(Float_t* values);
  Bool_t GetFillExchangeContainerWithQvectors() const  {return fProvideQnVectorsList;}
  Bool_t GetFillEventQA() const  {return fFillEventQA;}

private:
  Bool_t fCalibrateByRun;
  UInt_t fTriggerMask;
  TList* fListInputHistogramsQnCorrections;          //! List of input histograms for corrections
  TList* fEventQAList;
  QnCorrectionsCutsSet *fEventCuts;
  TString fLabel;
  TString fQAhistograms;
  Bool_t fFillEventQA;
  Bool_t fProvideQnVectorsList;
  Int_t fOutputSlotEventQA;
  Int_t fOutputSlotHistQA;
  Int_t fOutputSlotHistQn;
  Int_t fOutputSlotQnVectorsList;
  Int_t fOutputSlotTree;

  AnalysisTaskFlowVectorCorrections(const AnalysisTaskFlowVectorCorrections &c);
  AnalysisTaskFlowVectorCorrections& operator= (const AnalysisTaskFlowVectorCorrections &c);

  ClassDef(AnalysisTaskFlowVectorCorrections, 1);
};

#endif // ANALYSISTASKFLOWVECTORCORRECTION_H


