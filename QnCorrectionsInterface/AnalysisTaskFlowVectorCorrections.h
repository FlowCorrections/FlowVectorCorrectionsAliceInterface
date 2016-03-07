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

class AliAnalysis;
class QnCorrectionsManager;
class QnCorrectionsCutsSet;
class AliQnCorrectionsHistos;

class AnalysisTaskFlowVectorCorrections : public AliAnalysisTaskSE {
public:
  /* include the variable handling definitions */
#include "AnalysisTaskFlowVectorCorrectionsVar.h"

public:
  AnalysisTaskFlowVectorCorrections();
  AnalysisTaskFlowVectorCorrections(const char *name);
  virtual ~AnalysisTaskFlowVectorCorrections(){}


  virtual void UserExec(Option_t *);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();


  void SetUseTPCStandaloneTracks(Bool_t enable = kTRUE) { fUseTPCStandaloneTracks = enable; }
  void SetRunByRunCalibration(Bool_t enable) { fCalibrateByRun = enable; }
  void SetQnCorrectionsManager(QnCorrectionsManager* QnManager)  {fQnCorrectionsManager = QnManager;}
  void SetEventCuts(QnCorrectionsCutsSet *cuts)  {fEventCuts = cuts;}
  void SetProvideQnVectors(Bool_t enable = kTRUE) { fProvideQnVectorsList = enable; }
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
  Bool_t IsFillExchangeContainerWithQvectors() const  {return fProvideQnVectorsList;}
  Bool_t IsFillEventQA() const  {return fFillEventQA;}

private:
  /* Fill event data methods */
  void FillEventData();

  void FillDetectors();
  void FillTPC();
  void FillEsdTPC();
  void FillAodTPC();
  void FillVZERO();
  void FillTZERO();
  void FillZDC();
  void FillFMD();
  void FillRawFMD();
  void FillSPDTracklets();

  void FillEventInfo();
  void FillTrackInfo(AliESDtrack* p);
  void FillTrackInfo(AliVParticle* p);

  void SetDetectors();

private:
  /* variable handling methods */
  const Char_t* VarName(Int_t var) const;
  const Char_t* VarUnits(Int_t var) const;
  void SetDefaultVarNames();


private:
  Bool_t fCalibrateByRun;
  UInt_t fTriggerMask;
  TList* fListInputHistogramsQnCorrections;          //! List of input histograms for corrections
  TList* fEventQAList;
  QnCorrectionsManager *fQnCorrectionsManager;
  QnCorrectionsCutsSet *fEventCuts;
  AliQnCorrectionsHistos* fEventHistos;
  TString fLabel;
  TString fQAhistograms;
  Bool_t fFillEventQA;
  Bool_t fProvideQnVectorsList;
  Int_t fOutputSlotEventQA;
  Int_t fOutputSlotHistQA;
  Int_t fOutputSlotHistQn;
  Int_t fOutputSlotQnVectorsList;
  Int_t fOutputSlotTree;

  /* for filling event data */
  AliVEvent* fEvent;
  Float_t *fDataBank;
  Bool_t fUseTPCStandaloneTracks;
  Bool_t fFillVZERO;
  Bool_t fFillTPC;
  Bool_t fFillZDC;
  Bool_t fFillTZERO;
  Bool_t fFillFMD;
  Bool_t fFillRawFMD;
  Bool_t fFillSPD;
  Bool_t fIsAOD;
  Bool_t fIsESD;

  /* for variable handling */
  static const Char_t* fTrackingFlagNames[kNTrackingFlags];
  static const Char_t* fOfflineTriggerNames[64];
  Char_t* fVariableNames[kNVars][2];

  AnalysisTaskFlowVectorCorrections(const AnalysisTaskFlowVectorCorrections &c);
  AnalysisTaskFlowVectorCorrections& operator= (const AnalysisTaskFlowVectorCorrections &c);

  ClassDef(AnalysisTaskFlowVectorCorrections, 1);
};

inline const Char_t* AnalysisTaskFlowVectorCorrections::VarName(Int_t var) const {
  if ((!(var < 0)) && (var < kNVars))
    return fVariableNames[var][0];
  else
    return "";
}

inline const Char_t* AnalysisTaskFlowVectorCorrections::VarUnits(Int_t var) const {
  if ((!(var < 0)) && (var < kNVars))
    return fVariableNames[var][1];
  else
    return "";
}

#endif // ANALYSISTASKFLOWVECTORCORRECTION_H


