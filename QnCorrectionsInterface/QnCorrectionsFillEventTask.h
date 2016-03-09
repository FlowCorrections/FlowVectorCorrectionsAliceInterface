#ifndef QNCORRECTIONS_FILLEVENT_H
#define QNCORRECTIONS_FILLEVENT_H

/***************************************************************************
 * Package:       FlowVectorCorrections ALICE glue                         *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2014-2016                                                *
 ***************************************************************************/

#include "QnCorrectionsManager.h"
#include "AliQnCorrectionsHistos.h"
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisTaskSE.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>

class QnCorrectionsFillEvent : public TNamed {
public:

  QnCorrectionsFillEvent();
  ~QnCorrectionsFillEvent();


  void Process(AliAnalysisTaskSE* task, AliVEvent* event, Float_t* values);

  void FillDetectors(AliAnalysisTaskSE* task, Float_t* values);
  void FillTPC(Float_t* values);
  void FillEsdTPC(Float_t* values);
  void FillAodTPC(Float_t* values);
  void FillVZERO();
  void FillTZERO();
  void FillZDC();
  void FillFMD(AliAnalysisTaskSE* task);
  void FillRawFMD(Float_t* values);
  void FillSPDTracklets(Float_t* values);

  void FillEventInfo(Float_t* values);
  void FillTrackInfo(AliESDtrack* p, Float_t* values);
  void FillTrackInfo(AliVParticle* p, Float_t* values);

  void SetUseTPCStandaloneTracks(Bool_t b=kTRUE) { fUseTPCStandaloneTracks=b; }

  void SetEvent(AliVEvent* ev) {fEvent=ev;}
  void SetQnCorrectionsManager(QnCorrectionsManager *QnManager) { fQnCorrectionsManager = QnManager; SetDetectors(); }
  void SetQnCorrectionsQAHistos(AliQnCorrectionsHistos* QAHistos) { fQAhistos = QAHistos; }
  void SetDetectors();


private:

  QnCorrectionsFillEvent(const QnCorrectionsFillEvent &c);
  QnCorrectionsFillEvent& operator= (const QnCorrectionsFillEvent &c);

  AliVEvent* fEvent;
  QnCorrectionsManager *fQnCorrectionsManager;
  AliQnCorrectionsHistos* fQAhistos;

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

  ClassDef(QnCorrectionsFillEvent, 1);
};

#endif
