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

#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisTaskSE.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>

#include "QnCorrectionsManager.h"
#include "QnCorrectionsVarManagerTask.h"
#include "AliQnCorrectionsHistos.h"

class AliESDtrack;
class AliVParticle;

class QnCorrectionsFillEventTask : public QnCorrectionsVarManagerTask {
public:

  QnCorrectionsFillEventTask();
  QnCorrectionsFillEventTask(const char *name);
  ~QnCorrectionsFillEventTask();

public:

  virtual void UserExec(Option_t *) = 0;
  virtual void UserCreateOutputObjects() = 0;
  virtual void FinishTaskOutput() = 0;


  void SetUseTPCStandaloneTracks(Bool_t enable = kTRUE) { fUseTPCStandaloneTracks = enable; }

protected:
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

  QnCorrectionsFillEventTask(const QnCorrectionsFillEventTask &c);
  QnCorrectionsFillEventTask& operator= (const QnCorrectionsFillEventTask &c);

protected:
  AliVEvent* fEvent;
  QnCorrectionsManager *fQnCorrectionsManager;
  AliQnCorrectionsHistos* fEventHistos;
  Float_t *fDataBank;                             //!<! The event variables values data bank. Transient!
private:
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

  ClassDef(QnCorrectionsFillEventTask, 1);
};

#endif
