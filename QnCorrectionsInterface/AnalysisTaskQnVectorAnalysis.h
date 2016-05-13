/*
***********************************************************
  event plane corrections framework
  contact: jaap onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
  2014/12/10
  *********************************************************
*/

//#include "AliSysInfo.h"

#ifndef ANALYSISTASKQNANALYSIS_H
#define ANALYSISTASKQNANALYSIS_H

#include "TFile.h"
#include "TTree.h"
#include "AliAnalysisTaskSE.h"
#include "QnCorrectionsFillEventTask.h"

class AliAnalysis;
class QnCorrectionsManager;
class QnCorrectionsCutsSet;
class AliQnCorrectionsHistos;
class TList;
class TProfile;
class TGraphErrors;


//_________________________________________________________
class AnalysisTaskQnVectorAnalysis : public QnCorrectionsFillEventTask {

public:

  enum enumTrackDetectors{
    kTPC=0,
    kSPD,
    nTrackDetectors
  };

  enum enumEPdetectors{
    kVZEROA=0,
    kVZEROC,
    kTZEROA,
    kTZEROC,
    kFMDA,
    kFMDC,
    kRawFMDA,
    kRawFMDC,
    nEPDetectors
  };

  enum Constants{
    kNharmonics=4,
    kNresolutionComponents=3, /* valid for 3-(sub-event)detector method */
    kNxy=2
  };

  enum CorrelationConstants{
    kXX = 0,
    kXY,
    kYX,
    kYY,
    kNcorrelationComponents,
  };

  /* this is only valid if we are using 3-(sub-event)detector method
   * for evaluating the detector resolution
   */
  enum SubEventConstants {
    kABXX = 0,
    kABYY,
    kACXX,
    kACYY,
    kBCXX,
    kBCYY,
    nCorrelationPerDetector
  };

  AnalysisTaskQnVectorAnalysis();
  AnalysisTaskQnVectorAnalysis(const char *name);
  virtual ~AnalysisTaskQnVectorAnalysis();


  virtual void UserExec(Option_t *);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();

  AliQnCorrectionsHistos* GetHistograms() {return fEventPlaneHistos;}
  QnCorrectionsCutsSet* EventCuts()  const {return fEventCuts;}
  Bool_t IsEventSelected(Float_t* values);

  void SetEventCuts(QnCorrectionsCutsSet* cuts)  {fEventCuts = cuts;}
  void SetCentralityVariable(Int_t var) { fCentralityVariable = var; }

 private:
  TList* fEventQAList;
  QnCorrectionsCutsSet *fEventCuts;
  AliQnCorrectionsHistos* fEventPlaneHistos;

  AnalysisTaskQnVectorAnalysis(const AnalysisTaskQnVectorAnalysis &c);
  AnalysisTaskQnVectorAnalysis& operator= (const AnalysisTaskQnVectorAnalysis &c);

  TProfile* fVn[nTrackDetectors*nEPDetectors][kNharmonics][kNcorrelationComponents];

  Int_t fNDetectorResolutions;
  TGraphErrors* fDetectorResolution[nTrackDetectors*nEPDetectors][kNharmonics];
  Int_t fDetectorResolutionContributors[nTrackDetectors*nEPDetectors][kNresolutionComponents];
  TProfile *fDetectorResolutionCorrelations[nTrackDetectors*nEPDetectors][nCorrelationPerDetector][kNharmonics];

  TString fTrackDetectorNameInFile[nTrackDetectors];
  TString fEPDetectorNameInFile[nEPDetectors];

  Int_t fCentralityVariable;

  ClassDef(AnalysisTaskQnVectorAnalysis, 1);
};

#endif

