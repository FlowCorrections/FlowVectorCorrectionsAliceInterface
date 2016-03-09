/**************************************************************************************************
 *                                                                                                *
 * Package:       FlowVectorCorrections                                                           *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch                              *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com                             *
 *                Víctor González, UCM, victor.gonzalez@cern.ch                                   *
 *                Contributors are mentioned in the code where appropriate.                       *
 * Development:   2012-2016                                                                       *
 *                                                                                                *
 * This file is part of FlowVectorCorrections, a software package that corrects Q-vector          *
 * measurements for effects of nonuniform detector acceptance. The corrections in this package    *
 * are based on publication:                                                                      *
 *                                                                                                *
 *  [1] "Effects of non-uniform acceptance in anisotropic flow measurements"                      *
 *  Ilya Selyuzhenkov and Sergei Voloshin                                                         *
 *  Phys. Rev. C 77, 034904 (2008)                                                                *
 *                                                                                                *
 * The procedure proposed in [1] is extended with the following steps:                            *
 * (*) alignment correction between subevents                                                     *
 * (*) possibility to extract the twist and rescaling corrections                                 *
 *      for the case of three detector subevents                                                  *
 *      (currently limited to the case of two “hit-only” and one “tracking” detectors)            *
 * (*) (optional) channel equalization                                                            *
 * (*) flow vector width equalization                                                             *
 *                                                                                                *
 * FlowVectorCorrections is distributed under the terms of the GNU General Public License (GPL)   *
 * (https://en.wikipedia.org/wiki/GNU_General_Public_License)                                     *
 * either version 3 of the License, or (at your option) any later version.                        *
 *                                                                                                *
 **************************************************************************************************/
/***********************************************************
 Variable fill for Flow Qn vector corrections framework
 Based on work of Ionut-Cristian Arsene
 ***********************************************************/

#include <Riostream.h>

#include "QnCorrectionsFillEvent.h"
#include "QnCorrectionsVarManager.h"

#include "QnCorrectionsDataVector.h"
#include "QnCorrectionsDetector.h"
#include "QnCorrectionsManager.h"

#include "AliQnCorrectionsHistos.h"

#include <TChain.h>
#include <TH1D.h>
#include <TFile.h>

#include <AliCFContainer.h>
#include <AliInputEventHandler.h>
#include <AliAnalysisManager.h>
#include <AliCentrality.h>
#include <AliMultSelection.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisTaskSE.h>
#include <AliAODForwardMult.h>
#include <AliForwardUtil.h>

#include <AliVEvent.h>
#include <AliVZDC.h>
#include <AliVVZERO.h>
#include <AliVParticle.h>

#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliESDHeader.h>
#include <AliESDtrack.h>
#include <AliESDtrackCuts.h>
#include <AliESDFMD.h>

#include <AliAODInputHandler.h>
#include <AliAODEvent.h>
#include <AliAODHeader.h>
#include <AliAODTrack.h>

#include <AliLog.h>


ClassImp(QnCorrectionsFillEvent)

#ifdef QNCORRECTIONS_VARMANAGER_H
#define VAR QnCorrectionsVarManager
#endif


QnCorrectionsFillEvent::QnCorrectionsFillEvent() :
TNamed("AliQnCorrectionsFillEvent","Fill functions"),
fEvent(NULL),
fQnCorrectionsManager(NULL),
fQAhistos(NULL),
fUseTPCStandaloneTracks(kFALSE),
fFillVZERO(kFALSE),
fFillTPC(kFALSE),
fFillZDC(kFALSE),
fFillTZERO(kFALSE),
fFillFMD(kFALSE),
fFillRawFMD(kFALSE),
fFillSPD(kFALSE),
fIsAOD(kFALSE),
fIsESD(kFALSE)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
QnCorrectionsFillEvent::~QnCorrectionsFillEvent()
{
  //
  // Destructor
  //
}



//__________________________________________________________________
void QnCorrectionsFillEvent::SetDetectors() {
  //
  // determine which detectors are used (to call only the necessary detector fill functions)

  if(fQnCorrectionsManager->FindDetector(VAR::kTPC  ) != NULL)    fFillTPC = kTRUE;
  if(fQnCorrectionsManager->FindDetector(VAR::kVZERO) != NULL)  fFillVZERO = kTRUE;
  if(fQnCorrectionsManager->FindDetector(VAR::kTZERO) != NULL)  fFillTZERO = kTRUE;
  if(fQnCorrectionsManager->FindDetector(VAR::kZDC  ) != NULL)    fFillZDC = kTRUE;
  if(fQnCorrectionsManager->FindDetector(VAR::kFMD  ) != NULL)    fFillFMD = kTRUE;
  if(fQnCorrectionsManager->FindDetector(VAR::kFMDraw) != NULL)fFillRawFMD = kTRUE;
  if(fQnCorrectionsManager->FindDetector(VAR::kSPD  ) != NULL)    fFillSPD = kTRUE;


}


//__________________________________________________________________
void QnCorrectionsFillEvent::Process(AliAnalysisTaskSE* task, AliVEvent* event, Float_t* values) {

  fEvent=event;

  TString aod = "AOD";
  TString esd = "ESD";

  fIsAOD = ( aod.EqualTo(fEvent->Whoami()) ? kTRUE : kFALSE );
  fIsESD = ( esd.EqualTo(fEvent->Whoami()) ? kTRUE : kFALSE );

  FillEventInfo(values);
  FillDetectors(task, values);

}

//__________________________________________________________________
void QnCorrectionsFillEvent::FillEventInfo(Float_t* values) {
  //
  // fill event info
  //


  values[VAR::kRunNo]       = fEvent->GetRunNumber();
  values[VAR::kVtxX]        = -999.;
  values[VAR::kVtxY]        = -999.;
  values[VAR::kVtxZ]        = -999.;
  const AliVVertex *primVtx = fEvent->GetPrimaryVertex();
  if (primVtx){
    values[VAR::kVtxX]        = primVtx->GetX();
    values[VAR::kVtxY]        = primVtx->GetY();
    values[VAR::kVtxZ]        = primVtx->GetZ();
    values[VAR::kNVtxContributors]    = primVtx->GetNContributors();
  }

  AliMultSelection *MultSelection = (AliMultSelection * ) fEvent->FindListObject("MultSelection");
  values[VAR::kVZEROMultPercentile] = MultSelection->GetMultiplicityPercentile("V0M", kTRUE);

  AliESDEvent* esdEvent = static_cast<AliESDEvent*>(fEvent);
  AliCentrality* cent = esdEvent->GetCentrality();
  if(cent){
    values[VAR::kCentVZERO]   = cent->GetCentralityPercentile("V0M");
    values[VAR::kCentSPD]     = cent->GetCentralityPercentile("CL1");
    values[VAR::kCentTPC]     = cent->GetCentralityPercentile("TRK");
    values[VAR::kCentQuality] = cent->GetQuality();
  }


  AliVVZERO* vzero = fEvent->GetVZEROData();
  values[VAR::kVZEROATotalMult]     = vzero->GetMTotV0A();
  values[VAR::kVZEROCTotalMult]     = vzero->GetMTotV0C();
  values[VAR::kVZEROTotalMult]      = values[VAR::kVZEROATotalMult]+values[VAR::kVZEROCTotalMult];

  AliMultiplicity* spdmult = (AliMultiplicity*) fEvent->GetMultiplicity();
  values[VAR::kSPDntracklets]      = spdmult->GetNumberOfTracklets();
  values[VAR::kSPDnSingleClusters] = spdmult->GetNumberOfSingleClusters();
}



//__________________________________________________________________
void QnCorrectionsFillEvent::FillTrackInfo(AliVParticle* particle, Float_t* values) {

  Float_t dcaxy=0.0;
  Float_t dcaz=0.0;

  values[VAR::kPx]        = particle->Px();
  values[VAR::kPy]        = particle->Py();
  values[VAR::kPz]        = particle->Pz();
  values[VAR::kPt]        = particle->Pt();
  values[VAR::kP]         = particle->P();
  values[VAR::kPhi]       = particle->Phi();
  values[VAR::kTheta]     = particle->Theta();
  values[VAR::kEta]       = particle->Eta();
  values[VAR::kCharge]    = particle->Charge();
  values[VAR::kDcaXY]     = dcaxy;
  values[VAR::kDcaZ]      = dcaz;

  AliAODTrack* aodTrack=static_cast<AliAODTrack*>(particle);

  //values[VAR::kITSncls]       = particle->GetNcls(0); 
  values[VAR::kTPCncls]       = aodTrack->GetTPCNcls();
  values[VAR::kTPCchi2]       = aodTrack->Chi2perNDF();
  values[VAR::kTPCsignal]     = aodTrack->GetTPCsignal();
  for(Int_t ibit=0; ibit<9; ibit++) values[VAR::kFilterBit+ibit]     = aodTrack->TestFilterBit(BIT(ibit));

}


//__________________________________________________________________
void QnCorrectionsFillEvent::FillTrackInfo(AliESDtrack* particle, Float_t* values) {

  Float_t dcaxy=0.0;
  Float_t dcaz=0.0;
  particle->GetImpactParameters(dcaxy,dcaz);

  values[VAR::kPx]        = particle->Px();
  values[VAR::kPy]        = particle->Py();
  values[VAR::kPz]        = particle->Pz();
  values[VAR::kPt]        = particle->Pt();
  values[VAR::kP]         = particle->P();
  values[VAR::kPhi]       = particle->Phi();
  values[VAR::kTheta]     = particle->Theta();
  values[VAR::kEta]       = particle->Eta();
  values[VAR::kCharge]    = particle->Charge();
  values[VAR::kDcaXY]     = dcaxy;
  values[VAR::kDcaZ]      = dcaz;

  values[VAR::kTPCncls]       = particle->GetTPCNcls();
  values[VAR::kTPCnclsIter1]  = particle->GetTPCNclsIter1();
  values[VAR::kTPCchi2]       = values[VAR::kTPCncls]>0 ? particle->GetTPCchi2()/values[VAR::kTPCncls] : 0.0;
  values[VAR::kTPCchi2Iter1]  = values[VAR::kTPCnclsIter1]>0 ? particle->GetTPCchi2Iter1()/values[VAR::kTPCnclsIter1] : 0.0;
  values[VAR::kTPCsignal]     = particle->GetTPCsignal();



}

//_________________________________
void QnCorrectionsFillEvent::FillDetectors(AliAnalysisTaskSE* task, Float_t* values){

  if(fFillTPC)   FillTPC(values);
  if(fFillVZERO) FillVZERO();
  if(fFillZDC)   FillZDC();
  if(fFillTZERO) FillTZERO();
  if(fFillFMD)   FillFMD(task);
  if(fFillRawFMD)FillRawFMD(values);
  if(fFillSPD) FillSPDTracklets(values);

}


//_________________________________
void QnCorrectionsFillEvent::FillTPC(Float_t* values){
  //
  // fill TPC info
  //

  if(fIsAOD) FillAodTPC(values);
  if(fIsESD) FillEsdTPC(values);


}


//_________________________________
void QnCorrectionsFillEvent::FillAodTPC(Float_t* values){
  //
  // fill AOD TPC info
  //

  AliVParticle* vTrack;

  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); ++iTrack)
  {
    vTrack = fEvent->GetTrack(iTrack); //carefull do not modify it othwise  need to work with a copy 
    if (!vTrack) continue;

    FillTrackInfo(vTrack, values);
    fQAhistos->FillHistClass("TrackQA_NoCuts", values);

    Int_t nNoOfAcceptedConf = fQnCorrectionsManager->AddDataVector(VAR::kTPC, vTrack->Phi());

    for(Int_t conf=0; conf < nNoOfAcceptedConf; conf++){
      fQAhistos->FillHistClass(Form("TrackQA_%s",
          fQnCorrectionsManager->GetAcceptedDataDetectorConfigurationName(VAR::kTPC, conf)), values);
    }
  }
}




//_________________________________
void QnCorrectionsFillEvent::FillEsdTPC(Float_t* values){
  //
  // fill ESD TPC info
  //

  AliESDtrack* esdTrack;

  const AliESDEvent* esdEvent = static_cast<AliESDEvent*>(fEvent);
  const AliESDEvent& esd = *esdEvent;

  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); ++iTrack)
  {
    AliESDtrack* track = NULL;
    esdTrack = esd.GetTrack(iTrack); //carefull do not modify it othwise  need to work with a copy 
    if(fUseTPCStandaloneTracks) track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(&esd),esdTrack->GetID());
    else track = esdTrack;
    if (!track) continue;

    FillTrackInfo(track, values);
    fQAhistos->FillHistClass("TrackQA_NoCuts", values);

    Int_t nNoOfAcceptedConf = fQnCorrectionsManager->AddDataVector(VAR::kTPC, track->Phi());

    for(Int_t conf=0; conf < nNoOfAcceptedConf; conf++){
      fQAhistos->FillHistClass(Form("TrackQA_%s",
          fQnCorrectionsManager->GetAcceptedDataDetectorConfigurationName(VAR::kTPC, conf)), values);
    }

    if(fUseTPCStandaloneTracks) delete track;
  }
}



//_________________________________________________________________________________
void QnCorrectionsFillEvent::FillSPDTracklets(Float_t* values) {
  //
  // fill SPD info
  //

  Int_t nTracklets = 0;

  AliMultiplicity* mult = (AliMultiplicity*) fEvent->GetMultiplicity();
  nTracklets = mult->GetNumberOfTracklets();
  for(Int_t iTracklet=0; iTracklet<nTracklets; ++iTracklet) {
    values[VAR::kSPDtrackletEta]    = mult->GetEta(iTracklet);
    values[VAR::kSPDtrackletPhi]    = mult->GetPhi(iTracklet);

    Int_t nNoOfAcceptedConf = fQnCorrectionsManager->AddDataVector(VAR::kSPD, values[VAR::kSPDtrackletPhi]);

    for(Int_t conf=0; conf < nNoOfAcceptedConf; conf++){
      fQAhistos->FillHistClass(Form("TrackletQA_%s",
          fQnCorrectionsManager->GetAcceptedDataDetectorConfigurationName(VAR::kSPD, conf)), values);
    }
  }
}

void QnCorrectionsFillEvent::FillVZERO(){
  //
  // fill VZERO info
  //

  Double_t weight=0.;
  const Double_t kX[8] = {0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268, 0.38268, 0.92388};    // cosines of the angles of the VZERO sectors (8 per ring)
  const Double_t kY[8] = {0.38268, 0.92388, 0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268};    // sines     -- " --

  AliVVZERO* vzero = fEvent->GetVZEROData();

  for(Int_t ich=0; ich<64; ich++){
    weight=vzero->GetMultiplicity(ich);
    if(weight<0.01) weight=0.;

    fQnCorrectionsManager->AddDataVector(VAR::kVZERO, TMath::ATan2(kY[ich%8],kX[ich%8]), weight, ich);   // 1st ich is position in array, 2nd ich is channel id

  }
}



//_________________________________
void QnCorrectionsFillEvent::FillTZERO(){
  //
  // fill ESD TZERO info
  //

  Double_t weight=0.0;
  const Double_t kX[24] = {/* Cside */ 0.905348,0.571718,0.0848977,-0.424671,-0.82045,-0.99639,-0.905348,-0.571718,-0.0848977,0.424671,0.82045,0.99639, /* Aside */ 0.99995,0.870982,0.508635,0.00999978,-0.491315,-0.860982,-0.99995,-0.870982,-0.508635,-0.0100001,0.491315,0.860982};
  const Double_t kY[24] = {/* Cside */ 0.424671,0.82045,0.99639,0.905348,0.571718,0.0848976,-0.424671,-0.82045,-0.99639,-0.905348,-0.571719,-0.0848975, /* Aside */ -0.00999983,0.491315,0.860982,0.99995,0.870982,0.508635,0.00999974,-0.491315,-0.860982,-0.99995,-0.870982,-0.508635};

  const AliESDTZERO* tzero= ((AliESDEvent*)fEvent)->GetESDTZERO();


  for(Int_t ich=0; ich<24; ich++){
    weight=tzero->GetT0amplitude()[ich];
    if(weight<0.01) weight=0.;

    fQnCorrectionsManager->AddDataVector(VAR::kTZERO, TMath::ATan2(kY[ich%8],kX[ich%8]), weight, ich);   // 1st ich is position in array, 2nd ich is channel id

  }
}




//_________________________________
void QnCorrectionsFillEvent::FillZDC(){
  //
  // fill ZDC info
  //


  Double_t weight=0.0;
  const Double_t kX[10] = { /* Cside */ 0.0,  -1.75,  1.75, -1.75, 1.75, /* Aside */  0.0,  1.75, -1.75, 1.75, -1.75  };
  const Double_t kY[10] = { /* Cside */ 0.0,  -1.75, -1.75,  1.75, 1.75, /* Aside */  0.0, -1.75, -1.75, 1.75,  1.75  };


  AliVZDC* zdc = (AliVZDC*) fEvent->GetZDCData();

  Double_t ZDCenergy[10];
  for(Int_t i=0; i<5; ++i)    ZDCenergy[i]  = zdc->GetZNCTowerEnergy()[i];
  for(Int_t i=5; i<10; ++i)   ZDCenergy[i]  = zdc->GetZNATowerEnergy()[i-5];

  for(Int_t ich=1; ich<10; ich++){
    if(ich==5) continue;
    weight=ZDCenergy[ich];
    if(weight<100.) weight=0.;

    fQnCorrectionsManager->AddDataVector(VAR::kZDC, TMath::ATan2(kY[ich%8],kX[ich%8]), weight, ich);   // 1st ich is position in array, 2nd ich is channel id
  }
}


//_________________________________________________________________________________
void QnCorrectionsFillEvent::FillFMD(AliAnalysisTaskSE* task)
{
  //
  // fill ESD FMD info
  //

  Float_t m,eta,phi;

  AliAODEvent* aodEvent = AliForwardUtil::GetAODEvent(task);


  if (!aodEvent) {
    AliFatal("didn't get AOD\n");
    return;
  }


  TObject* obj = aodEvent->FindListObject("Forward");  
  if (!obj) return;

  AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(obj);

  const TH2D& d2Ndetadphi = aodForward->GetHistogram();

  Float_t FMDtotalmult=0.0;
  Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
  Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();


  // Loop over eta 
  Int_t nFMD=-1;
  for (Int_t iEta = 1; iEta <= nEta; iEta++) {
    Int_t valid = d2Ndetadphi.GetBinContent(iEta, 0);
    if (!valid) continue; // No data expected for this eta 

    eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
    // Loop over phi 
    for (Int_t iPhi = 1; iPhi <= nPhi; iPhi++) {
      phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);
      m     =  d2Ndetadphi.GetBinContent(iEta, iPhi);
      if(m<0.01) continue;
      nFMD++;

      fQnCorrectionsManager->AddDataVector(VAR::kFMD, phi, m, iEta*nPhi+iPhi);   // 1st ich is position in array, 2nd ich is channel id
    }
  }
}



//_________________________________
void QnCorrectionsFillEvent::FillRawFMD(Float_t* values)
{
  //
  // fill Raw FMD info
  //
  Bool_t isESD = (fEvent->IsA()==AliESDEvent::Class());
  if(!isESD) return;

  AliESDEvent* esdEvent = static_cast<AliESDEvent*>(fEvent);

  AliESDFMD* esdFmd = esdEvent->GetFMDData();

  Int_t id=-1;
  Int_t maxDet=3;
  Int_t maxRing=2;
  Int_t maxSector;
  Int_t maxStrip;
  Float_t m=0.0;
  Double_t phi,eta;
  Char_t ring;

  for(UShort_t det = 1; det <= maxDet; ++det) {
    (det == 1 ? maxRing=1 : maxRing=2);
    for(UShort_t ir = 0; ir < maxRing; ++ir) {
      ring = (ir == 0 ? 'I' : 'O');
      (ir == 0 ? maxSector=20 : maxSector=40);
      (ir == 0 ? maxStrip=512 : maxStrip=256);
      for(UShort_t sec = 0; sec < maxSector; ++sec) {
        phi  =  esdFmd->Phi(det, ring, sec, 0)/180.*TMath::Pi();
        for(UShort_t str = 0; str < maxStrip; ++str) {
          ++id;
          eta  =  esdFmd->Eta(det, ring, sec, str);
          m    =  esdFmd->Multiplicity(det, ring, sec, str);
          if(m ==  AliESDFMD::kInvalidMult) m=0;
          values[AliQnCorrectionsVarManager::kFMDEta] = eta;
          fQnCorrectionsManager->AddDataVector(VAR::kFMDraw, phi, m, id);   // 1st ich is position in array, 2nd ich is channel id
        }  // end loop over strips
      }  // end loop over sectors      
    }  // end loop over rings
  } // end loop over detectors
}


