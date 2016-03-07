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

#include "QnCorrectionsDataVector.h"
#include "QnCorrectionsDetector.h"
#include "QnCorrectionsManager.h"

#include "AliQnCorrectionsHistos.h"

#include "AnalysisTaskFlowVectorCorrections.h"

void AnalysisTaskFlowVectorCorrections::SetDetectors() {
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
void AnalysisTaskFlowVectorCorrections::FillEventData() {

  TString aod = "AOD";
  TString esd = "ESD";

  fIsAOD = ( aod.EqualTo(fEvent->Whoami()) ? kTRUE : kFALSE );
  fIsESD = ( esd.EqualTo(fEvent->Whoami()) ? kTRUE : kFALSE );

  FillEventInfo();
  FillDetectors();
}

//__________________________________________________________________
void AnalysisTaskFlowVectorCorrections::FillEventInfo() {
  //
  // fill event info
  //


  fDataBank[VAR::kRunNo]       = fEvent->GetRunNumber();
  fDataBank[VAR::kVtxX]        = -999.;
  fDataBank[VAR::kVtxY]        = -999.;
  fDataBank[VAR::kVtxZ]        = -999.;
  const AliVVertex *primVtx = fEvent->GetPrimaryVertex();
  if (primVtx){
    fDataBank[VAR::kVtxX]        = primVtx->GetX();
    fDataBank[VAR::kVtxY]        = primVtx->GetY();
    fDataBank[VAR::kVtxZ]        = primVtx->GetZ();
    fDataBank[VAR::kNVtxContributors]    = primVtx->GetNContributors();
  }

  AliMultSelection *MultSelection = (AliMultSelection * ) fEvent->FindListObject("MultSelection");
  fDataBank[VAR::kVZEROMultPercentile] = MultSelection->GetMultiplicityPercentile("V0M", kTRUE);

  AliESDEvent* esdEvent = static_cast<AliESDEvent*>(fEvent);
  AliCentrality* cent = esdEvent->GetCentrality();
  if(cent){
    fDataBank[VAR::kCentVZERO]   = cent->GetCentralityPercentile("V0M");
    fDataBank[VAR::kCentSPD]     = cent->GetCentralityPercentile("CL1");
    fDataBank[VAR::kCentTPC]     = cent->GetCentralityPercentile("TRK");
    fDataBank[VAR::kCentQuality] = cent->GetQuality();
  }


  AliVVZERO* vzero = fEvent->GetVZEROData();
  fDataBank[VAR::kVZEROATotalMult]     = vzero->GetMTotV0A();
  fDataBank[VAR::kVZEROCTotalMult]     = vzero->GetMTotV0C();
  fDataBank[VAR::kVZEROTotalMult]      = fDataBank[VAR::kVZEROATotalMult]+fDataBank[VAR::kVZEROCTotalMult];

  AliMultiplicity* spdmult = (AliMultiplicity*) fEvent->GetMultiplicity();
  fDataBank[VAR::kSPDntracklets]      = spdmult->GetNumberOfTracklets();
  fDataBank[VAR::kSPDnSingleClusters] = spdmult->GetNumberOfSingleClusters();
}



//__________________________________________________________________
void AnalysisTaskFlowVectorCorrections::FillTrackInfo(AliVParticle* particle) {

  Float_t dcaxy=0.0;
  Float_t dcaz=0.0;

  fDataBank[VAR::kPx]        = particle->Px();
  fDataBank[VAR::kPy]        = particle->Py();
  fDataBank[VAR::kPz]        = particle->Pz();
  fDataBank[VAR::kPt]        = particle->Pt();
  fDataBank[VAR::kP]         = particle->P();
  fDataBank[VAR::kPhi]       = particle->Phi();
  fDataBank[VAR::kTheta]     = particle->Theta();
  fDataBank[VAR::kEta]       = particle->Eta();
  fDataBank[VAR::kCharge]    = particle->Charge();
  fDataBank[VAR::kDcaXY]     = dcaxy;
  fDataBank[VAR::kDcaZ]      = dcaz;

  AliAODTrack* aodTrack=static_cast<AliAODTrack*>(particle);

  //fDataBank[VAR::kITSncls]       = particle->GetNcls(0);
  fDataBank[VAR::kTPCncls]       = aodTrack->GetTPCNcls();
  fDataBank[VAR::kTPCchi2]       = aodTrack->Chi2perNDF();
  fDataBank[VAR::kTPCsignal]     = aodTrack->GetTPCsignal();
  for(Int_t ibit=0; ibit<9; ibit++) fDataBank[VAR::kFilterBit+ibit]     = aodTrack->TestFilterBit(BIT(ibit));

}


//__________________________________________________________________
void AnalysisTaskFlowVectorCorrections::FillTrackInfo(AliESDtrack* particle) {

  Float_t dcaxy=0.0;
  Float_t dcaz=0.0;
  particle->GetImpactParameters(dcaxy,dcaz);

  fDataBank[VAR::kPx]        = particle->Px();
  fDataBank[VAR::kPy]        = particle->Py();
  fDataBank[VAR::kPz]        = particle->Pz();
  fDataBank[VAR::kPt]        = particle->Pt();
  fDataBank[VAR::kP]         = particle->P();
  fDataBank[VAR::kPhi]       = particle->Phi();
  fDataBank[VAR::kTheta]     = particle->Theta();
  fDataBank[VAR::kEta]       = particle->Eta();
  fDataBank[VAR::kCharge]    = particle->Charge();
  fDataBank[VAR::kDcaXY]     = dcaxy;
  fDataBank[VAR::kDcaZ]      = dcaz;

  fDataBank[VAR::kTPCncls]       = particle->GetTPCNcls();
  fDataBank[VAR::kTPCnclsIter1]  = particle->GetTPCNclsIter1();
  fDataBank[VAR::kTPCchi2]       = fDataBank[VAR::kTPCncls]>0 ? particle->GetTPCchi2()/fDataBank[VAR::kTPCncls] : 0.0;
  fDataBank[VAR::kTPCchi2Iter1]  = fDataBank[VAR::kTPCnclsIter1]>0 ? particle->GetTPCchi2Iter1()/fDataBank[VAR::kTPCnclsIter1] : 0.0;
  fDataBank[VAR::kTPCsignal]     = particle->GetTPCsignal();



}

//_________________________________
void AnalysisTaskFlowVectorCorrections::FillDetectors(){

  if(fFillTPC)   FillTPC();
  if(fFillVZERO) FillVZERO();
  if(fFillZDC)   FillZDC();
  if(fFillTZERO) FillTZERO();
  if(fFillFMD)   FillFMD();
  if(fFillRawFMD)FillRawFMD();
  if(fFillSPD) FillSPDTracklets();
}


//_________________________________
void AnalysisTaskFlowVectorCorrections::FillTPC(){
  //
  // fill TPC info
  //

  if(fIsAOD) FillAodTPC();
  if(fIsESD) FillEsdTPC();
}


//_________________________________
void AnalysisTaskFlowVectorCorrections::FillAodTPC(){
  //
  // fill AOD TPC info
  //

  AliVParticle* vTrack;

  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); ++iTrack)
  {
    vTrack = fEvent->GetTrack(iTrack); //carefull do not modify it othwise  need to work with a copy 
    if (!vTrack) continue;

    FillTrackInfo(vTrack);
    fEventHistos->FillHistClass("TrackQA_NoCuts", fDataBank);

    Int_t nNoOfAcceptedConf = fQnCorrectionsManager->AddDataVector(VAR::kTPC, vTrack->Phi());

    for(Int_t conf=0; conf < nNoOfAcceptedConf; conf++){
        fEventHistos->FillHistClass(Form("TrackQA_%s",
            fQnCorrectionsManager->GetAcceptedDataDetectorConfigurationName(VAR::kTPC, conf)),
            fDataBank);
    }
  }
}




//_________________________________
void AnalysisTaskFlowVectorCorrections::FillEsdTPC(){
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

    FillTrackInfo(track);
    fQAhistos->FillHistClass("TrackQA_NoCuts", fDataBank);

    Int_t nNoOfAcceptedConf = fQnCorrectionsManager->AddDataVector(VAR::kTPC, track->Phi());

    for(Int_t conf=0; conf < nNoOfAcceptedConf; conf++){
        fEventHistos->FillHistClass(Form("TrackQA_%s",
            fQnCorrectionsManager->GetAcceptedDataDetectorConfigurationName(VAR::kTPC, conf)),
            fDataBank);
    }

    if(fUseTPCStandaloneTracks) delete track;
  }
}



//_________________________________________________________________________________
void AnalysisTaskFlowVectorCorrections::FillSPDTracklets() {
  //
  // fill SPD info
  //

  Int_t nTracklets = 0;

  AliMultiplicity* mult = (AliMultiplicity*) fEvent->GetMultiplicity();
  nTracklets = mult->GetNumberOfTracklets();
  for(Int_t iTracklet=0; iTracklet<nTracklets; ++iTracklet) {
    fDataBank[VAR::kSPDtrackletEta]    = mult->GetEta(iTracklet);
    fDataBank[VAR::kSPDtrackletPhi]    = mult->GetPhi(iTracklet);

    Int_t nNoOfAcceptedConf = fQnCorrectionsManager->AddDataVector(VAR::kSPD, fDataBank[VAR::kSPDtrackletPhi]);

    for(Int_t conf=0; conf < nNoOfAcceptedConf; conf++){
      fEventHistos->FillHistClass(Form("TrackletQA_%s",
          fQnCorrectionsManager->GetAcceptedDataDetectorConfigurationName(VAR::kSPD, conf)),
          fDataBank);
    }
  }
}

void AnalysisTaskFlowVectorCorrections::FillVZERO(){
  //
  // fill VZERO info
  //

  Double_t weight=0.;
  const Double_t kX[8] = {0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268, 0.38268, 0.92388};    // cosines of the angles of the VZERO sectors (8 per ring)
  const Double_t kY[8] = {0.38268, 0.92388, 0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268};    // sines     -- " --
  const Double_t phi[8] = {TMath::ATan2(kY[0],kX[0]), TMath::ATan2(kY[1],kX[1]), TMath::ATan2(kY[2],kX[2]),
      TMath::ATan2(kY[3],kX[3]), TMath::ATan2(kY[4],kX[4]), TMath::ATan2(kY[5],kX[5]),
      TMath::ATan2(kY[6],kX[6]), TMath::ATan2(kY[7],kX[7]) };

  AliVVZERO* vzero = fEvent->GetVZEROData();

  for(Int_t ich=0; ich<64; ich++){
    weight=vzero->GetMultiplicity(ich);
    if(weight<0.01) weight=0.;

    fQnCorrectionsManager->AddDataVector(VAR::kVZERO, phi[ich%8], weight, ich);   // 1st ich is position in array, 2nd ich is channel id

  }
}



void AnalysisTaskFlowVectorCorrections::FillTZERO(){
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
void AnalysisTaskFlowVectorCorrections::FillZDC(){
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


void AnalysisTaskFlowVectorCorrections::FillFMD()
{
  //
  // fill ESD FMD info
  //

  Float_t m,eta,phi;

  AliAODEvent* aodEvent = AliForwardUtil::GetAODEvent(this);


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
void AnalysisTaskFlowVectorCorrections::FillRawFMD()
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
          fDataBank[AliQnCorrectionsVarManager::kFMDEta] = eta;
          fQnCorrectionsManager->AddDataVector(VAR::kFMDraw, phi, m, id);   // 1st ich is position in array, 2nd ich is channel id
        }  // end loop over strips
      }  // end loop over sectors      
    }  // end loop over rings
  } // end loop over detectors
}


