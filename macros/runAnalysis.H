/* Copyright(c) 2013-2014, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef RUNANALYSIS_H_
#define RUNANALYSIS_H_

#include "runs.H"

extern Bool_t bGRIDPlugin;
extern Bool_t bMC;

extern Double_t *globalEtaMin;
extern Double_t *globalEtaMax;
extern Double_t *globalPtMin;
extern Double_t *globalPtMax;
extern Double_t *centralityMin;
extern Double_t *centralityMax;
extern Double_t *zvertexMin;
extern Double_t *zvertexMax;
extern Double_t *sigmasTPCMin;
extern Double_t *sigmasTPCMax;
extern Double_t *sigmasTOFMin;
extern Double_t *sigmasTOFMax;
extern Int_t *nTrackSelectionCriteria;
extern Int_t *nPIDCriteria;
extern Float_t *fMixingPercentage;
extern Int_t *multiplicityMin;
extern Int_t *multiplicityMax;
extern Int_t *minSwarmSize;
extern Int_t nGlobalCuts;

extern Double_t *swarmPhiMin;
extern Double_t *swarmPhiMax;
extern Double_t *swarmEtaMin;
extern Double_t *swarmEtaMax;
extern Double_t *swarmPtMin;
extern Double_t *swarmPtMax;
extern Int_t nSwarmCuts;

/* Task level cuts */
extern Int_t nCentralityCuts;
extern Double_t *minCentrality;
extern Double_t *maxCentrality;

/* swarm level cuts */
extern Int_t nDeltasPt;
extern Double_t *deltaPt;
extern Int_t *noOfPtPoints;
extern Double_t **pointPt;

extern Int_t nDeltasEta;
extern Double_t *deltaEta;
extern Int_t *noOfEtaPoints;
extern Double_t **pointEta;

extern Int_t nDeltasPhi;
extern Double_t *deltaPhi;
extern Int_t *noOfPhiPoints;
extern Double_t **pointPhi;

extern TString szAliPhysicsVersion;
extern TString szDataPattern;
extern TString szDataDir;
extern TString szRunPrefix;
extern Int_t nNoOfInputFiles;
extern TObjArray listOfRuns;

extern TString szLocalFileList;

Bool_t loadRunOptions(Bool_t verb = kFALSE,const char *filename = "runoptions.txt");

#endif /* RUNANALYSIS_H_ */
