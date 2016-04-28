/* Copyright(c) 2013-2016, ALICE Experiment at CERN, All rights reserved. *
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

extern Bool_t bGRIDPlugin;
extern Bool_t bMC;

/* temporal flag to use multiplicity instead of centrality and to inhibit detectors for 2015 dataset */
extern Bool_t bUseMultiplicity;
extern Bool_t b2015DataSet;

/* Task level cuts */
extern Double_t centralityMin;
extern Double_t centralityMax;
extern Double_t zvertexMin;
extern Double_t zvertexMax;

/* Running conditions */
extern TString szAliPhysicsVersion;
extern TString szDataPattern;
extern TString szDataDir;
extern TString szRunPrefix;
extern Int_t nNoOfInputFiles;
extern TObjArray listOfRuns;
extern TObjArray listOfActiveRuns;
extern TString szCorrectionsSource;
extern TString szCorrectionsFilePath;
extern TString szCorrectionsFileName;

extern TString szLocalFileList;

Bool_t loadRunOptions(Bool_t verb = kFALSE,const char *filename = "runoptions.txt");

#endif /* RUNANALYSIS_H_ */
