/**************************************************************************
 * Copyright(c) 2013-2014, ALICE Experiment at CERN, All rights reserved. *
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

#include "Riostream.h"
#include "TSystem.h"
#include "TROOT.h"

#include "runAnalysisCriteria.H"

Bool_t bOptionsLoaded = kFALSE;
void CleanOptions();

TString szRootVersion;
TString szAliRootVersion;
TString szAliPhysicsVersion;
TString szDataPattern;
TString szDataDir;
TString szRunPrefix;
TString szLocalFileList;
Int_t nNoOfInputFiles;
TObjArray listOfRuns;

Bool_t bGRIDPlugin;
Bool_t bMC;
Double_t *globalEtaMin;
Double_t *globalEtaMax;
Double_t *globalPtMin;
Double_t *globalPtMax;
Double_t *centralityMin;
Double_t *centralityMax;
Double_t *zvertexMin;
Double_t *zvertexMax;
Double_t *sigmasTPCMin;
Double_t *sigmasTPCMax;
Double_t *sigmasTOFMin;
Double_t *sigmasTOFMax;
Int_t *nTrackSelectionCriteria;
Int_t *nPIDCriteria;
Float_t *fMixingPercentage;
Int_t *multiplicityMin;
Int_t *multiplicityMax;
Int_t *minSwarmSize;
Int_t nGlobalCuts;

Double_t *swarmPhiMin;
Double_t *swarmPhiMax;
Double_t *swarmEtaMin;
Double_t *swarmEtaMax;
Double_t *swarmPtMin;
Double_t *swarmPtMax;
Int_t nSwarmCuts;

/* Task level cuts */
Int_t nCentralityCuts;
Int_t nTaskPtCuts;
Int_t nTaskEtaCuts;
Int_t nMultiplicityCuts;
Int_t nSwarmSizeCuts;
Int_t nVertexCuts;
Int_t nTrackCriteriaCuts;
Int_t nPIDCriteriaCuts;
Int_t nSigmaTPCCuts;
Int_t nSigmaTOFCuts;
Int_t nEventMixingPercent;

Double_t *minCentrality;
Double_t *maxCentrality;
Double_t *minPt;
Double_t *maxPt;
Double_t *minEta;
Double_t *maxEta;
Double_t *minzVertex;
Double_t *maxzVertex;
Double_t *minsTPC;
Double_t *maxsTPC;
Double_t *minsTOF;
Double_t *maxsTOF;
Double_t *mixPerc;
Int_t *trackSel;
Int_t *pidSel;
Int_t *minMult;
Int_t *maxMult;
Int_t *swarmSize;

/* swarm level cuts */
#define MAXNOOFDELTACUTS 10
Int_t nMaxDeltas = MAXNOOFDELTACUTS;
Int_t nDeltasPt;
Double_t *deltaPt;
Int_t *noOfPtPoints;
Double_t **pointPt;

Int_t nDeltasEta;
Double_t *deltaEta;
Int_t *noOfEtaPoints;
Double_t **pointEta;

Int_t nDeltasPhi;
Double_t *deltaPhi;
Int_t *noOfPhiPoints;
Double_t **pointPhi;


Bool_t loadRunOptions(Bool_t verb,const char *filename) {
  ifstream optionsfile;
  TString currline;

  if (bOptionsLoaded) CleanOptions();

  nCentralityCuts = 0;
  nTaskPtCuts = 0;
  nTaskEtaCuts = 0;
  nMultiplicityCuts = 0;
  nSwarmSizeCuts = 0;
  nVertexCuts = 0;
  nTrackCriteriaCuts = 0;
  nPIDCriteriaCuts = 0;
  nSigmaTPCCuts = 0;
  nSigmaTOFCuts = 0;
  nEventMixingPercent = 0;
  nDeltasPt = 0;
  nDeltasEta = 0;
  nDeltasPhi = 0;

  listOfRuns.SetOwner(kTRUE);
  optionsfile.open(filename);
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);

  /* Load run options */
  if (!currline.EqualTo("Run options:")) { printf("ERROR: wrong run options file %s\n", filename); return kFALSE; }
  /* grid option */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.EqualTo("grid"))
    bGRIDPlugin = kTRUE;
  else if (currline.EqualTo("local"))
    bGRIDPlugin = kFALSE;
  else
    { printf("ERROR: wrong grid option in options file %s\n", filename); return kFALSE; }
  printf("  Running in %s\n", bGRIDPlugin ? "grid" : "local");

  /* MC option */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.EqualTo("real"))
    bMC = kFALSE;
  else if (currline.EqualTo("MC"))
    bMC = kTRUE;
  else
    { printf("ERROR: wrong MC option in options file %s\n", filename); return kFALSE; }
  printf("  with %s data\n", bMC ? "simulated MC" : "real");

  /* SW versions */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("RootVersion: ")) {
    szRootVersion = currline.Remove(0,strlen("RootVersion: "));
  }
  else
    { printf("ERROR: wrong root SW version in options file %s\n", filename); return kFALSE; }
  if (bGRIDPlugin)
    cout << "    Grid ROOT version: " << szRootVersion << endl;
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("AliRootVersion: ")) {
    szAliRootVersion = currline.Remove(0,strlen("AliRootVersion: "));
  }
  else
    { printf("ERROR: wrong Aliroot SW version in options file %s\n", filename); return kFALSE; }
  if (bGRIDPlugin)
    cout << "    Grid AliROOT version: " << szAliRootVersion << endl;
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("AliPhysicsVersion: ")) {
    szAliPhysicsVersion = currline.Remove(0,strlen("AliPhysicsVersion: "));
  }
  else
    { printf("ERROR: wrong Aliphysics SW version in options file %s\n", filename); return kFALSE; }
  if (bGRIDPlugin)
    cout << "    Grid AliPhysics version: " << szAliPhysicsVersion << endl;

  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.EqualTo("Task level:")) {
    currline.ReadLine(optionsfile);
    while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
    while(!currline.EqualTo("end")) {
      /* centrality cuts */
      if (currline.BeginsWith("Centrality: ")) {
        currline.Remove(0, strlen("Centrality: "));
        if (verb) cout << "    Centrality cuts string: " << currline << endl;
        TObjArray *pTokenArray = currline.Tokenize(",");
        TIter nextToken(pTokenArray);
        nCentralityCuts = pTokenArray->GetEntriesFast();
        minCentrality = new Double_t[nCentralityCuts];
        maxCentrality = new Double_t[nCentralityCuts];
        if (verb) printf("     %d centrality cuts\n", nCentralityCuts);
        TObjString *ptoken;
        for (Int_t i = 0; i < nCentralityCuts; i++) {
          ptoken = (TObjString *) nextToken();
          Double_t min;
          Double_t max;
          sscanf(ptoken->GetString().Data(), "%lf-%lf", &min, &max);
          minCentrality[i] = min;
          maxCentrality[i] = max;
        }
        delete pTokenArray;
        for (Int_t i = 0; i < nCentralityCuts; i++) {
          printf ("    centrality cut %d: %.1f-%.1f\n",i,minCentrality[i],maxCentrality[i]);
        }
        currline.ReadLine(optionsfile);
        while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
      } /* end centrality cuts */

      /* z Vertex cuts */
      if (currline.BeginsWith("zvertex: ")) {
        currline.Remove(0, strlen("zvertex: "));
        if (verb) cout << "    z vertex cuts string: " << currline << endl;
        TObjArray *pTokenArray = currline.Tokenize(",");
        TIter nextToken(pTokenArray);
        nVertexCuts = pTokenArray->GetEntriesFast();
        minzVertex = new Double_t[nVertexCuts];
        maxzVertex = new Double_t[nVertexCuts];
        if (verb) printf("     %d z vertex cuts\n", nVertexCuts);
        TObjString *ptoken;
        for (Int_t i = 0; i < nVertexCuts; i++) {
          ptoken = (TObjString *) nextToken();
          Double_t min;
          Double_t max;
          sscanf(ptoken->GetString().Data(), "%lf:%lf", &min, &max);
          minzVertex[i] = min;
          maxzVertex[i] = max;
        }
        delete pTokenArray;
        for (Int_t i = 0; i < nVertexCuts; i++) {
          printf ("    z vertex cut %d: %.1f:%.1f\n",i,minzVertex[i],maxzVertex[i]);
        }
        currline.ReadLine(optionsfile);
        while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
      } /* end centrality cuts */

      /* Load track selection criteria */
      if (currline.BeginsWith("Track selection: ")) {
        currline.Remove(0, strlen("Track selection: "));
        if (verb) cout << "    Track selection cuts string: " << currline << endl;
        TObjArray *pTokenArray = currline.Tokenize(", ");
        TIter nextToken(pTokenArray);
        nTrackCriteriaCuts = pTokenArray->GetEntriesFast();
        trackSel = new Int_t[nTrackCriteriaCuts];
        if (verb) printf("     %d track selection cuts\n", nTrackCriteriaCuts);
        TObjString *ptoken;
        for (Int_t i = 0; i < nTrackCriteriaCuts; i++) {
          ptoken = (TObjString *) nextToken();
          if (ptoken->GetString().EqualTo("clusters"))
            trackSel[i] = kClusters;
          else if (ptoken->GetString().EqualTo("ratiofindable"))
            trackSel[i] = kRatioFindable;
          else
            { printf("ERROR: wrong track selection criteria in options file %s\n", filename); return kFALSE; }
        }
        delete pTokenArray;
        for (Int_t i = 0; i < nTrackCriteriaCuts; i++) {
          switch(trackSel[i]) {
          case kClusters:
            printf ("    track selection cut %d: %s\n",i,"clusters");
            break;
          case kRatioFindable:
            printf ("    track selection cut %d: %s\n",i,"ratiofindable");
            break;
          default:
            printf("ERROR: wrong track selection criteria in options file %s\n", filename);
            return kFALSE;
          }
        }
        currline.ReadLine(optionsfile);
        while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
      } /* end track selection criteria cuts */

      /* Load PID selection criteria */
      if (currline.BeginsWith("PID: ")) {
        currline.Remove(0, strlen("PID: "));
        if (verb) cout << "    PID selection cuts string: " << currline << endl;
        TObjArray *pTokenArray = currline.Tokenize(", ");
        TIter nextToken(pTokenArray);
        nPIDCriteriaCuts = pTokenArray->GetEntriesFast();
        pidSel = new Int_t[nPIDCriteriaCuts];
        if (verb) printf("     %d PID selection cuts\n", nPIDCriteriaCuts);
        TObjString *ptoken;
        for (Int_t i = 0; i < nPIDCriteriaCuts; i++) {
          ptoken = (TObjString *) nextToken();
          if (ptoken->GetString().EqualTo("tpcplustof"))
            pidSel[i] = kTPCplusTOF;
          else if (ptoken->GetString().EqualTo("tpcandtof"))
            pidSel[i] = kTPCandTOF;
          else
            { printf("ERROR: wrong PID criteria in options file %s\n", filename); return kFALSE; }
        }
        delete pTokenArray;
        for (Int_t i = 0; i < nPIDCriteriaCuts; i++) {
          switch(pidSel[i]) {
          case kTPCplusTOF:
            printf ("    PID selection cut %d: %s\n",i,"tpcplustof");
            break;
          case kTPCandTOF:
            printf ("    PID selection cut %d: %s\n",i,"tpcandtof");
            break;
          default:
            printf("ERROR: wrong track selection criteria in options file %s\n", filename);
            return kFALSE;
          }
        }
        currline.ReadLine(optionsfile);
        while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
      } /* end PID criteria cuts */

      /* TOF n sigmas cuts */
      if (currline.BeginsWith("NsTOF: ")) {
        currline.Remove(0, strlen("NsTOF: "));
        if (verb) cout << "    n sigmas TOF cuts string: " << currline << endl;
        TObjArray *pTokenArray = currline.Tokenize(",");
        TIter nextToken(pTokenArray);
        nSigmaTOFCuts = pTokenArray->GetEntriesFast();
        minsTOF = new Double_t[nSigmaTOFCuts];
        maxsTOF = new Double_t[nSigmaTOFCuts];
        if (verb) printf("     %d TOF n sigmas cuts\n", nSigmaTOFCuts);
        TObjString *ptoken;
        for (Int_t i = 0; i < nSigmaTOFCuts; i++) {
          ptoken = (TObjString *) nextToken();
          Double_t min;
          Double_t max;
          sscanf(ptoken->GetString().Data(), "%lf:%lf", &min, &max);
          minsTOF[i] = min;
          maxsTOF[i] = max;
        }
        delete pTokenArray;
        for (Int_t i = 0; i < nSigmaTOFCuts; i++) {
          printf ("    TOF n sigmas cut %d: %.1f:%.1f\n",i,minsTOF[i],maxsTOF[i]);
        }
        currline.ReadLine(optionsfile);
        while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
      } /* end TOF n sigmas cuts */

      /* TPC n sigmas cuts */
      if (currline.BeginsWith("NsTPC: ")) {
        currline.Remove(0, strlen("NsTPC: "));
        if (verb) cout << "    n sigmas TPC cuts string: " << currline << endl;
        TObjArray *pTokenArray = currline.Tokenize(",");
        TIter nextToken(pTokenArray);
        nSigmaTPCCuts = pTokenArray->GetEntriesFast();
        minsTPC = new Double_t[nSigmaTPCCuts];
        maxsTPC = new Double_t[nSigmaTPCCuts];
        if (verb) printf("     %d TPC n sigmas cuts\n", nSigmaTPCCuts);
        TObjString *ptoken;
        for (Int_t i = 0; i < nSigmaTPCCuts; i++) {
          ptoken = (TObjString *) nextToken();
          Double_t min;
          Double_t max;
          sscanf(ptoken->GetString().Data(), "%lf:%lf", &min, &max);
          minsTPC[i] = min;
          maxsTPC[i] = max;
        }
        delete pTokenArray;
        for (Int_t i = 0; i < nSigmaTPCCuts; i++) {
          printf ("    TPC n sigmas cut %d: %.1f:%.1f\n",i,minsTPC[i],maxsTPC[i]);
        }
        currline.ReadLine(optionsfile);
        while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
      } /* end TPC n sigmas cuts */

      /* event mixing percentage */
      if (currline.BeginsWith("mixing: ")) {
        currline.Remove(0, strlen("mixing: "));
        if (verb) cout << "    mixing percentage string: " << currline << endl;
        TObjArray *pTokenArray = currline.Tokenize(",");
        TIter nextToken(pTokenArray);
        nEventMixingPercent = pTokenArray->GetEntriesFast();
        mixPerc = new Double_t[nEventMixingPercent];
        if (verb) printf("     %d event mixing percentages\n", nEventMixingPercent);
        TObjString *ptoken;
        for (Int_t i = 0; i < nEventMixingPercent; i++) {
          ptoken = (TObjString *) nextToken();
          Float_t percent;
          sscanf(ptoken->GetString().Data(), "%f", &percent);
          mixPerc[i] = percent;
        }
        delete pTokenArray;
        for (Int_t i = 0; i < nEventMixingPercent; i++) {
          printf ("    event mixing percentage %d: %.1f\n",i,mixPerc[i]);
        }
        currline.ReadLine(optionsfile);
        while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
      } /* end event mixing percentage */
    }
  }
  else
    { printf("ERROR: wrong Task cuts section in options file %s\n", filename); return kFALSE; }

  /* sanity checks to fill not included cuts */
  /* Pt and eta cuts are left after the swarm cuts */
  if (nCentralityCuts == 0) {
    nCentralityCuts = 1;
    minCentrality = new Double_t[16];
    maxCentrality = new Double_t[16];
    minCentrality[0] = 0.0;
    maxCentrality[0] = 100.0;
  }
  if (nMultiplicityCuts == 0) {
    nMultiplicityCuts = 1;
    minMult = new Int_t[16];
    maxMult = new Int_t[16];
    minMult[0] = 0.0;
    maxMult[0] = 999999.0;
  }
  if (nSwarmSizeCuts == 0) {
    nSwarmSizeCuts = 1;
    swarmSize = new Int_t[16];
    swarmSize[0] = 5;
  }

  if (nVertexCuts == 0) {
    nVertexCuts = 1;
    minzVertex = new Double_t[16];
    maxzVertex = new Double_t[16];
    minzVertex[nVertexCuts] = -10.0;
    maxzVertex[nVertexCuts] = 10.0;
  }

  if (nTrackCriteriaCuts == 0) {
    nTrackCriteriaCuts = 1;
    trackSel = new Int_t[16];
    trackSel[nTrackCriteriaCuts] = kClusters;
  }

  if (nPIDCriteriaCuts == 0) {
    nPIDCriteriaCuts = 1;
    pidSel = new Int_t[16];
    pidSel[nPIDCriteriaCuts] = kTPCplusTOF;
  }

  if (nSigmaTPCCuts == 0) {
    nSigmaTPCCuts = 1;
    minsTPC = new Double_t[16];
    maxsTPC = new Double_t[16];
    minsTPC[nSigmaTPCCuts] = -3.0;
    maxsTPC[nSigmaTPCCuts] = 3.0;
  }

  if (nSigmaTOFCuts == 0) {
    nSigmaTOFCuts = 1;
    minsTOF = new Double_t[16];
    maxsTOF = new Double_t[16];
    minsTOF[nSigmaTOFCuts] = -3.0;
    maxsTOF[nSigmaTOFCuts] = 3.0;
  }

  if (nEventMixingPercent == 0) {
    nEventMixingPercent = 1;
    mixPerc = new Double_t[16];
    mixPerc[nEventMixingPercent] = 100.0;
  }

  /* now the swarm cuts data */
  nDeltasPt = 0;
  deltaPt = new Double_t [MAXNOOFDELTACUTS];
  noOfPtPoints = new Int_t [MAXNOOFDELTACUTS];
  pointPt = new Double_t* [MAXNOOFDELTACUTS];

  nDeltasEta = 0;
  deltaEta = new Double_t [MAXNOOFDELTACUTS];
  noOfEtaPoints = new Int_t [MAXNOOFDELTACUTS];
  pointEta = new Double_t* [MAXNOOFDELTACUTS];

  nDeltasPhi = 0;
  deltaPhi = new Double_t [MAXNOOFDELTACUTS];
  noOfPhiPoints = new Int_t [MAXNOOFDELTACUTS];
  pointPhi = new Double_t* [MAXNOOFDELTACUTS];

  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.EqualTo("Swarm level:")) {
    currline.ReadLine(optionsfile);
    while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
    while(!currline.EqualTo("end")) {
      if (currline.EqualTo("Delta pT:")) {
        currline.ReadLine(optionsfile);
        while(!currline.EqualTo("Delta eta:")) {
          if (verb) cout << "    Delta Pt cuts string: " << currline << endl;
          TObjArray *pTokenArray = currline.Tokenize(">,");
          TIter nextToken(pTokenArray);
          TObjString *pToken;
          /* the first one is the delta value */
          noOfPtPoints[nDeltasPt] = pTokenArray->GetEntriesFast() - 1;
          pointPt[nDeltasPt] = new Double_t[noOfPtPoints[nDeltasPt]];
          pToken = (TObjString *) nextToken();
          deltaPt[nDeltasPt] = pToken->GetString().Atof();
          for (Int_t i = 0; i < noOfPtPoints[nDeltasPt]; i++){
            pToken = (TObjString *) nextToken();
            pointPt[nDeltasPt][i] = pToken->GetString().Atof();
          }
          delete pTokenArray;
          printf("    delta Pt = %.1f; %d cuts centred at: ", deltaPt[nDeltasPt],noOfPtPoints[nDeltasPt]);
          for (Int_t i = 0; i < noOfPtPoints[nDeltasPt]; i++) {
            printf ("%.1f, ", pointPt[nDeltasPt][i]);
          }
          printf("\n");
          currline.ReadLine(optionsfile);
          nDeltasPt++;
          if (!(nDeltasPt < nMaxDeltas))
            {printf("ERROR: too many pT deltas in options file %s\n", filename); return kFALSE;}
        }
      }
      if (currline.EqualTo("Delta eta:")) {
        currline.ReadLine(optionsfile);
        while(!currline.EqualTo("Delta phi:")) {
          if (verb) cout << "    Delta eta cuts string: " << currline << endl;
          TObjArray *pTokenArray = currline.Tokenize(">,");
          TIter nextToken(pTokenArray);
          TObjString *pToken;
          /* the first one is the delta value */
          noOfEtaPoints[nDeltasEta] = pTokenArray->GetEntriesFast() - 1;
          pointEta[nDeltasEta] = new Double_t[noOfEtaPoints[nDeltasEta]];
          pToken = (TObjString *) nextToken();
          deltaEta[nDeltasEta] = pToken->GetString().Atof();
          for (Int_t i = 0; i < noOfEtaPoints[nDeltasEta]; i++){
            pToken = (TObjString *) nextToken();
            pointEta[nDeltasEta][i] = pToken->GetString().Atof();
          }
          delete pTokenArray;
          printf("    delta eta = %.1f; %d cuts centred at: ", deltaEta[nDeltasEta],noOfEtaPoints[nDeltasEta]);
          for (Int_t i = 0; i < noOfEtaPoints[nDeltasEta]; i++) {
            printf ("%.2f, ", pointEta[nDeltasEta][i]);
          }
          printf("\n");
          currline.ReadLine(optionsfile);
          nDeltasEta++;
          if (!(nDeltasEta < nMaxDeltas))
            {printf("ERROR: too many eta deltas in options file %s\n", filename); return kFALSE;}
        }
      }
      if (currline.EqualTo("Delta phi:")) {
        currline.ReadLine(optionsfile);
        while(!currline.EqualTo("end")) {
          if (verb) cout << "    Delta phi cuts string: " << currline << endl;
          TObjArray *pTokenArray = currline.Tokenize(">,");
          TIter nextToken(pTokenArray);
          TObjString *pToken;
          /* the first one is the delta value */
          noOfPhiPoints[nDeltasPhi] = pTokenArray->GetEntriesFast() - 1;
          pointPhi[nDeltasPhi] = new Double_t[noOfPhiPoints[nDeltasPhi]];
          pToken = (TObjString *) nextToken();
          deltaPhi[nDeltasPhi] = pToken->GetString().Atof();
          for (Int_t i = 0; i < noOfPhiPoints[nDeltasPhi]; i++){
            pToken = (TObjString *) nextToken();
            pointPhi[nDeltasPhi][i] = pToken->GetString().Atof();
          }
          delete pTokenArray;
          printf("    delta phi = %.1f; %d cuts centred at: ", deltaPhi[nDeltasPhi],noOfPhiPoints[nDeltasPhi]);
          for (Int_t i = 0; i < noOfPhiPoints[nDeltasPhi]; i++) {
            printf ("%.1f, ", pointPhi[nDeltasPhi][i]);
          }
          printf("\n");
          currline.ReadLine(optionsfile);
          nDeltasPhi++;
          if (!(nDeltasPhi < nMaxDeltas))
            {printf("ERROR: too many phi deltas in options file %s\n", filename); return kFALSE;}
        }
      }
    }
  }
  else
    { printf("ERROR: wrong Swarm cuts section in options file %s\n", filename); return kFALSE; }

  /* so everything went ok! let's build the cuts table */
  /* First the swarm cuts */
  /* let's evaluate the number of swarms */
  Int_t nPtCuts = 0;
  Int_t nEtaCuts = 0;
  Int_t nPhiCuts = 0;
  for (Int_t i=0;i<nDeltasPt;i++)
    nPtCuts += noOfPtPoints[i];
  for (Int_t i=0;i<nDeltasEta;i++)
    nEtaCuts += noOfEtaPoints[i];
  for (Int_t i=0;i<nDeltasPhi;i++)
    nPhiCuts += noOfPhiPoints[i];
  nSwarmCuts = nPtCuts*nEtaCuts*nPhiCuts;
  printf("  Total swarms per task: %d\n",nSwarmCuts);

  /* now allocate the space for the cuts data */
  swarmPhiMin = new Double_t[nSwarmCuts];
  swarmPhiMax = new Double_t[nSwarmCuts];
  swarmEtaMin = new Double_t[nSwarmCuts];
  swarmEtaMax = new Double_t[nSwarmCuts];
  swarmPtMin = new Double_t[nSwarmCuts];
  swarmPtMax = new Double_t[nSwarmCuts];

  /* now fill it keeping track of max and min eta and Pt */
  Double_t __globalPtMin = 1e7;
  Double_t __globalPtMax = 0.;
  Double_t __globalEtaMin = 1e7;
  Double_t __globalEtaMax = -1e7;
  Int_t ix = 0;
  for (Int_t ixDeltaPt = 0; ixDeltaPt < nDeltasPt; ixDeltaPt++) {
    for (Int_t ixPtPoint = 0; ixPtPoint < noOfPtPoints[ixDeltaPt]; ixPtPoint++) {
      Double_t __minPt = pointPt[ixDeltaPt][ixPtPoint] - deltaPt[ixDeltaPt] / 2.;
      Double_t __maxPt = pointPt[ixDeltaPt][ixPtPoint] + deltaPt[ixDeltaPt] / 2.;
      if (__minPt < __globalPtMin) __globalPtMin = __minPt;
      if (__globalPtMax < __maxPt) __globalPtMax = __maxPt;
      for (Int_t ixDeltaEta = 0; ixDeltaEta < nDeltasEta; ixDeltaEta++) {
        for (Int_t ixEtaPoint = 0; ixEtaPoint < noOfEtaPoints[ixDeltaEta]; ixEtaPoint++) {
          Double_t __minEta = pointEta[ixDeltaEta][ixEtaPoint] - deltaEta[ixDeltaEta] / 2.;
          Double_t __maxEta = pointEta[ixDeltaEta][ixEtaPoint] + deltaEta[ixDeltaEta] / 2.;
          if (__minEta < __globalEtaMin) __globalEtaMin = __minEta;
          if (__globalEtaMax < __maxEta) __globalEtaMax = __maxEta;
          for (Int_t ixDeltaPhi = 0; ixDeltaPhi < nDeltasPhi; ixDeltaPhi++) {
            for (Int_t ixPhiPoint = 0; ixPhiPoint < noOfPhiPoints[ixDeltaPhi]; ixPhiPoint++) {
              Double_t __minPhi = pointPhi[ixDeltaPhi][ixPhiPoint] - deltaPhi[ixDeltaPhi] / 2.;
              Double_t __maxPhi = pointPhi[ixDeltaPhi][ixPhiPoint] + deltaPhi[ixDeltaPhi] / 2.;
              swarmPhiMin[ix] = __minPhi;
              swarmPhiMax[ix] = __maxPhi;
              swarmEtaMin[ix] = __minEta;
              swarmEtaMax[ix] = __maxEta;
              swarmPtMin[ix] = __minPt;
              swarmPtMax[ix] = __maxPt;
              if (verb) printf("swarm ix = %d\n", ix);
              ix++;
            }
          }
        }
      }
    }
  }

  /* now the task cuts */
  /* first complete the sanity checks */
  if (nTaskPtCuts == 0) {
    nTaskPtCuts = 1;
    minPt = new Double_t[16];
    maxPt = new Double_t[16];
    minPt[0] = __globalPtMin;
    maxPt[0] = __globalPtMax;
  }
  if (nTaskEtaCuts == 0) {
    nTaskEtaCuts = 1;
    minEta = new Double_t[16];
    maxEta = new Double_t[16];
    minEta[0] = __globalEtaMin;
    maxEta[0] = __globalEtaMax;
  }




  /* now evaluate the number of cuts */
  nGlobalCuts = nCentralityCuts * nTaskPtCuts * nTaskEtaCuts * nMultiplicityCuts * nSwarmSizeCuts
      * nVertexCuts * nTrackCriteriaCuts * nPIDCriteriaCuts * nSigmaTPCCuts * /* nSigmaTOFCuts */ nEventMixingPercent;
  printf("  Total tasks: %d\n", nGlobalCuts);
  /* allocate memory space */
  globalEtaMin = new Double_t [nGlobalCuts];
  globalEtaMax = new Double_t [nGlobalCuts];
  globalPtMin = new Double_t [nGlobalCuts];
  globalPtMax = new Double_t [nGlobalCuts];
  centralityMin = new Double_t [nGlobalCuts];
  centralityMax = new Double_t [nGlobalCuts];
  zvertexMin = new Double_t [nGlobalCuts];
  zvertexMax = new Double_t [nGlobalCuts];
  sigmasTPCMin = new Double_t [nGlobalCuts];
  sigmasTPCMax = new Double_t [nGlobalCuts];
  sigmasTOFMin = new Double_t [nGlobalCuts];
  sigmasTOFMax = new Double_t [nGlobalCuts];
  nTrackSelectionCriteria = new Int_t [nGlobalCuts];
  nPIDCriteria = new Int_t [nGlobalCuts];
  fMixingPercentage = new Float_t [nGlobalCuts];
  multiplicityMin = new Int_t [nGlobalCuts];
  multiplicityMax = new Int_t [nGlobalCuts];
  minSwarmSize = new Int_t [nGlobalCuts];

  if (verb) {
    printf(
        "  nCentralityCuts: %d\n"
        "  nTaskPtCuts: %d\n"
        "  nTaskEtaCuts: %d\n"
        "  nVertexCuts: %d\n"
        "  nTrackCriteriaCuts: %d\n"
        "  nPIDCriteriaCuts: %d\n"
        "  nSigmaTPCCuts: %d\n"
        "  nSigmaTOFCuts: %d\n"
        "  nEventMixingPercent: %d\n"
        "  nMultiplicityCuts: %d\n"
        "  nSwarmSizeCuts: %d\n",
        nCentralityCuts,
        nTaskPtCuts,
        nTaskEtaCuts,
        nVertexCuts,
        nTrackCriteriaCuts,
        nPIDCriteriaCuts,
        nSigmaTPCCuts,
        nSigmaTOFCuts,
        nEventMixingPercent,
        nMultiplicityCuts,
        nSwarmSizeCuts);
  }

  ix =0;
  for (Int_t ixCentrality = 0; ixCentrality < nCentralityCuts; ixCentrality++) {
    if (verb) printf("   loop A\n");
    for (Int_t ixPt = 0; ixPt < nTaskPtCuts; ixPt++) {
      if (verb) printf("   loop B\n");
      for (Int_t ixEta = 0; ixEta < nTaskEtaCuts; ixEta++) {
        if (verb) printf("   loop C\n");
        for (Int_t ixVertex = 0; ixVertex < nVertexCuts; ixVertex++) {
          if (verb) printf("   loop D\n");
          for (Int_t ixTrack = 0; ixTrack < nTrackCriteriaCuts; ixTrack++) {
            if (verb) printf("   loop E\n");
            for (Int_t ixPID = 0; ixPID < nPIDCriteriaCuts; ixPID++) {
              if (verb) printf("   loop F\n");
              for (Int_t ixsTPC = 0; ixsTPC < nSigmaTPCCuts; ixsTPC++) {
                if (verb) printf("   loop G\n");
/*                for (Int_t ixsTOF = 0; ixsTOF < nSigmaTOFCuts; ixsTOF++) {
                  if (verb) printf("   loop H\n"); */
                  for (Int_t ixMix = 0; ixMix < nEventMixingPercent; ixMix++) {
                    if (verb) printf("   loop I\n");
                    for (Int_t ixMult = 0; ixMult < nMultiplicityCuts; ixMult++) {
                      if (verb) printf("   loop J\n");
                      for (Int_t ixSwarm = 0; ixSwarm < nSwarmSizeCuts; ixSwarm++) {
                        globalEtaMin[ix] = minEta[ixEta];
                        globalEtaMax[ix] = maxEta[ixEta];
                        globalPtMin[ix] = minPt[ixPt];
                        globalPtMax[ix] = maxPt[ixPt];
                        centralityMin[ix] = minCentrality[ixCentrality];
                        centralityMax[ix] = maxCentrality[ixCentrality];
                        zvertexMin[ix] = minzVertex[ixVertex];
                        zvertexMax[ix] = maxzVertex[ixVertex];
                        sigmasTPCMin[ix] = minsTPC[ixsTPC];
                        sigmasTPCMax[ix] = maxsTPC[ixsTPC];
                        sigmasTOFMin[ix] = minsTOF[ixsTPC];
                        sigmasTOFMax[ix] = maxsTOF[ixsTPC];
                        nTrackSelectionCriteria[ix] = trackSel[ixTrack];
                        nPIDCriteria[ix] = pidSel[ixPID];
                        fMixingPercentage[ix] = mixPerc[ixMix];
                        multiplicityMin[ix] = minMult[ixMult];
                        multiplicityMax[ix] = maxMult[ixMult];
                        minSwarmSize[ix] = swarmSize[ixSwarm];
                        if (verb) printf("task ix = %d\n", ix);
                        ix++;
                      }
                    }
                  }
/*                } */
              }
            }
          }
        }
      }
    }
  }

  /* closing the options file */
  if (verb) printf(" Closing the options file\n");
  optionsfile.close();
  if (verb) printf(" Options file closed\n");

  /* now get the data files location */
  if (bGRIDPlugin) {
    TString szDataLocFile;
    ifstream datalocfile;

    if (bMC) { szDataLocFile = "GRIDMCdata.txt"; }
    else { szDataLocFile = "GRIDrealdata.txt"; }

    if (verb) printf(" Opening the data location file: %s\n", szDataLocFile.Data());
    datalocfile.open(szDataLocFile.Data());
    if (datalocfile.fail()) { cout << "ERROR: cannot open location file " << szDataLocFile << endl; return kFALSE; }

    /* the data pattern */
    currline.ReadLine(datalocfile);
    while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(datalocfile);
    if (currline.BeginsWith("Data pattern: ")) {
      szDataPattern = currline.Remove(0,strlen("Data pattern: "));
    }
    else
      { cout << "ERROR: wrong Data pattern in data location file " << szDataLocFile << endl; return kFALSE; }
    cout << "  GRID data pattern: " << szDataPattern << endl;
    /* the grid data dir */
    currline.ReadLine(datalocfile);
    while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(datalocfile);
    if (currline.BeginsWith("Grid data dir: ")) {
      szDataDir = currline.Remove(0,strlen("Grid data dir: "));
    }
    else
      { cout << "ERROR: wrong GRID Data dir in data location file " << szDataLocFile << endl; return kFALSE; }
    cout << "  GRID data dir: " << szDataDir << endl;
    /* the number of input files */
    currline.ReadLine(datalocfile);
    while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(datalocfile);
    if (currline.BeginsWith("Number of input files: ")) {
      currline.Remove(0,strlen("Number of input files: "));
      nNoOfInputFiles = currline.Atoi();
    }
    else
      { cout << "ERROR: wrong number of input files in data location file " << szDataLocFile << endl; return kFALSE; }
    cout << "  GRID number of input files: " << nNoOfInputFiles << endl;
    /* the run prefix */
    currline.ReadLine(datalocfile);
    while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(datalocfile);
    if (currline.BeginsWith("Run prefix: ")) {
      szRunPrefix = currline.Remove(0,strlen("Run prefix: "));
    }
    else
      { cout << "ERROR: wrong Run prefix in data location file " << szDataLocFile << endl; return kFALSE; }
    if (szRunPrefix.IsWhitespace())
      printf("  NO run prefix\n");
    else
      cout << "  Run prefix: " << szRunPrefix << endl;

    /* now the list of runs to process */
    TObjString *pszRunNo;
    currline.ReadLine(datalocfile);
    while (!currline.IsWhitespace()) {
      if (!currline.BeginsWith("#")) {
        pszRunNo = new TObjString(currline);
        listOfRuns.Add(pszRunNo);
      }
      currline.ReadLine(datalocfile);
    }
    printf("  List of runs: ");
    for (Int_t i=0;i<listOfRuns.GetEntriesFast();i++) {
      if (i%8 == 0) {
        printf("\n  ");
      }
      printf("%s, ", ((TObjString*) listOfRuns.At(i))->GetString().Data());
    }
    printf("\n");
    datalocfile.close();
  }
  else {
    if (bMC) {
      szLocalFileList = "filelist_mc.txt";
    }
    else
      szLocalFileList = "filelist.txt";
    printf("  Local data list file: %s\n", (const char*) szLocalFileList);
  }
  bOptionsLoaded = kTRUE;
  return kTRUE;
}

void CleanOptions() {
  listOfRuns.Clear();
  delete globalEtaMin;
  delete globalEtaMax;
  delete globalPtMin;
  delete globalPtMax;
  delete centralityMin;
  delete centralityMax;
  delete zvertexMin;
  delete zvertexMax;
  delete sigmasTPCMin;
  delete sigmasTPCMax;
  delete sigmasTOFMin;
  delete sigmasTOFMax;
  delete nTrackSelectionCriteria;
  delete nPIDCriteria;
  delete fMixingPercentage;
  delete multiplicityMin;
  delete multiplicityMax;
  delete minSwarmSize;
  delete minCentrality;
  delete maxCentrality;
  delete minPt;
  delete maxPt;
  delete minEta;
  delete maxEta;
  delete minzVertex;
  delete maxzVertex;
  delete minsTPC;
  delete maxsTPC;
  delete minsTOF;
  delete maxsTOF;
  delete mixPerc;
  delete trackSel;
  delete pidSel;
  delete minMult;
  delete maxMult;
  delete swarmSize;
  delete swarmPhiMin;
  delete swarmPhiMax;
  delete swarmEtaMin;
  delete swarmEtaMax;
  delete swarmPtMin;
  delete swarmPtMax;

  delete deltaPt;
  delete noOfPtPoints;
  for (Int_t i = 0; i < nDeltasPt;i++)
    delete pointPt[i];
  delete pointPt;

  delete deltaEta;
  delete noOfEtaPoints;
  for (Int_t i = 0; i < nDeltasEta;i++)
    delete pointEta[i];
  delete pointEta;

  delete deltaPhi;
  delete noOfPhiPoints;
  for (Int_t i = 0; i < nDeltasPhi;i++)
    delete pointPhi[i];
  delete pointPhi;
}
