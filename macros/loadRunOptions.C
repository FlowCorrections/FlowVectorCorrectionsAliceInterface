/**************************************************************************
 * Copyright(c) 2013-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "TSystem.h"
#include "TROOT.h"

using namespace std;

#include "runAnalysis.H"


Bool_t bOptionsLoaded = kFALSE;
void CleanOptions();

/// \brief Load the run options for the current task
/// The run options as present in the input file are taken over the global variables
/// \param verb flag for verbose output
/// \param filename the options filename
Bool_t loadRunOptions(Bool_t verb,const char *filename) {
  ifstream optionsfile;
  TString currline;

  if (bOptionsLoaded) CleanOptions();

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
  if (currline.BeginsWith("AliPhysicsVersion: ")) {
    szAliPhysicsVersion = currline.Remove(0,strlen("AliPhysicsVersion: "));
  }
  else
    { printf("ERROR: wrong Aliphysics SW version in options file %s\n", filename); return kFALSE; }
  if (bGRIDPlugin)
    cout << "    Grid AliPhysics version: " << szAliPhysicsVersion << endl;

  /* Execution conditions */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("Use multiplicity: ")) {
    currline.Remove(0,strlen("Use multiplicity: "));
    if (currline.Contains("yes"))
      bUseMultiplicity = kTRUE;
    else if (currline.Contains("no"))
      bUseMultiplicity = kFALSE;
    else
      { printf("ERROR: wrong Use multiplicity option in options file %s\n", filename); return kFALSE; }
  }
  else
    { printf("ERROR: wrong Use multiplicity option in options file %s\n", filename); return kFALSE; }
  printf("  Use multiplicity: %s\n", bUseMultiplicity ? "yes" : "no");

  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("Is 2015 dataset")) {
    currline.Remove(0,strlen("Is 2015 dataset"));
    if (currline.Contains("yes"))
      b2015DataSet = kTRUE;
    else if (currline.Contains("no"))
      b2015DataSet = kFALSE;
    else
      { printf("ERROR: wrong Is 2015 dataset option in options file %s\n", filename); return kFALSE; }
  }
  else
    { printf("ERROR: wrong Is 2015 dataset option in options file %s\n", filename); return kFALSE; }
  printf("  Is 2015 dataset: %s\n", b2015DataSet ? "yes" : "no");

  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.EqualTo("Task level:")) {
    currline.ReadLine(optionsfile);
    while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
    while(!currline.EqualTo("end")) {
      /* centrality cuts */
      if (currline.BeginsWith("Centrality: ")) {
        currline.Remove(0, strlen("Centrality: "));
        if (verb) cout << "    Centrality cut string: " << currline << endl;
        Double_t min;
        Double_t max;
        sscanf(currline.Data(), "%lf-%lf", &min, &max);
        centralityMin = min;
        centralityMax = max;
        printf ("      centrality cut: %.1f-%.1f\n",centralityMin,centralityMax);
        currline.ReadLine(optionsfile);
        while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
      } /* end centrality cuts */

      /* z Vertex cuts */
      if (currline.BeginsWith("zvertex: ")) {
        currline.Remove(0, strlen("zvertex: "));
        if (verb) cout << "    z vertex cut string: " << currline << endl;
        Double_t min;
        Double_t max;
        sscanf(currline.Data(), "%lf:%lf", &min, &max);
        zvertexMin = min;
        zvertexMax = max;
        printf ("      z vertex cut: %.1f:%.1f\n",zvertexMin,zvertexMax);
        currline.ReadLine(optionsfile);
        while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
      } /* end centrality cuts */
    }
  }
  else
    { printf("ERROR: wrong Task cuts section in options file %s\n", filename); return kFALSE; }

  /* now the corrections file */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (!currline.EqualTo("Corrections file:")) { printf("ERROR: wrong corrections file location in options file %s\n", filename); return kFALSE; }
  printf(" Corrections file:\n");
  /* source option */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("source: ")) {
    currline.Remove(0,strlen("source: "));
    if (currline.Contains("local"))
      szCorrectionsSource = "local";
    else if (currline.Contains("alien"))
      szCorrectionsSource = "alien";
    else if (currline.Contains("OADB"))
      szCorrectionsSource = "OADB";
    else if (currline.Contains("OCDB"))
      szCorrectionsSource = "OCDB";
    else
      { printf("ERROR: wrong corrections file source in options file %s\n", filename); return kFALSE; }
  }
  else
    { printf("ERROR: wrong corrections file source in options file %s\n", filename); return kFALSE; }
  printf("  source: %s\n", (const char *) szCorrectionsSource);
  /* path option */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("path: ")) {
    currline.Remove(0,strlen("path: "));
    szCorrectionsFilePath = currline;
  }
  else
    { printf("ERROR: wrong corrections file path in options file %s\n", filename); return kFALSE; }
  printf("  path: %s\n", (const char *) szCorrectionsFilePath);

  /* file name option */
  currline.ReadLine(optionsfile);
  while (currline.BeginsWith("#") || currline.IsWhitespace()) currline.ReadLine(optionsfile);
  if (currline.BeginsWith("filename: ")) {
    currline.Remove(0,strlen("filename: "));
    szCorrectionsFileName = currline;
  }
  else
    { printf("ERROR: wrong corrections file naem in options file %s\n", filename); return kFALSE; }
  printf("  filename: %s\n", (const char *) szCorrectionsFileName);

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
    TObjString *pszActiveRunNo;
    currline.ReadLine(datalocfile);
    while (!currline.IsWhitespace()) {
      if (!currline.BeginsWith("#")) {
        pszRunNo = new TObjString(currline);
        if (pszRunNo->GetString().BeginsWith(">")) {
          /* the run is included in the current analysis */
          pszActiveRunNo = new TObjString(pszRunNo->GetString().Remove(0,1));
          listOfActiveRuns.Add(pszActiveRunNo);
          listOfRuns.Add(pszActiveRunNo->Clone());
          delete pszRunNo;
        }
        else {
          /* the run is included as a potential concurrent run */
          listOfRuns.Add(pszRunNo);
        }
      }
      currline.ReadLine(datalocfile);
    }
    printf("  List of concurrent runs: ");
    for (Int_t i=0;i<listOfRuns.GetEntriesFast();i++) {
      if (i%8 == 0) {
        printf("\n    ");
      }
      printf("%s, ", ((TObjString*) listOfRuns.At(i))->GetString().Data());
    }
    printf("\n");
    printf("  List of runs to analyze: ");
    for (Int_t i=0;i<listOfActiveRuns.GetEntriesFast();i++) {
      if (i%8 == 0) {
        printf("\n    ");
      }
      printf("%s, ", ((TObjString*) listOfActiveRuns.At(i))->GetString().Data());
    }
    printf("\n");
    if (verb) printf(" Closing the data location file: %s\n", szDataLocFile.Data());
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
  listOfActiveRuns.Clear();
}
