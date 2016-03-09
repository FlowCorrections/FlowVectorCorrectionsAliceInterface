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

#include "runs.H"

#define NOOFTOTALRUNS 92
const char *sRunNo[NOOFTOTALRUNS] = {  "137161","137162","137231","137232","137235","137236","137243","137366","137431","137432","137434","137439",
    "137440","137441","137443","137530","137531","137539","137541","137544","137546","137549","137595","137608","137638","137639","137685","137686",
    "137691","137692","137693","137704","137718","137722","137724","137751","137752","137844","137848","138190","138192","138197","138201","138225",
    "138275","138364","138396","138438","138439","138442","138469","138534","138578","138579","138582","138583","138621","138624","138638","138652",
    "138653","138662","138666","138730","138732","138837","138870","138871","138872","139028","139029","139036","139037","139038","139105","139107",
    "139173","139309","139310","139314","139328","139329","139360","139437","139438","139465","139503","139505","139507","139510","169099","170040"};

TObjArray *runs = NULL;

void loadRuns() {
  runs = new TObjArray(NOOFTOTALRUNS);
  runs->SetOwner(kTRUE);

  for (Int_t i = 0; i < NOOFTOTALRUNS; i++) {
    runs->AddAt(new TObjString(sRunNo[i]),i);
  }
}

Double_t getRunIndex(const char *run) {
  return runs->IndexOf(runs->FindObject(run));
}

#define NOOFREASONS 8
const char *sReasons[NOOFREASONS] = { "norec", "diffrecwtruth", "notruth", "reason4", "reason5", "reason6", "reason7", "reason8" };

TObjArray *reasons = NULL;

void loadReasons() {
  reasons = new TObjArray(NOOFREASONS);
  reasons->SetOwner(kTRUE);

  for (Int_t i = 0; i < NOOFREASONS; i++) {
    reasons->AddAt(new TObjString(sReasons[i]),i);
  }
}

Double_t getReasonIndex(const char *reason) {
  return reasons->IndexOf(reasons->FindObject(reason));
}








