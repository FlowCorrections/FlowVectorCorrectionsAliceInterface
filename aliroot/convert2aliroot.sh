#!/bin/bash
# this script prepares the FlowVectorCorrections package for insertion in AliROOT
# it prepends all file and class names with Ali
# $1 points to the source directory
# $2 points to the target directory


replacenames ()
{
  file=$1
  for j in $listclasses; do
    sed -i "s/QnCorrections${j}/AliQnCorrections${j}/g" $file
  done
  sed -i "s/QNCORRECTIONS/ALIQNCORRECTIONS/g" $file
}

replacenames1 ()
{
  file=$1
  for j in $listfiles1; do
    sed -i "s/AnalysisTask${j}/AliAnalysisTask${j}/g" $file
  done
  sed -i "s/ANALYSISTASK/ALIANALYSISTASK/g" $file
}

replacenames2 ()
{
  file=$1
  for j in $listfiles2; do
    sed -i "s/QnCorrections${j}/AliQnCorrections${j}/g" $file
  done
  sed -i "s/QNCORRECTIONS/ALIQNCORRECTIONS/g" $file
}

replacetracingnames ()
{
  file=$1
  for j in $tracingnames; do
    sed -i "s/QnCorrections${j}/Ali${j}/g" $file
  done
}


inputfolder=$1
outputfolder=$2

if [ $# -lt 1 ]; then
  inputfolder=../QnCorrectionsInterface
fi
if [ $# -lt 2 ]; then
  outputfolder=AliQnCorrectionsInterface
fi

rsync -av $inputfolder/ $outputfolder
rsync -av $inputfolder/../macros/ $outputfolder/macros

listclasses="CorrectionOnInputData
CorrectionOnQvector
CorrectionsSetOnInputData
CorrectionsSetOnQvector
CorrectionStepBase
CutAbove
CutsBase
CutBelow
CutSetBit
CutOutside
CutsSet
CutValue
CutWithin
DataVector
Detector
EventClassVariable
HistogramBase
HistogramChannelized
InputGainEqualization
Manager
Profile
QnVector"

listfiles1="FlowVectorCorrections 
QnVectorAnalysis"

listfiles2="FillEventTask
VarManagerTask"

# handle now the tracing fucntions
tracingnames="Log Info Warning Error Fatal"

for j in $listfiles1; do
  mv $outputfolder/AnalysisTask${j}.cxx $outputfolder/AliAnalysisTask${j}.cxx
  mv $outputfolder/AnalysisTask${j}.h $outputfolder/AliAnalysisTask${j}.h
done

for j in $listfiles2; do
  mv $outputfolder/QnCorrections${j}.cxx $outputfolder/AliQnCorrections${j}.cxx
  mv $outputfolder/QnCorrections${j}.h $outputfolder/AliQnCorrections${j}.h
done

listC=`find $outputfolder/ -name '*.cxx'`
listH=`find $outputfolder/ -name '*.h'`
listTXT=`find $outputfolder/ -name '*.txt'`
listMacros="$outputfolder/macros/AddTaskFlowQnVectorCorrections.C 
$outputfolder/macros/AddTaskQnVectorAnalysis.C 
$outputfolder/macros/runAnalysis.C"

for k in $listC; do
  replacenames "$k"
  replacenames1 "$k"
  replacenames2 "$k"
  replacetracingnames "$k"
done
for k in $listH; do
  replacenames "$k"
  replacenames1 "$k"
  replacenames2 "$k"
  replacetracingnames "$k"
done
for k in $listTXT; do
  replacenames1 "$k"
  replacenames2 "$k"
done
for k in $listMacros; do
  replacenames "$k"
  replacenames1 "$k"
  replacenames2 "$k"
  replacetracingnames "$k"
done


