// read root file made by simulation program and tree with sorted data

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

inline Float_t GetArrayMin(const Float_t *arr, Int_t size)
{
  Double_t min = arr[0];
  for(Int_t i = 1;i < size;i++)
    min = min>arr[i]?arr[i]:min;
  return min;
}

inline Float_t GetArrayMax(const Float_t *arr, Int_t size)
{
  Double_t max = arr[0];
  for(Int_t i = 1;i < size;i++)
    max = max<arr[i]?arr[i]:max;
  return max;
}

inline Float_t GetArrayRms(const Float_t *arr, Int_t size)
{
  Double_t ss = 0;
  for(Int_t i = 0;i < size;i++)
    ss += arr[i]*arr[i];
  return TMath::Sqrt(ss/size);
}

inline Float_t GetArrayMean(const Float_t *arr, Int_t size)
{
  Double_t s = 0;
  for(Int_t i = 0;i < size;i++)
    s += arr[i];
  return s/size;
}

int main(int argc, char *argv[])
{
  auto file = TFile::Open(argv[1], "UPDATE");
  auto t0 = (TTree*)file->Get("t0");
  auto t1 = (TTree*)file->Get("t1");
  auto tree = new TTree("tree", "creator -> 0 : cherenkov, 1 : scintillation");
  Float_t wavelen[300], tSignal[300], tFirst, tLast, tRms, tMean;
  Int_t creator[300];
  Int_t evtID, Nscint, Ncheren, Ncol, Npe;

  t0->SetBranchAddress("evtID", &evtID);
  t0->SetBranchAddress("Ncheren", &Ncheren);
  t0->SetBranchAddress("Nscint", &Nscint);
  t0->SetBranchAddress("Ncol", &Ncol);
  t0->SetBranchAddress("Npe", &Npe);

  tree->Branch("evtID", &evtID, "evtID/I");
  tree->Branch("Ncheren", &Ncheren, "Ncheren/I");
  tree->Branch("Nscint", &Nscint, "Nscint/I");
  tree->Branch("Ncol", &Ncol, "Ncol/I");
  tree->Branch("Npe", &Npe, "Npe/I");
  tree->Branch("wavelen", wavelen, "wavelen[Npe]/F");
  tree->Branch("t_signal", tSignal, "t_signal[Npe]/F");
  tree->Branch("t_first", &tFirst, "t_first/F");
  tree->Branch("t_last", &tLast, "t_last/F");
  tree->Branch("t_rms", &tRms, "t_rms/F");
  tree->Branch("t_mean", &tMean, "t_mean/F");
  tree->Branch("creator", creator, "creator[Npe]/I");

  Float_t wavelength_old, tSignal_old;
  Int_t creator_old;

  t1->SetBranchAddress("wavelen", &wavelength_old);
  t1->SetBranchAddress("t_signal", &tSignal_old);
  t1->SetBranchAddress("creator", &creator_old);
  
  Int_t cnt = 0;  
  for(Long64_t i = 0;i < t0->GetEntries();i++)
  {
    t0->GetEntry(i);
    for(Int_t j = 0;j < Npe;j++)
    {
      t1->GetEntry(cnt++);
      wavelen[j] = wavelength_old;
      tSignal[j] = tSignal_old;
      creator[j] = creator_old;
    }
    tFirst = GetArrayMin(tSignal, Npe);
    tLast = GetArrayMax(tSignal, Npe);
    tRms = GetArrayRms(tSignal, Npe);
    tMean = GetArrayMean(tSignal, Npe);
    tree->Fill();
  }
  gDirectory->Delete("t1;1");
  gDirectory->Delete("t0;1");
  tree->Write();
  file->Close();
  return 0;
}