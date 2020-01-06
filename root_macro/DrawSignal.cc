#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"

class  SignalFunc {
  public:
  SignalFunc(Float_t *signalTime, const Int_t Npe):fNpe(Npe)
  {
    fSignalTime = new Double_t[Npe];
    for(Int_t i = 0;i < fNpe;i++)
    {
      fSignalTime[i] = (Double_t)signalTime[i];
      // cout << fSignalTime[i] << endl;
    }
    fSigma = TMath::Sqrt((-fA+TMath::Sqrt(fA*fA-fB*fB))/(fB*fB))*fK;
    cout << fSigma << endl;
  }

  ~SignalFunc(){
    // delete[] fSignalTime;
  }
  Double_t operator() (Double_t *x, Double_t *p);
  private:
  const Int_t fNpe;
  const Double_t fTrise = 2.9;
  const Double_t fA = TMath::Log(0.09);
  const Double_t fB = TMath::Log(9);
  const Double_t fK = TMath::Log(2.5);
  Double_t fSigma;
  Double_t *fSignalTime;
};

Double_t SignalFunc::operator() (Double_t *x, Double_t *p)
{
  Double_t xin = x[0];
  Double_t xout = 0;
  cout << TMath::Exp(-fSigma*fSigma) << endl;
  for(Int_t i = 0;i < 1;i++)
  {
    if(xin < fSignalTime[i] - TMath::Exp(-fSigma*fSigma))
      xout += 0;
    else
      xout += TMath::LogNormal(xin + fSignalTime[i] - TMath::Exp(-fSigma*fSigma), fSigma, 0., 1e9);
  }
  // Npe * electron charge * gain * impedence
  return fNpe*TMath::Qe()*1e7*1e9*50.*xout;
}

Int_t DrawSignal(Int_t evt = 0)
{
  auto file = TFile::Open("./build/simulateAC.root", "READ");
  auto tree = (TTree*)file->Get("tree");
  Int_t Npe;
  Float_t tSignal[100];
  tree->SetBranchAddress("Npe", &Npe);
  tree->SetBranchAddress("t_signal", tSignal);
  tree->GetEntry(evt);
  auto signal = new SignalFunc(tSignal, Npe);
  cout << Npe << endl;
  cout << "hello " << endl;
  auto f1 = new TF1("f1", *signal, 0., 300., 0);
  f1->Eval(20.);
  f1->Draw();
  return 0;
}