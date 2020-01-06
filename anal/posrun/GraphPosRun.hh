#ifndef GRAPHPOSRUN_HH
#define GRAPHPOSRUN_HH

const Int_t distMin__ = 0;
const Int_t distMax__ = 240;
const Int_t step__ = 10;
const Int_t nStep__ = (distMax__ - distMin__)/step__ + 1;
Double_t NormalizeADC(Double_t *arr1, Int_t size1, Double_t *arr2, Int_t size2)
{
  Double_t sum1 = 0;
  Double_t sum2 = 0;
  for(Int_t i = 0;i < size1;i++)
    sum1 += arr1[i];
  for(Int_t i = 0;i < size2;i++)
    sum2 += arr2[i];
  for(Int_t i = 0;i < size1;i++)
    arr1[i] = arr1[i]*sum2*size1/(sum1*size2);
  return sum2*size1/(sum1*size2);
}

void DivdeArrByArr(Double_t *nom, Double_t *den, Double_t *dst, Int_t size)
{
  for(Int_t i =0 ;i < size;i++)
    dst[i] = nom[i]/den[i];
}

void DivideArr(Double_t *nom, Double_t den, Double_t *dst, Int_t size)
{
  for(Int_t i =0 ;i < size;i++)
    dst[i] = nom[i]/den;
}
TGraphErrors* GraphPion1000(Int_t tNum, TString branch = "Npe")
{
  TString typeDir = TString::Format("../../result/posrun/TYPE%02d/", tNum);
  TString fileStr;
  TFile *file;
  TTree *t;
  TH1I *hist;
  Double_t pos[nStep__], mean[nStep__], error[nStep__];
  for(Int_t i = 0;i < nStep__;i++)
  {
    fileStr = TString::Format("posrun_pi+_%d_1000.root", step__*i);
    
    cout << fileStr << endl;
    file = TFile::Open(typeDir + fileStr, "READ");
    t = (TTree*)file->Get("tree");
    t->Draw(branch + ">>hist(100, 0, 100)", "", "goff");
    hist = (TH1I*)gDirectory->Get("hist");
    pos[i] = step__*i;
    mean[i] = hist->GetMean();
    error[i] = hist->GetStdDev()/TMath::Sqrt(t->GetEntries());
    file->Close();
  }
  auto gr = new TGraphErrors(nStep__, pos, mean, nullptr, error);
  return gr;
}

TGraphErrors* GraphKaon1000(Int_t tNum, TString branch = "Npe")
{
  TString typeDir = TString::Format("../../result/posrun/TYPE%02d/", tNum);
  TString fileStr;
  TFile *file;
  TTree *t;
  TH1I *hist;
  Double_t pos[nStep__], mean[nStep__], error[nStep__];
  for(Int_t i = 0;i < nStep__;i++)
  {
    fileStr = TString::Format("posrun_kaon+_%d_1000.root", step__*i);
    cout << fileStr << endl;
    file = TFile::Open(typeDir + fileStr, "READ");
    t = (TTree*)file->Get("tree");
    t->Draw(branch + ">>hist(100, 0, 100)", "", "goff");
    hist = (TH1I*)gDirectory->Get("hist");
    pos[i] = step__*i;
    mean[i] = hist->GetMean();
    error[i] = hist->GetStdDev()/TMath::Sqrt(t->GetEntries());
    file->Close();
  }
  auto gr = new TGraphErrors(nStep__, pos, mean, nullptr, error);
  return gr;
}

TGraphErrors* GraphElectron500(Int_t tNum, TString branch = "Npe")
{
  TString typeDir = TString::Format("../../result/posrun/TYPE%02d/", tNum);
  TString fileStr;
  TFile *file;
  TTree *t;
  TH1I *hist;
  Double_t pos[nStep__], mean[nStep__], error[nStep__];
  for(Int_t i = 0;i < nStep__;i++)
  {
    fileStr = TString::Format("posrun_e-_%d_500.root", step__*i);
    cout << fileStr << endl;
    file = TFile::Open(typeDir + fileStr, "READ");
    t = (TTree*)file->Get("tree");
    t->Draw(branch + ">>hist(100, 0, 100)", "", "goff");
    hist = (TH1I*)gDirectory->Get("hist");
    pos[i] = step__*i;
    mean[i] = hist->GetMean();
    error[i] = hist->GetStdDev()/TMath::Sqrt(t->GetEntries());
    file->Close();
  }
  auto gr = new TGraphErrors(nStep__, pos, mean, nullptr, error);
  return gr;
}

#endif
