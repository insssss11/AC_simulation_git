// written by Hyunmin Yang, Korea University, HANUL
// this calculated mean free path using simulation result

void anal_mfp()
{
  TString rootDir = "../../result/testbench/";
  TFile *file;
  TTree *tree;
  TH1I *hist;
  Double_t mean[4];
  
  for(Int_t i = 0;i < 4;i++)
  {
    file = TFile::Open(rootDir + TString::Format("setup%d.root", i + 1));
    tree = (TTree*)file->Get("tree");
    tree->Draw("Ncol>>hist(10000, 0, 10000)", "", "goff");
    hist = (TH1I*)gDirectory->Get("hist");
    mean[i] = hist->GetMean();
    file->Close();
  }
  cout << "------------------------------------------------------------" << endl;
  cout << "# of photons : " << mean[3] << endl;
  cout << "---------------------Mean Free Path-------------------------" << endl;
  cout << "Aerogel - Aerogel : " << -40./TMath::Log(mean[0]/mean[3]) << " mm" << endl;
  cout << "Aerogel - empty   : " << -20./TMath::Log(mean[1]/mean[3]) << " mm" << endl;
  cout << "empty   - Aerogel : " << -20./TMath::Log(mean[2]/mean[3]) << " mm" << endl;
  cout << "------------------------------------------------------------" << endl;
}