#include "GraphPosRun.hh"

int anal_kaon_abs(int rad, int emp)
{
  TGraphErrors *gkd[3];
  Int_t runno[] = {rad, emp};
  for(Int_t i = 0;i < 2;i++)
  {
    gkd[i] = GraphKaon1000(runno[i]);
    gkd[i]->SetMarkerStyle(kFullDotMedium);
    gkd[i]->SetLineColor(kBlack + i);gkd[i]->SetMarkerColor(kBlack + i);
  }
  gkd[2] = GraphElectron500(runno[1]);
  gkd[2]->SetMarkerStyle(kFullDotMedium);
  gkd[2]->SetLineColor(kBlack + 2);gkd[2]->SetMarkerColor(kBlack + 2);
  Double_t postest[] = {0, 50, 100, 150, 175, 200};
  Double_t meantest[] = {119.842, 112.384, 109.741, 131.642, 120.699, 107.633};

  printf("-----------------------------------------------------------------------------------\n");
  printf("(area under graph of test)/(area under graph of sim)\n");
  printf("TYPE A : \t%f \b\n", NormalizeADC(meantest, sizeof(postest)/sizeof(Double_t), gkd[2]->GetY(), gkd[2]->GetN()));
  printf("-----------------------------------------------------------------------------------\n");
  cout << meantest[0] << " " << meantest[1] << " " << meantest[2] << endl;
  gkd[3] = new TGraphErrors(sizeof(postest)/sizeof(Double_t), postest, meantest);
  gkd[3]->SetMarkerStyle(kFullDotMedium);
  gkd[3]->SetLineColor(kBlack + 3);gkd[3]->SetMarkerColor(kBlack + 3);
  TCanvas *c1 = new TCanvas("c1", "c1", 750, 750);c1->SetGrid();c1->SetLeftMargin(0.12);
  TMultiGraph *mg1 = new TMultiGraph();
  mg1->SetTitle("incident at 40 degrees(TYPEA);incident position(mm);Npe");
  mg1->Add(gkd[0]);mg1->Add(gkd[1]);mg1->Add(gkd[2]);mg1->Add(gkd[3]);
  auto lg1 = new TLegend(0.12, 0.78, 0.67, 0.9);
  mg1->GetHistogram()->GetYaxis()->SetRangeUser(0, 3.);
  lg1->AddEntry(gkd[0], "1000 MeV/c K+ with aerogel", "LP");
  lg1->AddEntry(gkd[1], "1000 MeV/c K+ without aerogel", "LP");
  lg1->AddEntry(gkd[2], "500 MeV/c e- without aerogel", "LP");
  lg1->AddEntry(gkd[3], "500 MeV/c e- without aerogel(measured, normalized)", "LP");
  mg1->Draw("APL");lg1->Draw();
  c1->Print("anal_kaon.svg", "svg");
  return 0;
}
