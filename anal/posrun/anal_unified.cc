#include "GraphPosRun.hh"

Int_t anal1()
{
  TGraph *ged[3]; // electron davis model
  TGraph *gpd[3]; // pion davis model
  TGraph *geu[3]; // electron unified model 
  TGraph *gpu[3]; // pion unified model
  Int_t dRun[] = {20, 21, 22};
  Int_t uRun[] = {34, 35, 36};

  for(Int_t i = 0;i < 3;i++)
  {
    ged[i] = GraphElectron1000(dRun[i]);
    ged[i]->SetMarkerStyle(kFullDotMedium);
    ged[i]->SetLineColor(kBlack + i);ged[i]->SetMarkerColor(kBlack + i);
    gpd[i] = GraphPion1000(dRun[i]);
    gpd[i]->SetMarkerStyle(kFullDotMedium);
    gpd[i]->SetLineColor(kBlack + i);gpd[i]->SetMarkerColor(kBlack + i);
    geu[i] = GraphElectron1000(uRun[i]);
    geu[i]->SetMarkerStyle(kFullDotMedium);
    geu[i]->SetLineColor(kBlack + i);geu[i]->SetMarkerColor(kBlack + i);
    gpu[i] = GraphPion1000(uRun[i]);
    gpu[i]->SetMarkerStyle(kFullDotMedium);
    gpu[i]->SetLineColor(kBlack + i);gpu[i]->SetMarkerColor(kBlack + i);
  }
  TMultiGraph *mged, *mgpd, *mgeu, *mgpu;
  mged->SetTitle("500 MeV/c electron incident at 40 degrees(DAVIS LUT model);incident position(mm);Npe");
  mgpd->SetTitle("1000 MeV/c pion incident at 40 degrees(DAVIS LUT model);incident position(mm);Npe");
  mgeu->SetTitle("500 MeV/c electron incident at 40 degrees(unified model);incident position(mm);Npe");
  mgpu->SetTitle("1000 MeV/c pion incident at 40 degrees(unified model);incident position(mm);Npe");
  mged->GetHistogram()->SetMaximum(60.);
  mgpd->GetHistogram()->SetMaximum(60.);
  mgeu->GetHistogram()->SetMaximum(60.);
  mgpu->GetHistogram()->SetMaximum(60.);
  mged->Add(ged[0]);mged->Add(ged[1]);mged->Add(ged[2]);
  mgpd->Add(gpd[0]);mgpd->Add(gpd[1]);mgpd->Add(gpd[2]);
  mgeu->Add(geu[0]);mgeu->Add(geu[1]);mgeu->Add(geu[2]);
  mgpu->Add(gpu[0]);mgpu->Add(gpu[1]);mgpu->Add(gpu[2]);
  TCanvas *c1 = new TCanvas("c1", "c1", 750, 750);c1->SetGrid();mged->Draw("APL");
  TCanvas *c2 = new TCanvas("c2", "c2", 750, 750);c2->SetGrid();mgpd->Draw("APL");
  TCanvas *c3 = new TCanvas("c3", "c3", 750, 750);c3->SetGrid();mgeu->Draw("APL");
  TCanvas *c4 = new TCanvas("c4", "c4", 750, 750);c4->SetGrid();mgpu->Draw("APL");
  return 0;
}

