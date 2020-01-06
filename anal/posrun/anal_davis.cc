// written by Hyunmin Yang, Korea University, HANUL
// this file draw position vs Npe
// each pair of measured light collection efficiency and simulated Npe(LUTDAVIS)
// and also compare simulation result with beam test data

#include "GraphPosRun.hh"

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

Int_t anal_davis()
{
  TGraphErrors *ged[3]; // electron davis model
  TGraphErrors *gpd[3]; // pion davis model

  const Int_t dRun[] = {20, 21, 22};

  for(Int_t i = 0;i < 3;i++)
  {
    ged[i] = GraphElectron500(dRun[i]);
    ged[i]->SetMarkerStyle(kFullDotMedium);
    ged[i]->SetLineColor(kBlack + i);ged[i]->SetMarkerColor(kBlack + i);
    gpd[i] = GraphPion1000(dRun[i]);
    gpd[i]->SetMarkerStyle(kFullDotMedium);
    gpd[i]->SetLineColor(kBlack + i);gpd[i]->SetMarkerColor(kBlack + i);
  }
  
  TMultiGraph *mged, *mgpd;
  mged = new TMultiGraph();
  mgpd = new TMultiGraph();
  mged->SetTitle("500 MeV/c electron incident at 40 degrees(DAVIS LUT model);incident position(mm);Npe");
  mgpd->SetTitle("1000 MeV/c pion incident at 40 degrees(DAVIS LUT model);incident position(mm);Npe");
  mged->Add(ged[0]);mged->Add(ged[1]);mged->Add(ged[2]);
  mgpd->Add(gpd[0]);mgpd->Add(gpd[1]);mgpd->Add(gpd[2]);
  
  auto lg1 = new TLegend(0.1, 0.78, 0.45, 0.9);
  lg1->AddEntry(ged[0], "TYPEA", "LP"); lg1->AddEntry(ged[1], "TYPEB", "LP");lg1->AddEntry(ged[2], "TYPEC", "LP");
  auto lg2 = new TLegend(0.1, 0.78, 0.45, 0.9);
  lg2->AddEntry(gpd[0], "TYPEA", "LP"); lg2->AddEntry(gpd[1], "TYPEB", "LP");lg2->AddEntry(gpd[2], "TYPEC", "LP");    
  
  TCanvas *c1 = new TCanvas("c1", "c1", 750, 750);c1->SetGrid();mged->Draw("APL");lg1->Draw();
  mged->GetHistogram()->GetYaxis()->SetRangeUser(0., 45.);
  c1->Print("anal_davis_c1.svg", "svg");
  
  TCanvas *c2 = new TCanvas("c2", "c2", 750, 750);c2->SetGrid();mgpd->Draw("APL");lg2->Draw();
  mgpd->GetHistogram()->GetYaxis()->SetRangeUser(0., 45.);
  c2->Print("anal_davis_c2.svg", "svg");
  
  // comparing with test performed at 2019, Jul
  TGraph *get[3]; // electron test model
  TGraph *ger[3]; // ratio of type A and type C
  
  // position and ADC value
  Double_t posetA[] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 210};
  Double_t meanetA[] = {321.972, 375.979, 370.348, 336.28, 368.619, 410.583, 423.365, 487.766, 504.752, 604.2,  572.406, 504.421};
  Double_t posetC[] = {0, 20, 50, 80, 100, 120, 150, 175, 200};
  Double_t meanetC[] = {291.394, 299.957, 323.364, 365.941, 398.873, 441.689, 412.698, 428.3, 397.247};

  printf("-----------------------------------------------------------------------------------\n");
  printf("(area under graph of test)/(area under graph of sim)\n");
  printf("TYPE A : \t%f \b\n", NormalizeADC(meanetA, sizeof(posetA)/sizeof(Double_t), ged[0]->GetY(), ged[0]->GetN()));
  printf("TYPE C : \t%f \b\n", NormalizeADC(meanetC, sizeof(posetC)/sizeof(Double_t), ged[2]->GetY(), ged[2]->GetN()));
  printf("-----------------------------------------------------------------------------------\n");
  get[0] = new TGraph(sizeof(posetA)/sizeof(Double_t), posetA, meanetA);
  get[2] = new TGraph(sizeof(posetC)/sizeof(Double_t), posetC, meanetC);

  ged[0]->SetLineColor(kBlack);ged[0]->SetMarkerColor(kBlack);
  ged[2]->SetLineColor(kBlack);ged[2]->SetMarkerColor(kBlack);

  get[0]->SetLineColor(kRed);get[0]->SetMarkerColor(kRed);ged[0]->SetMarkerStyle(kFullDotMedium);
  get[2]->SetLineColor(kRed);get[2]->SetMarkerColor(kRed);ged[2]->SetMarkerStyle(kFullDotMedium);
  TMultiGraph *mgedt[3];mgedt[0] = new TMultiGraph();mgedt[2] = new TMultiGraph();
  mgedt[0]->SetTitle("500 MeV/c electron incident at 40 degrees(TYPEA);incident position(mm);Npe");
  mgedt[2]->SetTitle("500 MeV/c electron incident at 40 degrees(TYPEC);incident position(mm);Npe");  
  mgedt[0]->Add(ged[0]);mgedt[0]->Add(get[0]);
  mgedt[2]->Add(ged[2]);mgedt[2]->Add(get[2]);
  auto lg3 = new TLegend(0.1, 0.78, 0.45, 0.9);
  lg3->AddEntry(ged[0], "sim(DAVIS)", "LP"); lg3->AddEntry(get[0], "beam test(normalized)", "LP");
  auto lg4 = new TLegend(0.1, 0.78, 0.45, 0.9);
  lg4->AddEntry(ged[2], "sim(DAVIS)", "LP"); lg4->AddEntry(get[2], "beam test(normalized)", "LP");
  
  TCanvas *c3 = new TCanvas("c3", "c3", 750, 750);c3->SetGrid();mgedt[0]->Draw("APL");lg3->Draw();
  mgedt[0]->GetHistogram()->GetYaxis()->SetRangeUser(0., 45.);
  c3->Print("anal_davis_c3.svg", "svg");
  
  TCanvas *c4 = new TCanvas("c4", "c4", 750, 750);c4->SetGrid();mgedt[2]->Draw("APL");lg4->Draw();
  mgedt[2]->GetHistogram()->GetYaxis()->SetRangeUser(0., 45.);
  c4->Print("anal_davis_c4.svg", "svg");
  return 0;
}

