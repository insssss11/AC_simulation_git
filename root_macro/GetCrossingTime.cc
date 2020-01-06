#ifndef _GETCROSSINGTIME_
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TMath.h"

#define _GETCROSSINGTIME_
TGraph *GetVoltageFile(const char *fileName, const char *cut = "0")
{
    const Double_t impedence = 50;
    const Double_t gain = 1e7;
    auto file = TFile::Open(fileName, "READ");
    auto t = (TTree*)file->Get("tree");
    t->Draw("t_signal>>h_tmp(30, 10, 40)", TString::Format("evtID == %s", cut), "goff");
    auto h = (TH1D*)gDirectory->Get("h_tmp");
    Int_t nbinx = h->GetNbinsX();
    auto g = new TGraph(nbinx);
    for(Int_t i = 0;i < nbinx;i++)
        g->SetPoint(i, h->GetBinCenter(i), impedence*TMath::Qe()*gain*h->GetBinContent(i)*1e9);
    return g;
}

Double_t GetTimeGraph(TGraph *g, Double_t thresh = 0.25)
{
    Int_t np = g->GetN();
    Double_t x1,y1, x2, y2;
    Double_t myTime = 0.;
    for(Int_t i = 0;i < np-1;i++)
    {
        g->GetPoint(i, x1, y1);
        g->GetPoint(i, x2, y2);
        if((y1-thresh)*(y2-thresh) < 0)
        {
            myTime = (thresh - y1)/(y2 - y1);
            break;
        }
    }
    return myTime;
}
#endif

void GetCrossingTime()
{
    GetVoltageFile("./build/simulateAC.root")->Draw();
}