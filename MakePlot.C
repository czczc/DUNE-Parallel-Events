#include "MakePlot.h"
#include "FastMC.h"

#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TString.h"
#include "TProfile.h"
#include "TCut.h"

#include <iostream>

using namespace std;

MakePlot::MakePlot()
{
    fm = new FastMC();
    SetStyle();
}

MakePlot::~MakePlot()
{
    delete fm;
}

void MakePlot::PlotEv_vs_Angle(const char* saveAs)
{
    int nTotal = fm->fChain->Draw("","cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    c1->Divide(3,1);
    c1->cd(1);
    fm->fChain->Draw("Ev:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCQE(100,0,50,200,0,20)","cc&&qel","colz");
    TH2F *h = (TH2F*)gDirectory->Get("EvsAngleCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Neutrino Energy (GeV)");

    c1->cd(2);
    fm->fChain->Draw("Ev:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCRes(100,0,50,200,0,20)","cc&&res","colz");
    h = (TH2F*)gDirectory->Get("EvsAngleCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Neutrino Energy (GeV)");

    c1->cd(3);
    fm->fChain->Draw("Ev:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCDis(100,0,50,200,0,20)","cc&&dis","colz");
    h = (TH2F*)gDirectory->Get("EvsAngleCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Neutrino Energy (GeV)");

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEl_vs_Angle(const char* saveAs)
{
    int nTotal = fm->fChain->Draw("","cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    c1->Divide(3,1);
    c1->cd(1);
    fm->fChain->Draw("El:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCQE(100,0,50,200,0,20)","cc&&qel","colz");
    TH2F *h = (TH2F*)gDirectory->Get("EvsAngleCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    c1->cd(2);
    fm->fChain->Draw("El:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCRes(100,0,50,200,0,20)","cc&&res","colz");
    h = (TH2F*)gDirectory->Get("EvsAngleCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    c1->cd(3);
    fm->fChain->Draw("El:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCDis(100,0,50,200,0,20)","cc&&dis","colz");
    h = (TH2F*)gDirectory->Get("EvsAngleCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEl_vs_Length(const char* saveAs)
{
    int nTotal = fm->fChain->Draw("","cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    c1->Divide(3,1);
    c1->cd(1);
    fm->fChain->Draw("El:0.4/sin(abs(asin(pxl/pl)))>>El_vs_LengthCCQE(100,0,20,200,0,20)","cc&&qel","colz");
    TH2F *h = (TH2F*)gDirectory->Get("El_vs_LengthCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    c1->cd(2);
    fm->fChain->Draw("El:0.4/sin(abs(asin(pxl/pl)))>>El_vs_LengthCCRes(100,0,20,200,0,20)","cc&&res","colz");
    h = (TH2F*)gDirectory->Get("El_vs_LengthCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    c1->cd(3);
    fm->fChain->Draw("El:0.4/sin(abs(asin(pxl/pl)))>>El_vs_LengthCCDis(100,0,20,200,0,20)","cc&&dis","colz");
    h = (TH2F*)gDirectory->Get("El_vs_LengthCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEv(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",800, 600);

    fm->fChain->Draw("Ev>>Ev(400,0.1,40)");
    TH1F *h = (TH1F*)gDirectory->Get("Ev");
    h->SetTitle("");
    h->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    gPad->SetLogy();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEl_vs_Ev(const char* saveAs)
{
    int nTotal = fm->fChain->Draw("","cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    fm->fChain->Draw("El:Ev>>EvElCCQE(200,0,20,200,0,20)","cc&&qel","colz");
    TH2F *h = (TH2F*)gDirectory->Get("EvElCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetYaxis()->SetTitle("Lepton energy (GeV)");
    h->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    c1->cd(4);
    TProfile *px = h->ProfileX("pxCCQE");
    TH1D *pxp = px->ProjectionX("pxpCCQE");
    for (int i=1; i<=200; i++) {
        pxp->SetBinContent(i, pxp->GetBinContent(i)/i*10) ;
        pxp->SetBinError(i, pxp->GetBinError(i)/i*10);
    }
    pxp->GetYaxis()->SetTitle("Lepton Energy / Neutrino Energy");
    pxp->Draw();

    c1->cd(2);
    fm->fChain->Draw("El:Ev>>EvElCCRes(200,0,20,200,0,20)","cc&&res","colz");
    h = (TH2F*)gDirectory->Get("EvElCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetYaxis()->SetTitle("Lepton energy (GeV)");
    h->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    c1->cd(5);
    px = h->ProfileX("pxCCRes");
    pxp = px->ProjectionX("pxpCCRes");
    for (int i=1; i<=200; i++) {
        pxp->SetBinContent(i, pxp->GetBinContent(i)/i*10) ;
        pxp->SetBinError(i, pxp->GetBinError(i)/i*10);
    }
    pxp->GetYaxis()->SetTitle("Lepton Energy / Neutrino Energy");
    pxp->Draw();

    c1->cd(3);
    fm->fChain->Draw("El:Ev>>EvElCCDis(200,0,20,200,0,20)","cc&&dis","colz");
    h = (TH2F*)gDirectory->Get("EvElCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetYaxis()->SetTitle("Lepton energy (GeV)");
    h->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    c1->cd(6);
    px = h->ProfileX("pxCCDis");
    pxp = px->ProjectionX("pxpCCDis");
    for (int i=1; i<=200; i++) {
        pxp->SetBinContent(i, pxp->GetBinContent(i)/i*10) ;
        pxp->SetBinError(i, pxp->GetBinError(i)/i*10);
    }
    pxp->GetYaxis()->SetTitle("Lepton Energy / Neutrino Energy");
    pxp->Draw();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotAngle_cut_E(const char* saveAs)
{
    TCut ecut("El<6");

    int nTotal = fm->fChain->Draw("", ecut && "cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCQE(100,0,50)", ecut && "cc&&qel");
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180>>AngleCCQE(180,0,90)", ecut && "cc&&qel");
    TH1F *h = (TH1F*)gDirectory->Get("AngleCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(4);
    TH1* hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;

    c1->cd(2);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCRes(100,0,50)", ecut && "cc&&res");
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180>>AngleCCRes(180,0,90)", ecut && "cc&&res");
    h = (TH1F*)gDirectory->Get("AngleCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;

    c1->cd(3);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCDis(100,0,50)", ecut && "cc&&dis");
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180>>AngleCCDis(180,0,90)", ecut && "cc&&dis");
    h = (TH1F*)gDirectory->Get("AngleCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotAngle_cut_E_Perp(const char* saveAs)
{
    TCut ecut("El<6");

    int nTotal = fm->fChain->Draw("", ecut && "cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCQE(100,0,50)", ecut && "cc&&qel");
    fm->fChain->Draw("abs(asin(pzl/pl))/3.14159*180>>AngleCCQE(180,0,90)", ecut && "cc&&qel");
    TH1F *h = (TH1F*)gDirectory->Get("AngleCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(4);
    TH1* hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;

    c1->cd(2);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCRes(100,0,50)", ecut && "cc&&res");
    fm->fChain->Draw("abs(asin(pzl/pl))/3.14159*180>>AngleCCRes(180,0,90)", ecut && "cc&&res");
    h = (TH1F*)gDirectory->Get("AngleCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;

    c1->cd(3);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCDis(100,0,50)", ecut && "cc&&dis");
    fm->fChain->Draw("abs(asin(pzl/pl))/3.14159*180>>AngleCCDis(180,0,90)", ecut && "cc&&dis");
    h = (TH1F*)gDirectory->Get("AngleCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotLength_cut_E(const char* saveAs)
{
    TCut ecut("El<6");
    TH1F *h = 0;
    TH1 *hc = 0;

    int nTotal = fm->fChain->Draw("", ecut && "cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    fm->fChain->Draw("0.4/sin(abs(asin(pxl/pl)))>>LengthCCQE(100,0,20)", ecut && "cc&&qel");
    h = (TH1F*)gDirectory->Get("LengthCCQE");
    cout << "average length: " << h->GetMean() << endl;
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    c1->cd(4);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,101));
    for (int i=1; i<=100; i++) {
        hc->SetBinContent(i, 1-hc->GetBinContent(i));
    }
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(3*5) << ", " << hc->GetBinContent(6*5) << ", " << hc->GetBinContent(10*5) << endl;

    c1->cd(2);
    fm->fChain->Draw("0.4/sin(abs(asin(pxl/pl)))>>LengthCCRes(100,0,20)", ecut && "cc&&res");
    h = (TH1F*)gDirectory->Get("LengthCCRes");
    cout << "average length: " << h->GetMean() << endl;
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,101));
    for (int i=1; i<=100; i++) {
        hc->SetBinContent(i, 1-hc->GetBinContent(i));
    }
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(3*5) << ", " << hc->GetBinContent(6*5) << ", " << hc->GetBinContent(10*5) << endl;

    c1->cd(3);
    fm->fChain->Draw("0.4/sin(abs(asin(pxl/pl)))>>LengthCCDis(100,0,20)", ecut && "cc&&dis");
    h = (TH1F*)gDirectory->Get("LengthCCDis");
    cout << "average length: " << h->GetMean() << endl;
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,101));
    for (int i=1; i<=100; i++) {
        hc->SetBinContent(i, 1-hc->GetBinContent(i));
    }
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(3*5) << ", " << hc->GetBinContent(6*5) << ", " << hc->GetBinContent(10*5) << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotLength_cut_E_Perp(const char* saveAs)
{
    TCut ecut("El<6");
    TH1F *h = 0;
    TH1 *hc = 0;

    int nTotal = fm->fChain->Draw("", ecut && "cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    fm->fChain->Draw("0.4/sin(abs(asin(pzl/pl)))>>LengthCCQE(100,0,20)", ecut && "cc&&qel");
    h = (TH1F*)gDirectory->Get("LengthCCQE");
    cout << "average length: " << h->GetMean() << endl;
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    c1->cd(4);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,101));
    for (int i=1; i<=100; i++) {
        hc->SetBinContent(i, 1-hc->GetBinContent(i));
    }
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(3*5) << ", " << hc->GetBinContent(6*5) << ", " << hc->GetBinContent(10*5) << endl;

    c1->cd(2);
    fm->fChain->Draw("0.4/sin(abs(asin(pzl/pl)))>>LengthCCRes(100,0,20)", ecut && "cc&&res");
    h = (TH1F*)gDirectory->Get("LengthCCRes");
    cout << "average length: " << h->GetMean() << endl;
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,101));
    for (int i=1; i<=100; i++) {
        hc->SetBinContent(i, 1-hc->GetBinContent(i));
    }
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(3*5) << ", " << hc->GetBinContent(6*5) << ", " << hc->GetBinContent(10*5) << endl;

    c1->cd(3);
    fm->fChain->Draw("0.4/sin(abs(asin(pzl/pl)))>>LengthCCDis(100,0,20)", ecut && "cc&&dis");
    h = (TH1F*)gDirectory->Get("LengthCCDis");
    cout << "average length: " << h->GetMean() << endl;
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,101));
    for (int i=1; i<=100; i++) {
        hc->SetBinContent(i, 1-hc->GetBinContent(i));
    }
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(3*5) << ", " << hc->GetBinContent(6*5) << ", " << hc->GetBinContent(10*5) << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEv_vs_Ef(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    fm->fChain->Draw("Ev:Ef>>EvEfCCQE(50,0,5,200,0,20)","cc&&qel&&pf>0","colz");
    TH2F *h = (TH2F*)gDirectory->Get("EvEfCCQE");
    h->SetTitle("CC-QE");
    h->GetYaxis()->SetTitle("Neutrino energy (GeV)");
    h->GetXaxis()->SetTitle("Hadron Energy (GeV)");
    gPad->SetLogz();
    c1->cd(4);
    h->ProjectionX()->Draw();


    c1->cd(2);
    fm->fChain->Draw("Ev:Ef>>EvEfCCRes(50,0,5,200,0,20)","cc&&res&&pf>0","colz");
    h = (TH2F*)gDirectory->Get("EvEfCCRes");
    h->SetTitle("CC-Res");
    h->GetYaxis()->SetTitle("Neutrino energy (GeV)");
    h->GetXaxis()->SetTitle("Hadron Energy (GeV)");
    gPad->SetLogz();
    c1->cd(5);
    h->ProjectionX()->Draw();

    c1->cd(3);
    fm->fChain->Draw("Ev:Ef>>EvEfCCDis(50,0,5,200,0,20)","cc&&dis&&pf>0","colz");
    h = (TH2F*)gDirectory->Get("EvEfCCDis");
    h->SetTitle("CC-Dis");
    h->GetYaxis()->SetTitle("Neutrino energy (GeV)");
    h->GetXaxis()->SetTitle("Hadron Energy (GeV)");
    gPad->SetLogz();
    c1->cd(6);
    h->ProjectionX()->Draw();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotNf_vs_Ratio(const char* particle, const char* title, const char* saveAs)
{
    TH2F *h = 0;
    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    c1->Divide(3,1);
    c1->cd(1);
    fm->fChain->Draw(TString::Format("nf:(%s)/nf>>Nf_RatioCCQE(51,0,1.02,51,0,51)", particle), "cc&&qel","colz");
    h = (TH2F*)gDirectory->Get("Nf_RatioCCQE");
    h->SetTitle("CC-QE");
    h->GetYaxis()->SetTitle("# Final state particles (hadron side)");
    h->GetXaxis()->SetTitle(title);
    gPad->SetLogz();

    c1->cd(2);
    fm->fChain->Draw(TString::Format("nf:(%s)/nf>>Nf_RatioCCRes(51,0,1.02,51,0,51)", particle),"cc&&res","colz");
    h = (TH2F*)gDirectory->Get("Nf_RatioCCRes");
    h->SetTitle("CC-Res");
    h->GetYaxis()->SetTitle("# Final state particles (hadron side)");
    h->GetXaxis()->SetTitle(title);
    gPad->SetLogz();

    c1->cd(3);
    fm->fChain->Draw(TString::Format("nf:(%s)/nf>>Nf_RatioCCDis(51,0,1.02,51,0,51)", particle),"cc&&res","colz");
    h = (TH2F*)gDirectory->Get("Nf_RatioCCDis");
    h->SetTitle("CC-Dis");
    h->GetYaxis()->SetTitle("# Final state particles (hadron side)");
    h->GetXaxis()->SetTitle(title);
    gPad->SetLogz();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotAnglel_vs_Angleh(const char* saveAs)
{
    TH2F *h = 0;
    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180:abs(asin(pxf/pf))/3.14159*180>>Angleh_vs_AnglelCCQE(180,0,90,180,0,90)","cc&&qel&&El<6&&pf>0&&pdgf!=2112","colz");
    h = (TH2F*)gDirectory->Get("Angleh_vs_AnglelCCQE");
    h->SetTitle("CC-QE");
    h->GetXaxis()->SetTitle("hadron angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("lepton angle to wire plane (deg)");
    gPad->SetLogz();
    c1->cd(4);
    h->ProjectionX()->Draw();

    c1->cd(2);
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180:abs(asin(pxf/pf))/3.14159*180>>Angleh_vs_AnglelCCRes(180,0,90,180,0,90)","cc&&res&&El<6&&pf>0&&pdgf!=2112","colz");
    h = (TH2F*)gDirectory->Get("Angleh_vs_AnglelCCRes");
    h->SetTitle("CC-Res");
    h->GetXaxis()->SetTitle("hadron angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("lepton angle to wire plane (deg)");
    gPad->SetLogz();
    c1->cd(5);
    h->ProjectionX()->Draw();

    c1->cd(3);
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180:abs(asin(pxf/pf))/3.14159*180>>Angleh_vs_AnglelCCDis(180,0,90,180,0,90)","cc&&dis&&El<6&&pf>0&&pdgf!=2112","colz");
    h = (TH2F*)gDirectory->Get("Angleh_vs_AnglelCCDis");
    h->SetTitle("CC-Dis");
    h->GetXaxis()->SetTitle("hadron angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("lepton angle to wire plane (deg)");
    gPad->SetLogz();
    c1->cd(6);
    h->ProjectionX()->Draw();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::CreateTree(const char* filename)
{
    TFile f(filename,"recreate");
    TTree *T = new TTree("T", "New Tree from Fast MC");
    float h_min_angle = 91;
    int h_min_angle_pdg = 0;
    float h_min_angle_Ef = 0;
    // float h_sec_min_angle = 91;
    float l_angle = 91;

    T->Branch("h_min_angle", &h_min_angle);
    T->Branch("h_min_angle_pdg", &h_min_angle_pdg);
    T->Branch("h_min_angle_Ef", &h_min_angle_Ef);
    // T->Branch("h_sec_min_angle", &h_sec_min_angle);
    T->Branch("l_angle", &l_angle);
    T->Branch("cc", &(fm->cc) );
    T->Branch("qel", &(fm->qel) );
    T->Branch("res", &(fm->res) );
    T->Branch("dis", &(fm->dis) );
    T->Branch("El", &(fm->El));
    T->Branch("nf", &(fm->nf));

    int nEntries = fm->fChain->GetEntries();
    for (int i=0; i<nEntries; i++) {
        fm->fChain->GetEntry(i);

        l_angle = abs(asin(fm->pxl/fm->pl))/3.14159*180;

        for (int j=0; j<fm->nf; j++) {
            if (fm->pf[j]<=0 || fm->pdgf[j] == 2112) continue;
            float a = abs(asin(fm->pxf[j]/fm->pf[j]))/3.14159*180;
            if (a<h_min_angle) {
                h_min_angle = a;
                h_min_angle_pdg = fm->pdgf[j];
                h_min_angle_Ef = fm->Ef[j];
            }
            // if (a < h_sec_min_angle_wp && a>h_min_angle_wp) h_sec_min_angle_wp = a;
        }
        T->Fill();

        h_min_angle = 91;
        h_min_angle_pdg = 0;
        h_min_angle_Ef = 0;
        l_angle = 91;
    }
    T->Write();
    f.Close();
}

void MakePlot::PlotMinAngle_hadron(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 800);

    TFile *f = new TFile("ah.root");
    TTree *t = (TTree*)f->Get("T");
    TCut ecut("El<6");

    int nTotal = fm->fChain->Draw("", ecut && "cc");

    c1->Divide(3,2);
    c1->cd(1);
    t->Draw("h_min_angle>>AngleCCQE(180,0,90)", ecut && "cc&&qel");
    TH1F *h = (TH1F*)gDirectory->Get("AngleCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Minimum hadron angle to wire plane (deg)");
    c1->cd(4);
    TH1* hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    cout << "CCQE:" << endl;
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;
    cout << hc->GetBinContent(5*2)*0.315 << ", " << hc->GetBinContent(7.5*2)*0.453 << ", " << hc->GetBinContent(10*2)*0.563 << endl;
    cout << "-----------" << endl;

    c1->cd(2);
    t->Draw("h_min_angle>>AngleCCRes(180,0,90)", ecut && "cc&&res");
    h = (TH1F*)gDirectory->Get("AngleCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Minimum hadron angle to wire plane (deg)");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    cout << "CCRes:" << endl;
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;
    cout << hc->GetBinContent(5*2)*0.280 << ", " << hc->GetBinContent(7.5*2)*0.401 << ", " << hc->GetBinContent(10*2)*0.507 << endl;
    cout << "-----------" << endl;

    c1->cd(3);
    t->Draw("h_min_angle>>AngleCCDis(180,0,90)", ecut && "cc&&dis");
    h = (TH1F*)gDirectory->Get("AngleCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Minimum hadron angle to wire plane (deg)");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    cout << "CCDis:" << endl;
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;
    cout << hc->GetBinContent(5*2)*0.176 << ", " << hc->GetBinContent(7.5*2)*0.270 << ", " << hc->GetBinContent(10*2)*0.364 << endl;
    cout << "-----------" << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::SetStyle()
{
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleFont(42, "x");
    gStyle->SetTitleFont(42, "y");
    gStyle->SetTitleFont(42, "z");
    gStyle->SetLabelFont(42, "x");
    gStyle->SetLabelFont(42, "y");
    gStyle->SetLabelFont(42, "z");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.25, "y");
    gStyle->SetHistLineWidth(2);
    gStyle->SetHistLineColor(kBlack);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetPalette(1);
    gROOT->ForceStyle();
    gROOT->UseCurrentStyle();
    fm->fChain->UseCurrentStyle();

}
