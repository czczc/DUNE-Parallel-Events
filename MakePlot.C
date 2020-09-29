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
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "THStack.h"
#include "TLegend.h"

#include <iostream>
#include <iomanip>

using namespace std;

MakePlot::MakePlot(const Double_t angle, const Bool_t r)
{
    DetMass = 40.0;
    exposure = 3.0;
    POTperYr = 1.1e21;

    fm = new FastMC();
    SetStyle();
    n_events = fm->fChain->GetEntriesFast();
    cutAngle = angle;
    rot = r;
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
    fm->fChain->Draw("Ev_reco:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCQE(100,0,50,200,0,20)","(cc&&qel)*OSCWeight","colz");
    TH2F *h = (TH2F*)gDirectory->Get("EvsAngleCCQE");
    h->SetTitle("CC-QE");
    h->GetXaxis()->SetTitle("Electron angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");

    c1->cd(2);
    fm->fChain->Draw("Ev_reco:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCRes(100,0,50,200,0,20)","(cc&&res)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvsAngleCCRes");
    h->SetTitle("CC-Res");
    h->GetXaxis()->SetTitle("Electron angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");

    c1->cd(3);
    fm->fChain->Draw("Ev_reco:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCDis(100,0,50,200,0,20)","(cc&&dis)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvsAngleCCDis");
    h->SetTitle("CC-Dis");
    h->GetXaxis()->SetTitle("Electron angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");

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
    fm->fChain->Draw("El:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCQE(100,0,50,200,0,20)","(cc&&qel)*OSCWeight","colz");
    TH2F *h = (TH2F*)gDirectory->Get("EvsAngleCCQE");
    h->SetTitle("CC-QE");
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    c1->cd(2);
    fm->fChain->Draw("El:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCRes(100,0,50,200,0,20)","(cc&&res)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvsAngleCCRes");
    h->SetTitle("CC-Res");
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    c1->cd(3);
    fm->fChain->Draw("El:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCDis(100,0,50,200,0,20)","(cc&&dis)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvsAngleCCDis");
    h->SetTitle("CC-Dis");
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
    fm->fChain->Draw("El:0.4/sin(abs(asin(pxl/pl)))>>El_vs_LengthCCQE(100,0,20,200,0,20)","(cc&&qel)*OSCWeight","colz");
    TH2F *h = (TH2F*)gDirectory->Get("El_vs_LengthCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    c1->cd(2);
    fm->fChain->Draw("El:0.4/sin(abs(asin(pxl/pl)))>>El_vs_LengthCCRes(100,0,20,200,0,20)","(cc&&res)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("El_vs_LengthCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    c1->cd(3);
    fm->fChain->Draw("El:0.4/sin(abs(asin(pxl/pl)))>>El_vs_LengthCCDis(100,0,20,200,0,20)","(cc&&dis)*OSCWeight","colz");
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

    // fm->fChain->Draw("Ev>>Ev(400,0.1,40)",
    //     Form("POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr)
    // );
    // TH1F *h = (TH1F*)gDirectory->Get("Ev");
    fm->fChain->Draw("Ev>>EvOsc(200,0,20)",
        Form("(cc)*OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr)
    );
    fm->fChain->Draw("Ev_reco>>EvOsc_reco(200,0,20)",
        Form("(cc)*OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr)
    );
    TH1F *h = (TH1F*)gDirectory->Get("EvOsc");
    TH1F *h2 = (TH1F*)gDirectory->Get("EvOsc_reco");
    h->Draw();
    h->SetTitle("");
    h->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    h->GetYaxis()->SetTitle("Events / 100 MeV");
    h2->SetLineColor(kRed);
    h2->Draw("same");
    cout << "totoal cc events: " << h->Integral() << endl;
    // hOsc->SetLineColor(kRed);
    // hOsc->Draw("same");
    // gPad->SetLogy();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotOscWeight_vs_Ev(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",800, 600);

    fm->fChain->Draw("OSCWeight:Ev>>Osc(200,0.1,20, 100, 0, 1)");
    TH2F *h = (TH2F*)gDirectory->Get("Osc");
    TProfile *px = h->ProfileX("pxOsc");
    px->SetTitle("");
    px->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    px->GetYaxis()->SetTitle("Oscillation Probability");
    px->Draw("");

    gPad->SetLogx();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEl_vs_Ev(const char* saveAs)
{   // don't apply osc weight here
    int nTotal = fm->fChain->Draw("","cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    fm->fChain->Draw("El:Ev>>EvElCCQE(200,0,20,200,0,20)","(cc&&qel)","colz");
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
    fm->fChain->Draw("El:Ev>>EvElCCRes(200,0,20,200,0,20)","(cc&&res)","colz");
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
    fm->fChain->Draw("El:Ev>>EvElCCDis(200,0,20,200,0,20)","(cc&&dis)","colz");
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
    TCut ecut("Ev_reco<6");

    int nTotal = fm->fChain->Draw("", (ecut && "cc")*"OSCWeight");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCQE(100,0,50)", ecut && "cc&&qel");
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180>>AngleCCQE(180,0,90)", (ecut && "cc&&qel")*"OSCWeight");
    TH1F *h = (TH1F*)gDirectory->Get("AngleCCQE");
    h->SetTitle("CC-QE");
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(4);
    TH1* hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(2.5*2) << ", " << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << endl;

    c1->cd(2);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCRes(100,0,50)", ecut && "cc&&res");
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180>>AngleCCRes(180,0,90)", (ecut && "cc&&res")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("AngleCCRes");
    h->SetTitle("CC-Res");
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(2.5*2) << ", " << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << endl;

    c1->cd(3);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCDis(100,0,50)", ecut && "cc&&dis");
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180>>AngleCCDis(180,0,90)", (ecut && "cc&&dis")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("AngleCCDis");
    h->SetTitle("CC-Dis");
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(2.5*2) << ", " << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEv_langle_cut(const char* saveAs)
{
    TCut ecut("Ev_reco<10");

    int nTotal = fm->fChain->Draw("", (ecut && "cc")*"OSCWeight");

    TString weight = Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr);
    TCut anglecut1, anglecut2, anglecut3;
    if (!rot) {
        anglecut1 = "abs(asin(pxl/pl))/3.14159*180<7.5";
        anglecut2 = "abs(asin(pxl/pl))/3.14159*180<5.0";
        anglecut3 = "abs(asin(pxl/pl))/3.14159*180<2.5";
    }
    else {
        anglecut1 = "abs(asin(pzl/pl))/3.14159*180<7.5";
        anglecut2 = "abs(asin(pzl/pl))/3.14159*180<5.0";
        anglecut3 = "abs(asin(pzl/pl))/3.14159*180<2.5";
    }

    TCanvas *c1 = new TCanvas("c1","",1200, 400);

    c1->Divide(3,1);
    c1->cd(1);

    fm->fChain->Draw("Ev_reco>>Ev_langle_cut1_qe(100,0,10)", (ecut && anglecut1 && "cc&&qel")*weight.Data());
    fm->fChain->Draw("Ev_reco>>Ev_langle_cut2_qe(100,0,10)", (ecut && anglecut2 && "cc&&qel")*weight.Data());
    fm->fChain->Draw("Ev_reco>>Ev_langle_cut3_qe(100,0,10)", (ecut && anglecut3 && "cc&&qel")*weight.Data());
    fm->fChain->Draw("Ev_reco>>Ev_qe(100,0,10)", (ecut && "cc&&qel")*weight.Data());

    TH1F *ht = (TH1F*)gDirectory->Get("Ev_qe");
    TH1F *h1 = (TH1F*)gDirectory->Get("Ev_langle_cut1_qe");
    TH1F *h2 = (TH1F*)gDirectory->Get("Ev_langle_cut2_qe");
    TH1F *h3 = (TH1F*)gDirectory->Get("Ev_langle_cut3_qe");
    ht->SetTitle("CC-QE");
    ht->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
    ht->GetYaxis()->SetTitle("Event / 100 MeV");
    ht->Draw();
    h1->SetLineColor(kRed);
    h1->Draw("same");
    h2->SetLineColor(kBlue);
    h2->Draw("same");
    h3->SetLineColor(kMagenta);
    h3->Draw("same");
    TLegend *leg1 = new TLegend(.52,.58,.89,.89);
    leg1->AddEntry(ht, " Total", "l");
    leg1->AddEntry(h1, " #theta_{e} < 7.5#circ", "l");
    leg1->AddEntry(h2, " #theta_{e} < 5.0#circ", "l");
    leg1->AddEntry(h3, " #theta_{e} < 2.5#circ", "l");
    leg1->Draw();

    c1->cd(2);
    fm->fChain->Draw("Ev_reco>>Ev_langle_cut1_res(100,0,10)", (ecut && anglecut1 && "cc&&res")*weight.Data());
    fm->fChain->Draw("Ev_reco>>Ev_langle_cut2_res(100,0,10)", (ecut && anglecut2 && "cc&&res")*weight.Data());
    fm->fChain->Draw("Ev_reco>>Ev_langle_cut3_res(100,0,10)", (ecut && anglecut3 && "cc&&res")*weight.Data());
    fm->fChain->Draw("Ev_reco>>Ev_res(100,0,10)", (ecut && "cc&&res")*weight.Data());

    ht = (TH1F*)gDirectory->Get("Ev_res");
    h1 = (TH1F*)gDirectory->Get("Ev_langle_cut1_res");
    h2 = (TH1F*)gDirectory->Get("Ev_langle_cut2_res");
    h3 = (TH1F*)gDirectory->Get("Ev_langle_cut3_res");
    ht->SetTitle("CC-RES");
    ht->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
    ht->GetYaxis()->SetTitle("Event / 100 MeV");
    ht->Draw();
    h1->SetLineColor(kRed);
    h1->Draw("same");
    h2->SetLineColor(kBlue);
    h2->Draw("same");
    h3->SetLineColor(kMagenta);
    h3->Draw("same");
    TLegend *leg2 = new TLegend(.52,.58,.89,.89);
    leg2->AddEntry(ht, " Total", "l");
    leg2->AddEntry(h1, " #theta_{e} < 7.5#circ", "l");
    leg2->AddEntry(h2, " #theta_{e} < 5.0#circ", "l");
    leg2->AddEntry(h3, " #theta_{e} < 2.5#circ", "l");
    leg2->Draw();

    c1->cd(3);
    fm->fChain->Draw("Ev_reco>>Ev_langle_cut1_dis(100,0,10)", (ecut && anglecut1 && "cc&&dis")*weight.Data());
    fm->fChain->Draw("Ev_reco>>Ev_langle_cut2_dis(100,0,10)", (ecut && anglecut2 && "cc&&dis")*weight.Data());
    fm->fChain->Draw("Ev_reco>>Ev_langle_cut3_dis(100,0,10)", (ecut && anglecut3 && "cc&&dis")*weight.Data());
    fm->fChain->Draw("Ev_reco>>Ev_dis(100,0,10)", (ecut && "cc&&dis")*weight.Data());

    ht = (TH1F*)gDirectory->Get("Ev_dis");
    h1 = (TH1F*)gDirectory->Get("Ev_langle_cut1_dis");
    h2 = (TH1F*)gDirectory->Get("Ev_langle_cut2_dis");
    h3 = (TH1F*)gDirectory->Get("Ev_langle_cut3_dis");
    ht->SetTitle("CC-DIS");
    ht->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
    ht->GetYaxis()->SetTitle("Event / 100 MeV");
    ht->Draw();
    h1->SetLineColor(kRed);
    h1->Draw("same");
    h2->SetLineColor(kBlue);
    h2->Draw("same");
    h3->SetLineColor(kMagenta);
    h3->Draw("same");
    TLegend *leg3 = new TLegend(.52,.58,.89,.89);
    leg3->AddEntry(ht, " Total", "l");
    leg3->AddEntry(h1, " #theta_{e} < 7.5#circ", "l");
    leg3->AddEntry(h2, " #theta_{e} < 5.0#circ", "l");
    leg3->AddEntry(h3, " #theta_{e} < 2.5#circ", "l");
    leg3->Draw();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}


void MakePlot::PlotAngle_cut_E_Perp(const char* saveAs)
{
    TCut ecut("Ev_reco<6");

    int nTotal = fm->fChain->Draw("", ecut && "cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCQE(100,0,50)", ecut && "cc&&qel");
    fm->fChain->Draw("abs(asin(pzl/pl))/3.14159*180>>AngleCCQE(180,0,90)", (ecut && "cc&&qel")*"OSCWeight");
    TH1F *h = (TH1F*)gDirectory->Get("AngleCCQE");
    h->SetTitle("CC-QE");
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(4);
    TH1* hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(2.5*2) << ", " << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << endl;

    c1->cd(2);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCRes(100,0,50)", ecut && "cc&&res");
    fm->fChain->Draw("abs(asin(pzl/pl))/3.14159*180>>AngleCCRes(180,0,90)", (ecut && "cc&&res")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("AngleCCRes");
    h->SetTitle("CC-Res");
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(2.5*2) << ", " << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << endl;

    c1->cd(3);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCDis(100,0,50)", ecut && "cc&&dis");
    fm->fChain->Draw("abs(asin(pzl/pl))/3.14159*180>>AngleCCDis(180,0,90)", (ecut && "cc&&dis")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("AngleCCDis");
    h->SetTitle("CC-Dis");
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(2.5*2) << ", " << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << endl;

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
    fm->fChain->Draw("0.4/sin(abs(asin(pxl/pl)))>>LengthCCQE(100,0,20)", (ecut && "cc&&qel")*"OSCWeight");
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
    fm->fChain->Draw("0.4/sin(abs(asin(pxl/pl)))>>LengthCCRes(100,0,20)", (ecut && "cc&&res")*"OSCWeight");
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
    fm->fChain->Draw("0.4/sin(abs(asin(pxl/pl)))>>LengthCCDis(100,0,20)", (ecut && "cc&&dis")*"OSCWeight");
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
    fm->fChain->Draw("0.4/sin(abs(asin(pzl/pl)))>>LengthCCQE(100,0,20)", (ecut && "cc&&qel")*"OSCWeight");
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
    fm->fChain->Draw("0.4/sin(abs(asin(pzl/pl)))>>LengthCCRes(100,0,20)", (ecut && "cc&&res")*"OSCWeight");
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
    fm->fChain->Draw("0.4/sin(abs(asin(pzl/pl)))>>LengthCCDis(100,0,20)", (ecut && "cc&&dis")*"OSCWeight");
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
    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    c1->Divide(3,1);
    c1->cd(1);
    fm->fChain->Draw("Ev:Ef>>EvEfCCQE(50,0,5,200,0,20)","(cc&&qel&&pf>0)*OSCWeight","colz");
    TH2F *h = (TH2F*)gDirectory->Get("EvEfCCQE");
    h->SetTitle("CC-QE");
    h->GetYaxis()->SetTitle("Neutrino energy (GeV)");
    h->GetXaxis()->SetTitle("Hadron Energy (GeV)");
    gPad->SetLogz();
    // c1->cd(4);
    // h->ProjectionX()->Draw();


    c1->cd(2);
    fm->fChain->Draw("Ev:Ef>>EvEfCCRes(50,0,5,200,0,20)","(cc&&res&&pf>0)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvEfCCRes");
    h->SetTitle("CC-Res");
    h->GetYaxis()->SetTitle("Neutrino energy (GeV)");
    h->GetXaxis()->SetTitle("Hadron Energy (GeV)");
    gPad->SetLogz();
    // c1->cd(5);
    // h->ProjectionX()->Draw();

    c1->cd(3);
    fm->fChain->Draw("Ev:Ef>>EvEfCCDis(50,0,5,200,0,20)","(cc&&dis&&pf>0)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvEfCCDis");
    h->SetTitle("CC-Dis");
    h->GetYaxis()->SetTitle("Neutrino energy (GeV)");
    h->GetXaxis()->SetTitle("Hadron Energy (GeV)");
    gPad->SetLogz();
    // c1->cd(6);
    // h->ProjectionX()->Draw();

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
    fm->fChain->Draw(TString::Format("nf:(%s)/nf>>Nf_RatioCCQE(51,0,1.02,42,0,42)", particle), "(cc&&qel)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("Nf_RatioCCQE");
    h->SetTitle("CC-QE");
    h->GetYaxis()->SetTitle("# Final state particles (hadron side)");
    h->GetXaxis()->SetTitle(title);
    gPad->SetLogz();

    c1->cd(2);
    fm->fChain->Draw(TString::Format("nf:(%s)/nf>>Nf_RatioCCRes(51,0,1.02,42,0,42)", particle),"(cc&&res)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("Nf_RatioCCRes");
    h->SetTitle("CC-Res");
    h->GetYaxis()->SetTitle("# Final state particles (hadron side)");
    h->GetXaxis()->SetTitle(title);
    gPad->SetLogz();

    c1->cd(3);
    fm->fChain->Draw(TString::Format("nf:(%s)/nf>>Nf_RatioCCDis(51,0,1.02,42,0,42)", particle),"(cc&&res)*OSCWeight","colz");
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
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180:abs(asin(pxf/pf))/3.14159*180>>Angleh_vs_AnglelCCQE(180,0,90,180,0,90)","(cc&&qel&&El<6&&pf>0&&pdgf!=2112&&pdgf!=22)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("Angleh_vs_AnglelCCQE");
    h->SetTitle("CC-QE");
    h->GetXaxis()->SetTitle("hadron angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("lepton angle to wire plane (deg)");
    gPad->SetLogz();
    c1->cd(4);
    h->ProjectionX()->Draw();

    c1->cd(2);
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180:abs(asin(pxf/pf))/3.14159*180>>Angleh_vs_AnglelCCRes(180,0,90,180,0,90)","(cc&&res&&El<6&&pf>0&&pdgf!=2112&&pdgf!=22)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("Angleh_vs_AnglelCCRes");
    h->SetTitle("CC-Res");
    h->GetXaxis()->SetTitle("hadron angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("lepton angle to wire plane (deg)");
    gPad->SetLogz();
    c1->cd(5);
    h->ProjectionX()->Draw();

    c1->cd(3);
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180:abs(asin(pxf/pf))/3.14159*180>>Angleh_vs_AnglelCCDis(180,0,90,180,0,90)","(cc&&dis&&El<6&&pf>0&&pdgf!=2112&&pdgf!=22)*OSCWeight","colz");
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

void MakePlot::Print(TString filename)
{
    TFile *f = new TFile(filename);
    TTree *t = (TTree*)f->Get("T");
    bool isnc = filename.Contains("-nc-");
    TCut ccCut = isnc ? "nc" : "cc";
    // TCut ecut = isnc ? "Ev_reco<6" : "El<6";
    TCut ecut = "Ev_reco<6";

    t->Draw("l_angle>>LAngleCC(180,0,90)", (ecut && ccCut)*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    TH1F *h = (TH1F*)gDirectory->Get("LAngleCC");
    float nTotal = h->Integral(0,181);
    cout << fixed << setprecision(1);
    cout << "total " << (isnc ? "nc" : "cc") << " events: " << nTotal << endl; // this include others like ccqel


    t->Draw("l_angle>>LAngleQE(180,0,90)", (ecut && ccCut && "qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    t->Draw("l_angle>>LAngleRes(180,0,90)", (ecut && ccCut && "res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    t->Draw("l_angle>>LAngleDis(180,0,90)", (ecut && ccCut && "dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("LAngleQE");
    float nQE = h->Integral(0,181);
    h = (TH1F*)gDirectory->Get("LAngleRes");
    float nRes = h->Integral(0,181);
    h = (TH1F*)gDirectory->Get("LAngleDis");
    float nDis = h->Integral(0,181);
    float n = nQE + nRes + nDis;
    cout << "QE, Res, Dis: " << nQE << ", " << nRes << ", " << nDis << " | " << n << endl;

    float nSignal = n;
    if (isnc) nSignal = 801.1;
    cout << "total CC signal: " << nSignal << endl;

    for (int i=1; i<=3; i++) {
        float angle = 2.5 * i;
        // TCut anglecut = Form("l_angle<%.1f && h_min_angle<%.1f", angle, angle);
        // TCut anglecut = Form("l_angle<%.1f", angle);
        // TCut anglecut = Form("h_min_angle<%.1f", angle);
        // TCut anglecut = Form("l_angle<%.1f && EM_min_angle<%.1f", angle, angle);
        // TCut anglecut = Form("h_min_angle<%.1f && EM_min_angle<%.1f", angle, angle);
        TCut anglecut = Form("h_min_angle<%.1f && EM_min_angle_gt_3cm<%.1f", angle, angle);

        t->Draw("l_angle>>LAngleQEAngle(180,0,90)" , (ecut && ccCut && anglecut && "qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
        t->Draw("l_angle>>LAngleResAngle(180,0,90)", (ecut && ccCut && anglecut && "res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
        t->Draw("l_angle>>LAngleDisAngle(180,0,90)", (ecut && ccCut && anglecut && "dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
        h = (TH1F*)gDirectory->Get("LAngleQEAngle");
        float nQEAngle = h->Integral(0,181);
        h = (TH1F*)gDirectory->Get("LAngleResAngle");
        float nResAngle = h->Integral(0,181);
        h = (TH1F*)gDirectory->Get("LAngleDisAngle");
        float nDisAngle = h->Integral(0,181);
        float nAll = (nQEAngle + nResAngle + nDisAngle);

        cout << "< " << angle << " deg : "
            << nQEAngle << " (" << nQEAngle*100./nQE << "%), "
            << nResAngle << " (" << nResAngle*100./nRes << "%), "
            << nDisAngle << " (" << nDisAngle*100./nDis << "%) | "
            << nAll << " (" << nAll*100./nSignal << "%)" << endl;
    }

}

void MakePlot::PlotMinAngle_lepton(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 800);

    TFile *f;
    if (!rot) f = new TFile("numu-nue-cc-hadron_angle.root");
    else f = new TFile("numu-nue-cc-hadron_angle_perp.root");
    TTree *t = (TTree*)f->Get("T");
    TCut ecut("El<6");

    int nTotal = fm->fChain->Draw("", (ecut && "cc")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    cout<<"-----------"<<endl;
    cout<<"Total simulated events: "<<fm->fChain->GetEntries()<<endl;
    cout<<"Number of CC events: "<<nTotal<<endl;
    cout<<"Number of long hadron events: "<<t->Draw("",Form("(El<6&&cc&&h_min_angle<91)*OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr))<<endl;
    c1->Divide(3,2);
    c1->cd(1);
    t->Draw("l_angle>>LAngleCCQE(180,0,90)", (ecut && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    TH1F *h = (TH1F*)gDirectory->Get("LAngleCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Minimum lepton angle to wire plane (deg)");
    c1->cd(4);
    TH1* hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    Double_t I = h->Integral();
    cout << "CCQE:" << endl;
    cout << I <<endl;
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;
    cout << hc->GetBinContent(5*2)*I << ", " << hc->GetBinContent(7.5*2)*I << ", " << hc->GetBinContent(10*2)*I << endl;
    cout << "-----------" << endl;

    c1->cd(2);
    t->Draw("l_angle>>LAngleCCRes(180,0,90)", (ecut && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("LAngleCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Minimum lepton angle to wire plane (deg)");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    I = h->Integral();
    cout << "CCRes:" << endl;
    cout << I <<endl;
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;
    cout << hc->GetBinContent(5*2)*I << ", " << hc->GetBinContent(7.5*2)*I << ", " << hc->GetBinContent(10*2)*I << endl;
    cout << "-----------" << endl;

    c1->cd(3);
    t->Draw("l_angle>>LAngleCCDis(180,0,90)", (ecut && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("LAngleCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Minimum lepton angle to wire plane (deg)");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    I = h->Integral();
    cout << "CCDis:" << endl;
    cout << I <<endl;
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;
    cout << hc->GetBinContent(5*2)*I << ", " << hc->GetBinContent(7.5*2)*I << ", " << hc->GetBinContent(10*2)*I << endl;
    cout << "-----------" << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotMinAngle_hadron(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 800);

    TFile *f;
    if (!rot) f = new TFile("numu-nue-cc-hadron_angle.root");
    else f = new TFile("numu-nue-cc-hadron_angle_perp.root");
    TTree *t = (TTree*)f->Get("T");
    TCut ecut("Ev_reco<6");

    int nTotal = fm->fChain->Draw("", (ecut && "cc")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    cout<<"-----------"<<endl;
    cout<<"Total simulated events: "<<fm->fChain->GetEntries()<<endl;
    cout<<"Number of CC events: "<<nTotal<<endl;
    cout<<"Number of long hadron events: "<<t->Draw("",Form("(El<6&&cc&&h_min_angle<91)*OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr))<<endl;
    c1->Divide(3,2);
    c1->cd(1);
    t->Draw("h_min_angle>>AngleCCQE(180,0,90)", (ecut && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    TH1F *h = (TH1F*)gDirectory->Get("AngleCCQE");
    h->SetTitle("CC-QE");
    h->GetXaxis()->SetTitle("Minimum hadron angle to wire plane (deg)");
    c1->cd(4);
    TH1* hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    Double_t I = h->Integral();
    cout << "CCQE:" << endl;
    cout << I <<endl;
    cout << hc->GetBinContent(2.5*2) << ", " << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << endl;
    // cout << hc->GetBinContent(5*2)*0.315 << ", " << hc->GetBinContent(7.5*2)*0.453 << ", " << hc->GetBinContent(10*2)*0.563 << endl;
    // cout << hc->GetBinContent(5*2)*I << ", " << hc->GetBinContent(7.5*2)*I << ", " << hc->GetBinContent(10*2)*I << endl;
    // cout << hc->GetBinContent(5*2)*0.315*I << ", " << hc->GetBinContent(7.5*2)*0.453*I << ", " << hc->GetBinContent(10*2)*0.563*I << endl;
    cout << "-----------" << endl;

    c1->cd(2);
    t->Draw("h_min_angle>>AngleCCRes(180,0,90)", (ecut && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("AngleCCRes");
    h->SetTitle("CC-Res");
    h->GetXaxis()->SetTitle("Minimum hadron angle to wire plane (deg)");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    I = h->Integral();
    cout << "CCRes:" << endl;
    cout << I <<endl;
    cout << hc->GetBinContent(2.5*2) << ", " << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << endl;
    // cout << hc->GetBinContent(5*2)*0.280 << ", " << hc->GetBinContent(7.5*2)*0.401 << ", " << hc->GetBinContent(10*2)*0.507 << endl;
    // cout << hc->GetBinContent(5*2)*I << ", " << hc->GetBinContent(7.5*2)*I << ", " << hc->GetBinContent(10*2)*I << endl;
    // cout << hc->GetBinContent(5*2)*0.280*I << ", " << hc->GetBinContent(7.5*2)*0.401*I << ", " << hc->GetBinContent(10*2)*0.507*I << endl;
    cout << "-----------" << endl;

    c1->cd(3);
    t->Draw("h_min_angle>>AngleCCDis(180,0,90)", (ecut && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("AngleCCDis");
    h->SetTitle("CC-Dis");
    h->GetXaxis()->SetTitle("Minimum hadron angle to wire plane (deg)");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    I = h->Integral();
    cout << "CCDis:" << endl;
    cout << I <<endl;
    cout << hc->GetBinContent(2.5*2) << ", " << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << endl;
    // cout << hc->GetBinContent(5*2)*0.176 << ", " << hc->GetBinContent(7.5*2)*0.270 << ", " << hc->GetBinContent(10*2)*0.364 << endl;
    // cout << hc->GetBinContent(5*2)*I << ", " << hc->GetBinContent(7.5*2)*I << ", " << hc->GetBinContent(10*2)*I << endl;
    // cout << hc->GetBinContent(5*2)*0.176*I << ", " << hc->GetBinContent(7.5*2)*0.270*I << ", " << hc->GetBinContent(10*2)*0.364*I << endl;
    cout << "-----------" << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEv_hangle_cut(const char* saveAs, bool cutBoth)
{
    TFile *f;
    if (!rot) f = new TFile("numu-nue-cc-hadron_angle.root");
    else f = new TFile("numu-nue-cc-hadron_angle_perp.root");
    TTree *t = (TTree*)f->Get("T");

    TCut ecut("Ev_reco<10");

    int nTotal = t->Draw("", (ecut && "cc")*"OSCWeight");

    TString weight = Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr);
    TCut anglecut1, anglecut2, anglecut3;
    if (!cutBoth) {
        anglecut1 = "h_min_angle<7.5";
        anglecut2 = "h_min_angle<5.0";
        anglecut3 = "h_min_angle<2.5";
    }
    else {
        anglecut1 = "l_angle<7.5 && h_min_angle<7.5";
        anglecut2 = "l_angle<5.0 && h_min_angle<5.0";
        anglecut3 = "l_angle<2.5 && h_min_angle<2.5";
    }
    TString prestr(" ");
    if (cutBoth) prestr = "#theta_{e},";

    TCanvas *c1 = new TCanvas("c1","",1200, 400);

    c1->Divide(3,1);
    c1->cd(1);

    t->Draw("Ev_reco>>Ev_langle_cut1_qe(100,0,10)", (ecut && anglecut1 && "cc&&qel")*weight.Data());
    t->Draw("Ev_reco>>Ev_langle_cut2_qe(100,0,10)", (ecut && anglecut2 && "cc&&qel")*weight.Data());
    t->Draw("Ev_reco>>Ev_langle_cut3_qe(100,0,10)", (ecut && anglecut3 && "cc&&qel")*weight.Data());
    t->Draw("Ev_reco>>Ev_qe(100,0,10)", (ecut && "cc&&qel")*weight.Data());

    TH1F *ht = (TH1F*)gDirectory->Get("Ev_qe");
    TH1F *h1 = (TH1F*)gDirectory->Get("Ev_langle_cut1_qe");
    TH1F *h2 = (TH1F*)gDirectory->Get("Ev_langle_cut2_qe");
    TH1F *h3 = (TH1F*)gDirectory->Get("Ev_langle_cut3_qe");
    ht->SetTitle("CC-QE");
    ht->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
    ht->GetYaxis()->SetTitle("Event / 100 MeV");
    ht->Draw();
    h1->SetLineColor(kRed);
    h1->Draw("same");
    h2->SetLineColor(kBlue);
    h2->Draw("same");
    h3->SetLineColor(kMagenta);
    h3->Draw("same");
    TLegend *leg1 = new TLegend(.52,.58,.89,.89);
    leg1->AddEntry(ht, " Total", "l");
    leg1->AddEntry(h1, prestr + "#theta_{h} < 7.5#circ", "l");
    leg1->AddEntry(h2, prestr + "#theta_{h} < 5.0#circ", "l");
    leg1->AddEntry(h3, prestr + "#theta_{h} < 2.5#circ", "l");
    leg1->Draw();

    c1->cd(2);
    t->Draw("Ev_reco>>Ev_langle_cut1_res(100,0,10)", (ecut && anglecut1 && "cc&&res")*weight.Data());
    t->Draw("Ev_reco>>Ev_langle_cut2_res(100,0,10)", (ecut && anglecut2 && "cc&&res")*weight.Data());
    t->Draw("Ev_reco>>Ev_langle_cut3_res(100,0,10)", (ecut && anglecut3 && "cc&&res")*weight.Data());
    t->Draw("Ev_reco>>Ev_res(100,0,10)", (ecut && "cc&&res")*weight.Data());

    ht = (TH1F*)gDirectory->Get("Ev_res");
    h1 = (TH1F*)gDirectory->Get("Ev_langle_cut1_res");
    h2 = (TH1F*)gDirectory->Get("Ev_langle_cut2_res");
    h3 = (TH1F*)gDirectory->Get("Ev_langle_cut3_res");
    ht->SetTitle("CC-RES");
    ht->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
    ht->GetYaxis()->SetTitle("Event / 100 MeV");
    ht->Draw();
    h1->SetLineColor(kRed);
    h1->Draw("same");
    h2->SetLineColor(kBlue);
    h2->Draw("same");
    h3->SetLineColor(kMagenta);
    h3->Draw("same");
    TLegend *leg2 = new TLegend(.52,.58,.89,.89);
    leg2->AddEntry(ht, " Total", "l");
    leg2->AddEntry(h1, prestr + "#theta_{h} < 7.5#circ", "l");
    leg2->AddEntry(h2, prestr + "#theta_{h} < 5.0#circ", "l");
    leg2->AddEntry(h3, prestr + "#theta_{h} < 2.5#circ", "l");
    leg2->Draw();

    c1->cd(3);
    t->Draw("Ev_reco>>Ev_langle_cut1_dis(100,0,10)", (ecut && anglecut1 && "cc&&dis")*weight.Data());
    t->Draw("Ev_reco>>Ev_langle_cut2_dis(100,0,10)", (ecut && anglecut2 && "cc&&dis")*weight.Data());
    t->Draw("Ev_reco>>Ev_langle_cut3_dis(100,0,10)", (ecut && anglecut3 && "cc&&dis")*weight.Data());
    t->Draw("Ev_reco>>Ev_dis(100,0,10)", (ecut && "cc&&dis")*weight.Data());

    ht = (TH1F*)gDirectory->Get("Ev_dis");
    h1 = (TH1F*)gDirectory->Get("Ev_langle_cut1_dis");
    h2 = (TH1F*)gDirectory->Get("Ev_langle_cut2_dis");
    h3 = (TH1F*)gDirectory->Get("Ev_langle_cut3_dis");
    ht->SetTitle("CC-DIS");
    ht->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
    ht->GetYaxis()->SetTitle("Event / 100 MeV");
    ht->Draw();
    h1->SetLineColor(kRed);
    h1->Draw("same");
    h2->SetLineColor(kBlue);
    h2->Draw("same");
    h3->SetLineColor(kMagenta);
    h3->Draw("same");
    TLegend *leg3 = new TLegend(.52,.58,.89,.89);
    leg3->AddEntry(ht, " Total", "l");
    leg3->AddEntry(h1, prestr+"#theta_{h} < 7.5#circ", "l");
    leg3->AddEntry(h2, prestr+"#theta_{h} < 5.0#circ", "l");
    leg3->AddEntry(h3, prestr+"#theta_{h} < 2.5#circ", "l");
    leg3->Draw();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEv_sb_overlapping(const char* saveAs, bool long_photon_pi)
{
    TFile *f;
    if (!rot) f = new TFile("numu-nue-cc-hadron_angle.root");
    else f = new TFile("numu-nue-cc-hadron_angle_perp.root");
    TTree *t = (TTree*)f->Get("T");

    TFile *f2;
    if (!rot) f2 = new TFile("numu-numu-nc-hadron_angle.root");
    else f2 = new TFile("numu-numu-nc-hadron_angle_perp.root");
    TTree *t2 = (TTree*)f2->Get("T");


    TCut ecut("Ev_reco<10");

    int nTotal = t->Draw("", (ecut && "cc")*"OSCWeight");

    TString weight = Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr);
    TCut anglecut1, anglecut2, anglecut3;
    TCut anglecut1b, anglecut2b, anglecut3b;


    TString prestr = "(NC#pi^{0}) #theta_{e},";

    TCanvas *c1 = new TCanvas("c1","",800, 800);

    anglecut1 = "l_angle<7.5 && h_min_angle<7.5";
    anglecut2 = "l_angle<5.0 && h_min_angle<5.0";
    anglecut3 = "l_angle<2.5 && h_min_angle<2.5";

    t->Draw("Ev_reco>>Ev_langle_cut1(100,0,10)", (ecut && anglecut1 && "cc")*weight.Data());
    t->Draw("Ev_reco>>Ev_langle_cut2(100,0,10)", (ecut && anglecut2 && "cc")*weight.Data());
    t->Draw("Ev_reco>>Ev_langle_cut3(100,0,10)", (ecut && anglecut3 && "cc")*weight.Data());
    t->Draw("Ev_reco>>Ev(100,0,10)", (ecut && "cc")*weight.Data());

    ecut = "Ev_reco<10";
    anglecut1b = "EM_min_angle<7.5 && h_min_angle<7.5";
    anglecut2b = "EM_min_angle<5.0 && h_min_angle<5.0";
    anglecut3b = "EM_min_angle<2.5 && h_min_angle<2.5";
    if(long_photon_pi) {
        anglecut1b = "EM_min_angle_gt_3cm<7.5 && h_min_angle<7.5";
        anglecut2b = "EM_min_angle_gt_4_6cm<5.0 && h_min_angle<5.0";
        anglecut3b = "EM_min_angle_gt_9_2cm<2.5 && h_min_angle<2.5";
    }
    t2->Draw("Ev_reco>>Ev_langle_cut1b(100,0,10)", (ecut && anglecut1b && "nc")*weight.Data());
    t2->Draw("Ev_reco>>Ev_langle_cut2b(100,0,10)", (ecut && anglecut2b && "nc")*weight.Data());
    t2->Draw("Ev_reco>>Ev_langle_cut3b(100,0,10)", (ecut && anglecut3b && "nc")*weight.Data());

    TH1F *ht = (TH1F*)gDirectory->Get("Ev");
    TH1F *h1 = (TH1F*)gDirectory->Get("Ev_langle_cut1");
    TH1F *h2 = (TH1F*)gDirectory->Get("Ev_langle_cut2");
    TH1F *h3 = (TH1F*)gDirectory->Get("Ev_langle_cut3");
    ht->SetTitle("");
    ht->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
    ht->GetYaxis()->SetTitle("Event / 100 MeV");
    ht->Draw();
    // h1->SetLineColor(kRed);
    // h1->Draw("same");
    // h2->SetLineColor(kBlue);
    // h2->Draw("same");
    // h3->SetLineColor(kMagenta);
    // h3->Draw("same");


    TH1F *h1b = (TH1F*)gDirectory->Get("Ev_langle_cut1b");
    TH1F *h2b = (TH1F*)gDirectory->Get("Ev_langle_cut2b");
    TH1F *h3b = (TH1F*)gDirectory->Get("Ev_langle_cut3b");
    // h1b->SetLineStyle(2);
    h1b->SetLineColor(kRed);
    h1b->Draw("same");
    // h2b->SetLineStyle(2);
    h2b->SetLineColor(kBlue);
    h2b->Draw("same");
    // h3b->SetLineStyle(2);
    h3b->SetLineColor(kMagenta);
    h3b->Draw("same");

    TLegend *leg1 = new TLegend(.52,.58,.89,.89);
    leg1->AddEntry(ht, " Total CC", "l");
    // leg1->AddEntry(h1, prestr + "#theta_{h} < 7.5#circ", "l");
    // leg1->AddEntry(h2, prestr + "#theta_{h} < 5.0#circ", "l");
    // leg1->AddEntry(h3, prestr + "#theta_{h} < 2.5#circ", "l");
    if (!long_photon_pi) {
        leg1->AddEntry(h1b, prestr + "#theta_{h} < 7.5#circ", "l");
        leg1->AddEntry(h2b, prestr + "#theta_{h} < 5.0#circ", "l");
        leg1->AddEntry(h3b, prestr + "#theta_{h} < 2.5#circ", "l");
    }
    else {
        leg1->AddEntry(h1b, "#gamma > 3.0 cm, #theta_{e},#theta_{#pi} < 7.5#circ", "l");
        leg1->AddEntry(h2b, "#gamma > 4.6 cm, #theta_{e},#theta_{#pi} < 5.0#circ", "l");
        leg1->AddEntry(h3b, "#gamma > 9.2 cm, #theta_{e},#theta_{#pi} < 2.5#circ", "l");
    }

    leg1->Draw();


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
    // float h_min_angle_eside = 91;
    // int h_min_angle_eside_pdg = 0;
    // float h_min_angle_eside_Ef = 0;
    // float h_sec_min_angle = 91;
    float l_angle = 91;
    bool is_lh_sameside = false;
    float a;

    T->Branch("h_min_angle", &h_min_angle);
    T->Branch("h_min_angle_pdg", &h_min_angle_pdg);
    T->Branch("h_min_angle_Ef", &h_min_angle_Ef);
    // T->Branch("h_min_angle_eside", &h_min_angle_eside);
    // T->Branch("h_min_angle_eside_pdg", &h_min_angle_eside_pdg);
    // T->Branch("h_min_angle_eside_Ef", &h_min_angle_eside_Ef);
    // T->Branch("h_sec_min_angle", &h_sec_min_angle);
    T->Branch("l_angle", &l_angle);
    T->Branch("is_lh_sameside", &is_lh_sameside);

    // photon info
    float ph_min_angle = 91;
    float l_ph_convert = 0;
    int nEMf = 0;
    float Eph = 0;
    T->Branch("EM_min_angle", &ph_min_angle);
    T->Branch("l_convert", &l_ph_convert);
    T->Branch("nEMf", &(fm->nEMf));
    TRandom3 r;

    // pi,photon overlap
    // h_min_angle_pdg == 211 or -211
    // h_minangle_Ef == 25, 30, 35 MeV
    float EM_min_angle_gt_3cm = 91;
    float EM_min_angle_gt_4_6cm = 91;
    float EM_min_angle_gt_9_2cm = 91;
    T->Branch("EM_min_angle_gt_3cm", &EM_min_angle_gt_3cm);
    T->Branch("EM_min_angle_gt_4_6cm", &EM_min_angle_gt_4_6cm);
    T->Branch("EM_min_angle_gt_9_2cm", &EM_min_angle_gt_9_2cm);


    T->Branch("cc", &(fm->cc) );
    T->Branch("nc", &(fm->nc) );
    T->Branch("qel", &(fm->qel) );
    T->Branch("res", &(fm->res) );
    T->Branch("dis", &(fm->dis) );
    T->Branch("El", &(fm->El));
    T->Branch("Ev", &(fm->Ev));
    T->Branch("Ev_reco", &(fm->Ev_reco));
    T->Branch("nf", &(fm->nf));
    T->Branch("POTWeight",&(fm->POTWeight));
    T->Branch("OSCWeight",&(fm->OSCWeight));
    T->Branch("n_events",&n_events);
    T->Branch("DetMass",&DetMass);
    T->Branch("exposure",&exposure);
    T->Branch("POTperYr",&POTperYr);

    int nEntries = fm->fChain->GetEntries();
    for (int i=0; i<nEntries; i++) {
        fm->fChain->GetEntry(i);

        if(!rot) l_angle = asin(fm->pxl/fm->pl)/3.14159*180;
	    else l_angle = asin(fm->pzl/fm->pl)/3.14159*180;

        for (int j=0; j<fm->nf; j++) {

            if (fm->pf[j]<=0 || fm->pdgf[j] == 2112 || fm->pdgf[j] == 111 || fm->pdgf[j] == 22) continue;
	        //Cut hadrons below threshold
	        if (fm->pdgf[j] == 2212 && fm->Ef[j] < 0.039) continue; //cut protons below 39 MeV
	        if ((fm->pdgf[j] == 321 || fm->pdgf[j] == -321) && fm->Ef[j] < 0.03) continue; //cut kaons below 30 MeV
	        if ((fm->pdgf[j] == 211 || fm->pdgf[j] == -211) && fm->Ef[j] < 0.0175) continue; //cut pions below 17.5 MeV
            if(!rot) a = asin(fm->pxf[j]/fm->pf[j])/3.14159*180;
	        else a = asin(fm->pzf[j]/fm->pf[j])/3.14159*180;
            // if (abs(a)<h_min_angle_eside && l_angle*a>0) {
            //     h_min_angle_eside = abs(a);
            //     h_min_angle_eside_pdg = fm->pdgf[j];
            //     h_min_angle_eside_Ef = fm->Ef[j];
            // }
    	    if (abs(a)<h_min_angle) {
                h_min_angle = abs(a);
                h_min_angle_pdg = fm->pdgf[j];
                h_min_angle_Ef = fm->Ef[j];
                if (l_angle*a>0) is_lh_sameside = true;
                else is_lh_sameside = false;
            }
            // if (a < h_sec_min_angle_wp && a>h_min_angle_wp) h_sec_min_angle_wp = a;
        }
        bool isE = false;
        for (int j=0; j<fm->nEMf; j++) {
            // if(fm->EMEf[j]>6) continue;
            if(!rot) a = abs(asin(fm->EMpxf[j]/fm->EMpf[j]))/3.14159*180;
            else a = abs(asin(fm->EMpzf[j]/fm->EMpf[j]))/3.14159*180;
            if(a < ph_min_angle) {
                if(fm->EMpdgf[j] == 22) isE = false;
                else isE = true;
                float l = r.Exp(18);
                if( (!isE) && l>3) continue; // photon must convert within 3 cm;
                l_ph_convert = l;
                ph_min_angle = a;
                Eph = fm->EMEf[j];
            }
        }
        if(isE) l_ph_convert = 0;


        // photon overlap with pion far away
        if (h_min_angle_pdg == 211 || h_min_angle_pdg == -211) {
            for (int j=0; j<fm->nEMf; j++) {
                if(fm->EMpdgf[j] != 22) continue;
                if(!rot) a = abs(asin(fm->EMpxf[j]/fm->EMpf[j]))/3.14159*180;
                else a = abs(asin(fm->EMpzf[j]/fm->EMpf[j]))/3.14159*180;
                float l = r.Exp(18);
                if (l>3 && h_min_angle_Ef>25e-3) {
                    if (a<EM_min_angle_gt_3cm) {
                        EM_min_angle_gt_3cm = a;
                    }
                }
                if (l>4.6 && h_min_angle_Ef>30e-3) {
                    if (a<EM_min_angle_gt_4_6cm) {
                        EM_min_angle_gt_4_6cm = a;
                    }
                }
                if (l>9.2 && h_min_angle_Ef>35e-3) {
                    if (a<EM_min_angle_gt_9_2cm) {
                        EM_min_angle_gt_9_2cm = a;
                    }
                }
            }
        }
	    l_angle = abs(l_angle);
        T->Fill();

        h_min_angle = 91;
        h_min_angle_pdg = 0;
        h_min_angle_Ef = 0;
	    // h_min_angle_eside = 91;
     //    h_min_angle_eside_pdg = 0;
     //    h_min_angle_eside_Ef = 0;
        l_angle = 91;
        is_lh_sameside = false;

        ph_min_angle = 91;
        l_ph_convert = 0;
        nEMf = 0;
        Eph = 0;

        EM_min_angle_gt_3cm = 91;
        EM_min_angle_gt_4_6cm = 91;
        EM_min_angle_gt_9_2cm = 91;

    }
    T->Write();
    f.Close();
}

void MakePlot::CreateCCEMTree(const char* filename)
{
	TFile f(filename,"recreate");
    	TTree *T = new TTree("T", "New Tree from Fast MC");

	float l_angle = 91;
	float ph_min_angle = 91;
	float l_ph_convert = 0;
	float Eph = 0;
	float a;
	TRandom3 r;

	bool isE = 0;

    	T->Branch("l_angle", &l_angle);
	T->Branch("EM_min_angle", &ph_min_angle);
	T->Branch("l_convert", &l_ph_convert);
    	T->Branch("cc", &(fm->cc) );
	T->Branch("nc", &(fm->nc) );
    	T->Branch("qel", &(fm->qel) );
    	T->Branch("res", &(fm->res) );
    	T->Branch("dis", &(fm->dis) );
    	T->Branch("El", &(fm->El));
    	T->Branch("EEM", &Eph);
	T->Branch("POTWeight",&(fm->POTWeight));
	T->Branch("OSCWeight",&(fm->OSCWeight));
	T->Branch("n_events",&n_events);
	T->Branch("DetMass",&DetMass);
	T->Branch("exposure",&exposure);
	T->Branch("POTperYr",&POTperYr);

	int nEntries = fm->fChain->GetEntries();

    	for (int i=0; i<nEntries; i++) {
        	fm->fChain->GetEntry(i);
		if(!rot) l_angle = abs(asin(fm->pxl/fm->pl))/3.14159*180;
		if(rot) l_angle = abs(asin(fm->pzl/fm->pl))/3.14159*180;
        	for (int j=0; j<fm->nEMf; j++) {
			if(fm->EMEf[j]>6) continue;
			if(!rot) a = abs(asin(fm->EMpxf[j]/fm->EMpf[j]))/3.14159*180;
			if(rot) a = abs(asin(fm->EMpzf[j]/fm->EMpf[j]))/3.14159*180;
			if(a < ph_min_angle) {
				if(fm->EMpdgf[j] == 22) isE = 0;
				else isE = 1;
				float l = r.Exp(18);
				if(isE == 0 && l>3) continue;
				l_ph_convert = l;
				ph_min_angle = a;
				Eph = fm->EMEf[j];
			}
        	}
		if(isE) l_ph_convert = 0;
		T->Fill();
		ph_min_angle = 91;
		l_angle = 91;
    	}
	T->Write();
	f.Close();
}

void MakePlot::CreateNCEMTree(const char* filename)
{
	TFile f(filename,"recreate");
    	TTree *T = new TTree("T", "New Tree from Fast MC");

	float l_angle = 91;
	float ph_min_angle = 91;
	float l_ph_convert = 0;
	float Eph = 0;
	float a, ap, ah;
	TRandom3 r;
	bool isE = 0, hasPi0 = 0;

    T->Branch("l_angle", &l_angle);
	T->Branch("EM_min_angle", &ph_min_angle);
	T->Branch("l_convert", &l_ph_convert);
    T->Branch("cc", &(fm->cc) );
	T->Branch("nc", &(fm->nc) );
	T->Branch("qel", &(fm->qel) );
	T->Branch("res", &(fm->res) );
	T->Branch("dis", &(fm->dis) );
    T->Branch("El", &(fm->El));
    T->Branch("Ev", &(fm->Ev));
	T->Branch("EEM", &Eph);
	T->Branch("isE", &isE);
	T->Branch("hasPi0", &hasPi0);
	T->Branch("POTWeight",&(fm->POTWeight));
	T->Branch("OSCWeight",&(fm->OSCWeight));
	T->Branch("n_events",&n_events);
	T->Branch("DetMass",&DetMass);
	T->Branch("exposure",&exposure);
	T->Branch("POTperYr",&POTperYr);

	int nEntries = fm->fChain->GetEntries();

    for (int i=0; i<nEntries; i++) {
        fm->fChain->GetEntry(i);
	    hasPi0 = 0;
	    for (int k=0; k<fm->nf; k++) {
       		if (fm->pf[k]>0 && fm->pdgf[k] == 111) hasPi0 = 1;
    	}
	    l_angle = asin(fm->pxl/fm->pl)/3.14159*180;
    	for (int j=0; j<fm->nEMf; j++) {
    		if(fm->EMEf[j]>6) continue;
    		if(!rot) a = asin(fm->EMpxf[j]/fm->EMpf[j])/3.14159*180;
    		else a = asin(fm->EMpzf[j]/fm->EMpf[j])/3.14159*180;
    		if(abs(a) < ph_min_angle) {
    			if(fm->EMpdgf[j] != 22) {
    				for (int k=0; k<fm->nf; k++) {
    					if (fm->pf[k]<=0 || fm->pdgf[k] == 2112 || fm->pdgf[k] == 111 || fm->pdgf[k] == 22) continue;
            					if(!rot) ah = asin(fm->pxf[k]/fm->pf[k])/3.14159*180;
    					if(rot) ah = asin(fm->pzf[k]/fm->pf[k])/3.14159*180;
            					if (abs(ah)>cutAngle) continue;
    					if (ah*a < 0) continue;
    					ph_min_angle = abs(a);
    					isE = 1;
    					Eph = fm->EMEf[j];
        				}
    			}
    			else {
        			for (int k=0; k<fm->nf; k++) {
           					if (fm->pf[k]<=0 || (fm->pdgf[k] != 211 && fm->pdgf[k] != -211)) continue;
            					if(!rot) ap = asin(fm->pxf[k]/fm->pf[k])/3.14159*180;
        				if(rot) ap = asin(fm->pzf[k]/fm->pf[k])/3.14159*180;
        				if (abs(ap)>cutAngle) continue;
        				if (ap*a<0) continue;
        				float l = r.Exp(18);
        				if (l <1 || l>5 ) continue;
        				l_ph_convert = l;
        				ph_min_angle = abs(a);
        				isE = 0;
        				Eph = fm->EMEf[j];
        			}
    			}
    		}
        }
	   if(isE) l_ph_convert = 0;
		T->Fill();
		ph_min_angle = 91;
		l_angle = 91;
    }
	T->Write();
	f.Close();
}

void MakePlot::PlotEvents_vs_El_l(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    gStyle->SetLegendBorderSize(0);

    TFile *fccs = new TFile(Form("numu-nue-cc-hadron_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tccs = (TTree*)fccs->Get("T");
    TFile *fccb = new TFile(Form("numu-nue-cc-EM_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tccb = (TTree*)fccb->Get("T");
    TFile *fncb = new TFile(Form("numu-numu-nc-EM_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tncb = (TTree*)fncb->Get("T");
    TCut phcut2("isE==0&&hasPi0");
    TCut ecut("isE==1&&hasPi0");
    TCut EMacut(Form("EM_min_angle<%f",cutAngle));
    TCut lgood(Form("l_angle>%f",cutAngle));
    TCut lbad(Form("l_angle<%f",cutAngle));

    cout<<endl<<"----------------------------------------------------------------------------"<<endl;
    cout<<"Lepton cut"<<endl;
    cout<<"Is rotated: "<<rot<<endl;
    cout<<"Angle cut: "<<cutAngle<<endl;
    cout<<"QE"<<endl;
    cout<<"----------------------------"<<endl;
    c1->Divide(3,1);
    THStack *hs = new THStack("hs","QE");
    TLegend *leg = new TLegend(.52,.68,.89,.89);
    c1->cd(1);
    tncb->Draw("EEM>>GNCbQE(40,0,6)", (phcut2 && EMacut && "nc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    TH1F *h = (TH1F*)gDirectory->Get("GNCbQE");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow#gamma+X NC [#theta_{#gamma}<%.1f#circ]",cutAngle),"f");
    tncb->Draw("EEM>>ENCbQE(40,0,6)", (ecut && EMacut && "nc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbQE");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow e+X NC [#theta_{e}<%.1f#circ]",cutAngle),"f");
    tccb->Draw("EEM>>ECCbQE(40,0,6)", (EMacut && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbQE");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow e/#gamma CC [#theta_{EM}<%.1f#circ]",cutAngle),"f");
    tccs->Draw("El>>ECCsQE(40,0,6)", (lbad && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsQE");
    h->SetFillColor(kRed);
    cout<<"CC Problem (lepton angle): "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#nu_{e} CC [#theta_{e}<%.1f#circ]",cutAngle),"f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsQEg(40,0,6)", (lgood && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsQEg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (lepton angle): "<<h->Integral()<<endl<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#nu_{e} CC [#theta_{e}>%.1f#circ]",cutAngle),"f");
    hs->SetTitle("QE;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();

    cout<<"Res"<<endl;
    cout<<"----------------------------"<<endl;
    hs = new THStack("hs","Res");
    //leg = new TLegend(.52,.7,.89,.89);
    c1->cd(2);
    tncb->Draw("EEM>>GNCbRes(40,0,6)", (phcut2 && EMacut && "nc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("GNCbRes");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{#mu} NC [#theta_{EM}<7.5#circ]","f");
    tncb->Draw("EEM>>ENCbRes(40,0,6)", (ecut && EMacut && "nc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbRes");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    tccb->Draw("EEM>>ECCbRes(40,0,6)", (EMacut && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbRes");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{e} CC [#theta_{EM}<7.5#circ]","f");
    tccs->Draw("El>>ECCsRes(40,0,6)", (lbad && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsRes");
    h->SetFillColor(kRed);
    cout<<"CC Problem (lepton angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}<7.5#circ]","f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsResg(40,0,6)", (lgood && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsResg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (lepton angle): "<<h->Integral()<<endl<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}>7.5#circ]","f");
    hs->SetTitle("Res;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();

    cout<<"Dis"<<endl;
    cout<<"----------------------------"<<endl;
    hs = new THStack("hs","Dis");
    //leg = new TLegend(.52,.7,.89,.89);
    c1->cd(3);
    tncb->Draw("EEM>>GNCbDis(40,0,6)", (phcut2 && EMacut && "nc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("GNCbDis");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{#mu} NC [#theta_{EM}<7.5#circ]","f");
    tncb->Draw("EEM>>ENCbDis(40,0,6)", (ecut && EMacut && "nc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbDis");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    tccb->Draw("EEM>>ECCbDis(40,0,6)", (EMacut && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbDis");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{e} CC [#theta_{EM}<7.5#circ]","f");
    tccs->Draw("El>>ECCsDis(40,0,6)", (lbad && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsDis");
    h->SetFillColor(kRed);
    cout<<"CC Problem (lepton angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}<7.5#circ]","f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsDisg(40,0,6)", (lgood && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsDisg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (lepton angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}>7.5#circ]","f");
    hs->SetTitle("Dis;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();
    cout<<"----------------------------------------------------------------------------"<<endl<<endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEvents_vs_El_h(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    gStyle->SetLegendBorderSize(0);

    TFile *fccs = new TFile(Form("numu-nue-cc-hadron_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tccs = (TTree*)fccs->Get("T");
    TFile *fccb = new TFile(Form("numu-nue-cc-EM_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tccb = (TTree*)fccb->Get("T");
    TFile *fncb = new TFile(Form("numu-numu-nc-EM_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tncb = (TTree*)fncb->Get("T");
    TCut phcut2("isE==0&&hasPi0");
    TCut ecut("isE==1&&hasPi0");
    TCut EMacut(Form("EM_min_angle<%f",cutAngle));
    TCut hgood(Form("h_min_angle>%f",cutAngle));
    TCut hbad(Form("h_min_angle<%f",cutAngle));

    cout<<endl<<"----------------------------------------------------------------------------"<<endl;
    cout<<"Hadron cut"<<endl;
    cout<<"Is rotated: "<<rot<<endl;
    cout<<"Angle cut: "<<cutAngle<<endl;
    cout<<"QE"<<endl;
    cout<<"----------------------------"<<endl;
    c1->Divide(3,1);
    THStack *hs = new THStack("hs","QE");
    TLegend *leg = new TLegend(.52,.68,.89,.89);
    c1->cd(1);
    tncb->Draw("EEM>>GNCbQE(40,0,6)", (phcut2 && EMacut && "nc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    TH1F *h = (TH1F*)gDirectory->Get("GNCbQE");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow#gamma+X NC [#theta_{#gamma}<%.1f#circ]",cutAngle),"f");
    tncb->Draw("EEM>>ENCbQE(40,0,6)", (ecut && EMacut && "nc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbQE");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow e+X NC [#theta_{e}<%.1f#circ]",cutAngle),"f");
    tccb->Draw("EEM>>ECCbQE(40,0,6)", (EMacut && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbQE");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow e/#gamma CC [#theta_{EM}<%.1f#circ]",cutAngle),"f");
    tccs->Draw("El>>ECCsQE(40,0,6)", (hbad && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsQE");
    h->SetFillColor(kRed);
    cout<<"CC Problem (hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#nu_{e} CC [#theta_{h}<%.1f#circ]",cutAngle),"f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsQEg(40,0,6)", (hgood && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsQEg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (hadron angle): "<<h->Integral()<<endl<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#nu_{e} CC [#theta_{h}>%.1f#circ]",cutAngle),"f");
    hs->SetTitle("QE;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();

    cout<<"Res"<<endl;
    cout<<"----------------------------"<<endl;
    hs = new THStack("hs","Res");
    //leg = new TLegend(.52,.7,.89,.89);
    c1->cd(2);
    tncb->Draw("EEM>>GNCbRes(40,0,6)", (phcut2 && EMacut && "nc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("GNCbRes");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{#mu} NC [#theta_{EM}<7.5#circ]","f");
    tncb->Draw("EEM>>ENCbRes(40,0,6)", (ecut && EMacut && "nc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbRes");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    tccb->Draw("EEM>>ECCbRes(40,0,6)", (EMacut && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbRes");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{e} CC [#theta_{EM}<7.5#circ]","f");
    tccs->Draw("El>>ECCsRes(40,0,6)", (hbad && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsRes");
    h->SetFillColor(kRed);
    cout<<"CC Problem (hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}<7.5#circ]","f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsResg(40,0,6)", (hgood && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsResg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (hadron angle): "<<h->Integral()<<endl<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}>7.5#circ]","f");
    hs->SetTitle("Res;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();

    cout<<"Dis"<<endl;
    cout<<"----------------------------"<<endl;
    hs = new THStack("hs","Dis");
    //leg = new TLegend(.52,.7,.89,.89);
    c1->cd(3);
    tncb->Draw("EEM>>GNCbDis(40,0,6)", (phcut2 && EMacut && "nc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("GNCbDis");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{#mu} NC [#theta_{EM}<7.5#circ]","f");
    tncb->Draw("EEM>>ENCbDis(40,0,6)", (ecut && EMacut && "nc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbDis");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    tccb->Draw("EEM>>ECCbDis(40,0,6)", (EMacut && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbDis");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{e} CC [#theta_{EM}<7.5#circ]","f");
    tccs->Draw("El>>ECCsDis(40,0,6)", (hbad && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsDis");
    h->SetFillColor(kRed);
    cout<<"CC Problem (hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}<7.5#circ]","f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsDisg(40,0,6)", (hgood && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsDisg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}>7.5#circ]","f");
    hs->SetTitle("Dis;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();
    cout<<"----------------------------------------------------------------------------"<<endl<<endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEvents_vs_El_lh(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    gStyle->SetLegendBorderSize(0);

    TFile *fccs = new TFile(Form("numu-nue-cc-hadron_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tccs = (TTree*)fccs->Get("T");
    TFile *fccb = new TFile(Form("numu-nue-cc-EM_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tccb = (TTree*)fccb->Get("T");
    TFile *fncb = new TFile(Form("numu-numu-nc-EM_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tncb = (TTree*)fncb->Get("T");
    TCut phcut2("isE==0&&hasPi0");
    TCut ecut("isE==1&&hasPi0");
    TCut EMacut(Form("EM_min_angle<%f",cutAngle));
    TCut lhgood(Form("l_angle>%f||h_min_angle>%f",cutAngle,cutAngle));
    TCut lhbad(Form("l_angle<%f&&h_min_angle<%f",cutAngle,cutAngle));

    cout<<endl<<"----------------------------------------------------------------------------"<<endl;
    cout<<"Lepton and Hadron cut"<<endl;
    cout<<"Is rotated: "<<rot<<endl;
    cout<<"Angle cut: "<<cutAngle<<endl;
    cout<<"QE"<<endl;
    cout<<"----------------------------"<<endl;
    c1->Divide(3,1);
    THStack *hs = new THStack("hs","QE");
    TLegend *leg = new TLegend(.52,.68,.89,.89);
    c1->cd(1);
    tncb->Draw("EEM>>GNCbQE(40,0,6)", (phcut2 && EMacut && "nc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    TH1F *h = (TH1F*)gDirectory->Get("GNCbQE");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow#gamma+X NC [#theta_{#gamma}<%.1f#circ]",cutAngle),"f");
    tncb->Draw("EEM>>ENCbQE(40,0,6)", (ecut && EMacut && "nc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbQE");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow e+X NC [#theta_{e}<%.1f#circ]",cutAngle),"f");
    tccb->Draw("EEM>>ECCbQE(40,0,6)", (EMacut && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbQE");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow e/#gamma CC [#theta_{EM}<%.1f#circ]",cutAngle),"f");
    tccs->Draw("El>>ECCsQE(40,0,6)", (lhbad && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsQE");
    h->SetFillColor(kRed);
    cout<<"CC Problem (lepton and hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#nu_{e} CC [#theta_{e},#theta_{h}<%.1f#circ]",cutAngle),"f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsQEg(40,0,6)", (lhgood && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsQEg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (lepton and hadron angle): "<<h->Integral()<<endl<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#nu_{e} CC [#theta_{e} or #theta_{h}>%.1f#circ]",cutAngle),"f");
    hs->SetTitle("QE;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();

    cout<<"Res"<<endl;
    cout<<"----------------------------"<<endl;
    hs = new THStack("hs","Res");
    //leg = new TLegend(.52,.7,.89,.89);
    c1->cd(2);
    tncb->Draw("EEM>>GNCbRes(40,0,6)", (phcut2 && EMacut && "nc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("GNCbRes");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{#mu} NC [#theta_{EM}<7.5#circ]","f");
    tncb->Draw("EEM>>ENCbRes(40,0,6)", (ecut && EMacut && "nc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbRes");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    tccb->Draw("EEM>>ECCbRes(40,0,6)", (EMacut && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbRes");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{e} CC [#theta_{EM}<7.5#circ]","f");
    tccs->Draw("El>>ECCsRes(40,0,6)", (lhbad && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsRes");
    h->SetFillColor(kRed);
    cout<<"CC Problem (lepton and hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}<7.5#circ]","f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsResg(40,0,6)", (lhgood && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsResg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (lepton and hadron angle): "<<h->Integral()<<endl<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}>7.5#circ]","f");
    hs->SetTitle("Res;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();

    cout<<"Dis"<<endl;
    cout<<"----------------------------"<<endl;
    hs = new THStack("hs","Dis");
    //leg = new TLegend(.52,.7,.89,.89);
    c1->cd(3);
    tncb->Draw("EEM>>GNCbDis(40,0,6)", (phcut2 && EMacut && "nc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("GNCbDis");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{#mu} NC [#theta_{EM}<7.5#circ]","f");
    tncb->Draw("EEM>>ENCbDis(40,0,6)", (ecut && EMacut && "nc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbDis");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    tccb->Draw("EEM>>ECCbDis(40,0,6)", (EMacut && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbDis");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{e} CC [#theta_{EM}<7.5#circ]","f");
    tccs->Draw("El>>ECCsDis(40,0,6)", (lhbad && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsDis");
    h->SetFillColor(kRed);
    cout<<"CC Problem (lepton and hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}<7.5#circ]","f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsDisg(40,0,6)", (lhgood && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsDisg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (lepton and hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}>7.5#circ]","f");
    hs->SetTitle("Dis;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();
    cout<<"----------------------------------------------------------------------------"<<endl<<endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::TabulateNCEvents()
{
	TFile f("NCcuts.root","recreate");
    	TTree *T = new TTree("T", "New Tree from Fast MC");

	float l_ph_convert = 0;
	float a, ap, ah;
	TRandom3 r;

	bool hasPi0 = 0, hasShallowE = 0, EhasH = 0, hasShallowG = 0, GdecR = 0, GhasPi = 0;
	int t = 1;

	T->Branch("t",&t);
    	T->Branch("cc", &(fm->cc) );
	T->Branch("nc", &(fm->nc) );
    	T->Branch("qel", &(fm->qel) );
    	T->Branch("res", &(fm->res) );
    	T->Branch("dis", &(fm->dis) );
	T->Branch("hasPi0", &hasPi0);
	T->Branch("hasShallowE", &hasShallowE);
	T->Branch("EhasH", &EhasH);
	T->Branch("hasShallowG", &hasShallowG);
	T->Branch("GdecR", &GdecR);
	T->Branch("GhasPi", &GhasPi);
	T->Branch("POTWeight",&(fm->POTWeight));
	T->Branch("OSCWeight",&(fm->OSCWeight));
	T->Branch("n_events",&n_events);
	T->Branch("DetMass",&DetMass);
	T->Branch("exposure",&exposure);
	T->Branch("POTperYr",&POTperYr);

	int nEntries = fm->fChain->GetEntries();

    	for (int i=0; i<nEntries; i++) {
        	fm->fChain->GetEntry(i);
		for (int k=0; k<fm->nf; k++) {
           		if (fm->pf[k]>0 && fm->pdgf[k] == 111) hasPi0 = 1;
        	}
        	for (int j=0; j<fm->nEMf; j++) {
			if(fm->EMEf[j]>6) continue;
			if(!rot) a = asin(fm->EMpxf[j]/fm->EMpf[j])/3.14159*180;
			if(rot) a = asin(fm->EMpzf[j]/fm->EMpf[j])/3.14159*180;
			if(fm->EMpdgf[j] != 22 && abs(a)<cutAngle)
			{
				hasShallowE = 1;
				for (int k=0; k<fm->nf; k++) {
					if (fm->pf[k]<=0 || fm->pdgf[k] == 2112 || fm->pdgf[k] == 111 || fm->pdgf[k] == 22) continue;
            				if(!rot) ah = asin(fm->pxf[k]/fm->pf[k])/3.14159*180;
					if(rot) ah = asin(fm->pzf[k]/fm->pf[k])/3.14159*180;
            				if (abs(ah)>cutAngle) continue;
					if (ah*a < 0) continue;
					EhasH = 1;
        				}
			}
			if(fm->EMpdgf[j] == 22 && abs(a)<cutAngle){
				hasShallowG = 1;
				float l = r.Exp(18);
				if (l >=1 && l<=5 ) GdecR = 1;
				for (int k=0; k<fm->nf; k++) {
            				if (fm->pf[k]<=0 || (fm->pdgf[k] != 211 && fm->pdgf[k] != -211)) continue;
            				if(!rot) ap = asin(fm->pxf[k]/fm->pf[k])/3.14159*180;
					if(rot) ap = asin(fm->pzf[k]/fm->pf[k])/3.14159*180;
					if (abs(ap)>cutAngle) continue;
					if (ap*a<0) continue;
					if (GdecR == 0 ) continue;
					GhasPi = 1;
				}
			}
        	}
		T->Fill();
		hasPi0 = 0; hasShallowE = 0; EhasH = 0; hasShallowG = 0; GdecR = 0; GhasPi = 0; l_ph_convert = 0;
    	}
	T->Write();

	TCanvas can("can","can",500,500);
	cout<<endl<<"----------------------------------------------------------------------------"<<endl;
    	cout<<"Is rotated: "<<rot<<endl;
    	cout<<"Angle cut: "<<cutAngle<<endl;
	T->Draw("t>>h1(40,0,6)","(dis||res||qel)*nc*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	TH1F *h = (TH1F*)gDirectory->Get("h1");
	cout<<"Total NC: "<<h->Integral()<<endl;
	T->Draw("t>>h2(40,0,6)","(dis||res||qel)*nc*hasPi0*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h2");
	cout<<"NC with pi0: "<<h->Integral()<<endl;
	T->Draw("t>>h3(40,0,6)","(dis||res||qel)*nc*hasPi0*hasShallowE*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h3");
	cout<<"NC with pi0 and shallow e: "<<h->Integral()<<endl;
	T->Draw("t>>h3p1(40,0,6)","(dis||res||qel)*nc*hasPi0*hasShallowE*EhasH*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h3p1");
	cout<<"NC with pi0 and shallow e and hadron: "<<h->Integral()<<endl;
	T->Draw("t>>h4(40,0,6)","(dis||res||qel)*nc*hasPi0*hasShallowG*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h4");
	cout<<"NC with pi0 and shallow gamma: "<<h->Integral()<<endl;
	T->Draw("t>>h4p1(40,0,6)","(dis||res||qel)*nc*hasPi0*hasShallowG*GdecR*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h4p1");
	cout<<"NC with pi0 and shallow gamma in decay range: "<<h->Integral()<<endl;
	T->Draw("t>>h4p2(40,0,6)","(dis||res||qel)*nc*hasPi0*hasShallowG*GdecR*GhasPi*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h4p2");
	cout<<"NC with pi0 and shallow gamma in decay range and pion: "<<h->Integral()<<endl;
	T->Draw("t>>h4p2a3p1(40,0,6)","(dis||res||qel)*nc*hasPi0*hasShallowG*GdecR*GhasPi*hasShallowE*EhasH*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h4p2a3p1");
	cout<<"NC with pi0 and shallow gamma in decay range and pion, and shallow e and hadron: "<<h->Integral()<<endl;
	f.Close();
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
    gStyle->SetTitleOffset(1.35, "y");
    gStyle->SetHistLineWidth(2);
    gStyle->SetHistLineColor(kBlack);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetPalette(1);
    gROOT->ForceStyle();
    gROOT->UseCurrentStyle();
    fm->fChain->UseCurrentStyle();

}
