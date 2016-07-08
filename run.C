{
    // gROOT->ProcessLine(".x loadClasses.C" );
    // FastMC x;
    // x.GetEntry(0);
    // cout << x.nf << endl;
    // // x.Loop();
    // cout << "done" << endl;

    MakePlot a;
    // a.fm->GetEntry(0);
    // cout << a.fm->nf << endl;
    // a.fm->fChain->Draw("Ev");
    // a.PlotEv_vs_Angle("figs/Ev_vs_Angle.png");
    // a.PlotEl_vs_Angle("figs/El_vs_Angle.png");
    // a.PlotEv("figs/Ev.png");
    // a.PlotEl_vs_Ev("figs/El_vs_Ev.png");
    // a.PlotAngle_cut_E("figs/Angle_cut_E.png");
    // a.PlotAngle_cut_E_Perp("figs/Angle_cut_E_Perp.png");
    // a.PlotEl_vs_Length("figs/El_vs_Length.png");
    // a.PlotLength_cut_E("figs/Length_cut_E.png");
    // a.PlotLength_cut_E_Perp("figs/Length_cut_E_Perp.png");

    // a.PlotEv_vs_Ef("figs/Ev_vs_Ef.png");
    // a.PlotNf_vs_Ratio("nfp", "# protons / total", "figs/Nf_vs_Ratio_proton.png");
    // a.PlotNf_vs_Ratio("nfn", "# neutrons / total", "figs/Nf_vs_Ratio_neutron.png");
    // a.PlotNf_vs_Ratio("nfn+nfp", "(neutron + proton) / total", "figs/Nf_vs_Ratio_np.png");
    // a.PlotNf_vs_Ratio("nfpip+nfpim", "(#pi^{+} + #pi^{-}) / total", "figs/Nf_vs_Ratio_pipm.png");
    // a.PlotNf_vs_Ratio("nfpi0+nfem", "(#pi^{0} + gamma) / total", "figs/Nf_vs_Ratio_pi0g.png");
    // a.PlotNf_vs_Ratio("nfkp+nfkm", "(K^{+} + K^{-}) / total", "figs/Nf_vs_Ratio_kpm.png");
    // a.PlotNf_vs_Ratio("nfk0", "(K^{0} + #barK^{0}) / total", "figs/Nf_vs_Ratio_k0.png");
    // a.PlotAnglel_vs_Angleh("figs/Anglel_vs_Angleh.png");
    // a.CreateTree("ah.root");
    a.PlotMinAngle_hadron("figs/MinAngle_hadron.png");
}
