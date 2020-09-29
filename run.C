{
    gROOT->ProcessLine(".x loadClasses.C" );
    // FastMC x;
    // x.GetEntry(0);
    // cout << x.nf << endl;
    // // x.Loop();
    // cout << "done" << endl;

    double angle = 7.5;
    Bool_t rot = false;
    MakePlot a(angle, rot);
    // a.CreateTree("numu-nue-cc-hadron_angle.root");
    // a.CreateTree("numu-nue-cc-hadron_angle_perp.root");
    // a.Print("numu-nue-cc-hadron_angle.root");
    // a.Print("numu-nue-cc-hadron_angle_perp.root");
    // a.CreateTree("numu-numu-nc-hadron_angle.root");
    // a.CreateTree("numu-numu-nc-hadron_angle_perp.root");
    // a.Print("numu-numu-nc-hadron_angle.root");
    // a.Print("numu-numu-nc-hadron_angle_perp.root");

    // a.PlotEv("figs/Ev.png");
    // a.PlotOscWeight_vs_Ev("figs/OSCWeight_vs_Ev.png");
    // a.PlotEl_vs_Ev("figs/El_vs_Ev.png");
    a.PlotEv_vs_Angle("figs/Ev_vs_Angle.png");
    // a.PlotEl_vs_Angle("figs/El_vs_Angle.png");
    // a.PlotAngle_cut_E("figs/Angle_cut_E.png");
    // a.PlotAngle_cut_E_Perp("figs/Angle_cut_E_Perp.png");
    // a.PlotEv_langle_cut("figs/Ev_langle_cut.png");
    // a.PlotEv_langle_cut("figs/Ev_langle_cut_perp.png");
    // a.PlotEv_hangle_cut("figs/Ev_hangle_cut.png");
    // a.PlotEv_sb_overlapping("figs/Ev_sb_overlapping.png");
    // a.PlotEv_sb_overlapping("figs/Ev_sb_overlapping_perp.png");
    // a.PlotEv_sb_overlapping("figs/Ev_sb_overlapping_long_photon.png", true);
    // a.PlotEv_sb_overlapping("figs/Ev_sb_overlapping_long_photon_perp.png", true);

    // a.rot = true;
    // a.PlotEv_langle_cut("figs/Ev_langle_cut_perp.png");
    // a.PlotEv_hangle_cut("figs/Ev_hangle_cut_perp.png");
    // a.PlotEv_hangle_cut("figs/Ev_bothangle_cut.png", true);
    // a.PlotEv_hangle_cut("figs/Ev_bothangle_cut_perp.png", true);



    //a.PlotEl_vs_Length("figs/El_vs_Length.png");
    //a.PlotLength_cut_E("figs/Length_cut_E.png");
    //a.PlotLength_cut_E_Perp("figs/Length_cut_E_Perp.png");

    // a.PlotEv_vs_Ef("figs/Ev_vs_Ef.png");
    // a.PlotNf_vs_Ratio("nfp", "# protons / total", "figs/Nf_vs_Ratio_proton.png");
    // a.PlotNf_vs_Ratio("nfn", "# neutrons / total", "figs/Nf_vs_Ratio_neutron.png");
    // a.PlotNf_vs_Ratio("nfn+nfp", "(neutron + proton) / total", "figs/Nf_vs_Ratio_np.png");
    // a.PlotNf_vs_Ratio("nfpip+nfpim", "(#pi^{+} + #pi^{-}) / total", "figs/Nf_vs_Ratio_pipm.png");
    // a.PlotNf_vs_Ratio("nfpi0+nfem", "(#pi^{0} + gamma) / total", "figs/Nf_vs_Ratio_pi0g.png");
    // a.PlotNf_vs_Ratio("nfkp+nfkm", "(K^{+} + K^{-}) / total", "figs/Nf_vs_Ratio_kpm.png");
    // a.PlotNf_vs_Ratio("nfk0", "(K^{0} + #barK^{0}) / total", "figs/Nf_vs_Ratio_k0.png");
    // a.PlotAnglel_vs_Angleh("figs/Anglel_vs_Angleh.png");
    // a.PlotMinAngle_hadron("figs/MinAngle_hadron.png");
    // a.PlotMinAngle_hadron("figs/MinAngle_hadron_perp.png");
    // a.PlotMinAngle_lepton("figs/MinAngle_lepton.png");
    // a.CreateTree(Form("numu-nue-cc-hadron_angle_%.1f_%i.root",angle,rot));
    //a.CreateCCEMTree(Form("numu-nue-cc-EM_angle_%.1f_%i.root",angle,rot));//Be sure to reference correct MC file in FastMC.C
    //a.CreateNCEMTree(Form("numu-numu-nc-EM_angle_%.1f_%i.root",angle,rot));//Be sure to reference correct MC file in FastMC.C
    // a.TabulateNCEvents();

    // a.CalcOverlap();

	/*MakePlot a;
  	Bool_t rot;
  	for(double angle = 7.5; angle >=2.5; angle-=2.5)
	{
		rot = 0;
		a = MakePlot(angle, rot);
		a.PlotEvents_vs_El_l(Form("figs/SignalBackgrounds_l_%.1f_%i.png",angle,rot));
    		a.PlotEvents_vs_El_h(Form("figs/SignalBackgrounds_h_%.1f_%i.png",angle,rot));
    		a.PlotEvents_vs_El_lh(Form("figs/SignalBackgrounds_lh_%.1f_%i.png",angle,rot));
    		//a.TabulateNCEvents();
		rot = 1;
		a = MakePlot(angle, rot);
		a.PlotEvents_vs_El_l(Form("figs/SignalBackgrounds_l_%.1f_%i.png",angle,rot));
    		a.PlotEvents_vs_El_h(Form("figs/SignalBackgrounds_h_%.1f_%i.png",angle,rot));
    		a.PlotEvents_vs_El_lh(Form("figs/SignalBackgrounds_lh_%.1f_%i.png",angle,rot));
    		//a.TabulateNCEvents();

	}*/
}
