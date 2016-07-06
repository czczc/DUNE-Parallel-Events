//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 20 17:23:24 2016 by ROOT version 5.34/32
// from TTree gst/GENIE Summary Event Tree
// found on file: fastmcNtp_20140711_lbne_g4lbnev3r2p4_nuflux_numuflux_nue_LAr_1_g280_Ar40_5000_Default.root
//////////////////////////////////////////////////////////

#ifndef FastMC_h
#define FastMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class FastMC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Char_t          FluxVersion;
   Int_t           iev;
   Int_t           neu;
   Int_t           fspl;
   Int_t           tgt;
   Int_t           Z;
   Int_t           A;
   Int_t           hitnuc;
   Int_t           hitqrk;
   Int_t           resid;
   Bool_t          sea;
   Bool_t          qel;
   Bool_t          mec;
   Bool_t          res;
   Bool_t          dis;
   Bool_t          coh;
   Bool_t          cohe;
   Bool_t          dfr;
   Bool_t          imd;
   Bool_t          imda;
   Bool_t          nuel;
   Bool_t          em;
   Bool_t          cc;
   Bool_t          nc;
   Bool_t          charm;
   Int_t           neut_code;
   Int_t           nuance_code;
   Double_t        wght;
   Double_t        POTWeight;
   Double_t        POTperYear;
   Double_t        nomPOTperYear;
   Double_t        OSCWeight;
   Double_t        xs;
   Double_t        ys;
   Double_t        ts;
   Double_t        Q2s;
   Double_t        Ws;
   Double_t        x;
   Double_t        y;
   Double_t        t;
   Double_t        Q2;
   Double_t        P;
   Double_t        ET;
   Double_t        W;
   Double_t        W2;
   Double_t        Ev;
   Double_t        Ehad;
   Double_t        Evtx;
   Double_t        p0v;
   Double_t        pxv;
   Double_t        pyv;
   Double_t        pzv;
   Double_t        x_reco;
   Double_t        y_reco;
   Double_t        t_reco;
   Double_t        Q2_reco;
   Double_t        QE_Q2_reco;
   Double_t        QE_Ev_reco;
   Double_t        W_reco;
   Double_t        W2_reco;
   Double_t        Ev_reco;
   Double_t        EvB_reco;
   Double_t        EvUB_reco;
   Double_t        Ehad_reco;
   Double_t        EhadB_reco;
   Double_t        EhadUB_reco;
   Double_t        Ehad_unbias_default;
   Double_t        Ehad_unbias_numu;
   Double_t        Ehad_unbias_numubar;
   Double_t        Ehad_unbias_nue;
   Double_t        Ehad_unbias_nuebar;
   Double_t        Ehad_unbias_nc;
   Double_t        Ehad_unbias_ncbar;
   Double_t        p0v_reco;
   Double_t        pxv_reco;
   Double_t        pyv_reco;
   Double_t        pzv_reco;
   Int_t           EvClass_reco;
   Int_t           isQE_reco;
   Bool_t          IsContv_reco;
   Int_t           hasElikeEM;
   Double_t        En;
   Double_t        p0n;
   Double_t        pxn;
   Double_t        pyn;
   Double_t        pzn;
   Int_t           Pdgl;
   Double_t        El;
   Double_t        p0l;
   Double_t        pxl;
   Double_t        pyl;
   Double_t        pzl;
   Double_t        pl;
   Double_t        cthl;
   Int_t           PdgL;
   Double_t        EL;
   Double_t        p0L;
   Double_t        pxL;
   Double_t        pyL;
   Double_t        pzL;
   Double_t        pL;
   Double_t        cthL;
   Int_t           ParentPdgL;
   Double_t        ParentEL;
   Double_t        p0LParent;
   Double_t        pxLParent;
   Double_t        pyLParent;
   Double_t        pzLParent;
   Double_t        pLParent;
   Double_t        cthLParent;
   Int_t           PdgL_reco;
   Double_t        EL_reco;
   Double_t        p0L_reco;
   Double_t        pxL_reco;
   Double_t        pyL_reco;
   Double_t        pzL_reco;
   Double_t        pL_reco;
   Double_t        cthL_reco;
   Double_t        TrkL_reco;
   Bool_t          IsContL_reco;
   Int_t           Pdgl_reco;
   Double_t        El_reco;
   Double_t        p0l_reco;
   Double_t        pxl_reco;
   Double_t        pyl_reco;
   Double_t        pzl_reco;
   Double_t        pl_reco;
   Double_t        Trkl_reco;
   Double_t        cthl_reco;
   Int_t           nfp;
   Int_t           nfn;
   Int_t           nfpip;
   Int_t           nfpim;
   Int_t           nfpi0;
   Int_t           nfkp;
   Int_t           nfkm;
   Int_t           nfk0;
   Int_t           nfem;
   Int_t           nfother;
   Int_t           nip;
   Int_t           nin;
   Int_t           nipip;
   Int_t           nipim;
   Int_t           nipi0;
   Int_t           nikp;
   Int_t           nikm;
   Int_t           nik0;
   Int_t           niem;
   Int_t           niother;
   Int_t           ni;
   Int_t           pdgi[20];   //[ni]
   Int_t           resc[20];   //[ni]
   Double_t        Ei[20];   //[ni]
   Double_t        p0i[20];   //[ni]
   Double_t        pxi[20];   //[ni]
   Double_t        pyi[20];   //[ni]
   Double_t        pzi[20];   //[ni]
   Int_t           nf;
   Int_t           pdgf[50];   //[nf]
   Double_t        Ef[50];   //[nf]
   Double_t        p0f[50];   //[nf]
   Double_t        pxf[50];   //[nf]
   Double_t        pyf[50];   //[nf]
   Double_t        pzf[50];   //[nf]
   Double_t        pf[50];   //[nf]
   Double_t        cthf[50];   //[nf]
   Double_t        Ef_reco[50];   //[nf]
   Double_t        p0f_reco[50];   //[nf]
   Double_t        pxf_reco[50];   //[nf]
   Double_t        pyf_reco[50];   //[nf]
   Double_t        pzf_reco[50];   //[nf]
   Double_t        Trkf_reco[50];   //[nf]
   Double_t        cthf_reco[50];   //[nf]
   Double_t        pT_Imbalance;
   Int_t           nEMf;
   Int_t           EMpdgf[22];   //[nEMf]
   Double_t        EMp0f[22];   //[nEMf]
   Double_t        EMpxf[22];   //[nEMf]
   Double_t        EMpyf[22];   //[nEMf]
   Double_t        EMpzf[22];   //[nEMf]
   Double_t        EMpf[22];   //[nEMf]
   Double_t        EMcthf[22];   //[nEMf]
   Double_t        EMEf[22];   //[nEMf]
   Double_t        EMEf_reco[22];   //[nEMf]
   Double_t        EMp0f_reco[22];   //[nEMf]
   Double_t        EMpxf_reco[22];   //[nEMf]
   Double_t        EMpyf_reco[22];   //[nEMf]
   Double_t        EMpzf_reco[22];   //[nEMf]
   Double_t        EMTrkf_reco[22];   //[nEMf]
   Double_t        EMcthf_reco[22];   //[nEMf]
   Double_t        EMPizMass_reco[22];   //[nEMf]
   Double_t        EMElPizM_reco[22];   //[nEMf]
   Double_t        vtxx;
   Double_t        vtxy;
   Double_t        vtxz;
   Double_t        vtxt;
   Double_t        calresp0;
   Double_t        caltotal;
   Double_t        calneutron;
   Double_t        calproton;
   Double_t        calpion;
   Double_t        calmeson;
   Double_t        calbaryon;
   Double_t        calhadron;
   Double_t        calem;
   Double_t        calother;
   Double_t        calbinding;
   Double_t        Tau_Prob_nue;
   Double_t        Tau_Prob_numu;
   Int_t           rw7;
   Double_t        rw7_qelCCvalenciaRPA[7];   //[rw7]
   Double_t        rw7_qelNCvalenciaRPA[7];   //[rw7]
   Double_t        rw7_reslowq2[7];   //[rw7]
   Int_t           frw7;
   Int_t           fNsyst;
   Double_t        frw7_TargetYTilt0_5mm[7];   //[frw7]
   Double_t        frw7_TargetYOffset0_5mm[7];   //[frw7]
   Double_t        frw7_TargetXTilt0_5mm[7];   //[frw7]
   Double_t        frw7_TargetXOffset0_5mm[7];   //[frw7]
   Double_t        frw7_TargetBe[7];   //[frw7]
   Double_t        frw7_ProtonP90GeV[7];   //[frw7]
   Double_t        frw7_ProtonP80GeV[7];   //[frw7]
   Double_t        frw7_ProtonP70GeV[7];   //[frw7]
   Double_t        frw7_ProtonP60GeV[7];   //[frw7]
   Double_t        frw7_ProtonP50GeV[7];   //[frw7]
   Double_t        frw7_ProtonP40GeV[7];   //[frw7]
   Double_t        frw7_ProtonP30GeV[7];   //[frw7]
   Double_t        frw7_ProtonP20GeV[7];   //[frw7]
   Double_t        frw7_ProtonP150GeV[7];   //[frw7]
   Double_t        frw7_ProtonP140GeV[7];   //[frw7]
   Double_t        frw7_ProtonP130GeV[7];   //[frw7]
   Double_t        frw7_ProtonP110GeV[7];   //[frw7]
   Double_t        frw7_ProtonP100GeV[7];   //[frw7]
   Double_t        frw7_Horn2YOffset0_5mm[7];   //[frw7]
   Double_t        frw7_Horn2XOffset0_5mm[7];   //[frw7]
   Double_t        frw7_Horn1YOffset0_5mm[7];   //[frw7]
   Double_t        frw7_Horn1XOffset0_5mm[7];   //[frw7]
   Double_t        frw7_DecayPipeRadius3m[7];   //[frw7]
   Double_t        frw7_DecayPipeRadius1_9m[7];   //[frw7]
   Double_t        frw7_DecayPipeLength250m[7];   //[frw7]
   Double_t        frw7_BeamSigmaY1_4mm[7];   //[frw7]
   Double_t        frw7_BeamSigmaX1_4mm[7];   //[frw7]
   Double_t        frw7_700kW[7];   //[frw7]
   Double_t        frw7_230kA[7];   //[frw7]

   // List of branches
   TBranch        *b_FluxVersion;   //!
   TBranch        *b_iev;   //!
   TBranch        *b_neu;   //!
   TBranch        *b_fspl;   //!
   TBranch        *b_tgt;   //!
   TBranch        *b_Z;   //!
   TBranch        *b_A;   //!
   TBranch        *b_hitnuc;   //!
   TBranch        *b_hitqrk;   //!
   TBranch        *b_resid;   //!
   TBranch        *b_sea;   //!
   TBranch        *b_qel;   //!
   TBranch        *b_mec;   //!
   TBranch        *b_res;   //!
   TBranch        *b_dis;   //!
   TBranch        *b_coh;   //!
   TBranch        *b_cohe;   //!
   TBranch        *b_dfr;   //!
   TBranch        *b_imd;   //!
   TBranch        *b_imda;   //!
   TBranch        *b_nuel;   //!
   TBranch        *b_em;   //!
   TBranch        *b_cc;   //!
   TBranch        *b_nc;   //!
   TBranch        *b_charm;   //!
   TBranch        *b_neut_code;   //!
   TBranch        *b_nuance_code;   //!
   TBranch        *b_wght;   //!
   TBranch        *b_POTWeight;   //!
   TBranch        *b_POTperYear;   //!
   TBranch        *b_nomPOTperYear;   //!
   TBranch        *b_OSCWeight;   //!
   TBranch        *b_xs;   //!
   TBranch        *b_ys;   //!
   TBranch        *b_ts;   //!
   TBranch        *b_Q2s;   //!
   TBranch        *b_Ws;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_t;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_P;   //!
   TBranch        *b_ET;   //!
   TBranch        *b_W;   //!
   TBranch        *b_W2;   //!
   TBranch        *b_Ev;   //!
   TBranch        *b_Ehad;   //!
   TBranch        *b_Evtx;   //!
   TBranch        *b_p0v;   //!
   TBranch        *b_pxv;   //!
   TBranch        *b_pyv;   //!
   TBranch        *b_pzv;   //!
   TBranch        *b_x_reco;   //!
   TBranch        *b_y_reco;   //!
   TBranch        *b_t_reco;   //!
   TBranch        *b_Q2_reco;   //!
   TBranch        *b_QE_Q2_reco;   //!
   TBranch        *b_QE_Ev_reco;   //!
   TBranch        *b_W_reco;   //!
   TBranch        *b_W2_reco;   //!
   TBranch        *b_Ev_reco;   //!
   TBranch        *b_EvB_reco;   //!
   TBranch        *b_EvUB_reco;   //!
   TBranch        *b_Ehad_reco;   //!
   TBranch        *b_EhadB_reco;   //!
   TBranch        *b_EhadUB_reco;   //!
   TBranch        *b_Ehad_unbias_default;   //!
   TBranch        *b_Ehad_unbias_numu;   //!
   TBranch        *b_Ehad_unbias_numubar;   //!
   TBranch        *b_Ehad_unbias_nue;   //!
   TBranch        *b_Ehad_unbias_nuebar;   //!
   TBranch        *b_Ehad_unbias_nc;   //!
   TBranch        *b_Ehad_unbias_ncbar;   //!
   TBranch        *b_p0v_reco;   //!
   TBranch        *b_pxv_reco;   //!
   TBranch        *b_pyv_reco;   //!
   TBranch        *b_pzv_reco;   //!
   TBranch        *b_EvClass_reco;   //!
   TBranch        *b_isQE_reco;   //!
   TBranch        *b_IsContv_reco;   //!
   TBranch        *b_hasElikeEM;   //!
   TBranch        *b_En;   //!
   TBranch        *b_p0n;   //!
   TBranch        *b_pxn;   //!
   TBranch        *b_pyn;   //!
   TBranch        *b_pzn;   //!
   TBranch        *b_Pdgl;   //!
   TBranch        *b_El;   //!
   TBranch        *b_p0l;   //!
   TBranch        *b_pxl;   //!
   TBranch        *b_pyl;   //!
   TBranch        *b_pzl;   //!
   TBranch        *b_pl;   //!
   TBranch        *b_cthl;   //!
   TBranch        *b_PdgL;   //!
   TBranch        *b_EL;   //!
   TBranch        *b_p0L;   //!
   TBranch        *b_pxL;   //!
   TBranch        *b_pyL;   //!
   TBranch        *b_pzL;   //!
   TBranch        *b_pL;   //!
   TBranch        *b_cthL;   //!
   TBranch        *b_ParentPdgL;   //!
   TBranch        *b_ParentEL;   //!
   TBranch        *b_p0LParent;   //!
   TBranch        *b_pxLParent;   //!
   TBranch        *b_pyLParent;   //!
   TBranch        *b_pzLParent;   //!
   TBranch        *b_pLParent;   //!
   TBranch        *b_cthLParent;   //!
   TBranch        *b_PdgL_reco;   //!
   TBranch        *b_EL_reco;   //!
   TBranch        *b_p0L_reco;   //!
   TBranch        *b_pxL_reco;   //!
   TBranch        *b_pyL_reco;   //!
   TBranch        *b_pzL_reco;   //!
   TBranch        *b_pL_reco;   //!
   TBranch        *b_cthL_reco;   //!
   TBranch        *b_TrkL_reco;   //!
   TBranch        *b_IsContL_reco;   //!
   TBranch        *b_Pdgl_reco;   //!
   TBranch        *b_El_reco;   //!
   TBranch        *b_p0l_reco;   //!
   TBranch        *b_pxl_reco;   //!
   TBranch        *b_pyl_reco;   //!
   TBranch        *b_pzl_reco;   //!
   TBranch        *b_pl_reco;   //!
   TBranch        *b_Trkl_reco;   //!
   TBranch        *b_cthl_reco;   //!
   TBranch        *b_nfp;   //!
   TBranch        *b_nfn;   //!
   TBranch        *b_nfpip;   //!
   TBranch        *b_nfpim;   //!
   TBranch        *b_nfpi0;   //!
   TBranch        *b_nfkp;   //!
   TBranch        *b_nfkm;   //!
   TBranch        *b_nfk0;   //!
   TBranch        *b_nfem;   //!
   TBranch        *b_nfother;   //!
   TBranch        *b_nip;   //!
   TBranch        *b_nin;   //!
   TBranch        *b_nipip;   //!
   TBranch        *b_nipim;   //!
   TBranch        *b_nipi0;   //!
   TBranch        *b_nikp;   //!
   TBranch        *b_nikm;   //!
   TBranch        *b_nik0;   //!
   TBranch        *b_niem;   //!
   TBranch        *b_niother;   //!
   TBranch        *b_ni;   //!
   TBranch        *b_pdgi;   //!
   TBranch        *b_resc;   //!
   TBranch        *b_Ei;   //!
   TBranch        *b_p0i;   //!
   TBranch        *b_pxi;   //!
   TBranch        *b_pyi;   //!
   TBranch        *b_pzi;   //!
   TBranch        *b_nf;   //!
   TBranch        *b_pdgf;   //!
   TBranch        *b_Ef;   //!
   TBranch        *b_p0f;   //!
   TBranch        *b_pxf;   //!
   TBranch        *b_pyf;   //!
   TBranch        *b_pzf;   //!
   TBranch        *b_pf;   //!
   TBranch        *b_cthf;   //!
   TBranch        *b_Ef_reco;   //!
   TBranch        *b_p0f_reco;   //!
   TBranch        *b_pxf_reco;   //!
   TBranch        *b_pyf_reco;   //!
   TBranch        *b_pzf_reco;   //!
   TBranch        *b_Trkf_reco;   //!
   TBranch        *b_cthf_reco;   //!
   TBranch        *b_pT_Imbalance;   //!
   TBranch        *b_nEMf;   //!
   TBranch        *b_EMpdgf;   //!
   TBranch        *b_EMp0f;   //!
   TBranch        *b_EMpxf;   //!
   TBranch        *b_EMpyf;   //!
   TBranch        *b_EMpzf;   //!
   TBranch        *b_EMpf;   //!
   TBranch        *b_EMcthf;   //!
   TBranch        *b_EMEf;   //!
   TBranch        *b_EMEf_reco;   //!
   TBranch        *b_EMp0f_reco;   //!
   TBranch        *b_EMpxf_reco;   //!
   TBranch        *b_EMpyf_reco;   //!
   TBranch        *b_EMpzf_reco;   //!
   TBranch        *b_EMTrkf_reco;   //!
   TBranch        *b_EMcthf_reco;   //!
   TBranch        *b_EMPizMass_reco;   //!
   TBranch        *b_EMElPizM_reco;   //!
   TBranch        *b_vtxx;   //!
   TBranch        *b_vtxy;   //!
   TBranch        *b_vtxz;   //!
   TBranch        *b_vtxt;   //!
   TBranch        *b_calresp0;   //!
   TBranch        *b_caltotal;   //!
   TBranch        *b_calneutron;   //!
   TBranch        *b_calproton;   //!
   TBranch        *b_calPion;   //!
   TBranch        *b_calmeson;   //!
   TBranch        *b_calbaryon;   //!
   TBranch        *b_calhadron;   //!
   TBranch        *b_calem;   //!
   TBranch        *b_calother;   //!
   TBranch        *b_calbinding;   //!
   TBranch        *b_Tau_Prob_nue;   //!
   TBranch        *b_Tau_Prob_numu;   //!
   TBranch        *b_rw7;   //!
   TBranch        *b_rw7_qelCCvalenciaRPA;   //!
   TBranch        *b_rw7_qelNCvalenciaRPA;   //!
   TBranch        *b_rw7_reslowq2;   //!
   TBranch        *b_frw7;   //!
   TBranch        *b_fNsyst;   //!
   TBranch        *b_frw7_TargetYTilt0_5mm;   //!
   TBranch        *b_frw7_TargetYOffset0_5mm;   //!
   TBranch        *b_frw7_TargetXTilt0_5mm;   //!
   TBranch        *b_frw7_TargetXOffset0_5mm;   //!
   TBranch        *b_frw7_TargetBe;   //!
   TBranch        *b_frw7_ProtonP90GeV;   //!
   TBranch        *b_frw7_ProtonP80GeV;   //!
   TBranch        *b_frw7_ProtonP70GeV;   //!
   TBranch        *b_frw7_ProtonP60GeV;   //!
   TBranch        *b_frw7_ProtonP50GeV;   //!
   TBranch        *b_frw7_ProtonP40GeV;   //!
   TBranch        *b_frw7_ProtonP30GeV;   //!
   TBranch        *b_frw7_ProtonP20GeV;   //!
   TBranch        *b_frw7_ProtonP150GeV;   //!
   TBranch        *b_frw7_ProtonP140GeV;   //!
   TBranch        *b_frw7_ProtonP130GeV;   //!
   TBranch        *b_frw7_ProtonP110GeV;   //!
   TBranch        *b_frw7_ProtonP100GeV;   //!
   TBranch        *b_frw7_Horn2YOffset0_5mm;   //!
   TBranch        *b_frw7_Horn2XOffset0_5mm;   //!
   TBranch        *b_frw7_Horn1YOffset0_5mm;   //!
   TBranch        *b_frw7_Horn1XOffset0_5mm;   //!
   TBranch        *b_frw7_DecayPipeRadius3m;   //!
   TBranch        *b_frw7_DecayPipeRadius1_9m;   //!
   TBranch        *b_frw7_DecayPipeLength250m;   //!
   TBranch        *b_frw7_BeamSigmaY1_4mm;   //!
   TBranch        *b_frw7_BeamSigmaX1_4mm;   //!
   TBranch        *b_frw7_700kW;   //!
   TBranch        *b_frw7_230kA;   //!

   FastMC(TTree *tree=0);
   virtual ~FastMC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef FastMC_cxx
FastMC::FastMC(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("fastmcNtp_20140711_lbne_g4lbnev3r2p4_nuflux_numuflux_nue_LAr_1_g280_Ar40_5000_Default.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data/fastmcNtp_20140711_lbne_g4lbnev3r2p4_nuflux_numuflux_nue_LAr_1_g280_Ar40_5000_Default.root");
      }
      f->GetObject("gst",tree);

   }
   Init(tree);
}

FastMC::~FastMC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FastMC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FastMC::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void FastMC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("FluxVersion", &FluxVersion, &b_FluxVersion);
   fChain->SetBranchAddress("iev", &iev, &b_iev);
   fChain->SetBranchAddress("neu", &neu, &b_neu);
   fChain->SetBranchAddress("fspl", &fspl, &b_fspl);
   fChain->SetBranchAddress("tgt", &tgt, &b_tgt);
   fChain->SetBranchAddress("Z", &Z, &b_Z);
   fChain->SetBranchAddress("A", &A, &b_A);
   fChain->SetBranchAddress("hitnuc", &hitnuc, &b_hitnuc);
   fChain->SetBranchAddress("hitqrk", &hitqrk, &b_hitqrk);
   fChain->SetBranchAddress("resid", &resid, &b_resid);
   fChain->SetBranchAddress("sea", &sea, &b_sea);
   fChain->SetBranchAddress("qel", &qel, &b_qel);
   fChain->SetBranchAddress("mec", &mec, &b_mec);
   fChain->SetBranchAddress("res", &res, &b_res);
   fChain->SetBranchAddress("dis", &dis, &b_dis);
   fChain->SetBranchAddress("coh", &coh, &b_coh);
   fChain->SetBranchAddress("cohe", &cohe, &b_cohe);
   fChain->SetBranchAddress("dfr", &dfr, &b_dfr);
   fChain->SetBranchAddress("imd", &imd, &b_imd);
   fChain->SetBranchAddress("imda", &imda, &b_imda);
   fChain->SetBranchAddress("nuel", &nuel, &b_nuel);
   fChain->SetBranchAddress("em", &em, &b_em);
   fChain->SetBranchAddress("cc", &cc, &b_cc);
   fChain->SetBranchAddress("nc", &nc, &b_nc);
   fChain->SetBranchAddress("charm", &charm, &b_charm);
   fChain->SetBranchAddress("neut_code", &neut_code, &b_neut_code);
   fChain->SetBranchAddress("nuance_code", &nuance_code, &b_nuance_code);
   fChain->SetBranchAddress("wght", &wght, &b_wght);
   fChain->SetBranchAddress("POTWeight", &POTWeight, &b_POTWeight);
   fChain->SetBranchAddress("POTperYear", &POTperYear, &b_POTperYear);
   fChain->SetBranchAddress("nomPOTperYear", &nomPOTperYear, &b_nomPOTperYear);
   fChain->SetBranchAddress("OSCWeight", &OSCWeight, &b_OSCWeight);
   fChain->SetBranchAddress("xs", &xs, &b_xs);
   fChain->SetBranchAddress("ys", &ys, &b_ys);
   fChain->SetBranchAddress("ts", &ts, &b_ts);
   fChain->SetBranchAddress("Q2s", &Q2s, &b_Q2s);
   fChain->SetBranchAddress("Ws", &Ws, &b_Ws);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("P", &P, &b_P);
   fChain->SetBranchAddress("ET", &ET, &b_ET);
   fChain->SetBranchAddress("W", &W, &b_W);
   fChain->SetBranchAddress("W2", &W2, &b_W2);
   fChain->SetBranchAddress("Ev", &Ev, &b_Ev);
   fChain->SetBranchAddress("Ehad", &Ehad, &b_Ehad);
   fChain->SetBranchAddress("Evtx", &Evtx, &b_Evtx);
   fChain->SetBranchAddress("p0v", &p0v, &b_p0v);
   fChain->SetBranchAddress("pxv", &pxv, &b_pxv);
   fChain->SetBranchAddress("pyv", &pyv, &b_pyv);
   fChain->SetBranchAddress("pzv", &pzv, &b_pzv);
   fChain->SetBranchAddress("x_reco", &x_reco, &b_x_reco);
   fChain->SetBranchAddress("y_reco", &y_reco, &b_y_reco);
   fChain->SetBranchAddress("t_reco", &t_reco, &b_t_reco);
   fChain->SetBranchAddress("Q2_reco", &Q2_reco, &b_Q2_reco);
   fChain->SetBranchAddress("QE_Q2_reco", &QE_Q2_reco, &b_QE_Q2_reco);
   fChain->SetBranchAddress("QE_Ev_reco", &QE_Ev_reco, &b_QE_Ev_reco);
   fChain->SetBranchAddress("W_reco", &W_reco, &b_W_reco);
   fChain->SetBranchAddress("W2_reco", &W2_reco, &b_W2_reco);
   fChain->SetBranchAddress("Ev_reco", &Ev_reco, &b_Ev_reco);
   fChain->SetBranchAddress("EvB_reco", &EvB_reco, &b_EvB_reco);
   fChain->SetBranchAddress("EvUB_reco", &EvUB_reco, &b_EvUB_reco);
   fChain->SetBranchAddress("Ehad_reco", &Ehad_reco, &b_Ehad_reco);
   fChain->SetBranchAddress("EhadB_reco", &EhadB_reco, &b_EhadB_reco);
   fChain->SetBranchAddress("EhadUB_reco", &EhadUB_reco, &b_EhadUB_reco);
   fChain->SetBranchAddress("Ehad_unbias_default", &Ehad_unbias_default, &b_Ehad_unbias_default);
   fChain->SetBranchAddress("Ehad_unbias_numu", &Ehad_unbias_numu, &b_Ehad_unbias_numu);
   fChain->SetBranchAddress("Ehad_unbias_numubar", &Ehad_unbias_numubar, &b_Ehad_unbias_numubar);
   fChain->SetBranchAddress("Ehad_unbias_nue", &Ehad_unbias_nue, &b_Ehad_unbias_nue);
   fChain->SetBranchAddress("Ehad_unbias_nuebar", &Ehad_unbias_nuebar, &b_Ehad_unbias_nuebar);
   fChain->SetBranchAddress("Ehad_unbias_nc", &Ehad_unbias_nc, &b_Ehad_unbias_nc);
   fChain->SetBranchAddress("Ehad_unbias_ncbar", &Ehad_unbias_ncbar, &b_Ehad_unbias_ncbar);
   fChain->SetBranchAddress("p0v_reco", &p0v_reco, &b_p0v_reco);
   fChain->SetBranchAddress("pxv_reco", &pxv_reco, &b_pxv_reco);
   fChain->SetBranchAddress("pyv_reco", &pyv_reco, &b_pyv_reco);
   fChain->SetBranchAddress("pzv_reco", &pzv_reco, &b_pzv_reco);
   fChain->SetBranchAddress("EvClass_reco", &EvClass_reco, &b_EvClass_reco);
   fChain->SetBranchAddress("isQE_reco", &isQE_reco, &b_isQE_reco);
   fChain->SetBranchAddress("IsContv_reco", &IsContv_reco, &b_IsContv_reco);
   fChain->SetBranchAddress("hasElikeEM", &hasElikeEM, &b_hasElikeEM);
   fChain->SetBranchAddress("En", &En, &b_En);
   fChain->SetBranchAddress("p0n", &p0n, &b_p0n);
   fChain->SetBranchAddress("pxn", &pxn, &b_pxn);
   fChain->SetBranchAddress("pyn", &pyn, &b_pyn);
   fChain->SetBranchAddress("pzn", &pzn, &b_pzn);
   fChain->SetBranchAddress("Pdgl", &Pdgl, &b_Pdgl);
   fChain->SetBranchAddress("El", &El, &b_El);
   fChain->SetBranchAddress("p0l", &p0l, &b_p0l);
   fChain->SetBranchAddress("pxl", &pxl, &b_pxl);
   fChain->SetBranchAddress("pyl", &pyl, &b_pyl);
   fChain->SetBranchAddress("pzl", &pzl, &b_pzl);
   fChain->SetBranchAddress("pl", &pl, &b_pl);
   fChain->SetBranchAddress("cthl", &cthl, &b_cthl);
   fChain->SetBranchAddress("PdgL", &PdgL, &b_PdgL);
   fChain->SetBranchAddress("EL", &EL, &b_EL);
   fChain->SetBranchAddress("p0L", &p0L, &b_p0L);
   fChain->SetBranchAddress("pxL", &pxL, &b_pxL);
   fChain->SetBranchAddress("pyL", &pyL, &b_pyL);
   fChain->SetBranchAddress("pzL", &pzL, &b_pzL);
   fChain->SetBranchAddress("pL", &pL, &b_pL);
   fChain->SetBranchAddress("cthL", &cthL, &b_cthL);
   fChain->SetBranchAddress("ParentPdgL", &ParentPdgL, &b_ParentPdgL);
   fChain->SetBranchAddress("ParentEL", &ParentEL, &b_ParentEL);
   fChain->SetBranchAddress("p0LParent", &p0LParent, &b_p0LParent);
   fChain->SetBranchAddress("pxLParent", &pxLParent, &b_pxLParent);
   fChain->SetBranchAddress("pyLParent", &pyLParent, &b_pyLParent);
   fChain->SetBranchAddress("pzLParent", &pzLParent, &b_pzLParent);
   fChain->SetBranchAddress("pLParent", &pLParent, &b_pLParent);
   fChain->SetBranchAddress("cthLParent", &cthLParent, &b_cthLParent);
   fChain->SetBranchAddress("PdgL_reco", &PdgL_reco, &b_PdgL_reco);
   fChain->SetBranchAddress("EL_reco", &EL_reco, &b_EL_reco);
   fChain->SetBranchAddress("p0L_reco", &p0L_reco, &b_p0L_reco);
   fChain->SetBranchAddress("pxL_reco", &pxL_reco, &b_pxL_reco);
   fChain->SetBranchAddress("pyL_reco", &pyL_reco, &b_pyL_reco);
   fChain->SetBranchAddress("pzL_reco", &pzL_reco, &b_pzL_reco);
   fChain->SetBranchAddress("pL_reco", &pL_reco, &b_pL_reco);
   fChain->SetBranchAddress("cthL_reco", &cthL_reco, &b_cthL_reco);
   fChain->SetBranchAddress("TrkL_reco", &TrkL_reco, &b_TrkL_reco);
   fChain->SetBranchAddress("IsContL_reco", &IsContL_reco, &b_IsContL_reco);
   fChain->SetBranchAddress("Pdgl_reco", &Pdgl_reco, &b_Pdgl_reco);
   fChain->SetBranchAddress("El_reco", &El_reco, &b_El_reco);
   fChain->SetBranchAddress("p0l_reco", &p0l_reco, &b_p0l_reco);
   fChain->SetBranchAddress("pxl_reco", &pxl_reco, &b_pxl_reco);
   fChain->SetBranchAddress("pyl_reco", &pyl_reco, &b_pyl_reco);
   fChain->SetBranchAddress("pzl_reco", &pzl_reco, &b_pzl_reco);
   fChain->SetBranchAddress("pl_reco", &pl_reco, &b_pl_reco);
   fChain->SetBranchAddress("Trkl_reco", &Trkl_reco, &b_Trkl_reco);
   fChain->SetBranchAddress("cthl_reco", &cthl_reco, &b_cthl_reco);
   fChain->SetBranchAddress("nfp", &nfp, &b_nfp);
   fChain->SetBranchAddress("nfn", &nfn, &b_nfn);
   fChain->SetBranchAddress("nfpip", &nfpip, &b_nfpip);
   fChain->SetBranchAddress("nfpim", &nfpim, &b_nfpim);
   fChain->SetBranchAddress("nfpi0", &nfpi0, &b_nfpi0);
   fChain->SetBranchAddress("nfkp", &nfkp, &b_nfkp);
   fChain->SetBranchAddress("nfkm", &nfkm, &b_nfkm);
   fChain->SetBranchAddress("nfk0", &nfk0, &b_nfk0);
   fChain->SetBranchAddress("nfem", &nfem, &b_nfem);
   fChain->SetBranchAddress("nfother", &nfother, &b_nfother);
   fChain->SetBranchAddress("nip", &nip, &b_nip);
   fChain->SetBranchAddress("nin", &nin, &b_nin);
   fChain->SetBranchAddress("nipip", &nipip, &b_nipip);
   fChain->SetBranchAddress("nipim", &nipim, &b_nipim);
   fChain->SetBranchAddress("nipi0", &nipi0, &b_nipi0);
   fChain->SetBranchAddress("nikp", &nikp, &b_nikp);
   fChain->SetBranchAddress("nikm", &nikm, &b_nikm);
   fChain->SetBranchAddress("nik0", &nik0, &b_nik0);
   fChain->SetBranchAddress("niem", &niem, &b_niem);
   fChain->SetBranchAddress("niother", &niother, &b_niother);
   fChain->SetBranchAddress("ni", &ni, &b_ni);
   fChain->SetBranchAddress("pdgi", pdgi, &b_pdgi);
   fChain->SetBranchAddress("resc", resc, &b_resc);
   fChain->SetBranchAddress("Ei", Ei, &b_Ei);
   fChain->SetBranchAddress("p0i", p0i, &b_p0i);
   fChain->SetBranchAddress("pxi", pxi, &b_pxi);
   fChain->SetBranchAddress("pyi", pyi, &b_pyi);
   fChain->SetBranchAddress("pzi", pzi, &b_pzi);
   fChain->SetBranchAddress("nf", &nf, &b_nf);
   fChain->SetBranchAddress("pdgf", pdgf, &b_pdgf);
   fChain->SetBranchAddress("Ef", Ef, &b_Ef);
   fChain->SetBranchAddress("p0f", p0f, &b_p0f);
   fChain->SetBranchAddress("pxf", pxf, &b_pxf);
   fChain->SetBranchAddress("pyf", pyf, &b_pyf);
   fChain->SetBranchAddress("pzf", pzf, &b_pzf);
   fChain->SetBranchAddress("pf", pf, &b_pf);
   fChain->SetBranchAddress("cthf", cthf, &b_cthf);
   fChain->SetBranchAddress("Ef_reco", Ef_reco, &b_Ef_reco);
   fChain->SetBranchAddress("p0f_reco", p0f_reco, &b_p0f_reco);
   fChain->SetBranchAddress("pxf_reco", pxf_reco, &b_pxf_reco);
   fChain->SetBranchAddress("pyf_reco", pyf_reco, &b_pyf_reco);
   fChain->SetBranchAddress("pzf_reco", pzf_reco, &b_pzf_reco);
   fChain->SetBranchAddress("Trkf_reco", Trkf_reco, &b_Trkf_reco);
   fChain->SetBranchAddress("cthf_reco", cthf_reco, &b_cthf_reco);
   fChain->SetBranchAddress("pT_Imbalance", &pT_Imbalance, &b_pT_Imbalance);
   fChain->SetBranchAddress("nEMf", &nEMf, &b_nEMf);
   fChain->SetBranchAddress("EMpdgf", EMpdgf, &b_EMpdgf);
   fChain->SetBranchAddress("EMp0f", EMp0f, &b_EMp0f);
   fChain->SetBranchAddress("EMpxf", EMpxf, &b_EMpxf);
   fChain->SetBranchAddress("EMpyf", EMpyf, &b_EMpyf);
   fChain->SetBranchAddress("EMpzf", EMpzf, &b_EMpzf);
   fChain->SetBranchAddress("EMpf", EMpf, &b_EMpf);
   fChain->SetBranchAddress("EMcthf", EMcthf, &b_EMcthf);
   fChain->SetBranchAddress("EMEf", EMEf, &b_EMEf);
   fChain->SetBranchAddress("EMEf_reco", EMEf_reco, &b_EMEf_reco);
   fChain->SetBranchAddress("EMp0f_reco", EMp0f_reco, &b_EMp0f_reco);
   fChain->SetBranchAddress("EMpxf_reco", EMpxf_reco, &b_EMpxf_reco);
   fChain->SetBranchAddress("EMpyf_reco", EMpyf_reco, &b_EMpyf_reco);
   fChain->SetBranchAddress("EMpzf_reco", EMpzf_reco, &b_EMpzf_reco);
   fChain->SetBranchAddress("EMTrkf_reco", EMTrkf_reco, &b_EMTrkf_reco);
   fChain->SetBranchAddress("EMcthf_reco", EMcthf_reco, &b_EMcthf_reco);
   fChain->SetBranchAddress("EMPizMass_reco", EMPizMass_reco, &b_EMPizMass_reco);
   fChain->SetBranchAddress("EMElPizM_reco", EMElPizM_reco, &b_EMElPizM_reco);
   fChain->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
   fChain->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
   fChain->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
   fChain->SetBranchAddress("vtxt", &vtxt, &b_vtxt);
   fChain->SetBranchAddress("calresp0", &calresp0, &b_calresp0);
   fChain->SetBranchAddress("caltotal", &caltotal, &b_caltotal);
   fChain->SetBranchAddress("calneutron", &calneutron, &b_calneutron);
   fChain->SetBranchAddress("calproton", &calproton, &b_calproton);
   fChain->SetBranchAddress("calpion", &calpion, &b_calPion);
   fChain->SetBranchAddress("calmeson", &calmeson, &b_calmeson);
   fChain->SetBranchAddress("calbaryon", &calbaryon, &b_calbaryon);
   fChain->SetBranchAddress("calhadron", &calhadron, &b_calhadron);
   fChain->SetBranchAddress("calem", &calem, &b_calem);
   fChain->SetBranchAddress("calother", &calother, &b_calother);
   fChain->SetBranchAddress("calbinding", &calbinding, &b_calbinding);
   fChain->SetBranchAddress("Tau_Prob_nue", &Tau_Prob_nue, &b_Tau_Prob_nue);
   fChain->SetBranchAddress("Tau_Prob_numu", &Tau_Prob_numu, &b_Tau_Prob_numu);
   fChain->SetBranchAddress("rw7", &rw7, &b_rw7);
   fChain->SetBranchAddress("rw7_qelCCvalenciaRPA", rw7_qelCCvalenciaRPA, &b_rw7_qelCCvalenciaRPA);
   fChain->SetBranchAddress("rw7_qelNCvalenciaRPA", rw7_qelNCvalenciaRPA, &b_rw7_qelNCvalenciaRPA);
   fChain->SetBranchAddress("rw7_reslowq2", rw7_reslowq2, &b_rw7_reslowq2);
   fChain->SetBranchAddress("frw7", &frw7, &b_frw7);
   fChain->SetBranchAddress("fNsyst", &fNsyst, &b_fNsyst);
   fChain->SetBranchAddress("frw7_TargetYTilt0.5mm", frw7_TargetYTilt0_5mm, &b_frw7_TargetYTilt0_5mm);
   fChain->SetBranchAddress("frw7_TargetYOffset0.5mm", frw7_TargetYOffset0_5mm, &b_frw7_TargetYOffset0_5mm);
   fChain->SetBranchAddress("frw7_TargetXTilt0.5mm", frw7_TargetXTilt0_5mm, &b_frw7_TargetXTilt0_5mm);
   fChain->SetBranchAddress("frw7_TargetXOffset0.5mm", frw7_TargetXOffset0_5mm, &b_frw7_TargetXOffset0_5mm);
   fChain->SetBranchAddress("frw7_TargetBe", frw7_TargetBe, &b_frw7_TargetBe);
   fChain->SetBranchAddress("frw7_ProtonP90GeV", frw7_ProtonP90GeV, &b_frw7_ProtonP90GeV);
   fChain->SetBranchAddress("frw7_ProtonP80GeV", frw7_ProtonP80GeV, &b_frw7_ProtonP80GeV);
   fChain->SetBranchAddress("frw7_ProtonP70GeV", frw7_ProtonP70GeV, &b_frw7_ProtonP70GeV);
   fChain->SetBranchAddress("frw7_ProtonP60GeV", frw7_ProtonP60GeV, &b_frw7_ProtonP60GeV);
   fChain->SetBranchAddress("frw7_ProtonP50GeV", frw7_ProtonP50GeV, &b_frw7_ProtonP50GeV);
   fChain->SetBranchAddress("frw7_ProtonP40GeV", frw7_ProtonP40GeV, &b_frw7_ProtonP40GeV);
   fChain->SetBranchAddress("frw7_ProtonP30GeV", frw7_ProtonP30GeV, &b_frw7_ProtonP30GeV);
   fChain->SetBranchAddress("frw7_ProtonP20GeV", frw7_ProtonP20GeV, &b_frw7_ProtonP20GeV);
   fChain->SetBranchAddress("frw7_ProtonP150GeV", frw7_ProtonP150GeV, &b_frw7_ProtonP150GeV);
   fChain->SetBranchAddress("frw7_ProtonP140GeV", frw7_ProtonP140GeV, &b_frw7_ProtonP140GeV);
   fChain->SetBranchAddress("frw7_ProtonP130GeV", frw7_ProtonP130GeV, &b_frw7_ProtonP130GeV);
   fChain->SetBranchAddress("frw7_ProtonP110GeV", frw7_ProtonP110GeV, &b_frw7_ProtonP110GeV);
   fChain->SetBranchAddress("frw7_ProtonP100GeV", frw7_ProtonP100GeV, &b_frw7_ProtonP100GeV);
   fChain->SetBranchAddress("frw7_Horn2YOffset0.5mm", frw7_Horn2YOffset0_5mm, &b_frw7_Horn2YOffset0_5mm);
   fChain->SetBranchAddress("frw7_Horn2XOffset0.5mm", frw7_Horn2XOffset0_5mm, &b_frw7_Horn2XOffset0_5mm);
   fChain->SetBranchAddress("frw7_Horn1YOffset0.5mm", frw7_Horn1YOffset0_5mm, &b_frw7_Horn1YOffset0_5mm);
   fChain->SetBranchAddress("frw7_Horn1XOffset0.5mm", frw7_Horn1XOffset0_5mm, &b_frw7_Horn1XOffset0_5mm);
   fChain->SetBranchAddress("frw7_DecayPipeRadius3m", frw7_DecayPipeRadius3m, &b_frw7_DecayPipeRadius3m);
   fChain->SetBranchAddress("frw7_DecayPipeRadius1.9m", frw7_DecayPipeRadius1_9m, &b_frw7_DecayPipeRadius1_9m);
   fChain->SetBranchAddress("frw7_DecayPipeLength250m", frw7_DecayPipeLength250m, &b_frw7_DecayPipeLength250m);
   fChain->SetBranchAddress("frw7_BeamSigmaY1.4mm", frw7_BeamSigmaY1_4mm, &b_frw7_BeamSigmaY1_4mm);
   fChain->SetBranchAddress("frw7_BeamSigmaX1.4mm", frw7_BeamSigmaX1_4mm, &b_frw7_BeamSigmaX1_4mm);
   fChain->SetBranchAddress("frw7_700kW", frw7_700kW, &b_frw7_700kW);
   fChain->SetBranchAddress("frw7_230kA", frw7_230kA, &b_frw7_230kA);
   Notify();
}

Bool_t FastMC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void FastMC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FastMC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef FastMC_cxx
