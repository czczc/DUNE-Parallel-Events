#ifndef MakePlot_h
#define MakePlot_h

class FastMC;
class TCanvas;

class MakePlot {
public :
   MakePlot();
   virtual ~MakePlot();

   FastMC *fm;
   TCanvas *c1;

   void PlotEv_vs_Angle(const char* saveAs=0);
   void PlotEl_vs_Angle(const char* saveAs=0);
   void PlotEl_vs_Length(const char* saveAs=0);


   void PlotEv(const char* saveAs=0);
   void PlotEl_vs_Ev(const char* saveAs=0);

   void PlotAngle_cut_E(const char* saveAs=0);
   void PlotAngle_cut_E_Perp(const char* saveAs=0);
   void PlotLength_cut_E(const char* saveAs=0);
   void PlotLength_cut_E_Perp(const char* saveAs=0);

   // hadron side
   void PlotEv_vs_Ef(const char* saveAs=0);
   void PlotNf_vs_Ratio(const char* particle, const char* title, const char* saveAs=0);
   void PlotAnglel_vs_Angleh(const char* saveAs=0);
   void PlotMinAngle_hadron(const char* saveAs=0);

   // create new tree
   void CreateTree(const char* filename);



   void SetStyle();

};

#endif