#include <cstdlib>
#include <cmath>
#include <string>
#include "DrawBase.h"
#include "QG/QGLikelihood/interface/QGLikelihoodCalculator.h"

bool Summer12=true;
std::string plotsdir = (Summer12) ? "plots_Summer12" : "plots";



void drawSinglePtBin( DrawBase* db, QGLikelihoodCalculator* qglc, TTree* tree, float ptMin, float ptMax );
void drawPlot( DrawBase* db, TH1D* h1_gluon, TH1D* h1_quark, std::string name, float ptMin, float ptMax );
void drawRoC( DrawBase* db, float ptMin, float ptMax, const std::string& flag, TH1D* h1_new_gluon, TH1D* h1_new_quark, TH1D* h1_old_gluon, TH1D* h1_old_quark, TH1D* h1_bdt_gluon=0, TH1D* h1_bdt_quark=0);



int main() {

  DrawBase* db = new DrawBase("checkQG");

  //QGLikelihoodCalculator* qglc = new QGLikelihoodCalculator("/afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src/UserCode/pandolf/QGDev/Fit/Output/Histos.root");
  QGLikelihoodCalculator* qglc;
  if( Summer12 ) qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/CMSSW_5_3_6/src/QG/QGLikelihood/test/Histos_2012.root");
  else           qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/CMSSW_5_3_6/src/QG/QGLikelihood/test/Histos.root");

  TFile* file;
  if( Summer12 ) file = TFile::Open("/afs/cern.ch/work/a/amarini/2ndLevel/QG/QG/QG_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_2_TREE.root");
  else           file = TFile::Open("/afs/cern.ch/work/a/amarini/2ndLevel/QG/QG/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_TREE.root");


  TTree* tree = (TTree*)file->Get("tree_passedEvents");

  std::string mkdircommand = "mkdir -p " + plotsdir;
  system(mkdircommand.c_str());

  drawSinglePtBin( db, qglc, tree, 20., 25. );
  drawSinglePtBin( db, qglc, tree, 25., 30. );
  drawSinglePtBin( db, qglc, tree, 30., 40. );
  drawSinglePtBin( db, qglc, tree, 40., 50. );
  drawSinglePtBin( db, qglc, tree, 80., 100. );
  drawSinglePtBin( db, qglc, tree, 150., 200. );
  drawSinglePtBin( db, qglc, tree, 200., 250. );
  drawSinglePtBin( db, qglc, tree, 300., 400. );
  drawSinglePtBin( db, qglc, tree, 500., 600. );
  drawSinglePtBin( db, qglc, tree, 800., 1000. );

  return 0;

}




void drawSinglePtBin( DrawBase* db, QGLikelihoodCalculator* qglc, TTree* tree, float ptMin, float ptMax ) {

  std::cout << "-> Processing pt bin: " << ptMin << "-" << ptMax << " GeV..." << std::endl;

  bool doFwd = (ptMin<100.);


  float pt;
  tree->SetBranchAddress("ptJet0", &pt);
  float eta;
  tree->SetBranchAddress("etaJet0", &eta);
  int pdgId;
  tree->SetBranchAddress("pdgIdPartJet0", &pdgId);
  float rho;
  tree->SetBranchAddress("rhoPF", &rho);
  int nCharged;
  tree->SetBranchAddress("nChargedJet0", &nCharged);
  int nNeutral;
  tree->SetBranchAddress("nNeutralJet0", &nNeutral);
  float ptD;
  tree->SetBranchAddress("ptDJet0", &ptD);
  float ptD_QC;
  tree->SetBranchAddress("ptD_QCJet0", &ptD_QC);
  float axis2_QC;
  tree->SetBranchAddress("axis2_QCJet0", &axis2_QC);
  int nCharged_QC;
  tree->SetBranchAddress("nChg_QCJet0", &nCharged_QC);
  int nNeutral_ptCut;
  tree->SetBranchAddress("nNeutral_ptCutJet0", &nNeutral_ptCut);
  float qglPaoloJet0;
  tree->SetBranchAddress("qglPaoloJet0", &qglPaoloJet0);


  TH1D* h1_qgl_old_gluon = new TH1D("qgl_old_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_old_quark = new TH1D("qgl_old_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgl_new_gluon = new TH1D("qgl_new_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_new_quark = new TH1D("qgl_new_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgbdt_gluon = new TH1D("qgbdt_gluon", "", 100, -0.5, 0.5);
  TH1D* h1_qgbdt_quark = new TH1D("qgbdt_quark", "", 100, -0.5, 0.5);


  // in the forward:
  TH1D* h1_qgl_new_F_gluon = new TH1D("qgl_new_F_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_new_F_quark = new TH1D("qgl_new_F_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgbdt_F_gluon = new TH1D("qgbdt_F_gluon", "", 100, -0.5, 0.5);
  TH1D* h1_qgbdt_F_quark = new TH1D("qgbdt_F_quark", "", 100, -0.5, 0.5);


  int nentries = tree->GetEntries();

  for( unsigned int ientry=0; ientry<nentries; ++ientry ) {

    tree->GetEntry(ientry);

    if( pt<ptMin || pt>ptMax ) continue;
    //if( rho>22. ) continue;

    float qgl_new = qglc->computeQGLikelihood2012( pt, eta, rho, nCharged_QC+nNeutral_ptCut, ptD_QC, axis2_QC);

    if( fabs(eta)<2.5 && h1_qgl_old_gluon->GetEntries()<10000 && h1_qgl_old_quark->GetEntries()<10000 ) { //save time

      float qgl_old = qglc->computeQGLikelihoodPU( pt, rho, nCharged, nNeutral, ptD);

      if( fabs(pdgId)<5 ) {
        h1_qgl_old_quark->Fill( qgl_old );
        h1_qgl_new_quark->Fill( qgl_new );
        h1_qgbdt_quark->Fill( qglPaoloJet0 );
      }
      if( pdgId==21 ) {
        h1_qgl_old_gluon->Fill( qgl_old );
        h1_qgl_new_gluon->Fill( qgl_new );
        h1_qgbdt_gluon->Fill( qglPaoloJet0 );
      }

    } else if( fabs(eta)>2.5 ) {

      if( fabs(pdgId)<5 ) {
        h1_qgl_new_F_quark->Fill( qgl_new );
        h1_qgbdt_F_quark->Fill( qglPaoloJet0 );
      }
      if( pdgId==21 ) {
        h1_qgl_new_F_gluon->Fill( qgl_new );
        h1_qgbdt_F_gluon->Fill( qglPaoloJet0 );
      }

    }
    
    if( h1_qgl_old_gluon->GetEntries()>10000 
     && h1_qgl_old_quark->GetEntries()>10000 
     && ( !doFwd || (h1_qgl_new_F_quark->GetEntries()>10000
     && h1_qgl_new_F_gluon->GetEntries()>10000) ) ) break;

  }


  drawPlot( db, h1_qgl_old_gluon, h1_qgl_old_quark, "old", ptMin, ptMax );
  drawPlot( db, h1_qgl_new_gluon, h1_qgl_new_quark, "new", ptMin, ptMax );
  drawPlot( db, h1_qgbdt_gluon, h1_qgbdt_quark, "bdt", ptMin, ptMax );

  drawPlot( db, h1_qgl_new_F_gluon, h1_qgl_new_F_quark, "new_F", ptMin, ptMax );
  drawPlot( db, h1_qgbdt_F_gluon, h1_qgbdt_F_quark, "bdt_F", ptMin, ptMax );

  drawRoC(db, ptMin, ptMax, "", h1_qgl_new_gluon, h1_qgl_new_quark, h1_qgl_old_gluon, h1_qgl_old_quark, 0, 0);
  drawRoC(db, ptMin, ptMax, "_withBDT", h1_qgl_new_gluon, h1_qgl_new_quark, h1_qgl_old_gluon, h1_qgl_old_quark, h1_qgbdt_gluon, h1_qgbdt_quark);
  drawRoC(db, ptMin, ptMax, "_F", h1_qgl_new_gluon, h1_qgl_new_quark, 0, 0, h1_qgbdt_F_gluon, h1_qgbdt_F_quark);

  delete h1_qgl_old_gluon;
  delete h1_qgl_old_quark;

  delete h1_qgl_new_gluon;
  delete h1_qgl_new_quark;

  delete h1_qgbdt_gluon;
  delete h1_qgbdt_quark;

  delete h1_qgl_new_F_gluon;
  delete h1_qgl_new_F_quark;

  delete h1_qgbdt_F_gluon;
  delete h1_qgbdt_F_quark;


}





void drawPlot( DrawBase* db, TH1D* h1_gluon, TH1D* h1_quark, std::string name, float ptMin, float ptMax ) {


  h1_quark->Rebin(2);
  h1_gluon->Rebin(2);

  float norm_max_g = h1_gluon->GetMaximum()/h1_gluon->Integral();
  float norm_max_q = h1_quark->GetMaximum()/h1_quark->Integral();
  float hmax = (norm_max_q>norm_max_g) ? norm_max_q : norm_max_g;

  float ymax = hmax*1.2;

  TH2D* h2_axes = new TH2D("axes", "", 10, h1_gluon->GetXaxis()->GetXmin(), h1_gluon->GetXaxis()->GetXmax(), 10, 0., ymax);
  h2_axes->SetXTitle("Quark-Gluon Likelihood Discriminator");
  h2_axes->SetYTitle("Normalized To Unity");

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h1_quark->SetLineColor(38);
  h1_quark->SetLineWidth(2);
  h1_quark->SetFillColor(38);
  h1_quark->SetFillStyle(3005);

  h1_gluon->SetLineColor(46);
  h1_gluon->SetLineWidth(2);
  h1_gluon->SetFillColor(46);
  h1_gluon->SetFillStyle(3004);


  h2_axes->Draw();

  h1_quark->DrawNormalized("same");
  h1_gluon->DrawNormalized("same");

  char legendTitle[300];
  sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax );
  TLegend* legend = new TLegend( 0.55, 0.7, 0.8, 0.9, legendTitle );
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry( h1_quark, "Quark Jets", "F");
  legend->AddEntry( h1_gluon, "Gluon Jets", "F");
  legend->Draw("same");
  
  h1_quark->DrawNormalized("same");
  h1_gluon->DrawNormalized("same");



  TPaveText* label = db->get_labelTop();
  label->Draw("same");

  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "%s/qgl_%s_pt%.0f_%.0f.eps", plotsdir.c_str(), name.c_str(), ptMin, ptMax);
  c1->SaveAs(canvasName);
  sprintf( canvasName, "%s/qgl_%s_pt%.0f_%.0f.png", plotsdir.c_str(), name.c_str(), ptMin, ptMax);
  c1->SaveAs(canvasName);

  delete c1;
  delete h2_axes;
  delete legend;

}


void drawRoC( DrawBase* db, float ptMin, float ptMax, const std::string& flag, TH1D* h1_new_gluon, TH1D* h1_new_quark, TH1D* h1_old_gluon, TH1D* h1_old_quark, TH1D* h1_bdt_gluon, TH1D* h1_bdt_quark) {


  TGraph* gr_RoC_old = new TGraph(0);
  TGraph* gr_RoC_new = new TGraph(0);
  TGraph* gr_RoC_bdt = new TGraph(0);

  int nbins = h1_new_quark->GetNbinsX();

  for( unsigned int ibin=1; ibin<nbins+1; ++ibin ) {

    float eff_q_old = -1.;
    float eff_g_old = -1.;
  
    if( h1_old_quark!=0 && h1_old_gluon!=0 ) {
      eff_q_old = h1_old_quark->Integral( nbins-ibin, nbins )/h1_old_quark->Integral( 1, nbins );
      eff_g_old = h1_old_gluon->Integral( nbins-ibin, nbins )/h1_old_gluon->Integral( 1, nbins );
    }
  
    float eff_q_bdt = -1.;
    float eff_g_bdt = -1.;
  
    if( h1_bdt_quark!=0 && h1_bdt_gluon!=0 ) { //opposite convention:
      eff_q_bdt = h1_bdt_quark->Integral( 1, ibin )/h1_bdt_quark->Integral( 1, nbins );
      eff_g_bdt = h1_bdt_gluon->Integral( 1, ibin )/h1_bdt_gluon->Integral( 1, nbins );
    }
  
    float eff_q_new = h1_new_quark->Integral( nbins-ibin, nbins )/h1_new_quark->Integral( 1, nbins );
    float eff_g_new = h1_new_gluon->Integral( nbins-ibin, nbins )/h1_new_gluon->Integral( 1, nbins );
  
    gr_RoC_new->SetPoint( ibin-1, 1.-eff_g_new, eff_q_new );

    if( h1_old_quark!=0 && h1_old_gluon!=0 ) 
      gr_RoC_old->SetPoint( ibin-1, 1.-eff_g_old, eff_q_old );

    if( h1_bdt_quark!=0 && h1_bdt_gluon!=0 ) 
      gr_RoC_bdt->SetPoint( ibin-1, 1.-eff_g_bdt, eff_q_bdt );

  }


  gr_RoC_new->SetMarkerSize(1.3);
  gr_RoC_new->SetMarkerStyle(24);
  gr_RoC_new->SetMarkerColor(kRed+3);

  if( h1_old_quark!=0 && h1_old_gluon!=0 ) {
    gr_RoC_old->SetMarkerSize(1.3);
    gr_RoC_old->SetMarkerStyle(20);
    gr_RoC_old->SetMarkerColor(kOrange+1);
  }

  if( h1_bdt_quark!=0 && h1_bdt_gluon!=0 ) {
    gr_RoC_bdt->SetMarkerSize(1.3);
    gr_RoC_bdt->SetMarkerStyle(21);
    gr_RoC_bdt->SetMarkerColor(38);
  }

  TCanvas* c1 = new TCanvas("c1_roc", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes_roc", "", 10, 0., 1.0001, 10, 0., 1.0001);
  h2_axes->SetXTitle( "Gluon Jet Rejection" );
  h2_axes->SetYTitle( "Quark Jet Efficiency" );

  h2_axes->Draw();

  TLine* diag = new TLine(0., 1., 1., 0.);
  diag->Draw("same");


  char legendTitle[300];
  sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax );
  TLegend* legend = new TLegend( 0.2, 0.2, 0.45, 0.45, legendTitle );
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  if( h1_old_quark!=0 && h1_old_gluon!=0 )
    legend->AddEntry( gr_RoC_old, "Old LD", "P");
  legend->AddEntry( gr_RoC_new, "New LD", "P");
  if( h1_bdt_quark!=0 && h1_bdt_gluon!=0 )
    legend->AddEntry( gr_RoC_bdt, "BDT", "P");
  legend->Draw("same");

  TPaveText* label = db->get_labelTop();
  label->Draw("same");
  
  if( h1_bdt_quark!=0 && h1_bdt_gluon!=0 ) 
    gr_RoC_bdt->Draw("p same");
  if( h1_old_quark!=0 && h1_old_gluon!=0 ) 
    gr_RoC_old->Draw("p same");
  gr_RoC_new->Draw("p same");

  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "%s/RoC_pt%.0f_%.0f%s.eps", plotsdir.c_str(), ptMin, ptMax, flag.c_str());
  c1->SaveAs(canvasName);
  sprintf( canvasName, "%s/RoC_pt%.0f_%.0f%s.png", plotsdir.c_str(), ptMin, ptMax, flag.c_str());
  c1->SaveAs(canvasName);

  delete c1;
  delete h2_axes;
  delete legend;
}
