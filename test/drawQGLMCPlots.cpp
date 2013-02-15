#include <iostream>
#include <cmath>
#include <cstdlib>
#include "DrawBase.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "QG/QGLikelihood/interface/QGLikelihoodCalculator.h"


bool Summer12=true;


void drawOneVariable( DrawBase* db, TTree* tree, const std::string& varName, const std::string& axisName, int nbins, float xmin, float xmax, std::string treeVar="" );

void drawSinglePtBin( DrawBase* db, QGLikelihoodCalculator* qglc, TTree* tree, float ptMin, float ptMax );
void drawPlot( DrawBase* db, TH1D* h1_gluon, TH1D* h1_quark, std::string name, float ptMin, float ptMax, const std::string& labelText="" );
void drawRoC( DrawBase* db, float ptMin, float ptMax, const std::string& flag, TH1D* h1_new_gluon, TH1D* h1_new_quark, TH1D* h1_old_gluon, TH1D* h1_old_quark, TH1D* h1_MLP_gluon=0, TH1D* h1_MLP_quark=0, const std::string& labelText="" );


int main() {


  TChain* tree = new TChain("reducedTree");
  tree->Add("/cmsrm/pc25_2/pandolf/MC/Summer12/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_testQG_finali_withCHS/QG_2nd*.root/reducedTree");

  DrawBase* db = new DrawBase("prova");

  db->set_outputdir("NormalizedPlotsMC");
  

  //drawOneVariable( db, tree, "ptD_QCJet", "p_{T}D", 50, 0., 1.0001);
  //drawOneVariable( db, tree, "rmsCandJet", "RMSCand", 50, 0., 0.3);
  //drawOneVariable( db, tree, "axis1_QCJet", "Axis_{1}", 50, 0., 0.3);
  //drawOneVariable( db, tree, "axis2_QCJet", "Axis_{2}", 50, 0., 0.3);
  //drawOneVariable( db, tree, "pullJet",     "Pull", 50, 0., 0.03);
  //drawOneVariable( db, tree, "RJet",        "R", 50, 0., 1.0001);
  //drawOneVariable( db, tree, "nPFCandJet", "Total Multiplicity", 50, 0., 100., "nChargedJet[0]+nNeutralJet[0]");
  //drawOneVariable( db, tree, "nChargedJet", "Charged Multiplicity", 50, 0., 100.);
  //drawOneVariable( db, tree, "nNeutralJet", "Neutral Multiplicity", 50, 0., 100.);

  QGLikelihoodCalculator* qglc;
  //if( Summer12 ) qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/CMSSW_5_3_6/src/QG/QGLikelihood/test/Histos_2012.root");
  if( Summer12 ) qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/public/Histos_2012_NEW.root");
  else           qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/CMSSW_5_3_6/src/QG/QGLikelihood/test/Histos.root");


  //drawSinglePtBin( db, qglc, tree, 20., 25. );
  //drawSinglePtBin( db, qglc, tree, 25., 30. );
  drawSinglePtBin( db, qglc, tree, 30., 40. );
  //drawSinglePtBin( db, qglc, tree, 40., 50. );
  drawSinglePtBin( db, qglc, tree, 50., 65. );
  //drawSinglePtBin( db, qglc, tree, 65., 80. );
  drawSinglePtBin( db, qglc, tree, 80., 100. );
  //drawSinglePtBin( db, qglc, tree, 150., 200. );
  drawSinglePtBin( db, qglc, tree, 200., 250. );
  //drawSinglePtBin( db, qglc, tree, 300., 400. );
  drawSinglePtBin( db, qglc, tree, 500., 600. );
  //drawSinglePtBin( db, qglc, tree, 800., 1000. );

  return 0;

}




void drawOneVariable( DrawBase* db, TTree* tree, const std::string& varName, const std::string& axisName, int nbins, float xmin, float xmax, std::string treeVar ) {

  if( treeVar=="" ) treeVar = varName + "[0]";

  std::string treeSelection = "abs(etaJet[0])<2. && ptJet[0]>80. && ptJet[0]<120.";
  std::string quarkSelection = treeSelection + "&& abs(pdgIdJet[0])<4";
  std::string gluonSelection = treeSelection + "&& pdgIdJet[0]==21";

  std::string legendTitle = "80 < p_{T} < 120 GeV";

  TH1D* h1_quark = new TH1D("quark", "", nbins, xmin, xmax );
  TH1D* h1_gluon = new TH1D("gluon", "", nbins, xmin, xmax );

  tree->Project( "quark", treeVar.c_str(), quarkSelection.c_str() );
  tree->Project( "gluon", treeVar.c_str(), gluonSelection.c_str() );

  h1_quark->SetLineWidth(2);
  h1_quark->SetLineColor(38);
  h1_quark->SetFillColor(38);
  h1_quark->SetFillStyle(3004);

  h1_gluon->SetLineWidth(2);
  h1_gluon->SetLineColor(46);
  h1_gluon->SetFillColor(46);
  h1_gluon->SetFillStyle(3005);

  float yMax_quark = h1_quark->GetMaximum()/h1_quark->Integral();
  float yMax_gluon = h1_gluon->GetMaximum()/h1_gluon->Integral();

  float yMax = (yMax_quark>yMax_gluon) ? yMax_quark : yMax_gluon;
  yMax *= 1.3;

  TPaveText* label_eta = new TPaveText(0.2, 0.85, 0.3, 0.9, "brNDC" );
  label_eta->SetTextSize(0.038);
  label_eta->SetFillColor(0);
  label_eta->AddText("|#eta| < 2");

  TPaveText* label_top = db->get_labelTop();

  TLegend* legend = new TLegend(0.58, 0.65, 0.9, 0.9, legendTitle.c_str());
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->AddEntry( h1_quark, "Quark Jets", "F" );
  legend->AddEntry( h1_gluon, "Gluon Jets", "F" );

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  
  TH2D* h2_axes = new TH2D("axes", "", 10, xmin, xmax, 10, 0., yMax );
  h2_axes->SetXTitle(axisName.c_str());
  h2_axes->SetYTitle("Normalized to Unity");

  h2_axes->Draw();

  h1_gluon->DrawNormalized("same");
  h1_quark->DrawNormalized("same");

  label_top->Draw("same");
  label_eta->Draw("same");
  legend->Draw("same");

  gPad->RedrawAxis();

  std::string canvasName = db->get_outputdir() + "/" + varName.c_str() + ".eps";
  c1->SaveAs(canvasName.c_str());
  std::string command = "epstopdf " + canvasName;
  system( command.c_str() );

  delete c1;
  delete h2_axes;
  delete h1_quark;
  delete h1_gluon;
  delete legend;
  delete label_eta;

}



void drawSinglePtBin( DrawBase* db, QGLikelihoodCalculator* qglc, TTree* tree, float ptMin, float ptMax ) {

  std::cout << "-> Processing pt bin: " << ptMin << "-" << ptMax << " GeV..." << std::endl;

  bool doFwd = (ptMin<100.);

  int njet;
  tree->SetBranchAddress("nJet", &njet);
  float pt[20];
  tree->SetBranchAddress("ptJet", pt);
  float eta[20];
  tree->SetBranchAddress("etaJet", eta);
  int pdgId[20];
  tree->SetBranchAddress("pdgIdJet", pdgId);
  float rho;
  tree->SetBranchAddress("rhoPF", &rho);
  int nCharged[20];
  tree->SetBranchAddress("nChargedJet", nCharged);
  int nNeutral[20];
  tree->SetBranchAddress("nNeutralJet", nNeutral);
  float ptD[20];
  tree->SetBranchAddress("ptDJet", ptD);
  float ptD_QC[20];
  tree->SetBranchAddress("ptD_QCJet", ptD_QC);
  float axis2_QC[20];
  tree->SetBranchAddress("axis2_QCJet", axis2_QC);
  int nCharged_QC[20];
  tree->SetBranchAddress("nChg_QCJet", nCharged_QC);
  int nNeutral_ptCut[20];
  tree->SetBranchAddress("nNeutral_ptCutJet", nNeutral_ptCut);
  float qglMLPJet[20];
  tree->SetBranchAddress("qgMLPJet", qglMLPJet);
  float qglJet[20];
  tree->SetBranchAddress("qglJet", qglJet);


  TH1D* h1_qgl_old_gluon = new TH1D("qgl_old_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_old_quark = new TH1D("qgl_old_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgl_new_gluon = new TH1D("qgl_new_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_new_quark = new TH1D("qgl_new_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgl_newHisto_gluon = new TH1D("qgl_newHisto_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_newHisto_quark = new TH1D("qgl_newHisto_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgMLP_gluon = new TH1D("qgMLP_gluon", "", 100, 0., 1.);
  TH1D* h1_qgMLP_quark = new TH1D("qgMLP_quark", "", 100, 0., 1.);


  // in the transition:
  TH1D* h1_qgl_old_T_gluon = new TH1D("qgl_old_T_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_old_T_quark = new TH1D("qgl_old_T_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgl_new_T_gluon = new TH1D("qgl_new_T_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_new_T_quark = new TH1D("qgl_new_T_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgl_newHisto_T_gluon = new TH1D("qgl_newHisto_T_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_newHisto_T_quark = new TH1D("qgl_newHisto_T_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgMLP_T_gluon = new TH1D("qgMLP_T_gluon", "", 100, 0., 1.);
  TH1D* h1_qgMLP_T_quark = new TH1D("qgMLP_T_quark", "", 100, 0., 1.);

  // in the forward:
  TH1D* h1_qgl_new_F_gluon = new TH1D("qgl_new_F_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_new_F_quark = new TH1D("qgl_new_F_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgl_newHisto_F_gluon = new TH1D("qgl_newHisto_F_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_newHisto_F_quark = new TH1D("qgl_newHisto_F_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgMLP_F_gluon = new TH1D("qgMLP_F_gluon", "", 100, 0., 1.);
  TH1D* h1_qgMLP_F_quark = new TH1D("qgMLP_F_quark", "", 100, 0., 1.);


  int nentries = tree->GetEntries();

  for( unsigned int ientry=0; ientry<nentries; ++ientry ) {

    tree->GetEntry(ientry);

    if( njet==0 ) continue;

    if( pt[0]<ptMin || pt[0]>ptMax ) continue;
    if( ptD_QC[0]>0.9 ) continue; //this is to cut out anomalous (~single particle) jets

    float qgl_new = qglJet[0];

    if( fabs(eta[0])<2. && h1_qgl_old_gluon->GetEntries()<10000 && h1_qgl_old_quark->GetEntries()<10000 ) { //save time

      float qgl_old = qglc->computeQGLikelihoodPU( pt[0], rho, nCharged[0], nNeutral[0], ptD[0]);
      float qgl_newHisto = qglc->computeQGLikelihood2012( pt[0], eta[0], rho, nCharged_QC[0]+nNeutral_ptCut[0], ptD_QC[0], axis2_QC[0]);

      if( fabs(pdgId[0])<5 ) {
        h1_qgl_old_quark->Fill( qgl_old );
        h1_qgl_new_quark->Fill( qgl_new );
        h1_qgl_newHisto_quark->Fill( qgl_newHisto );
        h1_qgMLP_quark->Fill( qglMLPJet[0] );
      }
      if( pdgId[0]==21 ) {
        h1_qgl_old_gluon->Fill( qgl_old );
        h1_qgl_new_gluon->Fill( qgl_new );
        h1_qgl_newHisto_gluon->Fill( qgl_newHisto );
        h1_qgMLP_gluon->Fill( qglMLPJet[0] );
      }


    } else if( fabs(eta[0])<2.5 ) {

      float qgl_old = qglc->computeQGLikelihoodPU( pt[0], rho, nCharged[0], nNeutral[0], ptD[0]);
      float qgl_newHisto = qglc->computeQGLikelihood2012( pt[0], eta[0], rho, nCharged_QC[0]+nNeutral_ptCut[0], ptD_QC[0], axis2_QC[0]);

      if( fabs(pdgId[0])<5 ) {
        h1_qgl_old_T_quark->Fill( qgl_old );
        h1_qgl_new_T_quark->Fill( qgl_new );
        h1_qgl_newHisto_T_quark->Fill( qgl_newHisto );
        h1_qgMLP_T_quark->Fill( qglMLPJet[0] );
      }
      if( pdgId[0]==21 ) {
        h1_qgl_old_T_gluon->Fill( qgl_old );
        h1_qgl_new_T_gluon->Fill( qgl_new );
        h1_qgl_newHisto_T_gluon->Fill( qgl_newHisto );
        h1_qgMLP_T_gluon->Fill( qglMLPJet[0] );
      }

    } else if( fabs(eta[0])>3. ) {

      float qgl_newHisto = qglc->computeQGLikelihood2012( pt[0], eta[0], rho, nCharged_QC[0]+nNeutral_ptCut[0], ptD_QC[0], axis2_QC[0]);

      if( fabs(pdgId[0])<5 ) {
        h1_qgl_new_F_quark->Fill( qgl_new );
        h1_qgl_newHisto_F_quark->Fill( qgl_newHisto );
        h1_qgMLP_F_quark->Fill( qglMLPJet[0] );
      }
      if( pdgId[0]==21 ) {
        h1_qgl_new_F_gluon->Fill( qgl_new );
        h1_qgl_newHisto_F_gluon->Fill( qgl_newHisto );
        h1_qgMLP_F_gluon->Fill( qglMLPJet[0] );
      }

    }
    
    if( h1_qgl_old_gluon->GetEntries()>10000 
     && h1_qgl_old_quark->GetEntries()>10000 
     && h1_qgl_old_T_gluon->GetEntries()>10000 
     && h1_qgl_old_T_quark->GetEntries()>10000 
     && ( !doFwd || (h1_qgl_new_F_quark->GetEntries()>10000
     && h1_qgl_new_F_gluon->GetEntries()>10000) ) ) break;


  }


  drawPlot( db, h1_qgl_old_gluon, h1_qgl_old_quark, "old", ptMin, ptMax, "|#eta| < 2" );
  drawPlot( db, h1_qgl_new_gluon, h1_qgl_new_quark, "new", ptMin, ptMax, "|#eta| < 2" );
  drawPlot( db, h1_qgl_newHisto_gluon, h1_qgl_newHisto_quark, "newHisto", ptMin, ptMax, "|#eta| < 2" );
  drawPlot( db, h1_qgMLP_gluon, h1_qgMLP_quark, "MLP", ptMin, ptMax, "|#eta| < 2" );

  drawPlot( db, h1_qgl_old_T_gluon, h1_qgl_old_T_quark, "old_T", ptMin, ptMax, "2 < |#eta| < 2.5" );
  drawPlot( db, h1_qgl_new_T_gluon, h1_qgl_new_T_quark, "new_T", ptMin, ptMax, "2 < |#eta| < 2.5" );
  drawPlot( db, h1_qgl_newHisto_T_gluon, h1_qgl_newHisto_T_quark, "newHisto_T", ptMin, ptMax, "2 < |#eta| < 2.5" );
  drawPlot( db, h1_qgMLP_T_gluon, h1_qgMLP_T_quark, "MLP_T", ptMin, ptMax, "2 < |#eta| < 2.5" );

  drawPlot( db, h1_qgl_new_F_gluon, h1_qgl_new_F_quark, "new_F", ptMin, ptMax, "3 < |#eta| < 5" );
  drawPlot( db, h1_qgl_newHisto_F_gluon, h1_qgl_newHisto_F_quark, "newHisto_F", ptMin, ptMax, "3 < |#eta| < 5" );
  drawPlot( db, h1_qgMLP_F_gluon, h1_qgMLP_F_quark, "MLP_F", ptMin, ptMax, "3 < |#eta| < 5" );

  drawRoC(db, ptMin, ptMax, "", h1_qgl_new_gluon, h1_qgl_new_quark, h1_qgl_old_gluon, h1_qgl_old_quark, 0, 0, "|#eta| < 2");
  drawRoC(db, ptMin, ptMax, "_withMLP", h1_qgl_new_gluon, h1_qgl_new_quark, h1_qgl_old_gluon, h1_qgl_old_quark, h1_qgMLP_gluon, h1_qgMLP_quark, "|#eta| < 2");
  drawRoC(db, ptMin, ptMax, "_vsHisto", h1_qgl_new_gluon, h1_qgl_new_quark, h1_qgl_newHisto_gluon, h1_qgl_newHisto_quark, 0, 0, "|#eta| < 2");
  drawRoC(db, ptMin, ptMax, "_T", h1_qgl_new_T_gluon, h1_qgl_new_T_quark, h1_qgl_old_T_gluon, h1_qgl_old_T_quark, h1_qgMLP_T_gluon, h1_qgMLP_T_quark, "2 < |#eta| < 2.5");
  drawRoC(db, ptMin, ptMax, "_T_vsHisto", h1_qgl_new_T_gluon, h1_qgl_new_T_quark, h1_qgl_newHisto_T_gluon, h1_qgl_newHisto_T_quark, 0, 0, "2 < |#eta| < 2.5");
  drawRoC(db, ptMin, ptMax, "_F", h1_qgl_new_F_gluon, h1_qgl_new_F_quark, 0, 0, h1_qgMLP_F_gluon, h1_qgMLP_F_quark, "3 < |#eta| < 5");
  drawRoC(db, ptMin, ptMax, "_F_vsHisto", h1_qgl_new_F_gluon, h1_qgl_new_F_quark, h1_qgl_newHisto_F_gluon, h1_qgl_newHisto_F_quark, 0, 0, "3 < |#eta| < 5");

  delete h1_qgl_old_gluon;
  delete h1_qgl_old_quark;

  delete h1_qgl_new_gluon;
  delete h1_qgl_new_quark;

  delete h1_qgl_newHisto_gluon;
  delete h1_qgl_newHisto_quark;

  delete h1_qgMLP_gluon;
  delete h1_qgMLP_quark;

  delete h1_qgl_old_T_gluon;
  delete h1_qgl_old_T_quark;

  delete h1_qgl_new_T_gluon;
  delete h1_qgl_new_T_quark;

  delete h1_qgl_newHisto_T_gluon;
  delete h1_qgl_newHisto_T_quark;

  delete h1_qgMLP_T_gluon;
  delete h1_qgMLP_T_quark;

  delete h1_qgl_new_F_gluon;
  delete h1_qgl_new_F_quark;

  delete h1_qgl_newHisto_F_gluon;
  delete h1_qgl_newHisto_F_quark;

  delete h1_qgMLP_F_gluon;
  delete h1_qgMLP_F_quark;


}





void drawPlot( DrawBase* db, TH1D* h1_gluon, TH1D* h1_quark, std::string name, float ptMin, float ptMax, const std::string& labelText ) {


  h1_quark->Rebin(2);
  h1_gluon->Rebin(2);

  float norm_max_g = h1_gluon->GetMaximum()/h1_gluon->Integral();
  float norm_max_q = h1_quark->GetMaximum()/h1_quark->Integral();
  float hmax = (norm_max_q>norm_max_g) ? norm_max_q : norm_max_g;

  float ymax = hmax*1.2;

  bool isMLP = (name=="MLP" || name=="MLP_F");

  TH2D* h2_axes = new TH2D("axes", "", 10, h1_gluon->GetXaxis()->GetXmin(), h1_gluon->GetXaxis()->GetXmax(), 10, 0., ymax);
  if( isMLP )
    h2_axes->SetXTitle("Quark-Gluon MLP Discriminator");
  else
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



  float xMin_legend = (isMLP) ? 0.2 : 0.55;
  float xMax_legend = (isMLP) ? 0.5 : 0.8;

  char legendTitle[300];
  sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax );
  TLegend* legend = new TLegend( xMin_legend, 0.7, xMax_legend, 0.9, legendTitle );
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry( h1_quark, "Quark Jets", "F");
  legend->AddEntry( h1_gluon, "Gluon Jets", "F");
  legend->Draw("same");
  
  h1_quark->DrawNormalized("same");
  h1_gluon->DrawNormalized("same");



  TPaveText* labelTop = db->get_labelTop();
  labelTop->Draw("same");

  float xMin_label = isMLP ? 0.7 : 0.2;
  float xMax_label = isMLP ? 0.9 : 0.4;

  TPaveText* label = new TPaveText( xMin_label, 0.83, xMax_label, 0.9, "brNDC" );
  label->SetTextSize(0.04);
  label->SetFillColor(0);
  label->AddText(labelText.c_str());
  if( labelText!="" )
    label->Draw("same");


  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "%s/qgl_%s_pt%.0f_%.0f.eps", db->get_outputdir().c_str(), name.c_str(), ptMin, ptMax);
  c1->SaveAs(canvasName);
  sprintf( canvasName, "%s/qgl_%s_pt%.0f_%.0f.png", db->get_outputdir().c_str(), name.c_str(), ptMin, ptMax);
  c1->SaveAs(canvasName);

  delete c1;
  delete h2_axes;
  delete legend;

}


void drawRoC( DrawBase* db, float ptMin, float ptMax, const std::string& flag, TH1D* h1_new_gluon, TH1D* h1_new_quark, TH1D* h1_old_gluon, TH1D* h1_old_quark, TH1D* h1_MLP_gluon, TH1D* h1_MLP_quark, const std::string& labelText ) {


  TGraph* gr_RoC_old = new TGraph(0);
  TGraph* gr_RoC_new = new TGraph(0);
  TGraph* gr_RoC_MLP = new TGraph(0);

  int nbins = h1_new_quark->GetNbinsX();

  for( unsigned int ibin=1; ibin<nbins+1; ++ibin ) {

    float eff_q_old = -1.;
    float eff_g_old = -1.;
  
    if( h1_old_quark!=0 && h1_old_gluon!=0 ) {
      eff_q_old = h1_old_quark->Integral( nbins-ibin, nbins )/h1_old_quark->Integral( 1, nbins );
      eff_g_old = h1_old_gluon->Integral( nbins-ibin, nbins )/h1_old_gluon->Integral( 1, nbins );
    }
  
    float eff_q_MLP = -1.;
    float eff_g_MLP = -1.;
  
    if( h1_MLP_quark!=0 && h1_MLP_gluon!=0 ) { //opposite convention:
      eff_q_MLP = h1_MLP_quark->Integral( 1, ibin )/h1_MLP_quark->Integral( 1, nbins );
      eff_g_MLP = h1_MLP_gluon->Integral( 1, ibin )/h1_MLP_gluon->Integral( 1, nbins );
    }
  
    float eff_q_new = h1_new_quark->Integral( nbins-ibin, nbins )/h1_new_quark->Integral( 1, nbins );
    float eff_g_new = h1_new_gluon->Integral( nbins-ibin, nbins )/h1_new_gluon->Integral( 1, nbins );
  
    gr_RoC_new->SetPoint( ibin-1, 1.-eff_g_new, eff_q_new );

    if( h1_old_quark!=0 && h1_old_gluon!=0 ) 
      gr_RoC_old->SetPoint( ibin-1, 1.-eff_g_old, eff_q_old );

    if( h1_MLP_quark!=0 && h1_MLP_gluon!=0 ) 
      gr_RoC_MLP->SetPoint( ibin-1, 1.-eff_g_MLP, eff_q_MLP );

  }


  gr_RoC_new->SetMarkerSize(1.3);
  gr_RoC_new->SetMarkerStyle(24);
  gr_RoC_new->SetMarkerColor(kRed+3);

  if( h1_old_quark!=0 && h1_old_gluon!=0 ) {
    gr_RoC_old->SetMarkerSize(1.3);
    gr_RoC_old->SetMarkerStyle(20);
    gr_RoC_old->SetMarkerColor(kOrange+1);
  }

  if( h1_MLP_quark!=0 && h1_MLP_gluon!=0 ) {
    gr_RoC_MLP->SetMarkerSize(1.3);
    gr_RoC_MLP->SetMarkerStyle(21);
    gr_RoC_MLP->SetMarkerColor(29);
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
  if( h1_MLP_quark!=0 && h1_MLP_gluon!=0 )
    legend->AddEntry( gr_RoC_MLP, "MLP", "P");
  legend->Draw("same");

  TPaveText* labelTop = db->get_labelTop();
  labelTop->Draw("same");

  TPaveText* label = new TPaveText( 0.7, 0.83, 0.9, 0.9, "brNDC" );
  label->SetTextSize(0.04);
  label->SetFillColor(0);
  label->AddText(labelText.c_str());
  if( labelText!="" )
    label->Draw("same");

  
  if( h1_MLP_quark!=0 && h1_MLP_gluon!=0 ) 
    gr_RoC_MLP->Draw("p same");
  if( h1_old_quark!=0 && h1_old_gluon!=0 ) 
    gr_RoC_old->Draw("p same");
  gr_RoC_new->Draw("p same");

  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "%s/RoC_pt%.0f_%.0f%s.eps", db->get_outputdir().c_str(), ptMin, ptMax, flag.c_str());
  c1->SaveAs(canvasName);
  sprintf( canvasName, "%s/RoC_pt%.0f_%.0f%s.png", db->get_outputdir().c_str(), ptMin, ptMax, flag.c_str());
  c1->SaveAs(canvasName);

  delete c1;
  delete h2_axes;
  delete legend;
}

