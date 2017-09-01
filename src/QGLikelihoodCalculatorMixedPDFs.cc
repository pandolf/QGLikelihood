#include "../interface/QGLikelihoodCalculatorMixedPDFs.h"
#include "../interface/Bins.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TMath.h"

#include <map>
using namespace std;

//#define DEBUG



//void getBins_int( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog=true);
//void getBins( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog=true);




// constructor:

QGLikelihoodCalculatorMixedPDFs::QGLikelihoodCalculatorMixedPDFs( const std::string& fileName ) {

  histoFile_ = TFile::Open(fileName.c_str());

  //nPtBins_ = nPtBins;
  //nRhoBins_ = nRhoBins;

}

// ADD map destructor
QGLikelihoodCalculatorMixedPDFs::~QGLikelihoodCalculatorMixedPDFs()
{
}





float QGLikelihoodCalculatorMixedPDFs::computeQGLikelihood( float pt, float rhoPF, int mult, float ptd, float axis2 ) {


  double ptBins[100];
  Bins::getBins_int( Bins::nPtBins+1, ptBins, Bins::Pt0,Bins::Pt1,true);
  ptBins[Bins::nPtBins+1]=Bins::PtLastExtend;

  double rhoBins[100];
  Bins::getBins_int(Bins::nRhoBins+1,rhoBins,Bins::Rho0,Bins::Rho1,false);

  int i_pt  = Bins::getBin(Bins::nPtBins ,ptBins ,pt   );
  int i_rho = Bins::getBin(Bins::nRhoBins,rhoBins,rhoPF);

  if(i_pt  <0 ) return -1;
  if(i_rho <0 ) return -1;



  char histoName[300];
  sprintf( histoName, "pdf_dy_mult_pt%d_rho%d"  , i_pt, i_rho );
  if(plots_[histoName]==NULL)
    plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  TH1F* h1_mult_dy = plots_[histoName];

  sprintf( histoName, "pdf_dy_ptd_pt%d_rho%d"  , i_pt, i_rho );
  if(plots_[histoName]==NULL)
    plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  TH1F* h1_ptd_dy = plots_[histoName];

  sprintf( histoName, "pdf_dy_axis2_pt%d_rho%d"  , i_pt, i_rho );
  if(plots_[histoName]==NULL)
    plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  TH1F* h1_axis2_dy = plots_[histoName];


  sprintf( histoName, "pdf_qcd_mult_pt%d_rho%d"  , i_pt, i_rho );
  if(plots_[histoName]==NULL)
    plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  TH1F* h1_mult_qcd = plots_[histoName];

  sprintf( histoName, "pdf_qcd_ptd_pt%d_rho%d"  , i_pt, i_rho );
  if(plots_[histoName]==NULL)
    plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  TH1F* h1_ptd_qcd = plots_[histoName];

  sprintf( histoName, "pdf_qcd_axis2_pt%d_rho%d"  , i_pt, i_rho );
  if(plots_[histoName]==NULL)
    plots_[histoName]=(TH1F*)histoFile_->Get(histoName)->Clone();
  TH1F* h1_axis2_qcd = plots_[histoName];


  float qcdP = likelihoodProduct( mult, ptd, axis2, h1_mult_qcd, h1_ptd_qcd, h1_axis2_qcd );
  float dyP  = likelihoodProduct( mult, ptd, axis2, h1_mult_dy , h1_ptd_dy , h1_axis2_dy  );

  float QGLikelihood = dyP / ( qcdP + dyP );

 // if(h1_nCharged_gluon) delete h1_nCharged_gluon;
 // if(h1_nCharged_quark) delete h1_nCharged_quark;
 // if(h1_nNeutral_gluon) delete h1_nNeutral_gluon;
 // if(h1_nNeutral_quark) delete h1_nNeutral_quark;
 // if(h1_ptd_gluon) delete h1_ptd_gluon;
 // if(h1_ptd_quark) delete h1_ptd_quark;
 // if(h1_rmsCand_gluon) delete h1_rmsCand_gluon;
 // if(h1_rmsCand_quark) delete h1_rmsCand_quark;

  return QGLikelihood;

}



