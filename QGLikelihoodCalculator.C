#include "QGLikelihoodCalculator.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "fitTools.h"





// constructor:

QGLikelihoodCalculator::QGLikelihoodCalculator( const std::string& fileName, int nPtBins ) {

  histoFile_ = TFile::Open(fileName.c_str());

  nPtBins_ = nPtBins;

}





float QGLikelihoodCalculator::computeQGLikelihood( float pt, int nCharged, int nNeutral, float ptD, float rmsCand ) {


  float ptMin = 0.;
  float ptMax = 0.;

  const int nBinsPlusOne(nPtBins_+1);
  Double_t ptBins[nBinsPlusOne];
  fitTools::getBins_int( nBinsPlusOne, ptBins, 15., 1000. );


  if( pt>ptBins[nPtBins_] ) {
    ptMin = ptBins[nPtBins_-1];
    ptMax = ptBins[nPtBins_];
  } else {
    for( unsigned int iBin=0; iBin<nPtBins_; ++iBin ) {
      if( pt>ptBins[iBin] && pt<ptBins[iBin+1] ) {
        ptMin = ptBins[iBin];
        ptMax = ptBins[iBin+1];
      } //if
    } //for
  } //else
  
  if( ptMax==0. ) return -1.;


  char histoName[200];
  sprintf( histoName, "nCharged_gluon_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_nCharged_gluon = (TH1F*)histoFile_->Get(histoName);
  sprintf( histoName, "nCharged_quark_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_nCharged_quark = (TH1F*)histoFile_->Get(histoName);

  sprintf( histoName, "nNeutral_gluon_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_nNeutral_gluon = (TH1F*)histoFile_->Get(histoName);
  sprintf( histoName, "nNeutral_quark_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_nNeutral_quark = (TH1F*)histoFile_->Get(histoName);

  sprintf( histoName, "ptD_gluon_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_ptD_gluon = (ptD>=0.) ? (TH1F*)histoFile_->Get(histoName) : 0;
  sprintf( histoName, "ptD_quark_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_ptD_quark = (ptD>=0.) ? (TH1F*)histoFile_->Get(histoName) : 0;

  sprintf( histoName, "rmsCand_gluon_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_rmsCand_gluon = (rmsCand>=0.) ? (TH1F*)histoFile_->Get(histoName) : 0;
  sprintf( histoName, "rmsCand_quark_pt%.0f_%.0f", ptMin, ptMax);
  TH1F* h1_rmsCand_quark = (rmsCand>=0.) ? (TH1F*)histoFile_->Get(histoName) : 0;


  float gluonP = likelihoodProduct( nCharged, nNeutral, ptD, rmsCand, h1_nCharged_gluon, h1_nNeutral_gluon, h1_ptD_gluon, h1_rmsCand_gluon );
  float quarkP = likelihoodProduct( nCharged, nNeutral, ptD, rmsCand, h1_nCharged_quark, h1_nNeutral_quark, h1_ptD_quark, h1_rmsCand_quark );

  float QGLikelihood = gluonP / (gluonP + quarkP );

  return QGLikelihood;

}


float QGLikelihoodCalculator::likelihoodProduct( float nCharged, float nNeutral, float ptD, float rmsCand, TH1F* h1_nCharged, TH1F* h1_nNeutral, TH1F* h1_ptD, TH1F* h1_rmsCand) {

  h1_nCharged->Scale(1./h1_nCharged->Integral("width"));
  h1_nNeutral->Scale(1./h1_nNeutral->Integral("width"));
  if( h1_ptD!=0 )
    h1_ptD->Scale(1./h1_ptD->Integral("width"));
  if( h1_rmsCand!=0 )
    h1_rmsCand->Scale(1./h1_rmsCand->Integral("width"));

  float likeliProd =  h1_nCharged->GetBinContent(h1_nCharged->FindBin(nCharged))*
                      h1_nNeutral->GetBinContent(h1_nNeutral->FindBin(nNeutral));

  if( h1_ptD!=0 )
    likeliProd*=h1_ptD->GetBinContent(h1_ptD->FindBin(ptD));

  if( h1_rmsCand!=0 )
    likeliProd*=h1_rmsCand->GetBinContent(h1_rmsCand->FindBin(rmsCand));

  return likeliProd;

}
