// ------------------------------------------------------------
//  
//    QGLikelihoodCalculator - Class
//    for the computation of the QG likelihood.
//    Needs files provided by having run the
//    Ntp1Finalizer_QG on QCD samples.
//
// ------------------------------------------------------------

#include <string>

#include "TFile.h"
#include "TH1F.h"



class QGLikelihoodCalculator {

 public:

  QGLikelihoodCalculator( const std::string& fileName="QG_QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10.root", int nPtBins=20 );
  virtual ~QGLikelihoodCalculator() {};

  float computeQGLikelihood( float pt, int nCharged, int nNeutral, float ptD, float rmsCand );
  float likelihoodProduct( float nCharged, float nNeutral, float ptD, float rmsCand, TH1F* h1_nCharged, TH1F* h1_nNeutral, TH1F* h1_ptD, TH1F* h1_rmsCand);



 private:

  TFile* histoFile_;

  int nPtBins_;

};

