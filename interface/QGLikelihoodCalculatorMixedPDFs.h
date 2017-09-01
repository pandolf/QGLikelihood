// ------------------------------------------------------------
//  
//    QGLikelihoodCalculatorMixedPDFs - Class
//    for the computation of the QG likelihood.
//    Needs files provided by having run the
//    Ntp1Finalizer_QG on QCD samples.
//
// ------------------------------------------------------------

#ifndef QGLikelihoodCalculatorMixedPDFs_h
#define QGLikelihoodCalculatorMixedPDFs_h

#include <string>

#include "QGLikelihood/interface/Bins.h"

#include "TFile.h"
#include "TH1F.h"
#include <map>



class QGLikelihoodCalculatorMixedPDFs {

 public:

  QGLikelihoodCalculatorMixedPDFs( const std::string& fileName="/afs/cern.ch/work/p/pandolf/public/ReducedHisto_2012.root");//, unsigned nPtBins=Bins::nBinsPt, unsigned int nRhoBins=Bins::nBinsRho );
   ~QGLikelihoodCalculatorMixedPDFs();

  float computeQGLikelihood( float pt, float rho, int mult, float ptD, float axis2 );

  float likelihoodProduct( int mult, float ptD, float axis2, TH1F* h1_mult, TH1F* h1_ptD, TH1F* h1_axis2 );



 private:

  TFile* histoFile_;
  std::map<std::string,TH1F*> plots_;
  //unsigned int nPtBins_;
  //unsigned int nRhoBins_;

};


#endif
