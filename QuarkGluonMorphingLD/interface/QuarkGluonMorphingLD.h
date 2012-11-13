// ------------------------------------------------------------
//  
//    QuarkGluonMorphingLD - Class
//    for the computation of the QG likelihood.
//    Needs files provided by having run the
//    double fit procedure. output.txt?
//
// ------------------------------------------------------------

#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include <map>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>
#include "TF1.h"

#ifndef QGM_LD_H
#define QGM_LD_H

#define MAX_STR_LENGTH 1023

class QuarkGluonMorphingLD {

 public:

  QuarkGluonMorphingLD( const char * configName="data/config.ini");
  ~QuarkGluonMorphingLD();

  float computeQuarkGluonMorphingLD( float pt, float rhoPF, float *vars );



  // functions needed to read config files:
  char GetChar(FILE *fr);
  int ReadBinTxt(const char*fileName,int parameter,TH2F*quark,TH2F *gluon);
  int ReadParTxt( const char *fileName,std::map< std::pair<int,int>, double*> *parq,std::map< std::pair<int,int>, double*> *parg ,int NPar=2 );
  char* ReadParameterFromFile(const char*fileName,const char * parName);

  //inteded for debug purpose  - onyl new
  int ComputePars(float pt , float rho,const char varName[], char type, double*par);


 private:
  

  FILE *fr_;
  TH2F* ptD0_q,*ptD0_g;
  TH2F* ptD1_q,*ptD1_g;
  TH2F* ptD2_q,*ptD2_g;
  std::map< std::pair<int,int>,double* > *parqNC,*parqNN,*pargNC,*pargNN,*parqPTD,*pargPTD;


  int nVars;
  std::vector<std::string> varName;
  std::vector<std::string> varFunc;
  std::map< std::pair< std::string, char>, std::map< std::pair<int,int> ,double* >* >  AllPar;
  bool isOldStyle;
  double gammadistr_(double* x, double* par);
  double gammadistr2_(double* x, double* par);
  double functionPtD_(double * x ,double*par);

};
#endif
