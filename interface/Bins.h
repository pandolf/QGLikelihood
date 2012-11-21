#ifndef BINS_H
#define BINS_H

#include "TMath.h"


class Bins{
public:
const static int nRhoBins=25; //18
const static int nPtBins=20; //20
constexpr static double Pt0=20;  //20
constexpr static double Pt1=1000; //1000
constexpr static double Rho0=0;  //2
constexpr static double Rho1=25; //20
constexpr static double PtLastExtend=1000; //3500

static int getBins(double  *bins,int nBins,double MinBin=15.0,double MaxBin=1000.,bool log=false);
static int getBin(int nBins,double  *bins,double value,double*x0=0,double*x1=0);
static void getBins_int( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog=true); 

};



#endif
