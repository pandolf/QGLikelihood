#include "../interface/QuarkGluonMorphingLD.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include "TMath.h"
#include "TF1.h"
#include <map>
#include <stdio.h>
#include "../interface/PtBins.h"

using namespace std;

//#define DEBUG


// constructor:

QuarkGluonMorphingLD::QuarkGluonMorphingLD( const char * configName)
{
	const char * vars=ReadParameterFromFile(configName,"QGFIT4VARS");
	const char * funcs=ReadParameterFromFile(configName,"QGFIT4FUNCS");
	const char * dir=ReadParameterFromFile(configName,"QGFIT4TXTDIR");
	string dirName(dir);
	
	char str[1023];int n; 
	nVars=0;
	while(sscanf(vars,"%s%n",str,&n)==1)
		{
		vars+=n;
		varName.push_back( string(str) );
		nVars++;
		}
	while(sscanf(funcs,"%s%n",str,&n)==1)
		{
		funcs+=n;
		varFunc.push_back( string(str) );
		}
	if( varFunc.size() != varName.size() ) fprintf(stderr,"ERROR NUMBER OF VARS DIFFERS\n");
	for(int i=0; i<nVars;i++){
	AllPar[ pair<string,char>( varName[i], 'Q') ] = new map< pair<int,int>,double* >;
	AllPar[ pair<string,char>( varName[i], 'G') ] = new map< pair<int,int>,double* >;
	
	//printf("DEBUG Map %s'Q'\n",varName[i].c_str());
	//printf("DEBUG Map %s'G'\n",varName[i].c_str());
	if( varFunc[i] == string("gamma") ){
		//printf("DEBUG gamma\n");
		ReadParTxt( (dirName+varName[i]+string(".txt")).c_str(),AllPar[ pair<string,char>( varName[i], 'Q') ],AllPar[ pair<string,char>( varName[i], 'G') ] );
		}
	else if( varFunc[i] == string("gamma2") ){
		//printf("DEBUG gamma\n");
		ReadParTxt( (dirName+varName[i]+string(".txt")).c_str(),AllPar[ pair<string,char>( varName[i], 'Q') ],AllPar[ pair<string,char>( varName[i], 'G') ] );
		}
	else if (varFunc[i] == string("functionPtD") ){
		//printf("DEBUG functionPtD\n");
		ReadParTxt( ( dirName+varName[i]+string(".txt")).c_str(),AllPar[ pair<string,char>( varName[i], 'Q') ],AllPar[ pair<string,char>( varName[i], 'G') ],3);
		}
	else printf("DEBUG function ERROR ---%s---\n",varFunc[i].c_str());
	
	}	
	isOldStyle=false;
	
}

// ADD map destructor
QuarkGluonMorphingLD::~QuarkGluonMorphingLD()
{
}

double QuarkGluonMorphingLD::gammadistr_(double* x, double* par)
{
        return TMath::Exp( - x[0] *par[0]/par[1] ) * TMath::Power(x[0],par[0]-1) * TMath::Power(par[1]/par[0],-par[0])/TMath::Gamma(par[0]) ;
}
double QuarkGluonMorphingLD::gammadistr2_(double* x, double* par)
{
	double alpha=par[1] * par[1]/ (par[0]*par[0]); //par 0 = sigma;  par 1= mean
	double beta=par[1];
	return TMath::Exp( - x[0] *alpha/beta ) * TMath::Power(x[0],alpha-1) * TMath::Power(beta/alpha,-alpha)/TMath::Gamma(alpha) ;		
}

//half gamma+ offset
double QuarkGluonMorphingLD::functionPtD_(double * x ,double*par)
{
        if((x[0]-par[0])<0)return 0;
        return TMath::Exp( - (x[0]-par[0]) *par[1]/par[2] ) * TMath::Power((x[0]-par[0]),par[1]-1) * TMath::Power(par[2]/par[1],-par[1])/TMath::Gamma(par[1]) ;
}


float QuarkGluonMorphingLD::computeQuarkGluonMorphingLD( float pt, float rhoPF, float*vars ) {
double Q=1;
double G=1;


//TF1 *pol3=new TF1("pol3","[0]+[1]*TMath::Log(x)+[2]*TMath::Log(x)*TMath::Log(x)+[3]*TMath::Log(x)*TMath::Log(x)*TMath::Log(x)",20,3500);//LOG! ->PT
//TF1 *pol1=new TF1("pol1","[0]+[1]*x",0,20); //NOT LOG -> RHO

double *par=new double[5];
double *x=new double[5];
//double a,b;
for(unsigned int i=0;i<varName.size();++i){
	//printf("DEBUG %s\n",varName[i].c_str());
	//printf("DEBUG %s\n",varFunc[i].c_str());
	int R;	
	if( varFunc[i] == string("gamma") ){
	x[0]=vars[i];
	//NC Q
	R=ComputePars(pt,rhoPF,varName[i].c_str(),'Q',par);if(R!=2)fprintf(stderr,"ERROR nPar=%d instead of 2\n",R);
	Q*=gammadistr_(x,par);
	//NC G
	R=ComputePars(pt,rhoPF,varName[i].c_str(),'G',par);if(R!=2)fprintf(stderr,"ERROR nPar=%d instead of 2\n",R);
	G*=gammadistr_(x,par);
	}
	else if( varFunc[i] == string("gamma2") ){
		x[0]=vars[i];
		R=ComputePars(pt,rhoPF,varName[i].c_str(),'Q',par);if(R!=2)fprintf(stderr,"ERROR nPar=%d instead of 2\n",R);
		Q*=gammadistr2_(x,par);
		R=ComputePars(pt,rhoPF,varName[i].c_str(),'G',par);if(R!=2)fprintf(stderr,"ERROR nPar=%d instead of 2\n",R);
		G*=gammadistr2_(x,par);
	}
	else if( varFunc[i] == string("functionPtD") ){ //select the right function
	x[0]=vars[i]; 
	//PtD Q
	R=ComputePars(pt,rhoPF,varName[i].c_str(),'Q',par);if(R!=3)fprintf(stderr,"ERROR nPar=%d instead of 3\n",R);
	Q*=functionPtD_(x,par);
	//PtD G
	R=ComputePars(pt,rhoPF,varName[i].c_str(),'G',par);if(R!=3)fprintf(stderr,"ERROR nPar=%d instead of 3\n",R);
	G*=functionPtD_(x,par);
	}
 }

delete[] par;
delete[] x;
//delete pol3;
//delete pol1;
if(Q==0)return 0;
return float(Q/(Q+G));
}

int QuarkGluonMorphingLD::ComputePars(float pt , float rhoPF,const char varName[], char type, double*par)
{
int R=0;
TF1 *pol3=new TF1("pol3","[0]+[1]*TMath::Log(x)+[2]*TMath::Log(x)*TMath::Log(x)+[3]*TMath::Log(x)*TMath::Log(x)*TMath::Log(x)",20,3500);//LOG! ->PT
TF1 *pol1=new TF1("pol1","[0]+[1]*x",0,20); //NOT LOG -> RHO
float a,b;
//	printf("PAR1\n");
	if((*AllPar[pair<string,char>(varName,type)])[pair<int,int>(0,0)] !=NULL){
	pol3->SetParameters( (*AllPar[pair<string,char>(varName,type)])[pair<int,int>(0,0)]);//par0 a
		b=pol3->Eval(pt);
	pol3->SetParameters( (*AllPar[pair<string,char>(varName,type)])[pair<int,int>(0,1)]);//par0 a
		a=pol3->Eval(pt);
	//printf("a=%f b=%f\n",a,b);
	pol1->SetParameter(0,b);pol1->SetParameter(1,a);
		par[0]=pol1->Eval(rhoPF);
	R++;
	}

	//printf("PAR2\n");
	if((*AllPar[pair<string,char>(varName,type)])[pair<int,int>(1,0)] !=NULL){
	pol3->SetParameters( (*AllPar[pair<string,char>(varName,type)])[pair<int,int>(1,0)]);//par0 a
		b=pol3->Eval(pt);
	pol3->SetParameters( (*AllPar[pair<string,char>(varName,type)])[pair<int,int>(1,1)]);//par0 a
		a=pol3->Eval(pt);
	//printf("a=%f b=%f\n",a,b);
	pol1->SetParameter(0,b);pol1->SetParameter(1,a);
		par[1]=pol1->Eval(rhoPF);
	R++;
	}
	
	//to this part only if is for PtD style ... ?
	//printf("PAR3\n");
	if((*AllPar[pair<string,char>(varName,type)])[pair<int,int>(2,0)] !=NULL){
	pol3->SetParameters( (*AllPar[pair<string,char>(varName,type)])[pair<int,int>(2,0)]);//par0 a
		b=pol3->Eval(pt);
	pol3->SetParameters( (*AllPar[pair<string,char>(varName,type)])[pair<int,int>(2,1)]);//par0 a
		a=pol3->Eval(pt);
	pol1->SetParameter(0,b);pol1->SetParameter(1,a);
		par[2]=pol1->Eval(rhoPF);
	R++;
	}
	//for(int i=0;i<R;i++)printf("par[%d]=%f ",i,par[i]);
	//printf("\n");
	delete pol3;
	delete pol1;
	return R;
}



char QuarkGluonMorphingLD::GetChar(FILE *fr)
{
fpos_t pos;//position in the file
char c;
//get position of the line to be analized;
fgetpos(fr,&pos);
c=fgetc(fr);
fsetpos(fr,&pos); //moving back to the beginning of the line
return c;
}

int QuarkGluonMorphingLD::ReadBinTxt(const char*fileName,int parameter,TH2F*quark,TH2F *gluon)
	{
	FILE *fr=fopen(fileName,"r");
	if(fr==NULL) return -1;
	//getting binnig
	double RhoBins[25];int nRhoBins=20;
	double PtBins[25];int nPtBins=18;
	getBins_int(18,PtBins,20,1000,true);
	PtBins[18]=3500;
	getBins_int(21,RhoBins,0,20,false);

	double PtBinsMean[25];for(int i=0;i<nPtBins;++i){PtBinsMean[i]=(PtBins[i]+PtBins[i+1])/2.;}  
	double RhoBinsMean[25];for(int i=0;i<nRhoBins;++i){RhoBinsMean[i]=(RhoBins[i]+RhoBins[i+1])/2.;}  
		//---
		char c;	
		//fscanf(fr,"[quark]\n");//move into the file
		while (( GetChar(fr) == '[') || (GetChar(fr)== '{'))
			{
			c='\0';while(c!='\n'){fscanf(fr,"%c",&c);fprintf(stderr,"%c",c);}
			} //[quark] {bla}
		
		while (( GetChar(fr) != '[') && (GetChar(fr)!= EOF) )
		{
			//read a formatted line
			float ptmin, ptmax, rhomin, rhomax,par[10];int nPar;
			fscanf(fr,"%f %f %f %f %d %*f %*f",&ptmin,&ptmax,&rhomin,&rhomax,&nPar);nPar-=2;
			fprintf(stderr,"%f %f %f %f %d %d\n",ptmin,ptmax,rhomin,rhomax,nPar,parameter);	
		
			for(int i=0;i<nPar;++i)fscanf(fr,"%f",&par[i]);
			c='\0';while(c!='\n'){fscanf(fr,"%c",&c);} //go to the end of line
			//check
			if(nPar<=parameter){perror("I do not have that parameter\n");break;}
			//filling the histogram
			quark->SetBinContent(quark->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parameter] );
		}//end of while: loop on the lines
		
		//skip lines that begin with [ or  -> useless
		while (( GetChar(fr) == '[') || (GetChar(fr)== '{'))
			{
			c='\0';while(c!='\n'){fscanf(fr,"%c",&c);fprintf(stderr,"%c",c);}
			} //[gluon] {bla}
		
		while (( GetChar(fr) != '[') && (GetChar(fr)!= EOF))
		{
			//read a formatted line
			float ptmin, ptmax, rhomin, rhomax,par[10];int nPar;
			fscanf(fr,"%f %f %f %f %d %*f %*f",&ptmin,&ptmax,&rhomin,&rhomax,&nPar);nPar-=2;
			fprintf(stderr,"%f %f %f %f %d %d\n",ptmin,ptmax,rhomin,rhomax,nPar,parameter);	
		
			for(int i=0;i<nPar;++i)fscanf(fr,"%f",&par[i]);
			c='\0';while(c!='\n'){fscanf(fr,"%c",&c);} //go to the end of line
			//check
			if(nPar<=parameter){perror("I do not have that parameter\n");break;}
			//filling the histogram
			gluon->SetBinContent(gluon->FindBin((ptmin+ptmax)/2.,(rhomin+rhomax)/2.), par[parameter] );
		}//end of while: loop on the lines
		fclose(fr);
	return 0;
}

int QuarkGluonMorphingLD::ReadParTxt( const char *fileName,std::map< pair<int,int>, double*> *parq,std::map< pair<int,int>, double*> *parg ,int NPar )
	{
	FILE *fr=fopen(fileName,"r");
	if(fr==NULL) return -1;
	for(int P=0;P<NPar;P++)//NOT TRUE FOR PTD
		{
		//(*parq)[pair<int,int>(P,0)]=new double[8];
		//(*parg)[pair<int,int>(P,1)]=new double[8];
		double *x=new double[5];
		fscanf(fr,"%*d aq %lf %lf %lf %lf\n",&x[0],&x[1],&x[2],&x[3]);
			(*parq)[pair<int,int>(P,1)]=x;
		x=new double[5];
		fscanf(fr,"%*d ag %lf %lf %lf %lf\n",&x[0],&x[1],&x[2],&x[3]);
			(*parg)[pair<int,int>(P,1)]=x;
		x=new double[5];
		fscanf(fr,"%*d bq %lf %lf %lf %lf\n",&x[0],&x[1],&x[2],&x[3]);
			(*parq)[pair<int,int>(P,0)]=x;
		x=new double[5];
		fscanf(fr,"%*d bg %lf %lf %lf %lf\n",&x[0],&x[1],&x[2],&x[3]);
			(*parg)[pair<int,int>(P,0)]=x;

		}
	return 0;
	}



char* QuarkGluonMorphingLD::ReadParameterFromFile(const char*fileName,const char * parName)
{
FILE *fr=fopen(fileName,"r");
if(fr==NULL) return NULL;
char *R=new char[MAX_STR_LENGTH]; //must be a pointer - should survive outside
char P[MAX_STR_LENGTH];
char S[MAX_STR_LENGTH];

//leggi una linea su S
int STATUS=1;
while(STATUS!=EOF)
{
	char c='\0';
	int i=0;
	while ( (STATUS=fscanf(fr,"%c",&c))!=EOF ){
		if(c=='\0') break;
		if(c=='\n') break;
		S[i]=c;i++;
		}
	S[i]='\0';
	i=0;
	while(S[i]!='=' && S[i]!='\0'){P[i]=S[i];i++;}
	P[i]='\0';int j=0;
	if(S[i]=='=')i++;
	while(S[i]!='\0'){R[j]=S[i];i++;j++;}
	R[j]='\0';
	bool isPar=true;
	for(i=0;;i++)
		{
		if(P[i]=='\0' && parName[i]=='\0')break;
		if(P[i]=='\0' || parName[i]=='\0'){isPar=false; break;} // but not together
		if(parName[i]!=P[i]){isPar=false;break;}
		}
	
	if(isPar){fclose(fr);return R;}
}
fclose(fr);
return NULL;
}
