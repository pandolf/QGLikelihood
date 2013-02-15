CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/CondFormats/JetMETObjects/interface

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs` 

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -L$(CMSSW_BASE)/lib/slc5_amd64_gcc462/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit



checkQGL: DrawBase.o fitTools.o Bins.o QGLikelihoodCalculator.o checkQGL.cpp
	$(CC) -Wall $(INCLUDES) -o checkQGL checkQGL.cpp DrawBase.o fitTools.o Bins.o QGLikelihoodCalculator.o $(ROOTFLAG) $(EXTRALIBS)

drawQGLMCPlots: DrawBase.o fitTools.o Bins.o QGLikelihoodCalculator.o drawQGLMCPlots.cpp
	$(CC) -Wall $(INCLUDES) -o drawQGLMCPlots drawQGLMCPlots.cpp DrawBase.o fitTools.o Bins.o QGLikelihoodCalculator.o $(ROOTFLAG) $(EXTRALIBS)



TFitParticleEtEtaPhi.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitParticleEtEtaPhi.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitParticleEtEtaPhi.cc $(ROOTFLAG)

TAbsFitParticle.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TAbsFitParticle.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TAbsFitParticle.cc $(ROOTFLAG)

TFitConstraintM.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintM.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintM.cc $(ROOTFLAG)

TFitConstraintMGaus.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintMGaus.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintMGaus.cc $(ROOTFLAG)

TFitConstraintEp.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintEp.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintEp.cc $(ROOTFLAG)

TAbsFitConstraint.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TAbsFitConstraint.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TAbsFitConstraint.cc $(ROOTFLAG)

TKinFitter.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TKinFitter.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TKinFitter.cc $(ROOTFLAG)



DrawBase.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C fitTools.o
	$(CC) $(CFLAGS) $(INCLUDES) fitTools.o $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C $(ROOTFLAG) $(EXTRALIBS) 

fitTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C $(ROOTFLAG) $(EXTRALIBS)

StatTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/StatTools.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/StatTools.cc $(ROOTFLAG) $(EXTRALIBS)

QGLikelihoodCalculator.o: $(CMSSW_BASE)/src/QG/QGLikelihood/src/QGLikelihoodCalculator.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/QG/QGLikelihood/interface $(CMSSW_BASE)/src/QG/QGLikelihood/src/QGLikelihoodCalculator.cc $(ROOTFLAG)

Bins.o: $(CMSSW_BASE)/src/QG/QGLikelihood/src/Bins.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/QG/QGLikelihood/interface $(CMSSW_BASE)/src/QG/QGLikelihood/src/Bins.cc $(ROOTFLAG)





BTagSFUtil.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/BTagSFUtil.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/BTagSFUtil.cc $(ROOTFLAG)

SFlightFuncs.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/SFlightFuncs.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/SFlightFuncs.cc $(ROOTFLAG)

MistagFuncs.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/MistagFuncs.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/MistagFuncs.cc $(ROOTFLAG)

PUWeight.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/PUWeight.C
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/PUWeight.C $(ROOTFLAG)


AnalysisJet.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisJet.cc 
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisJet.cc $(ROOTFLAG)

AnalysisElectron.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisElectron.cc AnalysisLepton.o
	$(CC) $(CFLAGS) $(INCLUDES) -IUserCode/emanuele/CommonTools/include/ $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisElectron.cc $(ROOTFLAG)

AnalysisMuon.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisMuon.cc AnalysisLepton.o
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisMuon.cc $(ROOTFLAG)

AnalysisLepton.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisLepton.cc
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisLepton.cc $(ROOTFLAG)

AnalysisNeutrino.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisNeutrino.h
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisNeutrino.h $(ROOTFLAG)



clean:
	rm *.o
