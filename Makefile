CC = gcc
CXX = g++
CCFLAGS = -std=c++11 -g -W -Wall -Wno-unused-parameter -Wno-ignored-qualifiers #-pedantic -fPIC

##### ROOT
# REQUIRE ROOT
# (Assuming a defined ROOTSYS guarantees all ROOT-related stuff infra works...) 
ifeq ($(ROOTSYS),)
  $(error Environment ROOTSYS not set)
endif
ROOTFLAGS = `root-config --cflags --glibs`
ROOTVERSION = -D ROOT5   # Probably useless

##### PHAST
ifeq ($(PHAST),)
  $(error Environment PHAST not set)
endif
PHAST_LIBS =
PHAST_INCL =
LHAPDF_INCL := -I$(shell lhapdf-config --incdir)
LHAPDF_LIBS := $(shell lhapdf-config --libs)

ifeq ($(SITE),BW)
CCFLAGS = -DBW
endif

ifeq ($(DEBUG),1)
CCFLAGS += -DDEBUG
endif

all : analySIDIS acceptance comparison extractor plotter dvm
analySIDIS : analySIDIS_split analySIDIS_collect
acceptance : accsplit accfuse acccollect
comparison : compRDRD compMCRD compMCMC compMult compAcc
extractor : FFExtractor FFPlotter
plotter : plotMult
dvm : DVM


%.o: %.cc %.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -c -o $@ $<

analySIDIS_split: analySIDIS_split.cc analySIDIS_split.h
	@echo 'Building SIDIS analysis package..'
	@$(CXX) $(CCFLAGS) -Wno-ignored-qualifiers $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< $(PHAST_LIBS) $(PHAST_INCL)

analySIDIS_collect: analySIDIS_collect.cc analySIDIS_collect.h
	@$(CXX) $(CCFLAGS) -Wno-ignored-qualifiers $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<
	@echo 'SIDIS analysis package built !'

accsplit: acceptance_split.cc acceptance_split.h
	@echo 'Building acceptance package..'
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

accfuse: acceptance_fuse.cc acceptance_fuse.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

acccollect: acceptance_collect.cc acceptance_collect.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<
	@echo 'Acceptance package built !'

compRDRD: compRDRD.cc compRDRD.h
	@echo 'Building RD/RD.RD/MC.MC/MC.Mult/Mult package..'
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

compMCRD: compMCRD.cc compMCRD.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

compMCMC: compMCMC.cc compMCMC.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

compAcc: compAcc.cc compAcc.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<

compMult: compMult.cc compMult.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<
	@echo 'RD/RD.RD/MC.MC/MC.Mult/Mult package built !'

FFExtractor: FFExtractor.cc FFExtractor.h
	@echo 'Building FF extraction package..'
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< $(LHAPDF_LIBS) $(LHAPDF_INCL)

FFPlotter: FFPlotter.cc FFPlotter.h
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< $(LHAPDF_LIBS) $(LHAPDF_INCL)
	@echo 'FF extraction package built !'

plotMult: plotMult.cc plotMult.h
	@echo 'Building plotting device package..'
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $< $(LHAPDF_LIBS) $(LHAPDF_INCL)
	@echo 'Plotting device package built !'

DVM: DVM.cc DVM.h
	@echo 'Building DVM package..'
	@$(CXX) $(CCFLAGS) $(ROOTFLAGS) $(ROOTVERSION) -o $@ $<
	@echo 'DVM package built !'

clean :
	@rm -rf *.o accsplit accfuse acccollect analySIDIS_split analySIDIS_collect compRDRD compMCRD compMCMC compAcc compMult FFExtractor FFPlotter plotMult DVM
