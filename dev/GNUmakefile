# Makefile

NAME := VYP

STATIC := -static

EXTRALIBS += -LEventDisplay/lib -lEventDisplay$(STATIC) -lEve -lGeom -lRGL -lFTGL -lGui -lGed -lTreePlayer -lEG

#include ./config/global_lkr.gmk

CFLAGS = $(shell root-config --cflags)
LIBS   = $(shell root-config --libs)
GLIBS  = $(shell root-config --glibs)


all: VYP

VYP: VYP.cc pixel.o
	@echo "Compiling $< ..."
	g++ -o VYP VYP.cc Pixel.o $(CFLAGS) $(LIBS) $(GLIBS)

pixel.o: Pixel.cc
	g++ -c -o Pixel.o Pixel.cc $(CFLAGS) $(LIBS) $(GLIBS)

try:
	@echo CFLAGS
	@echo $(CFLAGS)
	@echo LIBS
	@echo $(LIBS)
	@echo GLIBS
	@echo $(GLIBS)
	@echo EXTRALIBS
	@echo $(EXTRALIBS)
	@echo LD_LIBRARY_PATH
	@echo $(LD_LIBRARY_PATH)

