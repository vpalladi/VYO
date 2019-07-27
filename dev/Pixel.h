#ifndef Pixel_H
#define PIXEL_H

#include <iostream>

#include <TBox.h>
#include <TH2Poly.h>

using namespace std;

class Pixel{

public:

  Pixel();
  Pixel(Double_t Xmin, Double_t Ymin, Double_t Xmax, Double_t Ymax, Double_t C0);
  ~Pixel();

  void SetT0C( Double_t C0 )          { fC0 = C0; fC = C0; }
  void SetDiffusivity( Double_t D )   { fD = D; } // nm2/s
  void SetLimits(Double_t Xmin, Double_t Ymin, Double_t Xmax, Double_t Ymax);

  void AddBin(TH2Poly* hpoly);

  void SetPixelPosition(Int_t X, Int_t Y) { fXid = X; fYid = Y; }

  Double_t GetT0C()          { return fC0; }
  Double_t GetC()            { return fC;  }
  Double_t GetDiffusivity()  { return fD; }

  Int_t    GetPixelID()      { return fPixelID; }
  Double_t GetXmin()         { return fXmin; }
  Double_t GetYmin()         { return fYmin; }
  Double_t GetXmax()         { return fXmax; }
  Double_t GetYmax()         { return fYmax; }

  Double_t GetXcenter()      { return fXcenter; }
  Double_t GetYcenter()      { return fYcenter; }


  void Evolve(Double_t dt, Double_t CXm, Double_t CXp, Double_t CYm, Double_t CYp, Double_t dX, Double_t dY, Double_t R0);

private:

 Bool_t isPixelEroded;
 Double_t fXmin, fXmax, fYmin, fYmax;
 Double_t fXcenter, fYcenter;

 Int_t fXid, fYid;
 
 TBox* fPixel;
 Double_t fC0, fD;
 Double_t fC;
 
 Int_t fPixelID;
 
};


#endif
