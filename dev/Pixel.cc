

#include "Pixel.h"

Pixel::Pixel(){

  isPixelEroded = 0; // at t0 pixel is not eroded //
  
}

Pixel::Pixel(Double_t Xmin, Double_t Ymin, Double_t Xmax, Double_t Ymax, Double_t C0){
  
  Pixel::SetLimits( Xmin, Ymin, Xmax, Ymax );
  
  Pixel::SetT0C( C0 );

  Pixel::SetDiffusivity( 670 ); // nm2/s

  isPixelEroded = 0; // at t0 pixel is not eroded //
  
}


void Pixel::SetLimits(Double_t Xmin, Double_t Ymin, Double_t Xmax, Double_t Ymax){
  
//  fXmin = Xmin;
//  fYmin = Ymin;
//  fXmax = Xmax;
//  fYmax = Ymax;

  // coordinates are choosen in order to have unifor dx and dy

  fXmin = Xmin*Xmin/(500*500);
  fYmin = Ymin;
  fXmax = Xmax*Xmax/(500*500);
  fYmax = Ymax;
  

  fXcenter = (fXmax+fXmin)/2;
  fYcenter = (fYmax+fYmin)/2;

  //  cout << "fXmin     "<< fXmin      << endl;;
  //  cout << "fYmin     "<< fYmin      << endl;;
  //  cout << "fXmax     "<< fXmax      << endl;;
  //  cout << "fYmax     "<< fYmax      << endl;;
  //  cout << "fXcenter "<< fXcenter  << endl;
  //  cout << "fYcenter "<< fYcenter  << endl;
  
  fPixel = new TBox(fXmin, fYmin, fXmax, fYmax);

}


Pixel::~Pixel() { 

  delete fPixel; 

}


void Pixel::AddBin(TH2Poly* hpoly){

  fPixelID = hpoly->AddBin( fXmin, fYmin, fXmax, fYmax );
  
  hpoly->SetBinContent(fPixelID, fC0);

}


void Pixel::Evolve(Double_t dt, Double_t CXm, Double_t CXp, Double_t CYm, Double_t CYp, Double_t dX, Double_t dY, Double_t X0){

  dX = fXmax - fXmin;
  
  Double_t dEta = dX*dX/(X0*X0);
  Double_t Eta = fXcenter*fXcenter/(X0*X0);

  Double_t A =  4. * Eta/(X0*X0*dEta*dEta);
  Double_t B = 1/(dY*dY);

  fC = fD * dt * (- 2 * fC * (A+B) + A * ( CXp + CXm ) + B * ( CYp + CYm ) );

  cout << "A    " << A    << endl;
  cout << "B    " << B    << endl;
  cout << "dY   " << dY   << endl;
  cout << "dEta " << dEta << endl;
  cout << "dX   " << dX   << endl;
  cout << "X0   " << X0   << endl;
  cout << "eta  " << Eta  << endl;
  cout << "- 2 * fC * (A+B) " << - 2 * fC * (A+B) << endl;
  cout << ">>>>>>>>>>>>> " << fC << endl;

  //  Double_t dX = fXcenter - Pixels[fXid][Ym]->GetXcenter();
  //  Double_t dY = fYcenter - Pixels[Xm][fYid]->GetYcenter();

  //  cout << "---------------------" << endl;
  //  cout << "fC " << fC << endl;
  //  cout << "(dY * dY) " << (dY * dY)  << endl;
  //  cout << "(dX * dX) " << (dX * dX) << endl;
  //  cout << "dt " << dt << endl;
  //  cout << "1/fXcenter " << 1/fXcenter << endl;
  //  cout << "(CXp - fC) " << (CXp - fC) << endl;
  //  cout << "( CXp - 2 * fC + CXm ) / (dX * dX) " << "(" << CXp << "-" << 2 << "*" << fC << "+" << CXm << ")" <<  "/" << "(" << dX << "*" << dX << ")" << " = " 
  //  << ( CXp - 2 * fC + CXm ) / (dX * dX) << endl;
  //  cout << "(1/fXcenter * ( ( CXp - fC ) / dX )) " << (1/fXcenter * ( ( CXp - fC ) / dX )) << endl;
  //  cout << "( CYp - 2 * fC + CYm ) / (dY * dY) " << ( CYp - 2 * fC + CYm ) / (dY * dY) << endl;

  // fC = fC + fD * dt * ( (1/fXcenter * ( ( CXp - fC ) / dX ))  + ( CXp - 2 * fC + CXm ) / (dX * dX) + ( CYp - 2 * fC + CYm ) / (dY * dY) );
  //  fC = fD * dt * ( fC * (8*TMath::Sqrt(fXcenter*fXcenter*fXcenter/(R0*dX*dX)) +) + () + () )
 
  // original // fC = fC + fD * dt * ( 1/fXcenter * ( ( CXp - fC ) / dX )  + ( CXp - 2 * fC + CXm ) / (dX * dX) + ( CYp - 2 * fC + CYm ) / (dY * dY) );

  //  fC = fC + fD * dt;
  // pay attention on the boundary condition for the inner part of the sphere 
  //cout <<  fD << " "  << dt << endl;
  //cout << "dopo " << fC << endl;

}
