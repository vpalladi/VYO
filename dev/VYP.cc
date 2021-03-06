#include "Riostream.h"
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include <signal.h>
#include <fcntl.h>
#include <tgmath.h>

#include <Rtypes.h>

#include <TCanvas.h>
#include <TApplication.h>
#include <TEllipse.h>
#include <TH2Poly.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TGraph.h>

#include "Pixel.h"

// UNITS
// distances = nm; time=seconds; mass=ug; volume=ul

// polimer density 1400 ug/ml
// drug consentatrion 15+-2 ug/mg => 21+-2.8 ug/ul

// diffusivity 6.7 × 10−12 cm2/s


int main(Int_t argc, char **argv){

  TApplication app("app", &argc, argv);

  Double_t Rmax = 10000;
  Int_t nZ=50;
  Int_t nR=0;
  vector<Double_t> Z;
  vector<Double_t> R;

  for(int i_z=0; i_z<nZ; i_z++)
    Z.push_back( i_z*Rmax/nZ );

  Double_t R0 = 500.;
  R.push_back( 0 );

  while( R.back()<Rmax )
    R.push_back( R0*TMath::Sqrt( ++nR ) );
  

  cout << "Max number of pixels along R: " << R.size() << endl
       << "Max number of pixels along Z: " << Z.size() << endl;


  // '+1' takes into account an extra row and column of with 0 concentration pixels (boundary conditions) 
  Pixel*** Pixels = new Pixel**[R.size()+1];
  Bool_t** PixelMask = new Bool_t*[R.size()+1];
  for(int i_col=0; i_col<(R.size()+1); i_col++){
    Pixels[i_col] = new Pixel*[Z.size()+1];
    PixelMask[i_col] = new Bool_t[Z.size()+1];
  }

  //  TH2Poly* Particle =  new TH2Poly("Particle", "Particle", 0, 11000, 0, 11000);
  TH2Poly* Particle =  new TH2Poly("Particle", "Particle", 0, 500, 0, 11000);

  Particle->SetBit(TH1::kCanRebin);

  TRandom3 rand;

  for(int i_z=0; i_z<(Z.size()+1); i_z++)
    for(int i_r=0; i_r<(R.size()+1); i_r++)
      PixelMask[i_r][i_z] = 0;

  for(int i_z=0; i_z<(Z.size()-1); i_z++){ // pay attention on the '-1' => the real sphere is 1 step smaller
    for(int i_r=0; i_r<(R.size()-1); i_r++){
     
      if( (R.at(i_r)*R.at(i_r) + Z.at(i_z)*Z.at(i_z)) < (Rmax*Rmax) ){
  
  	PixelMask[i_r][i_z] = 1;
  
	Double_t concentration = rand.Gaus(21, 0.0001);
	
	Pixels[i_r][i_z] = new Pixel(R.at(i_r), Z.at(i_z), R.at(i_r+1), Z.at(i_z+1), concentration);
	Pixels[i_r][i_z]->AddBin(Particle);
	Pixels[i_r][i_z]->SetPixelPosition(i_r, i_z);
	
      }      
      else
	Pixels[i_r][i_z] = new Pixel(R.at(i_r), Z.at(i_z), R.at(i_r+1), Z.at(i_z+1), 0);
      
    }
  }
  
  TGraph g;
  
  // time evolution
  for(int i_step=0; i_step<10; i_step++){
    for(int i_z=0; i_z<(Z.size()); i_z++){
      for(int i_r=0; i_r<(R.size()); i_r++){
  
  	if( PixelMask[i_r][i_z] == 1 ){

	  Double_t CRm, CRp, CZm, CZp, C;

	  Int_t i_Rm = i_r-1;
	  Int_t i_Rp = i_r+1;
	  Int_t i_Zm = i_z-1;
	  Int_t i_Zp = i_z+1;

	  C = Pixels[i_r][i_z]->GetC();

	  // Zp
	  if( PixelMask[i_r][i_Zp] == 0 ){
	    CZp = 0;
	  }
	  else{ 
	    CZp = Pixels[i_r][i_Zp]->GetC();
	  }

	  // Zm
	  if( i_Zm<0 ){
	    CZm = C;
	  }
	  else{
	    if( PixelMask[i_r][i_Zm] == 0 )
	      CZm = 0;
	    else 
	      CZm = Pixels[i_r][i_Zm]->GetC();
	  }    
	  
	  // Rp
	  if( PixelMask[i_Rp][i_z] == 0 )
	    CRp = 0;
	  else 
	    CRp = Pixels[i_Rp][i_z]->GetC();
	  
	  // Rm
	  if( i_Rm<0 ){
	    CRm = C;
	  }
	  else{
	    if( PixelMask[i_Rm][i_z] == 0 )
	      CRm = 0;
	    else 
	      CRm = Pixels[i_Rm][i_z]->GetC();
	  }
	  
	  
	  Double_t dR;
	  if(i_Zm<0){
	    dR = 2*Pixels[i_r][i_z]->GetXcenter();
	  }
	  else{
	    dR = Pixels[i_r][i_z]->GetXcenter() - Pixels[i_r][i_Zm]->GetXcenter();
	    cout << ".....dR " << Pixels[i_r][i_z]->GetXcenter() << endl;
	  }

	  Double_t dZ;
	  if(i_Rm<0) {
	    dZ = Pixels[i_r][i_z]->GetXcenter();
	  }
	  else{
	    dZ = Pixels[i_r][i_z]->GetXcenter() - Pixels[i_Rm][i_z]->GetXcenter();
	    cout << Pixels[i_r][i_z]->GetXcenter() << " " << Pixels[i_Rm][i_z]->GetXcenter() << endl;
	  }

	  Pixels[i_r][i_z]->Evolve(5., CRm, CRp, CZm, CZp, dR, dZ, R0);
	  

  	}
  	
  	
      }
    } 
  }

  for(int i_r=0; i_r<(R.size()-1); i_r++){
    //g.SetPoint( i_r+1, Pixels[i_r][0]->GetXcenter(), Pixels[i_r][0]->GetC() );
    if( PixelMask[i_r][0] == 1 )
      g.SetPoint( i_r, Pixels[i_r][0]->GetXcenter(), Pixels[i_r][0]->GetC() );
  }
  
  TCanvas c("c", "c", 1000, 1000);
  
  Particle->Draw("colz");
  
  c.Update();

  TCanvas c2("c2", "c2", 1000, 1000);
  c2.Update();
  //  g.Draw("A*");

  app.Run();

  return 0;

}
