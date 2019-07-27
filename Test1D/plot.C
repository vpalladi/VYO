{

  TFile f("Out.root", "read");

f.Get("Ct_0").Draw("AL");

for(int i=1; i<1e4; i++){

  f.Get(Form("Ct_%d",i)).Draw("sameL");

 }
 
}
