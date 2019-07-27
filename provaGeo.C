
{

  gSystem->Load("libGeom");

  TGeoManager *GeoMan = new TGeoManager("GeoMan", "prima prova");
  TGeoMedium* medium=0;
  TGeoVolume* TopBox = GeoMan->MakeBox("TopBox", medium, 100, 100, 100); // top volume
  GeoMan->SetTopVolume(TopBox);                           // set top volume


  TGeoTube* t1 = new TGeoTube("T1", 0, 10, 2);
  TGeoTube* t2 = new TGeoTube("T2", 10, 20, 2);
    
  TGeoTranslation *Tr = new TGeoTranslation("Tr", 0, 0, 10);

  Tr->RegisterYourself();

  TGeoCompositeShape* comp = new TGeoCompositeShape("fullShape", "(T1)+(T2:Tr)");

  TGeoVolume* full = new TGeoVolume("full", comp);
  TopBox->SetLineColor(5);
  TopBox->AddNode(full, 1);
  GeoMan->CloseGeometry();
  
  TopBox->Raytrace();
  TopBox->Draw();

  TAxis3D a;
  a.Draw();

}

