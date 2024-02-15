void test1()
{
  TFile *f = TFile::Open("pytree.root");
  TH2F* hdphi_deta = new TH2F("hdphi_deta","; #Delta#phi; #Delta#eta", 100, -5, 5, 100, -6.283185307,6.283185307);
  T -> Draw("eta:phi>>hdphi_deta","","goff");
  TFile* fout = new TFile("pythia_histos.root");
  hdphi_deta->Write();
  fout->Close();
}

