void test1()
{
  TFile *f = TFile::Open("pytree.root");
  TH3F *T; f->GetObject("T",T);

  T -> Draw("eta:phi>>h(100,-5,5,-3.1415,3.1415)","","colz");
  
}

