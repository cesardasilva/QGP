#include "TProfile.h"
#include "TRandom.h"

TProfile *hp, *hpy;

void make_prof()
{
  hp  = new TProfile("hp", "x-profile of shifted Gaussian",160,-4.,4.);
  hpy = new TProfile("hpy","y-profile of shifted Gaussian",160,-4.,4.);
  for( int i=0; i<100000; ++i ) {
    double x = gRandom->Gaus(0,1);
    double y = gRandom->Gaus(1,2);
    hp->Fill(x,y);
    hpy->Fill(y,x);
  }
}
