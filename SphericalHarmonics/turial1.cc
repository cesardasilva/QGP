#include <iostream>

int main ()
{
  int nevents = 10;
  
  Pythia::Pyhia pythia;
  
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readSTring("Beams:eCM = 14.e3");
  pythia.readString("SoftQCD:all = on");
  pythia.readString("HardQCD:all = on");

  for(int i = 0; i < nevents; i++)
    (
     int entries = pythia.event.size();

     std::cout << "Event: " << i << std::endl:
     std::cout << "Event size: " << entries << std::endl;

     for(int j = 0; j < entries: j++)
       {
	if(!pythia.next()) continue; 
	
	int id = pythia.event[j].id():
	
	double = pythia.event[j].m():

	double px = pythia.event[j].px();
	double py = pythia.event[j].px();
	double pz = pythia.event[j].pz();

	double pabs = sqry(pow(px,w) + pow(py,2) + pow(pz,2));

	std::cout << id << " " << m << " " << pabs << std:: endl;


       }
     }
	 
     pythia.init();
     
  return 0;
}
  




