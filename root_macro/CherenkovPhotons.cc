#include "TMath.h"

using namespace TMath;

const Double_t RINDEX = 1.05;
const Double_t WIDTH  = 60e-3; // 60 mm
const Double_t ANGLE  = 40.*DegToRad(); // radian
const Double_t L_MIN  = 200e-9;
const Double_t L_MAX  = 650e-9;

Double_t GetBeta(Double_t p, Double_t m) // unit is MeV/c
{
  return p/Sqrt(p*p+m*m);
}

Double_t GetParticleMass(const char* particle)
{
	TString str(particle);
	if(str=="pi+")
		return 139.57018;
	else if(str=="e-")
		return 0.5109989461;
	else if(str=="p+")
		return 938.2720813;
	else if(str=="K+")
		return 493.677;
	else
		return -1.;
}


Int_t CherenkovPhotons()
{
	Double_t p, m;
	std::cout << "Enter the momentum in MeV/c : ";
	std::cin >> p;
	char buffer[64];
	std::cout << "Enter the particle(pi+, e-, p+, K+) : ";
  std::cin >> buffer;
	while((m = GetParticleMass(buffer)) < 0){
		std::cout << "Please enter the proper particle name in list." << std::endl
			<< "Candidate : pi+, e-, p+, K+ : ";
		std::cin >> buffer;
	}
  Double_t beta = GetBeta(p, m);
  Double_t dz   = WIDTH/Cos(ANGLE);

  Int_t nPhotons = (Int_t)2*Pi()*(1. - 1./(RINDEX*RINDEX*beta*beta))*(1./L_MIN - 1./L_MAX)*dz/137;
  std::cout << "The unit of input momentum is MeV/c." << std::endl;
  
  std::cout << "beta \t\t: " << beta << std::endl;
  std::cout << "dz \t\t: " << dz*1e3 << " mm" << std::endl;
  std::cout << "The # of photons : " << nPhotons << std::endl;

  return 1;
}

