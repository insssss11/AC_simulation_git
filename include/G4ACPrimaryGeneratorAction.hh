//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4ACPrimaryGeneratorAction.hh 90623 2015-06-05 09:24:30Z gcosmo $
//
/// \file G4ACPrimaryGeneratorAction.hh
/// \brief Definition of the G4ACPrimaryGeneratorAction class

#ifndef G4ACPrimaryGeneratorAction_h
#define G4ACPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ACPrimaryGeneratorActionMessenger.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class G4ACPrimaryGeneratorActionMessenger;
class G4ACEventGenerator;
/// The primary generator action class with particle gun.
///
enum GUNMODE {
    SINGLE_DETECTOR = 0,
    LEPS2
};

class G4ACPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    G4ACPrimaryGeneratorAction();    
    virtual ~G4ACPrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);         
  
    // method to access particle gun
    G4ParticleGun* GetParticleGun1() const { return fParticleGun1;}
    G4ParticleGun* GetParticleGun2() const { return fParticleGun2;}
    G4double GetIncidentAngle() const {return fIncidentAngle;}
    GUNMODE GetGunMode() const {return fGunMode;}

    void SetIncidentAngle(G4double incidentAngle);
    void SetIncidentPosition(G4double incidentPosition);
    void SetMomentumAmp(G4double momentumAmp);

    void SetParticle(const char* particleName);

    void SetGunMode(const char* gunModeName);
    void PrintStat();

    G4ACEventGenerator *GetEventGenerator()
    {
      return fEventGenerator;
    }

  private:
    // pointer a to G4 gun class
    G4ParticleGun  *fParticleGun1, // for SINGLE_DETECTOR mode
    *fParticleGun2; // for LEPS2 mode
    
    G4double fIncidentAngle;
    G4double fIncidentPosition;
    G4double fMomentumAmp;
    G4ParticleDefinition *fParticle;

    G4ACEventGenerator *fEventGenerator;

    enum GUNMODE fGunMode;

    G4ACPrimaryGeneratorActionMessenger* fMessenger;
    void PreparePrimaryVertex();
};

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
// gamma + proton -> ??
// custom made event generator for LEPS2/SPring8 experiment
// for now, the amplitude is not considered

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
enum DECAY_CHANNEL {
  PROBABILITY_SUM_ERROR = -1, // used when probabilities failed to sum up to an identity
  PION_CH1, // gamma + proton -> pi+ pi- proton
  PION_CH2, // gamma + proton -> pi+ delta -> pi+ pi- proton
  PION_CH3, // gamma + proton -> pi- delta++ -> pi+ pi- proton
  PION_CH4, // gamma + proton -> rho proton -> pi+ pi- proton
  KAON_CH1  // gamma + proton -> Kaon+ lambda
};

class G4ACEventGenerator
{
  public:
    G4ACEventGenerator(G4double gammaE = 2.0*GeV);
    ~G4ACEventGenerator();

  public:
    void SetGammaE(G4double gammaE)
    {
      fGammaE = gammaE/GeV;
      (*fPinitial)[3] = fGammaE;
      (*fPinitial)[4] = fGammaE + PROTON_MASS;
    }
    
    DECAY_CHANNEL GetDecayChannel(void){return fDecayChannel;}

    void GenerateEvent(G4ParticleGun *particleGun, G4Event *anEvent);
  private:
    // in GeV/c^2
    const Double_t PION_CHARGED_MASS = 0.13957018,
      PROTON_MASS = 0.938272981,
      DELTA_MASS = 1.232,
      RHO_MASS = 0.77526,
      KAON_CHARGED_MASS = 0.493677,
      LAMBDA_MASS = 1.115683;
    
    Double_t fGammaE;    
    // Root class to calculate kinematics
    TGenPhaseSpace *fGenPhaseSpace1, *fGenPhaseSpace2;

    TLorentzVector *fPinitial, *fPpip, *fPpim, *fPproton, *fPkap, *fPlambda;
  private:
    void GenerateForPionCh1(G4ParticleGun *particleGun, G4Event *anEvent); // pi+ pi- proton decay channel
    void GenerateForPionCh2(G4ParticleGun *particleGun, G4Event *anEvent); // pi+ delta0 decay channel
    void GenerateForPionCh3(G4ParticleGun *particleGun, G4Event *anEvent); // pi- delta++ decay channel 
    void GenerateForPionCh4(G4ParticleGun *particleGun, G4Event *anEvent); // rho proton decay channel
    void GenerateForKaonCh1(G4ParticleGun *particleGun, G4Event *anEvent); // kaon+ lambda decay channel

    G4bool IsTouchDetector(const TLorentzVector *momentum, G4double dz);

    DECAY_CHANNEL fDecayChannel;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
