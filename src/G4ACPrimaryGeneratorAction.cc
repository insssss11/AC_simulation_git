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
// $Id: G4ACPrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACPrimaryGeneratorAction.cc
/// \brief Implementation of the G4ACPrimaryGeneratorAction class

#include "G4ACPrimaryGeneratorAction.hh"
#include "G4ACPrimaryGeneratorActionMessenger.hh"
#include "G4ACBaseDetectorConstruction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "G4UImanager.hh"

#include <math.h>
#include "Randomize.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

extern const G4double TARGET_LENGTH, TARGET_CENTER_POSITION,
  AC_DETECTOR_R, TPC_EDGE_R, TPC_DOWNSTREAM_Z;

G4ACPrimaryGeneratorAction::G4ACPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun1(0), fParticleGun2(0), fIncidentAngle(0), fIncidentPosition(0), fMomentumAmp(0), fParticle(0), fEventGenerator(0), fGunMode(SINGLE_DETECTOR), fMessenger(0)
{
  fMessenger = new G4ACPrimaryGeneratorActionMessenger(this);
  // default configuration
  G4int n_particle = 1;
  fParticleGun1  = new G4ParticleGun(n_particle);
  fEventGenerator = new G4ACEventGenerator(2.*GeV);
  SetIncidentAngle(40.*deg);
  SetIncidentPosition(100.*mm);
  SetMomentumAmp(1000.*MeV);
  SetParticle("pi+");

  // this is for leps2 mode
  fParticleGun2  = new G4ParticleGun();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACPrimaryGeneratorAction::~G4ACPrimaryGeneratorAction()
{
  delete fParticleGun1;
  delete fParticleGun2;
  delete fMessenger;
  delete fEventGenerator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(fGunMode == SINGLE_DETECTOR)
    fParticleGun1->GeneratePrimaryVertex(anEvent);
  else if(fGunMode == LEPS2)
  {
    fEventGenerator->GenerateEvent(fParticleGun2, anEvent);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACPrimaryGeneratorAction::SetIncidentAngle(G4double incidentAngle)
{
  fIncidentAngle = incidentAngle;
  fParticleGun1->SetParticleMomentumDirection(G4ThreeVector(0., 1.*cos(fIncidentAngle), 1.*sin(fIncidentAngle)));
  SetIncidentPosition(fIncidentPosition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACPrimaryGeneratorAction::SetIncidentPosition(G4double incidentPosition)
{
  fIncidentPosition = incidentPosition;
  fParticleGun1->SetParticlePosition(G4ThreeVector(0.*mm, -50.*mm, (-256.3 + fIncidentPosition - 50*tan(fIncidentAngle))*mm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACPrimaryGeneratorAction::SetMomentumAmp(G4double momentumAmp)
{
  fMomentumAmp = momentumAmp;
  fParticleGun1->SetParticleMomentum(momentumAmp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACPrimaryGeneratorAction::SetParticle(const char *particleName)
{
  fParticle = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  if(fParticle != nullptr)
    { fParticleGun1->SetParticleDefinition(fParticle); }
  else{
    G4cerr << "Cannot find particle " << particleName << "!" << G4endl
    << "Set particle as Geantino." << G4endl;
    fParticleGun1->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("geantino"));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACPrimaryGeneratorAction::SetGunMode(const char *gunModeName)
{
  G4String gunMode(gunModeName);
  if(gunMode == "single")
  {
    fGunMode = SINGLE_DETECTOR;  
  }
  else if(gunMode == "leps2")
  {
    fGunMode = LEPS2;
  }
  else
    G4cerr << "Invalid gun mode name." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACPrimaryGeneratorAction::PrintStat()
{
  G4cout << "Particle : " << fParticle->GetParticleName() << G4endl
  << "Incident Position : " << fIncidentPosition/mm << "mm" << G4endl
  << "Incident Angle : " << fIncidentAngle/deg << "deg" << G4endl
  << "Momentum Amp : " << fMomentumAmp/MeV << "MeV/c" << G4endl
  << "Gunmode : " << fGunMode << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACEventGenerator::G4ACEventGenerator(G4double gammaE):
  fGammaE((Double_t)gammaE/GeV), fGenPhaseSpace1(0), fGenPhaseSpace2(0),
  fPinitial(0), fPpip(0), fPpim(0), fPproton(0), fPkap(0), fPlambda(0)
{
  fGenPhaseSpace1 = new TGenPhaseSpace();
  fGenPhaseSpace2 = new TGenPhaseSpace();
  fPinitial = new TLorentzVector(0., 0., fGammaE, fGammaE + PROTON_MASS);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACEventGenerator::~G4ACEventGenerator(){
  delete fGenPhaseSpace1;
  delete fGenPhaseSpace2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACEventGenerator::GenerateEvent(G4ParticleGun *particleGun, G4Event *anEvent)
{
  G4double probability[] = {0.2, 0.2, 0.2, 0.2, 0.2};
  G4double sumNum = 0.0;
  G4double randNum = G4UniformRand();
  if(randNum < (sumNum +=  probability[0]))
    GenerateForPionCh1(particleGun, anEvent);
  else if(randNum < (sumNum +=  probability[1]))
    GenerateForPionCh2(particleGun, anEvent);
  else if(randNum < (sumNum +=  probability[2]))
    GenerateForPionCh3(particleGun, anEvent);
  else if(randNum < (sumNum +=  probability[3]))
    GenerateForPionCh4(particleGun, anEvent);
  else if(randNum < (sumNum +=  probability[4]))
    GenerateForKaonCh1(particleGun, anEvent);
  else
  {
    G4cerr
    << "----------------------------------------------------------------------------------------------" << G4endl
    << "Warning! The sum of decay probabilities in the event generator failed to be an identity." << G4endl
    << "Generating the event listed at the lastest order in the code." << G4endl
    << "----------------------------------------------------------------------------------------------" << G4endl;
    GenerateForKaonCh1(particleGun, anEvent);
    fDecayChannel = PROBABILITY_SUM_ERROR;
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACEventGenerator::GenerateForPionCh1(G4ParticleGun *particleGun, G4Event *anEvent)
{
  Double_t masses[3] = {PION_CHARGED_MASS, PION_CHARGED_MASS, PROTON_MASS};
  G4ThreeVector position;
  fGenPhaseSpace1->SetDecay(*fPinitial, 3, masses);
  while(true)
  {
    fGenPhaseSpace1->Generate();

    fPpip = fGenPhaseSpace1->GetDecay(0);
    fPpim = fGenPhaseSpace1->GetDecay(1);
    fPproton = fGenPhaseSpace1->GetDecay(2);
  
    position = G4ThreeVector(0.*mm, 0.*mm, TARGET_CENTER_POSITION + (G4UniformRand() - 0.5)*TARGET_LENGTH);
    if(IsTouchDetector(fPpip, TPC_DOWNSTREAM_Z - position[2]) || IsTouchDetector(fPpim, TPC_DOWNSTREAM_Z - position[2]))
      break;
  }

  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("pi+"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPpip->Px(), fPpip->Py(), fPpip->Pz())/fPpip->P());
  particleGun->SetParticleEnergy(fPpip->E()*GeV);
  particleGun->SetParticlePosition(position);
  /*
  G4cout << "Particle : " << particleGun->GetParticleDefinition()->GetParticleName() << G4endl
    << "Mass : " << sqrt(particleGun->GetParticleEnergy()*particleGun->GetParticleEnergy() - particleGun->GetParticleMomentum()*particleGun->GetParticleMomentum())<< G4endl
    << "Momentum : [" << particleGun->GetParticleMomentum()*particleGun->GetParticleMomentumDirection().x() << ", "
    << particleGun->GetParticleMomentum()*particleGun->GetParticleMomentumDirection().y() << ", "
    << particleGun->GetParticleMomentum()*particleGun->GetParticleMomentumDirection().z() << "]" << G4endl;  
  */
  particleGun->GeneratePrimaryVertex(anEvent);
  
  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("pi-"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPpim->Px(), fPpim->Py(), fPpim->Pz())/fPpim->P());
  particleGun->SetParticleEnergy(fPpim->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);

  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("proton"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPproton->Px(), fPproton->Py(), fPproton->Pz())/fPproton->P());
  particleGun->SetParticleEnergy(fPproton->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);

  fDecayChannel = PION_CH1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACEventGenerator::GenerateForPionCh2(G4ParticleGun *particleGun, G4Event *anEvent)
{
  Double_t masses1[2] = {PION_CHARGED_MASS, DELTA_MASS};
  Double_t masses2[2] = {PION_CHARGED_MASS, PROTON_MASS};
  G4ThreeVector position;
  fGenPhaseSpace1->SetDecay(*fPinitial, 2, masses1);
  while(true)
  { 
    fGenPhaseSpace1->Generate();

    fPpip = fGenPhaseSpace1->GetDecay(0);
    fPproton = fGenPhaseSpace1->GetDecay(1); // this is for delta
  
    fGenPhaseSpace2->SetDecay(*fPproton, 2, masses2);
    fGenPhaseSpace2->Generate();
  
    fPpim = fGenPhaseSpace2->GetDecay(0);
    fPproton = fGenPhaseSpace2->GetDecay(1);

    position = G4ThreeVector(0.*mm, 0.*mm, TARGET_CENTER_POSITION + (G4UniformRand() - 0.5)*TARGET_LENGTH);
    if(IsTouchDetector(fPpip, TPC_DOWNSTREAM_Z - position[2]) || IsTouchDetector(fPpim, TPC_DOWNSTREAM_Z - position[2]))
      break;
  }
  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("pi+"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPpip->Px(), fPpip->Py(), fPpip->Pz())/fPpip->P());
  particleGun->SetParticleEnergy(fPpip->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("pi-"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPpim->Px(), fPpim->Py(), fPpim->Pz())/fPpim->P());
  particleGun->SetParticleEnergy(fPpim->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);

  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("proton"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPproton->Px(), fPproton->Py(), fPproton->Pz())/fPproton->P());
  particleGun->SetParticleEnergy(fPproton->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);

  fDecayChannel = PION_CH2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACEventGenerator::GenerateForPionCh3(G4ParticleGun *particleGun, G4Event *anEvent)
{
  Double_t masses1[2] = {PION_CHARGED_MASS, DELTA_MASS};
  Double_t masses2[2] = {PION_CHARGED_MASS, PROTON_MASS};
  G4ThreeVector position;
  fGenPhaseSpace1->SetDecay(*fPinitial, 2, masses1);
  while(true)
  {
    fGenPhaseSpace1->Generate();

    fPpim = fGenPhaseSpace1->GetDecay(0);
    fPproton = fGenPhaseSpace1->GetDecay(1); // this is for delta++
    
    fGenPhaseSpace2->SetDecay(*fPproton, 2, masses2);
    fGenPhaseSpace2->Generate();
    
    fPpip = fGenPhaseSpace2->GetDecay(0);
    fPproton = fGenPhaseSpace2->GetDecay(1);

    position = G4ThreeVector(0.*mm, 0.*mm, TARGET_CENTER_POSITION + (G4UniformRand() - 0.5)*TARGET_LENGTH);
    if(IsTouchDetector(fPpip, TPC_DOWNSTREAM_Z - position[2]) || IsTouchDetector(fPpim, TPC_DOWNSTREAM_Z - position[2]))
      break;
  }
  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("pi+"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPpip->Px(), fPpip->Py(), fPpip->Pz())/fPpip->P());
  particleGun->SetParticleEnergy(fPpip->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("pi-"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPpim->Px(), fPpim->Py(), fPpim->Pz())/fPpim->P());
  particleGun->SetParticleEnergy(fPpim->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);

  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("proton"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPproton->Px(), fPproton->Py(), fPproton->Pz())/fPproton->P());
  particleGun->SetParticleEnergy(fPproton->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);

  fDecayChannel = PION_CH3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACEventGenerator::GenerateForPionCh4(G4ParticleGun *particleGun, G4Event *anEvent)
{
  Double_t masses1[2] = {RHO_MASS, PROTON_MASS};
  Double_t masses2[2] = {PION_CHARGED_MASS, PION_CHARGED_MASS};
  G4ThreeVector position;
  fGenPhaseSpace1->SetDecay(*fPinitial, 2, masses1);
  while(true)
  {
    fGenPhaseSpace1->Generate();

    fPpip = fGenPhaseSpace1->GetDecay(0); // this is for rho
    fPproton = fGenPhaseSpace1->GetDecay(1);
    
    fGenPhaseSpace2->SetDecay(*fPpip, 2, masses2);
    fGenPhaseSpace2->Generate();
    
    fPpip = fGenPhaseSpace2->GetDecay(0);
    fPpim = fGenPhaseSpace2->GetDecay(1);

    position = G4ThreeVector(0.*mm, 0.*mm, TARGET_CENTER_POSITION + (G4UniformRand() - 0.5)*TARGET_LENGTH);
    if(IsTouchDetector(fPpip, TPC_DOWNSTREAM_Z - position[2]) || IsTouchDetector(fPpim, TPC_DOWNSTREAM_Z - position[2]))
      break;
  }
  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("pi+"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPpip->Px(), fPpip->Py(), fPpip->Pz())/fPpip->P());
  particleGun->SetParticleEnergy(fPpip->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("pi-"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPpim->Px(), fPpim->Py(), fPpim->Pz())/fPpim->P());
  particleGun->SetParticleEnergy(fPpim->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);

  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("proton"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPproton->Px(), fPproton->Py(), fPproton->Pz())/fPproton->P());
  particleGun->SetParticleEnergy(fPproton->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);

  fDecayChannel = PION_CH4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACEventGenerator::GenerateForKaonCh1(G4ParticleGun *particleGun, G4Event *anEvent)
{
  Double_t masses[2] = {KAON_CHARGED_MASS, LAMBDA_MASS};
  G4ThreeVector position;
  fGenPhaseSpace1->SetDecay(*fPinitial, 2, masses);
  while(true)
  {
    fGenPhaseSpace1->Generate();

    fPkap = fGenPhaseSpace1->GetDecay(0);
    fPlambda = fGenPhaseSpace1->GetDecay(1);
  
    position = G4ThreeVector(0.*mm, 0.*mm, TARGET_CENTER_POSITION + (G4UniformRand() - 0.5)*TARGET_LENGTH);
    if(IsTouchDetector(fPkap, TPC_DOWNSTREAM_Z - position[2]))
      break;
  }
  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("kaon+"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPkap->Px(), fPkap->Py(), fPkap->Pz())/fPkap->P());
  particleGun->SetParticleEnergy(fPkap->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);
  
  particleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("lambda"));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(fPlambda->Px(), fPlambda->Py(), fPlambda->Pz())/fPlambda->P());
  particleGun->SetParticleEnergy(fPlambda->E()*GeV);
  particleGun->SetParticlePosition(position);
  particleGun->GeneratePrimaryVertex(anEvent);

  fDecayChannel = KAON_CH1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4ACEventGenerator::IsTouchDetector(const TLorentzVector *momentum, G4double dz)
{
  G4double rho = tan(momentum->Theta())*dz - TPC_EDGE_R;
  return (rho < 0) && (rho > -AC_DETECTOR_R); 
}
