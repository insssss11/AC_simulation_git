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
// $Id: G4ACPrimaryGeneratorActionMessenger.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACPrimaryGeneratorActionMessenger.cc
/// \brief Implementation of the G4ACPrimaryGeneratorActionMessenger class

#include "G4ACPrimaryGeneratorActionMessenger.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include <math.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACPrimaryGeneratorActionMessenger::G4ACPrimaryGeneratorActionMessenger(G4ACPrimaryGeneratorAction* primaryAction)
: G4UImessenger(),
  fPrimaryAction(primaryAction), fGunDir(0), fGunDirSingle(0), fGunDirLeps2(0),
  fIncidentPositionCmd(0), fIncidentAngleCmd(0), fMomentumAmpCmd(0), fBeamAmpCmd(0),
  fParticleCmd(0), fGunModeCmd(0),
  fPrintStatCmd(0)
{
  fGunDir = new G4UIdirectory("/simulateAC/gun/");
  fGunDir->SetGuidance("PrimaryGenerator control");

  fGunDirSingle = new G4UIdirectory("/simulateAC/gun/single/");
  fGunDirSingle->SetGuidance("Gun setting for only one detector simulation");  

  fGunDirLeps2 = new G4UIdirectory("/simulateAC/gun/leps2/");
  fGunDirLeps2->SetGuidance("Gun setting for entire AC1 detectors simulation(please change gun mode to change decay channels)");

  fIncidentAngleCmd = new G4UIcmdWithADoubleAndUnit("/simulateAC/gun/single/incidentAngle", this);
  fIncidentAngleCmd->SetGuidance("Set incident angle w.r.t. the detector plane(for single gunmode)");
  fIncidentAngleCmd->SetParameterName("IncidentAngle", true);
  fIncidentAngleCmd->SetDefaultUnit("deg");
  fIncidentAngleCmd->SetRange("IncidentAngle>=0");
  fIncidentAngleCmd->SetDefaultValue(40);

  fIncidentPositionCmd = new G4UIcmdWithADoubleAndUnit("/simulateAC/gun/single/incidentPosition", this);
  fIncidentPositionCmd->SetGuidance("Set incident position from the bottom edge(for single gunmode)");
  fIncidentPositionCmd->SetParameterName("IncidentPosition", true);
  fIncidentPositionCmd->SetDefaultUnit("mm");
  fIncidentPositionCmd->SetUnitCandidates("mm cm m");
  fIncidentPositionCmd->SetDefaultValue(100);
  fIncidentPositionCmd->SetRange("IncidentPosition>=0");

  fMomentumAmpCmd = new G4UIcmdWithADoubleAndUnit("/simulateAC/gun/single/momentumAmp", this);
  fMomentumAmpCmd->SetGuidance("Set incident position from w.r.t. the bottom(for single gunmode)");
  fMomentumAmpCmd->SetParameterName("MomentumAmp", false);
  fMomentumAmpCmd->SetDefaultUnit("MeV");  
  fMomentumAmpCmd->SetUnitCandidates("eV keV MeV GeV");
  fMomentumAmpCmd->SetRange("MomentumAmp>0");

  fBeamAmpCmd = new G4UIcmdWithADoubleAndUnit("/simulateAC/gun/leps2/beamAmp", this);
  fBeamAmpCmd->SetGuidance("Set gamma energy incident to target(for leps2 gunmode)");
  fBeamAmpCmd->SetParameterName("BeamAmp", false);
  fBeamAmpCmd->SetDefaultUnit("GeV");
  fBeamAmpCmd->SetUnitCandidates("eV keV MeV GeV");
  fBeamAmpCmd->SetRange("BeamAmp>0");

  fParticleCmd = new G4UIcmdWithAString("/simulateAC/gun/single/particle", this);
  fParticleCmd->SetGuidance("Define primary particle");
  fParticleCmd->SetParameterName("Particle", false);

  fGunModeCmd = new G4UIcmdWithAString("/simulateAC/gun/mode", this);
  fGunModeCmd->SetGuidance("Set gun Shooting modes(to simulate of single detector sim or LEPS2)");
  fGunModeCmd->SetCandidates("single leps2");
  fGunModeCmd->SetParameterName("GunMode", false);
  fGunModeCmd->SetDefaultValue("single");

  fPrintStatCmd = new G4UIcmdWithoutParameter("/simulateAC/gun/printStat", this);
  fPrintStatCmd->SetGuidance("Print Particle Gun Status");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACPrimaryGeneratorActionMessenger::~G4ACPrimaryGeneratorActionMessenger()
{
  delete fGunDir;
  delete fGunDirSingle;  
  delete fGunDirLeps2;
  delete fIncidentAngleCmd;
  delete fIncidentPositionCmd;
  delete fMomentumAmpCmd;
  delete fBeamAmpCmd;
  delete fParticleCmd;
  delete fGunModeCmd;
  delete fPrintStatCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACPrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
  if(command == fIncidentAngleCmd)
      fPrimaryAction->SetIncidentAngle(fIncidentAngleCmd->GetNewDoubleValue(newValue));
  else if(command == fIncidentPositionCmd)
      fPrimaryAction->SetIncidentPosition(fIncidentPositionCmd->GetNewDoubleValue(newValue));
  else if(command == fMomentumAmpCmd)
      fPrimaryAction->SetMomentumAmp(fMomentumAmpCmd->GetNewDoubleValue(newValue));
  else if(command == fBeamAmpCmd)
      fPrimaryAction->GetEventGenerator()->SetGammaE(fBeamAmpCmd->GetNewDoubleValue(newValue));
  else if(command == fParticleCmd)
      fPrimaryAction->SetParticle(newValue.data());
  else if(command == fGunModeCmd)
      fPrimaryAction->SetGunMode(newValue.data());
  else if(command == fPrintStatCmd)
      fPrimaryAction->PrintStat();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
