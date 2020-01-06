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
/// \file G4ACPrimaryGeneratorActionMessenger.hh
/// \brief Definition of the G4ACPrimaryGeneratorActionMessenger class

#ifndef G4ACPrimaryGeneratorActionMessenger_h
#define G4ACPrimaryGeneratorActionMessenger_h 1

#include "G4ACPrimaryGeneratorAction.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"

class G4ACPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

class G4ACPrimaryGeneratorActionMessenger : public G4UImessenger
{
  public:
    G4ACPrimaryGeneratorActionMessenger(G4ACPrimaryGeneratorAction*);
    virtual ~G4ACPrimaryGeneratorActionMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);
  private:
    G4ACPrimaryGeneratorAction *fPrimaryAction;
    G4UIdirectory              *fGunDir;
    G4UIdirectory              *fGunDirSingle; // for single gun
    G4UIdirectory              *fGunDirLeps2; // for leps2 gun mode

    G4UIcmdWithADoubleAndUnit  *fIncidentPositionCmd;
    G4UIcmdWithADoubleAndUnit  *fIncidentAngleCmd;
    G4UIcmdWithADoubleAndUnit  *fMomentumAmpCmd;
    G4UIcmdWithADoubleAndUnit  *fBeamAmpCmd;
    G4UIcmdWithAString *fParticleCmd;
    G4UIcmdWithAString *fGunModeCmd;
    G4UIcmdWithoutParameter *fPrintStatCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
