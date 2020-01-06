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
// $Id: G4ACRunActionMessenger.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file G4ACRunActionMessenger.cc
/// \brief Implementation of the G4ACRunActionMessenger class

#include "G4ACRunActionMessenger.hh"
#include "G4ACRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4ACRunActionMessenger::G4ACRunActionMessenger(G4ACRunAction* runAction)
: G4UImessenger(), fRunAction(runAction), fACdir(0), fFileCmd(0)
{
  fACdir = new G4UIdirectory("/simulateAC/save/");
  fACdir->SetGuidance("Custom UI commands of AC simulation");

  fFileCmd = new G4UIcmdWithAString("/simulateAC/save/fileName", this);
  fFileCmd->SetParameterName("fileName",true);
  fFileCmd->SetDefaultValue("simulateAC");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACRunActionMessenger::~G4ACRunActionMessenger()
{
  delete fACdir;
  delete fFileCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACRunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fFileCmd)
  {
    fRunAction->SetFileName(newValue);
  }
}