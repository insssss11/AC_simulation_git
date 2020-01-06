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
// $Id: G4ACActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file G4ACActionInitialization.cc
/// \brief Implementation of the G4ACActionInitialization class

#include "G4ACEventAction.hh"
#include "G4ACActionInitialization.hh"
#include "G4ACPrimaryGeneratorAction.hh"
#include "G4ACRunAction.hh"
#include "G4ACSteppingAction.hh"
#include "G4ACStackingAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACActionInitialization::G4ACActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACActionInitialization::~G4ACActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACActionInitialization::BuildForMaster() const
{
  SetUserAction(new G4ACRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACActionInitialization::Build() const
{
  SetUserAction(new G4ACPrimaryGeneratorAction());

  G4ACRunAction* runAction = new G4ACRunAction();
  SetUserAction(runAction);
  
  G4ACEventAction* eventAction = new G4ACEventAction(runAction);
  SetUserAction(eventAction);
  
  SetUserAction(new G4ACSteppingAction(eventAction));
  SetUserAction(new G4ACStackingAction());
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
