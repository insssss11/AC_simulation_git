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
// $Id: G4ACRunAction.hh 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file G4ACRunActionMessenger.hh
/// \brief Definition of the G4ACRunActionMessenger class

#ifndef G4ACRunActionMessenger_h
#define G4ACRunActionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"
class G4ACRunAction;
class G4UIdirectory;
class G4UIcmdWithAString;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class G4ACRunActionMessenger : public G4UImessenger
{
  public:
    G4ACRunActionMessenger(G4ACRunAction *runAction);
    virtual ~G4ACRunActionMessenger();

    // virtual G4Run* GenerateRun();
    virtual void   SetNewValue(G4UIcommand*, G4String);

  private:
    G4ACRunAction *fRunAction;
    G4UIdirectory *fACdir;
    G4UIcmdWithAString *fFileCmd;
};

#endif

