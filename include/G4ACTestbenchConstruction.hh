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
// $Id: G4ACTestbenchConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file G4ACTestbenchConstruction.hh
/// \brief Definition of the G4ACTestbenchConstruction class
// Design TESTBENCH  : The aerogel led test performed at 2019.08.30

#ifndef G4ACTestbenchConstruction_h
#define G4ACTestbenchConstruction_h 1

#include "G4PVPlacement.hh"
#include "G4ACBaseDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

/// Detector construction class to define materials and geometry.

class G4ACTestbenchConstruction : public G4ACBaseDetectorConstruction
{
  public:
    G4ACTestbenchConstruction();
    virtual ~G4ACTestbenchConstruction();

    virtual G4VPhysicalVolume *Construct();
  
  protected:
    void ConstructAerogel();
    void ConstructAerogelBoundary();
    void ConstructInnerAir();
    void ConstructGuide();
    void ConstructGuideBoundary();
    void ConstructScorer();

  protected:
    G4double fAerogel11X = 116.*mm/2, fAerogel12X = 91.5*mm/2, fAerogel1Y = 118.*mm/2, fAerogel1Z = fAerogelThickness; // 11-5a
    G4double fAerogel21X = 120.*mm/2, fAerogel22X = 95.1*mm/2, fAerogel2Y = 120.*mm/2, fAerogel2Z = fAerogelThickness; // 11-5b

  private:
    void ConstructAerogel115a();
    void ConstructAerogel115b();
};


namespace ac_testbench
{
  typedef G4ACTestbenchConstruction G4ACDetectorConstruction;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

