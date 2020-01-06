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
// $Id: G4ACDetectorConstruction4.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file G4ACDetectorConstruction4.hh
/// \brief Definition of the G4ACDetectorConstruction4 class
// Design TYPE4 : with aerogel alignments shifted and the light guide rotated by 180 degrees

#ifndef G4ACDetectorConstruction4_h
#define G4ACDetectorConstruction4_h 1

#include "G4PVPlacement.hh"
#include "G4ACBaseDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"



class G4LogicalVolume;
class G4VPhysicalVolume;

/// Detector construction class to define materials and geometry.
class G4ACDetectorConstruction4 : public G4ACBaseDetectorConstruction
{
  public:
    G4ACDetectorConstruction4();
    virtual ~G4ACDetectorConstruction4();
    virtual G4VPhysicalVolume* Construct();
  protected:
    void ConstructHolder();
    void ConstructAerogel();
    void ConstructAerogelBoundary();
    
    G4ThreeVector
    *fAerogelPosVec4, *fAerogelPosVec5, *fAerogelPosVec6;

    G4VSolid
    *fSolidAerogel4, *fSolidAerogel5, *fSolidAerogel6;

    G4LogicalVolume
    *fLogicAerogel4, *fLogicAerogel5, *fLogicAerogel6;
};

namespace ac_detector_type4
{
  typedef G4ACDetectorConstruction4 G4ACDetectorConstruction;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

