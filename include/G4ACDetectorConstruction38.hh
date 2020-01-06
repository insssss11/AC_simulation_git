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
// $Id: G4ACDetectorConstruction38.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file G4ACDetectorConstruction38.hh
/// \brief Definition of the G4ACDetectorConstruction38 class
// Design TYPE38  : Same with TYPE20 but reflector on Holder is diffusive material(tyvek, TYPEA)

#ifndef G4ACDetectorConstruction38_h
#define G4ACDetectorConstruction38_h 1

#include "G4PVPlacement.hh"
#include "G4ACDetectorConstruction21.hh"
#include "globals.hh"
#include "G4Material.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

/// Detector construction class to define materials and geometry.

class G4ACDetectorConstruction38 : public G4ACDetectorConstruction21
{
  public:
    G4ACDetectorConstruction38();
    virtual ~G4ACDetectorConstruction38();

    virtual G4VPhysicalVolume *Construct();
  
  protected:
    void ConstructGuideBoundary();
    void ConstructHolderBoundary();
};


namespace ac_detector_type38
{
  typedef G4ACDetectorConstruction38 G4ACDetectorConstruction;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

