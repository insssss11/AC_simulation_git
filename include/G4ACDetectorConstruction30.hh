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
// $Id: G4ACDetectorConstruction30.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file G4ACDetectorConstruction30.hh
/// \brief Definition of the G4ACDetectorConstruction30 class
// Design TYPE30  : The reflectors on both the holder and light guide are diffuse materials(TYPEC, tyvek, LUT model)

#ifndef G4ACDetectorConstruction30_h
#define G4ACDetectorConstruction30_h 1

#include "G4PVPlacement.hh"
#include "G4ACDetectorConstruction22.hh"
#include "globals.hh"
#include "G4Material.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

/// Detector construction class to define materials and geometry.

class G4ACDetectorConstruction30 : public G4ACDetectorConstruction22
{
  public:
    G4ACDetectorConstruction30();
    virtual ~G4ACDetectorConstruction30();

    virtual G4VPhysicalVolume *Construct();
  
  protected:
    void ConstructHolderBoundary();
    void ConstructGuideBoundary();
    void ConstructAerogel();    
    void ConstructMaterial();
};


namespace ac_detector_type30
{
  typedef G4ACDetectorConstruction30 G4ACDetectorConstruction;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

