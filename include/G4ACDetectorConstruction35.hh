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
// $Id: G4ACDetectorConstruction35.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file G4ACDetectorConstruction35.hh
/// \brief Definition of the G4ACDetectorConstruction35 class
// Design TYPE35  : The 30-segment-design to be located at the TPC edge(ESR, TYPEB, unified model)

#ifndef G4ACDetectorConstruction35_h
#define G4ACDetectorConstruction35_h 1

#include "G4PVPlacement.hh"
#include "G4ACDetectorConstruction21.hh"
#include "globals.hh"
#include "G4Material.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

/// Detector construction class to define materials and geometry.

class G4ACDetectorConstruction35 : public G4ACDetectorConstruction21
{
  public:
    G4ACDetectorConstruction35();
    virtual ~G4ACDetectorConstruction35();

    virtual G4VPhysicalVolume *Construct();
  protected:
    void ConstructGuideBoundary();
    void ConstructHolderBoundary();

  protected:
};


namespace ac_detector_type35
{
  typedef G4ACDetectorConstruction35 G4ACDetectorConstruction;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

