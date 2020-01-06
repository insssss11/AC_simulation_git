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
// $Id: G4ACDetectorConstructionLEPS2.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file G4ACDetectorConstructionLEPS2.hh
/// \brief Definition of the G4ACDetectorConstructionLEPS2 class
// Design LEPS2  : Full AC construction for LEPS2/SPring8 experiment (TYPE20, 21, 22 used)

#ifndef G4ACDetectorConstructionLEPS2_h
#define G4ACDetectorConstructionLEPS2_h 1

#include "G4PVPlacement.hh"
#include "G4ACBaseDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

/// Detector construction class to define materials and geometry.

class G4ACDetectorConstructionLEPS2 : public G4ACBaseDetectorConstruction
{
  public:
    G4ACDetectorConstructionLEPS2();
    virtual ~G4ACDetectorConstructionLEPS2();

    virtual G4VPhysicalVolume *Construct();
  
  private:

    const G4double fTPCradius; // shortest distance
    const G4double fOffsetAngle; // the staring angle of iteration when repeting AC1 volumes
    const G4double fOffsetTPCrocation; // the whole position will move along z axis as this quantity
    const G4int fNelem;

    G4VSolid *fSolidGuideCaseA, *fSolidGuideCaseB, *fSolidGuideCaseC,
      *fSolidGuideInA, *fSolidGuideInB, *fSolidGuideInC,
      *fSolidGuideOutA, *fSolidGuideOutB, *fSolidGuideOutC;
    G4LogicalVolume *fLogicGuideCaseA, *fLogicGuideCaseB, *fLogicGuideCaseC;
    
    G4PVPlacement **fPhysHolderCaseArr, **fPhysGuideCaseArr, **fPhysInArr, 
     **fPhysScorerArr, **fPhysGlassArr, **fPhysCathodeArr,
     *fPhysAerogel1Arr, *fPhysAerogel2Arr, *fPhysAerogel3Arr;
    
    G4ThreeVector *fGuideInPosVecA, *fGuideInPosVecB, *fGuideInPosVecC,
      *fGuideOutPosVecA, *fGuideOutPosVecB, *fGuideOutPosVecC,
      *fScorerPosVecA, *fScorerPosVecB, *fScorerPosVecC;

    G4RotationMatrix **fRotationElement, *fRotationScorerA, *fRotationScorerB, *fRotationScorerC,
    *fRotationGuideA, *fRotationGuideB, *fRotationGuideC;

  private:
    void ConstructHolders();
    /*
    void ConstructGuides();
    void ConstructScorers();
    */
    void ConstructInnerAirs();
    void ConstructAerogels();
    /*
    void ConstructGuideA();
    void constrcutGuideB();
    void constrcutGuideB();

    void ConstructHolderBoundaries(); // construct boundary processes
    void ConstructGuideBoundaries();
    void ConstructAerogelBoundaries();
    void ConstructGlassBoundaries();
    */

   virtual void ConstructSDandField();
};


namespace ac_detector_type_leps2
{
  typedef G4ACDetectorConstructionLEPS2 G4ACDetectorConstruction;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

