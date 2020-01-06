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
// $Id: G4ACBaseDetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file G4ACBaseDetectorConstruction.hh
/// \brief Definition of the G4ACBaseDetectorConstruction class
// Design DEFAULT : Default design with new aerogel shape

#ifndef G4ACBaseDetectorConstruction_h
#define G4ACBaseDetectorConstruction_h 1

#include "G4PVPlacement.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "TSpline.h"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;

const G4double TARGET_LENGTH = 150.*mm, TARGET_CENTER_POSITION = -529.4*mm,
  AC_DETECTOR_R = 256.3*mm, TPC_EDGE_R = 800.*cos(30.*deg)*mm, TPC_DOWNSTREAM_Z = -4.4*mm;

/// Detector construction class to define materials and geometry.
class G4ACBaseDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    G4ACBaseDetectorConstruction();
    G4ACBaseDetectorConstruction(G4double crossX, G4double crossY, G4double aerogelThickness = 20.*mm/2, G4double aerogelGap =  0.1*mm/2, G4double angleTilt = 30.*deg);
    G4ACBaseDetectorConstruction(G4double BoxThickness);
    virtual ~G4ACBaseDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

    G4LogicalVolume* GetScorerVolume() const { return fLogicScorer;} 
    G4double GetBoxThickness() const {return fBoxThickness;} 
    G4double GetAngleTilt() const {return fAngleTilt;} 
    G4double GetAerogelGap() const {return fAerogelGap;}
  protected:
  // these methode construct solid volume of Detector components
    void ConstructHolder();
    void ConstructGuide();
    void ConstructScorer();
    void ConstructInnerAir();
    void ConstructAerogel();

    void ConstructMaterial();
    void ConstructSDandField();
    
    void ConstructGuideBoundary();
    void ConstructHolderBoundary();
    void ConstructAerogelBoundary();
    void ConstructGlassBoundary();
    // construct data splines using TSpline3 class
    void ConstructSpline();
    // return in unit of cm
    G4double GetAerogelAbs(G4double PeV);
    G4double GetAerogelScat(G4double PeV);

    // information of cross section of Guide and Holder
    const G4double fCrossX, fCrossY; // the thickness is not included

    // aerogel dimension(piece 1,2,3)
    const G4double fAerogelThickness;
    const G4double fAerogel1Z;
    const G4double fAerogel11X, fAerogel11Y;
    const G4double fAerogel12X, fAerogel12Y;
    const G4double fAerogel2Z;
    const G4double fAerogel21X, fAerogel21Y;
    const G4double fAerogel22X, fAerogel22Y;
    const G4double fAerogel3Z;
    const G4double fAerogel31X, fAerogel31Y;
    const G4double fAerogel32X, fAerogel32Y;

    // parameters for geometry
    const G4double fBoxThickness; // thickness of detector case
    const G4double fReflectorThickness;  
    const G4double fAngleTilt; // angle of tiling board
    const G4double fAerogelGap; // gap distance between aerogels
    const G4double fPMTradius; // diameter of PMT
    const G4double fCathodeRadius; // the radius of Cathode of PMT
    const G4double fGlassThickness;
    const G4double fScorerZ; // half height of PMT(= scorer)
    const G4double fCathodeZ;
    // the logical and physical volumes
    const G4bool fCheckOverlaps;


    G4VSolid
    *fSolidWorld,
    *fSolidHolderCase, *fSolidHolderOut, *fSolidHolderIn,
    *fSolidGuideCase, *fSolidGuideOut, *fSolidGuideIn,
    *fSolidOut, *fSolidIn,
    *fSolidScorer, *fSolidGlass, *fSolidCathode,
    *fSolidAerogel1, *fSolidAerogel2, *fSolidAerogel3,
    *fSolidHolderReflector, *fSolidGuideReflector;
    
    G4LogicalVolume
    *fLogicWorld,
    *fLogicHolderCase, // *fLogicHolderOut, *fLogicHolderIn,
    *fLogicGuideCase, // *fLogicGuideOut, *fLogicGuideIn,
    // *fLogicOut,
    *fLogicIn,
    *fLogicScorer, *fLogicGlass, *fLogicCathode,
    *fLogicAerogel1, *fLogicAerogel2, *fLogicAerogel3,
    *fLogicHolderReflector, *fLogicGuideReflector;

    G4VPhysicalVolume
    *fPhysWorld,
    *fPhysHolderCase, //  *fPhysHolderOut, *fPhysHolderIn,
    *fPhysGuideCase, //  *fPhysGuideOut, *fPhysGuideIn,
    // *fPhysOut,
    *fPhysIn,
    *fPhysScorer, *fPhysGlass, *fPhysCathode,
    *fPhysAerogel1, *fPhysAerogel2, *fPhysAerogel3,
    *fPhysHolderReflector, *fPhysGuideReflector;

    G4ThreeVector
    *fHolderOutPosVec, *fHolderInPosVec,
    *fGuideOutPosVec, *fGuideInPosVec,
    *fScorerPosVec,
    *fAerogelPosVec1, *fAerogelPosVec2, *fAerogelPosVec3,
    *fHolderReflectorPosVec, *fGuideReflectorPosVec;
    
    G4RotationMatrix
    *fRotationScorer;

    // material definetion
    G4Material *fMatAir, *fMatBox, *fMatGlass, *fMatScorer, *fMatAerogel, *fMatReflector;
    TSpline3 *fSplineAeroAbs, *fSplineESRref, *fSplineAlMylarRef, *fSplineTyvekRef;
};

namespace ac_detector_base{
  typedef G4ACBaseDetectorConstruction G4ACDetectorConstruction;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

