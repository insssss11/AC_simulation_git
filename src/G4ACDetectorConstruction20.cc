//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
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
// $Id: G4ACDetectorConstruction20.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACDetectorConstruction20.cc
/// \brief Definition of the G4ACDetectorConstruction20 class

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4ExtrudedSolid.hh"
#include "G4GenericTrap.hh"
#include "G4Trap.hh"
#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "globals.hh"
#include <math.h>

#include "G4ACDetectorConstruction20.hh"
#include "G4ACFineMeshSD.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4ACDetectorConstruction20::G4ACDetectorConstruction20()
: G4ACBaseDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction20::~G4ACDetectorConstruction20()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACDetectorConstruction20::Construct()
{
  //.......oooOO0OOooo........oooOO0OOooo......
  

  /////////////////////////////////////////////////////////////////
  // world
  /////////////////////////////////////////////////////////////////
  G4double worldHalfSizeX = 200.*mm;
  G4double worldHalfSizeY = 300.*mm;
  G4double worldHalfSizeZ = 400.*mm;
  fSolidWorld   = new G4Box("World", worldHalfSizeX, worldHalfSizeY, worldHalfSizeZ);
  fLogicWorld   = new G4LogicalVolume(fSolidWorld, fMatAir, "World");
  fPhysWorld    = new G4PVPlacement(0,    G4ThreeVector(),  fLogicWorld,  "World",  0, false, 0, fCheckOverlaps);

  /////////////////////////////////////////////////////////////////
  // construct Aerogel Holder and Light Guide, Scorer, Aerogel
  /////////////////////////////////////////////////////////////////
  ConstructHolder();
  ConstructGuide();
  ConstructScorer();
  ConstructInnerAir();
  ConstructAerogel();

  ConstructHolderBoundary(); // construct boundary processes
  ConstructGuideBoundary();
  ConstructAerogelBoundary();
  ConstructGlassBoundary();
  
  return fPhysWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction20::ConstructGuide()
{
  /////////////////////////////////////////////////////////////////
  // Light Guide Part 1 In
  /////////////////////////////////////////////////////////////////
  G4double guide1InZ = 70.*mm/2;
  // bottop face
  G4double guide1In1X = fCrossX;
  G4double guide1In1Y = fCrossY;

  // top fae
  G4double guide1In2X = 80.*mm/2;
  G4double guide1In2Y = 20.*mm/2;

  std::vector<G4TwoVector> guide1InPtVec(8);
  guide1InPtVec[0] = G4TwoVector(-guide1In1X, fBoxThickness);
  guide1InPtVec[1] = G4TwoVector(-guide1In1X, 2*guide1In1Y + fBoxThickness);
  guide1InPtVec[2] = G4TwoVector(guide1In1X, 2*guide1In1Y + fBoxThickness);
  guide1InPtVec[3] = G4TwoVector(guide1In1X, fBoxThickness);
  guide1InPtVec[4] = G4TwoVector(-guide1In2X, fBoxThickness);
  guide1InPtVec[5] = G4TwoVector(-guide1In2X, 2*guide1In2Y + fBoxThickness);
  guide1InPtVec[6] = G4TwoVector(guide1In2X, 2*guide1In2Y + fBoxThickness);
  guide1InPtVec[7] = G4TwoVector(guide1In2X, fBoxThickness);

  G4double guide1InBoxX = guide1In2X, guide1InBoxY = 0.5*fBoxThickness, guide1InBoxZ = guide1In2X*cos(fAngleTilt);

  G4VSolid *solidGuide1In1 = new G4GenericTrap("Guide1In1", guide1InZ, guide1InPtVec);
  G4VSolid *solidGuide1In2 = new G4Box("Guide1In2", guide1InBoxX, guide1InBoxY, guide1InBoxZ);
  
  fGuideInPosVec = new G4ThreeVector(0.*mm, 0.*mm, guide1InZ);
  G4ThreeVector guide1In2PosVec(0.*mm, 0.5*fBoxThickness, 2*guide1InZ - guide1InBoxZ);
  G4UnionSolid *solidGuide1In = new G4UnionSolid("Guide1In", solidGuide1In1, solidGuide1In2, 0, guide1In2PosVec - *fGuideInPosVec);
  /////////////////////////////////////////////////////////////////
  // Light Guide Part 1 Out
  /////////////////////////////////////////////////////////////////  
  G4double guide1OutZ = guide1InZ + fBoxThickness/2;

  // bottom face half length
  G4double guide1Out1X = fCrossX + fBoxThickness;
  G4double guide1Out1Y = fCrossY + fBoxThickness;

  // top face
  G4double guide1Out2X = guide1In2X + fBoxThickness;
  G4double guide1Out2Y = guide1In2Y + fBoxThickness;
  std::vector<G4TwoVector> guide1OutPtVec(8);
  guide1OutPtVec[0] = G4TwoVector(-guide1Out1X, 0.*mm);
  guide1OutPtVec[1] = G4TwoVector(-guide1Out1X, 2*guide1Out1Y);
  guide1OutPtVec[2] = G4TwoVector(guide1Out1X, 2*guide1Out1Y);
  guide1OutPtVec[3] = G4TwoVector(guide1Out1X, 0.*mm);
  guide1OutPtVec[4] = G4TwoVector(-guide1Out2X, 0.*mm);
  guide1OutPtVec[5] = G4TwoVector(-guide1Out2X, 2*guide1Out2Y);
  guide1OutPtVec[6] = G4TwoVector(guide1Out2X, 2*guide1Out2Y);
  guide1OutPtVec[7] = G4TwoVector(guide1Out2X, 0.*mm);

  // to make the centers of masses locate at Origin
  G4VSolid *solidGuide1Out = new G4GenericTrap("Guide1Out", guide1OutZ, guide1OutPtVec);
  fGuideOutPosVec = new G4ThreeVector(0.*mm, 0.*mm, guide1OutZ);

  /////////////////////////////////////////////////////////////////
  // Light Guide Part 1 Reflector
  /////////////////////////////////////////////////////////////////  
  G4double guideReflector1OutZ = guide1InZ + fReflectorThickness/2;
  // top face
  std::vector<G4TwoVector> guideReflector1OutPtVec(8);
  guideReflector1OutPtVec[0] = guide1InPtVec[0] + G4TwoVector(-fReflectorThickness,  -fReflectorThickness);
  guideReflector1OutPtVec[1] = guide1InPtVec[1] + G4TwoVector(-fReflectorThickness,  fReflectorThickness);
  guideReflector1OutPtVec[2] = guide1InPtVec[2] + G4TwoVector(fReflectorThickness,  fReflectorThickness);
  guideReflector1OutPtVec[3] = guide1InPtVec[3] + G4TwoVector(fReflectorThickness,  -fReflectorThickness);
  guideReflector1OutPtVec[4] = guide1InPtVec[4] + G4TwoVector(-fReflectorThickness,  -fReflectorThickness);
  guideReflector1OutPtVec[5] = guide1InPtVec[5] + G4TwoVector(-fReflectorThickness,  fReflectorThickness);
  guideReflector1OutPtVec[6] = guide1InPtVec[6] + G4TwoVector(fReflectorThickness,  fReflectorThickness);
  guideReflector1OutPtVec[7] = guide1InPtVec[7] + G4TwoVector(fReflectorThickness,  -fReflectorThickness);

  // to make the centers of masses locate at Origin
  G4VSolid *solidGuideReflector1Out1 = new G4GenericTrap("GuideReflector1Out1", guideReflector1OutZ, guideReflector1OutPtVec);
  fGuideReflectorPosVec = new G4ThreeVector(0.*mm, 0.*mm, guideReflector1OutZ);
  G4VSolid *solidGuideReflector1Out2 = new G4Box("GuideReflector1Out2", guide1InBoxX + fReflectorThickness, 0.501*(fBoxThickness - fReflectorThickness), guide1InBoxZ + fReflectorThickness);
  G4ThreeVector guideReflector1Out2PosVec(0.*mm, 0.5*(fBoxThickness - fReflectorThickness), 2*guide1InZ - guide1InBoxZ);
  G4UnionSolid *solidGuideReflector1Out = new G4UnionSolid("GuideReflector1Out", solidGuideReflector1Out1, solidGuideReflector1Out2, 0, guideReflector1Out2PosVec - *fGuideReflectorPosVec);

  /////////////////////////////////////////////////////////////////
  // Light Guide Part2(tilt) In
  /////////////////////////////////////////////////////////////////
  // sakak-kidung
  G4double tiltInX = guide1InBoxZ, tiltInY = tiltInX*tan(fAngleTilt), tiltInZ = guide1InBoxX;
  std::vector <G4TwoVector> tiltInPolygon(3);
  tiltInPolygon[0] = G4TwoVector(-2.*tiltInX - fBoxThickness, 0.*mm);
  tiltInPolygon[1] = G4TwoVector(-fBoxThickness, 0.*mm);
  tiltInPolygon[2] = G4TwoVector(-fBoxThickness, -2.*tiltInY);
  G4ExtrudedSolid *solidTiltIn = new G4ExtrudedSolid("TiltIn", tiltInPolygon, tiltInZ, 
    G4TwoVector(), 1., G4TwoVector(), 1.);

  /////////////////////////////////////////////////////////////////
  // Light Guide Part2(tilt) Out
  /////////////////////////////////////////////////////////////////
  // samkak-kidung
  std::vector <G4TwoVector> tiltOutPolygon(3);
  G4double tiltOutX = tiltInX + fBoxThickness, tiltOutY = tiltOutX*tan(fAngleTilt), tiltOutZ = tiltInZ + fBoxThickness;
  G4ThreeVector tiltOutPosVec(0.*mm, 0.*mm, 2.*(guide1InZ - guide1InBoxZ) - fBoxThickness);
  G4RotationMatrix *rotationTilt = new G4RotationMatrix();
  rotationTilt->rotateY(-90.*deg);
  // define bottom face before rotation 90 deg with X axis
  tiltOutPolygon[0] = G4TwoVector(-2*tiltOutX, 0.*mm);
  tiltOutPolygon[1] = G4TwoVector(0.*mm,0.*mm);
  tiltOutPolygon[2] = G4TwoVector(0.*mm, -2.*tiltOutY);

  // before lotation with x axis
  G4ExtrudedSolid *solidTiltOut = new G4ExtrudedSolid("TiltOut", tiltOutPolygon, tiltOutZ, 
    G4TwoVector(), 1., G4TwoVector(), 1.);

  /////////////////////////////////////////////////////////////////
  // Light Guide Part3(tilt) Reflector
  /////////////////////////////////////////////////////////////////
  // sakak-kidung
  std::vector <G4TwoVector> guideReflectorOutPolygon(3);
  guideReflectorOutPolygon[0] = tiltInPolygon[0] + G4TwoVector(-fReflectorThickness, 0.*mm);
  guideReflectorOutPolygon[1] = tiltInPolygon[1] + G4TwoVector(fReflectorThickness, 0.*mm);
  guideReflectorOutPolygon[2] = tiltInPolygon[2] + G4TwoVector(fReflectorThickness, -2.*fReflectorThickness*tan(fAngleTilt));
  G4ExtrudedSolid *solidGuideReflector2Out = new G4ExtrudedSolid("GuideReflector2Out", guideReflectorOutPolygon, tiltOutZ - fBoxThickness + fReflectorThickness, 
    G4TwoVector(), 1., G4TwoVector(), 1.);

  /////////////////////////////////////////////////////////////////
  // Tubs : PMT hole to sutract with box and (tiltInHole = In =  PMT hole)
  /////////////////////////////////////////////////////////////////
  G4double tiltHoleZ = 0.5*fBoxThickness*sin(fAngleTilt);
  G4ThreeVector tiltHolePosVec(-tiltOutX + tiltHoleZ*sin(fAngleTilt), -tiltOutY + tiltHoleZ*cos(fAngleTilt), 0.*mm);
  fRotationScorer = new G4RotationMatrix();
  fRotationScorer->rotate(-90.*deg, G4ThreeVector(cos(180.*deg - fAngleTilt), sin(180.*deg - fAngleTilt)));
  G4Tubs *solidTiltHole = new G4Tubs("TiltHole", 0.*mm, fCathodeRadius, tiltHoleZ, 0.*deg, 360.*deg);
  
  G4UnionSolid *solidTiltInHole = new G4UnionSolid("TiltInHole", solidTiltIn, solidTiltHole, fRotationScorer, tiltHolePosVec);

  // G4UnionSolid *solidGuideTubOut     = new G4UnionSolid("GuideTubOut", solidGuide1Out, solidGuide2Out, rotationGuide2, guide2PosVec - *fGuideOutPosVec);
  // G4UnionSolid *solidGuideTubIn      = new G4UnionSolid("GuideTubIn", solidGuide1In, solidGuide2In, rotationGuide2, guide2PosVec - *fGuideInPosVec);

  fSolidGuideOut  = new G4UnionSolid("GuideOut", solidGuide1Out, solidTiltOut, rotationTilt, tiltOutPosVec - *fGuideOutPosVec);
  fSolidGuideIn   = new G4UnionSolid("GuideIn", solidGuide1In, solidTiltInHole, rotationTilt, tiltOutPosVec- *fGuideInPosVec);
  G4UnionSolid *solidGuideReflectorOutWithoutHole = new G4UnionSolid("GuideReflectorOutWithoutHole", solidGuideReflector1Out, solidGuideReflector2Out, rotationTilt,  tiltOutPosVec - *fGuideReflectorPosVec);
  G4SubtractionSolid *solidGuideCaseWithoutHole = new G4SubtractionSolid("GuideCaseWithoutHole", fSolidGuideOut, solidGuideReflectorOutWithoutHole, 0, *fGuideReflectorPosVec - *fGuideOutPosVec);
  fSolidGuideCase = new G4SubtractionSolid("GuideCase", solidGuideCaseWithoutHole, fSolidGuideIn, 0, *fGuideInPosVec - *fGuideOutPosVec);
  fSolidGuideReflector = new G4SubtractionSolid("GuideReflector", solidGuideReflectorOutWithoutHole, fSolidGuideIn, 0, *fGuideInPosVec - *fGuideReflectorPosVec);

  fLogicGuideCase = new G4LogicalVolume(fSolidGuideCase, fMatBox, "GuideCase");
  fLogicGuideReflector = new G4LogicalVolume(fSolidGuideReflector, fMatReflector, "GuideReflector");

  fPhysGuideCase  = new G4PVPlacement(0,    *fGuideOutPosVec,  fLogicGuideCase,  "GuideCase",  fLogicWorld, false, 0, fCheckOverlaps);
  fPhysGuideReflector  = new G4PVPlacement(0,    *fGuideReflectorPosVec,  fLogicGuideReflector,  "GuideReflector",  fLogicWorld, false, 0, fCheckOverlaps);
  
  fScorerPosVec = new G4ThreeVector(tiltOutPosVec + G4ThreeVector(0.*mm, -tiltOutY, tiltOutX) + fScorerZ*G4ThreeVector(0.*mm, -cos(fAngleTilt), sin(fAngleTilt)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction20::ConstructScorer()
{
  G4RotationMatrix *rotationScorer = new G4RotationMatrix();
  rotationScorer->rotateX(fAngleTilt + 90.*deg);
  /////////////////////////////////////////////////////////////////
  // PMT = Scorer volume (this has Glass, cathode volume as daughter volumes)
  /////////////////////////////////////////////////////////////////
  fSolidScorer = new G4Tubs("Scorer", 0.*mm, fPMTradius, fScorerZ, 0.*deg, 360.*deg);
  // the photon go to the cathode of FMPMT through UV silicon cookie(glass)
  fSolidGlass  = new G4Tubs("Glass", 0.*mm, fCathodeRadius, fGlassThickness, 0.*deg, 360.*deg);
  // cathode
  fSolidCathode  = new G4Tubs("Cathode", 0.*mm, fCathodeRadius, fCathodeZ, 0.*deg, 360.*deg);
  fLogicScorer  = new G4LogicalVolume(fSolidScorer, fMatScorer, "Scorer");
  fLogicGlass   = new G4LogicalVolume(fSolidGlass, fMatGlass, "Glass");
  fLogicCathode = new G4LogicalVolume(fSolidCathode, fMatAir, "Cathode");
  fPhysScorer     = new G4PVPlacement(rotationScorer, *fScorerPosVec, fLogicScorer, "Scorer", fLogicWorld, false, 0, fCheckOverlaps);
  fPhysGlass      = new G4PVPlacement(0,    G4ThreeVector(0.*mm, 0.*mm, fScorerZ - fGlassThickness) ,    fLogicGlass,    "Glass", fLogicScorer, false, 0, fCheckOverlaps);
  fPhysCathode    = new G4PVPlacement(0,    G4ThreeVector(0.*mm, 0.*mm, fScorerZ - 2*fGlassThickness - fCathodeZ)  ,  fLogicCathode,  "Cathode", fLogicScorer, false, 0, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
