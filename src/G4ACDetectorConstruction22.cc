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
// $Id: G4ACDetectorConstruction22.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACDetectorConstruction22.cc
/// \brief Definition of the G4ACDetectorConstruction22 class

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

#include "G4ACDetectorConstruction22.hh"
#include "G4ACFineMeshSD.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction22::G4ACDetectorConstruction22()
: G4ACDetectorConstruction21(24.*deg, 100.*mm)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction22::~G4ACDetectorConstruction22()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACDetectorConstruction22::Construct()
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

void G4ACDetectorConstruction22::ConstructGuide()
{
  /////////////////////////////////////////////////////////////////
  // Light Guide Part 1 In
  /////////////////////////////////////////////////////////////////
  G4double guide1InZ1 = 80.*mm/2; // distance between top face and TPC edge
  G4double guide1InZ2 = 2.*fCrossX*sin(fAngleWithTPCedge)/2; // distance between bottom face and TPC edge
  G4double guide1InZ3 = 2.*fCrossX*sin(fAngleWithTPCedge/2)/2;

  // l is leftward distance, r is rightward distance
  // top face 
  G4double guide1In2Xl = 70.*mm/2;
  G4double guide1In2Xr = 70.*mm/2;
  G4double guide1In2Yl = 26.21563480929*mm/2; // = 16.507901*mm/2
  G4double guide1In2Yr = 0.*mm;

  // cross-section cut by TPC edge
  // G4double guide1In3Xl = fDistantFromCenterPMT - fCrossX;

  // bottop face
  G4double guide1In1Xl = guide1In2Xl;
  G4double guide1In1Xr = 2*fCrossX*cos(fAngleWithTPCedge) - guide1In2Xl;
  G4double guide1In1Yl = (guide1InZ1 + guide1InZ2 + guide1InZ3)*(45.*mm-guide1In2Yl)/(guide1InZ1 + guide1InZ3) + guide1In2Yl;
  G4double guide1In1Yr = fCrossY;


  std::vector<G4TwoVector> guide1InPtVec(8);
  guide1InPtVec[0] = G4TwoVector(-guide1In1Xr, fBoxThickness);
  guide1InPtVec[1] = G4TwoVector(-guide1In1Xr, fBoxThickness + 2*guide1In1Yr);
  guide1InPtVec[2] = G4TwoVector(guide1In1Xl, fBoxThickness + 2*guide1In1Yl);
  guide1InPtVec[3] = G4TwoVector(guide1In1Xl, fBoxThickness);
  guide1InPtVec[4] = G4TwoVector(-guide1In2Xr, fBoxThickness);
  guide1InPtVec[5] = G4TwoVector(-guide1In2Xr, fBoxThickness + 2*guide1In2Yr);
  guide1InPtVec[6] = G4TwoVector(guide1In2Xl, fBoxThickness + 2*guide1In2Yl);
  guide1InPtVec[7] = G4TwoVector(guide1In2Xl, fBoxThickness);

  G4VSolid *solidGuide1In1 = new G4GenericTrap("Guide1In1", guide1InZ1 + guide1InZ2 + guide1InZ3, guide1InPtVec);
  G4double guide1InBoxX = guide1In2Xl, guide1InBoxY = 0.51*fBoxThickness, guide1InBoxZ = guide1InZ1;
  G4Box *solidGuide1In2 = new G4Box("Guide1In2", guide1InBoxX, guide1InBoxY, guide1InBoxZ);
  G4UnionSolid *solidGuide1In = new G4UnionSolid("Guide1In", solidGuide1In1, solidGuide1In2, 0, G4ThreeVector(0*mm, 0.5*fBoxThickness, guide1InZ2 + guide1InZ3));
  fGuideInPosVec =
    new G4ThreeVector(G4ThreeVector(0.*mm, 0.*mm, guide1InZ1 + guide1InZ3 - guide1InZ2) - G4ThreeVector(guide1In2Xl - cos(fAngleWithTPCedge)*fCrossX, 0.*mm, -sin(fAngleWithTPCedge)*fCrossX));
  fGuideInPosVec->rotateY(fAngleWithTPCedge);
  /////////////////////////////////////////////////////////////////
  // Light Guide Part 1 Out
  /////////////////////////////////////////////////////////////////  
  G4double guide1OutZ = guide1InZ1 + guide1InZ2 + guide1InZ3 + fBoxThickness;

  // l is leftward distance, r is rightward distance
  // top face 
  G4double guide1Out2Xl = guide1In2Xl + fBoxThickness;
  G4double guide1Out2Xr = guide1In2Xr + fBoxThickness;
  G4double guide1Out2Yl = guide1In2Yl + fBoxThickness;
  G4double guide1Out2Yr = guide1In2Yr + fBoxThickness;
  // cross-section cut by TPC edge

  // bottop face
  G4double guide1Out1Xl = guide1In1Xl + 0.9*fBoxThickness;// the vertex of box of guide(GuideCase) are adjusted to meet that of box of holder(HolderCase)
  G4double guide1Out1Xr = guide1In1Xr + 1.13*fBoxThickness;
  G4double guide1Out1Yl = guide1In1Yl + 1.35*fBoxThickness;
  G4double guide1Out1Yr = guide1In1Yr + 1.2*fBoxThickness;

  std::vector<G4TwoVector> guide1OutPtVec(8);
  guide1OutPtVec[0] = G4TwoVector(-guide1Out1Xr, 0.*mm);
  guide1OutPtVec[1] = G4TwoVector(-guide1Out1Xr, 2*guide1Out1Yr);
  guide1OutPtVec[2] = G4TwoVector(guide1Out1Xl, 2*guide1Out1Yl);
  guide1OutPtVec[3] = G4TwoVector(guide1Out1Xl, 0.*mm);
  guide1OutPtVec[4] = G4TwoVector(-guide1Out2Xr, 0.*mm);
  guide1OutPtVec[5] = G4TwoVector(-guide1Out2Xr, 2*guide1Out2Yr);
  guide1OutPtVec[6] = G4TwoVector(guide1Out2Xl, 2*guide1Out2Yl);
  guide1OutPtVec[7] = G4TwoVector(guide1Out2Xl, 0.*mm);

  G4VSolid *solidGuide1Out = new G4GenericTrap("Guide1Out", guide1OutZ, guide1OutPtVec);
  fGuideOutPosVec = new G4ThreeVector(*fGuideInPosVec);
  
  /////////////////////////////////////////////////////////////////
  // Light Guide Part 1 Reflector
  /////////////////////////////////////////////////////////////////  
  G4double guideReflector1OutZ = guide1InZ1 + guide1InZ2 + guide1InZ3 + fReflectorThickness;
  // top face
  std::vector<G4TwoVector> guideReflector1OutPtVec(8);
  guideReflector1OutPtVec[0] = guide1InPtVec[0] + G4TwoVector(-1.2*fReflectorThickness, -fReflectorThickness);
  guideReflector1OutPtVec[1] = guide1InPtVec[1] + G4TwoVector(-1.2*fReflectorThickness, 1.3*fReflectorThickness); // the vertex of reflector of guide(GuideReflector) are adjusted to meet that of reflector of holder(HolderReflector)
  guideReflector1OutPtVec[2] = guide1InPtVec[2] + G4TwoVector(0.9*fReflectorThickness, 1.8*fReflectorThickness);
  guideReflector1OutPtVec[3] = guide1InPtVec[3] + G4TwoVector(0.9*fReflectorThickness, -fReflectorThickness);
  guideReflector1OutPtVec[4] = guide1InPtVec[4] + G4TwoVector(-fReflectorThickness, -fReflectorThickness);
  guideReflector1OutPtVec[5] = guide1InPtVec[5] + G4TwoVector(-fReflectorThickness, fReflectorThickness);
  guideReflector1OutPtVec[6] = guide1InPtVec[6] + G4TwoVector(fReflectorThickness, fReflectorThickness);
  guideReflector1OutPtVec[7] = guide1InPtVec[7] + G4TwoVector(fReflectorThickness, -fReflectorThickness);

  // to make the centers of masses locate at Origin
  G4VSolid *solidGuideReflector1Out1 = new G4GenericTrap("GuideReflector1Out1", guideReflector1OutZ, guideReflector1OutPtVec);
  fGuideReflectorPosVec = new G4ThreeVector(*fGuideInPosVec);
  G4VSolid *solidGuideReflector1Out2 = new G4Box("GuideReflector1Out2", guide1InBoxX + fReflectorThickness, 0.501*fBoxThickness, guide1InBoxZ + fReflectorThickness);
  G4UnionSolid *solidGuideReflector1Out = new G4UnionSolid("GuideReflector1Out", solidGuideReflector1Out1, solidGuideReflector1Out2, 0, G4ThreeVector(0*mm, 0.5*(fBoxThickness - fReflectorThickness), guide1InZ2 + guide1InZ3));
// 180.53610900, 201.72577092
  // to cut the light guide volume to fit in TPC edge
  std::vector <G4TwoVector> tpcEdgeCutPolygon(3);
  G4double tpcEdgeCutXr = 2.*guide1Out1Xr, tpcEdgeCutXl = 2.*guide1Out1Xl, tpcEdgeCutY = (tpcEdgeCutXl + tpcEdgeCutXr)*tan(fAngleWithTPCedge), tpcEdgeCutZ = 3.*fCrossY;
  tpcEdgeCutPolygon[0] = G4TwoVector(-tpcEdgeCutXr, 0.*mm);
  tpcEdgeCutPolygon[1] = G4TwoVector(tpcEdgeCutXl,0.*mm);
  tpcEdgeCutPolygon[2] = G4TwoVector(tpcEdgeCutXl, -tpcEdgeCutY);

  G4ExtrudedSolid *solidTpcEdgeCut = new G4ExtrudedSolid("SolidTpcEdgeCut", tpcEdgeCutPolygon, tpcEdgeCutZ,
    G4TwoVector(), 1., G4TwoVector(), 1.);
  G4RotationMatrix *rotationTpcEdgeCut = new G4RotationMatrix();
  rotationTpcEdgeCut->rotateX(90*deg);
  G4ThreeVector tpcEdgeCutPosVec(0.*mm, 0.*mm, guide1InZ2 - guide1InZ3 - guide1InZ1 - tan(fAngleWithTPCedge)*(guide1In2Xl + tpcEdgeCutXr));

  /////////////////////////////////////////////////////////////////
  // Light Guide Part2(tilt) In
  /////////////////////////////////////////////////////////////////
  // sakak-kidung
  G4double tiltInX = guide1InBoxX, tiltInY = tiltInX*tan(fAngleTilt), tiltInZ = guide1InBoxZ;
  std::vector <G4TwoVector> tiltInPolygon(3);
  tiltInPolygon[0] = G4TwoVector(-2.*tiltInX - fBoxThickness, -2.*tiltInY);
  tiltInPolygon[1] = G4TwoVector(- fBoxThickness, 0.*mm);
  tiltInPolygon[2] = G4TwoVector(-2.*tiltInX - fBoxThickness);
  G4ExtrudedSolid *solidTiltIn = new G4ExtrudedSolid("TiltIn", tiltInPolygon, tiltInZ, 
    G4TwoVector(), 1., G4TwoVector(), 1.);

  /////////////////////////////////////////////////////////////////
  // Light Guide Part2(tilt) Out
  /////////////////////////////////////////////////////////////////
  // samkak-kidung
  std::vector <G4TwoVector> tiltOutPolygon(3);
  G4double tiltOutX = tiltInX + fBoxThickness, tiltOutY = tiltOutX*tan(fAngleTilt), tiltOutZ = tiltInZ + fBoxThickness;
  G4ThreeVector tiltOutPosVec(tiltOutX, 0*mm, guide1InZ2 + guide1InZ3);
  tiltOutPolygon[0] = G4TwoVector(-2*tiltOutX, -2.*tiltOutY);
  tiltOutPolygon[1] = G4TwoVector(0.*mm, 0.*mm);
  tiltOutPolygon[2] = G4TwoVector(-2*tiltOutX, 0.*mm);

  G4ExtrudedSolid *solidTiltOut = new G4ExtrudedSolid("TiltOut", tiltOutPolygon, tiltOutZ, 
    G4TwoVector(), 1., G4TwoVector(), 1.);

  /////////////////////////////////////////////////////////////////
  // Light Guide Part2(tilt) Reflector
  /////////////////////////////////////////////////////////////////
  // sakak-kidung
  std::vector <G4TwoVector> guideReflectorOutPolygon(3);
  guideReflectorOutPolygon[0] = tiltInPolygon[0] + G4TwoVector(-fReflectorThickness, -2.*fReflectorThickness*tan(fAngleTilt));
  guideReflectorOutPolygon[1] = tiltInPolygon[1] + G4TwoVector(fReflectorThickness, 0.*mm);
  guideReflectorOutPolygon[2] = tiltInPolygon[2] + G4TwoVector(-fReflectorThickness, 0.*mm);
  G4ExtrudedSolid *solidGuideReflector2Out = new G4ExtrudedSolid("GuideReflector2Out", guideReflectorOutPolygon, tiltInZ + fReflectorThickness, 
    G4TwoVector(), 1., G4TwoVector(), 1.);

  ///////////////////////////////////////////-//////////////////////
  // Tubs : PMT hole to sutract with box and (tiltInHole = In =  PMT hole)
  /////////////////////////////////////////////////////////////////
  G4double tiltHoleZ = fBoxThickness*sin(fAngleTilt);
  G4ThreeVector tiltHolePosVec(-tiltOutX - tiltHoleZ*sin(fAngleTilt), -tiltOutY + tiltHoleZ*cos(fAngleTilt), 0.*mm);

  G4RotationMatrix *rotationTiltHole = new G4RotationMatrix();
  rotationTiltHole->rotate(-90.*deg, G4ThreeVector(cos(fAngleTilt), sin(fAngleTilt)));

  G4Tubs *solidTiltHole = new G4Tubs("TiltHole", 0.*mm, fCathodeRadius, tiltHoleZ, 0.*deg, 360.*deg);
  G4UnionSolid  *solidTiltInHole = new G4UnionSolid("TiltInHole", solidTiltIn, solidTiltHole, rotationTiltHole, tiltHolePosVec);

  G4UnionSolid *solidGuideOutBeforeCut = new G4UnionSolid("GuideOutBeforeCut", solidGuide1Out, solidTiltOut, 0, tiltOutPosVec);
  G4UnionSolid *solidGuideInBeforeCut = new G4UnionSolid("GuideInBeforeCut", solidGuide1In, solidTiltInHole, 0, tiltOutPosVec);
  G4UnionSolid *solidGuideReflectorOutBeforeCut = new G4UnionSolid("GuideReflectorOutWithoutHole", solidGuideReflector1Out, solidGuideReflector2Out, 0, tiltOutPosVec);
  
  G4SubtractionSolid *solidGuideCaseWithoutHole = new G4SubtractionSolid("GuideCaseWithoutHole", solidGuideOutBeforeCut, solidGuideReflectorOutBeforeCut, 0, G4ThreeVector());
  G4SubtractionSolid *solidGuideCaseBeforeCut = new G4SubtractionSolid("GuideCaseBeforeCut", solidGuideCaseWithoutHole, solidGuideInBeforeCut, 0, *fGuideInPosVec - *fGuideOutPosVec);
  G4SubtractionSolid *solidGuideReflectorBeforeCut = new G4SubtractionSolid("GuideReflectorBeforeCut", solidGuideReflectorOutBeforeCut, solidGuideInBeforeCut, 0, *fGuideInPosVec - *fGuideReflectorPosVec);

  fSolidGuideOut  = new G4SubtractionSolid("GuideOut", solidGuideOutBeforeCut, solidTpcEdgeCut, rotationTpcEdgeCut, tpcEdgeCutPosVec);
  fSolidGuideIn   = new G4SubtractionSolid("GuideIn", solidGuideInBeforeCut, solidTpcEdgeCut, rotationTpcEdgeCut, tpcEdgeCutPosVec);
  fSolidGuideCase = new G4SubtractionSolid("GuideCase", solidGuideCaseBeforeCut, solidTpcEdgeCut, rotationTpcEdgeCut, tpcEdgeCutPosVec);
  fSolidGuideReflector = new G4SubtractionSolid("GuideReflector", solidGuideReflectorBeforeCut, solidTpcEdgeCut, rotationTpcEdgeCut, tpcEdgeCutPosVec);

  fLogicGuideCase = new G4LogicalVolume(fSolidGuideCase, fMatBox, "GuideCase");
  fLogicGuideReflector = new G4LogicalVolume(fSolidGuideReflector, fMatReflector, "GuideReflector");

  fRotationGuide = new G4RotationMatrix();
  fRotationGuide->rotateY(-fAngleWithTPCedge);
  fPhysGuideCase  = new G4PVPlacement(fRotationGuide, *fGuideOutPosVec,  fLogicGuideCase,  "GuideCase",  fLogicWorld, false, 0, fCheckOverlaps);
  fPhysGuideReflector = new G4PVPlacement(fRotationGuide, *fGuideReflectorPosVec, fLogicGuideReflector, "GuideReflector", fLogicWorld, false, 0, fCheckOverlaps);

  fRotationScorer = new G4RotationMatrix();
  fRotationScorer->rotateY(-fAngleWithTPCedge);
  fRotationScorer->rotate(90.*deg, G4ThreeVector(cos(fAngleTilt), sin(fAngleTilt)));

  fScorerPosVec = 
    new G4ThreeVector(G4ThreeVector(0.*mm, -tiltOutY, guide1InZ2 + guide1InZ3) + fScorerZ*G4ThreeVector(sin(fAngleTilt), -cos(fAngleTilt), 0.*mm));
  fScorerPosVec->rotateY(fAngleWithTPCedge);
  *fScorerPosVec = *fScorerPosVec + *fGuideOutPosVec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction22::ConstructScorer()
{
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
  fPhysScorer     = new G4PVPlacement(fRotationScorer, *fScorerPosVec, fLogicScorer, "Scorer", fLogicWorld, false, 0, fCheckOverlaps);
  fPhysGlass      = new G4PVPlacement(0,    G4ThreeVector(0.*mm, 0.*mm, fScorerZ - fGlassThickness) ,    fLogicGlass,    "Glass", fLogicScorer, false, 0, fCheckOverlaps);
  fPhysCathode    = new G4PVPlacement(0,    G4ThreeVector(0.*mm, 0.*mm, fScorerZ - 2*fGlassThickness - fCathodeZ)  ,  fLogicCathode,  "Cathode", fLogicScorer, false, 0, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
