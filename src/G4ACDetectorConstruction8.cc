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
// $Id: G4ACDetectorConstruction8.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACDetectorConstruction8.cc
/// \brief Implementation of the G4ACDetectorConstruction8 class

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

#include "G4ACDetectorConstruction8.hh"
#include "G4ACFineMeshSD.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction8::G4ACDetectorConstruction8()
: G4ACBaseDetectorConstruction(100.*mm/2, 100.*mm/2, 20.*mm/2, 0.05*mm, 30.*deg),
  fAerogelPosVec4(0), fAerogelPosVec5(0), fAerogelPosVec6(0),
  fSolidAerogel4(0), fSolidAerogel5(0), fSolidAerogel6(0),
  fLogicAerogel4(0), fLogicAerogel5(0), fLogicAerogel6(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction8::~G4ACDetectorConstruction8()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACDetectorConstruction8::Construct()
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

  // ConstructAerogel();

  ConstructHolderBoundary(); // construct boundary processes
  ConstructGuideBoundary();
  ConstructAerogelBoundary();
  ConstructGlassBoundary();
  
  return fPhysWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction8::ConstructHolder()
{
  /////////////////////////////////////////////////////////////////
  // Aerogel Holder part 1 (inside airs)
  /////////////////////////////////////////////////////////////////  
  // bottom face
  G4double holderIn1X = 100.*mm/2;
  G4double holderIn1Y = 60.*mm/2;
  
  // top face
  G4double holderIn2X = fCrossX;
  G4double holderIn2Y = fCrossY;
  
  G4double holderInZ = 200*mm/2;

  std::vector<G4TwoVector> holderInPtVec(8);

  holderInPtVec[0] = G4TwoVector(-holderIn1X,  fBoxThickness);
  holderInPtVec[1] = G4TwoVector(-holderIn1X,  fBoxThickness + 2*holderIn1Y);
  holderInPtVec[2] = G4TwoVector(holderIn1X,  fBoxThickness + 2*holderIn1Y);
  holderInPtVec[3] = G4TwoVector(holderIn1X, fBoxThickness);
  holderInPtVec[4] = G4TwoVector(-holderIn2X,  fBoxThickness);
  holderInPtVec[5] = G4TwoVector(-holderIn2X,  fBoxThickness + 2*holderIn2Y);  
  holderInPtVec[6] = G4TwoVector(holderIn2X,  fBoxThickness + 2*holderIn2Y);
  holderInPtVec[7] = G4TwoVector(holderIn2X, fBoxThickness);

  /////////////////////////////////////////////////////////////////
  // Aerogel Holder part 2(volume including case)
  /////////////////////////////////////////////////////////////////
  G4double holderOutZ = holderInZ + fBoxThickness/2;
  
  // bottom face
  G4double segmentAngle = 5.*deg; // total 30 segments
  G4double holderOut1X = holderIn1X + fBoxThickness*(1 - sin(segmentAngle));
  G4double holderOut1Y = holderIn1Y + fBoxThickness*(3 - sin(segmentAngle))/2;
  // top face
  G4double holderOut2X = holderIn2X + fBoxThickness, holderOut2Y = holderIn2Y + fBoxThickness;

  std::vector<G4TwoVector> holderOutPtVec(8);

  holderOutPtVec[0] = G4TwoVector(-holderOut1X, 0.*mm);
  holderOutPtVec[1] = G4TwoVector(-holderOut1X, 2*holderOut1Y);
  holderOutPtVec[2] = G4TwoVector(holderOut1X, 2*holderOut1Y);
  holderOutPtVec[3] = G4TwoVector(holderOut1X, 0.*mm);
  
  holderOutPtVec[4] = G4TwoVector(-holderOut2X, 0.*mm);
  holderOutPtVec[5] = G4TwoVector(-holderOut2X, 2*holderOut2Y);
  holderOutPtVec[6] = G4TwoVector(holderOut2X, 2*holderOut2Y);
  holderOutPtVec[7] = G4TwoVector(holderOut2X, 0.*mm);


  // position vector of the centers of masses(center of mass)
  fSolidHolderIn = new G4GenericTrap("HolderIn", holderInZ, holderInPtVec);
  fSolidHolderOut = new G4GenericTrap("HolderOut", holderOutZ, holderOutPtVec);
  
  fHolderInPosVec = new G4ThreeVector(0.*mm, 0.*mm, -holderInZ);
  fHolderOutPosVec = new G4ThreeVector(0.*mm, 0.*mm, -holderOutZ);
  
  fSolidHolderCase  = new G4SubtractionSolid("HolderCase", fSolidHolderOut, fSolidHolderIn, 0, *fHolderInPosVec - *fHolderOutPosVec);
  fLogicHolderCase = new G4LogicalVolume(fSolidHolderCase, fMatBox, "HolderCase");
  fPhysHolderCase = new G4PVPlacement(0,    *fHolderOutPosVec,  fLogicHolderCase,  "HolderCase",  fLogicWorld, false, 0, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction8::ConstructGuide()
{
  /////////////////////////////////////////////////////////////////
  // Light Guide1 (inside air)
  /////////////////////////////////////////////////////////////////
  G4double guide1InZ = 100.*mm/2;
  // bottop face
  G4double guide1In1X = fCrossX;
  G4double guide1In1Y = fCrossY;

  // top fae
  G4double guide1In2X = 80.*mm/2;
  G4double guide1In2Y = 70.*mm/2;

  std::vector<G4TwoVector> guide1InPtVec(8);
  guide1InPtVec[0] = G4TwoVector(-guide1In1X, fBoxThickness);
  guide1InPtVec[1] = G4TwoVector(-guide1In1X, fBoxThickness + 2*guide1In1Y);
  guide1InPtVec[2] = G4TwoVector(guide1In1X, fBoxThickness + 2*guide1In1Y);
  guide1InPtVec[3] = G4TwoVector(guide1In1X, fBoxThickness);
  guide1InPtVec[4] = G4TwoVector(-guide1In2X, fBoxThickness);
  guide1InPtVec[5] = G4TwoVector(-guide1In2X, fBoxThickness + 2*guide1In2Y);
  guide1InPtVec[6] = G4TwoVector(guide1In2X, fBoxThickness + 2*guide1In2Y);
  guide1InPtVec[7] = G4TwoVector(guide1In2X, fBoxThickness);


  // to make the centers of masses locate at Origin
  G4VSolid *solidGuide1In = new G4GenericTrap("Guide1In", guide1InZ, guide1InPtVec);
  fGuideInPosVec = new G4ThreeVector(0.*mm, 0.*mm, guide1InZ);

  /////////////////////////////////////////////////////////////////
  // Light Guide Outside part 1
  /////////////////////////////////////////////////////////////////  
  G4double guide1OutZ = guide1InZ;

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
  // Light Guide  Outside (tub) part2
  /////////////////////////////////////////////////////////////////
  G4RotationMatrix *rotationGuide2 = new G4RotationMatrix();
  rotationGuide2->rotateY(90.*deg);
  G4ThreeVector guide2PosVec(0.*mm, 0.*mm, 2.*guide1OutZ);

  G4double guide2OutRmin = 0, guide2OutRmax = 2.*guide1Out2Y;
  G4double guide2OutZ = guide1Out2X;
  G4double guide2OutSphi = 0.*deg, guide2OutDphi = 90.*deg;

  G4Tubs *solidGuide2Out = new G4Tubs("Guide2Out", guide2OutRmin, guide2OutRmax, guide2OutZ, guide2OutSphi, guide2OutDphi);
  
  /////////////////////////////////////////////////////////////////
  // Light Guide  Inside (tub) part2
  /////////////////////////////////////////////////////////////////
  G4double guide2InRmin = 0.*mm, guide2InRmax = guide2OutRmax-fBoxThickness;
  G4double guide2InZ = guide2OutZ - fBoxThickness;
  G4double guide2InSphi = 0.*deg, guide2InDphi = 90.*deg;

  G4Tubs *solidGuide2In = new G4Tubs("Guide2In", guide2InRmin, guide2InRmax, guide2InZ, guide2InSphi, guide2InDphi);

  /////////////////////////////////////////////////////////////////
  // Tilt board Out
  /////////////////////////////////////////////////////////////////
  // samkak-kidung

  std::vector <G4TwoVector> tiltOutPolygon(3);
  G4double tiltOutX = guide1Out2Y, tiltOutY = tiltOutX*tan(fAngleTilt), tiltOutZ = guide1Out2X;
  G4ThreeVector tiltOutPosVec(0.*mm, 0.*mm, 2.*guide1OutZ);
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
  // Tilt board In
  /////////////////////////////////////////////////////////////////
  // sakak-kidung
  std::vector <G4TwoVector> tiltInPolygon(3);
  tiltInPolygon[0] = G4TwoVector(-2.*tiltOutX + fBoxThickness, 0.*mm);
  tiltInPolygon[1] = G4TwoVector(-fBoxThickness, 0.*mm);
  tiltInPolygon[2] = G4TwoVector(-fBoxThickness, 2.*(-tiltOutX + fBoxThickness)*tan(fAngleTilt));
  G4ExtrudedSolid *solidTiltIn = new G4ExtrudedSolid("TiltIn", tiltInPolygon, tiltOutZ - fBoxThickness, 
    G4TwoVector(), 1., G4TwoVector(), 1.);

  ///////////////////////////////////////////W//////////////////////
  // Tubs : PMT hole to sutract with box and (tiltInHole = In =  PMT hole)
  /////////////////////////////////////////////////////////////////
  G4double tiltHoleZ = 0.5*fBoxThickness*sin(fAngleTilt);
  G4ThreeVector tiltHolePosVec(-tiltOutX + tiltHoleZ*sin(fAngleTilt), -tiltOutY + tiltHoleZ*cos(fAngleTilt), 0.*mm);
  fRotationScorer = new G4RotationMatrix();
  fRotationScorer->rotate(-90.*deg, G4ThreeVector(cos(180.*deg - fAngleTilt), sin(180.*deg - fAngleTilt)));

  G4Tubs *solidTiltHole = new G4Tubs("TiltHole", 0.*mm, fCathodeRadius, tiltHoleZ, 0.*deg, 360.*deg);
  G4UnionSolid  *solidTiltInHole = new G4UnionSolid("TiltInHole", solidTiltIn, solidTiltHole, fRotationScorer, tiltHolePosVec);

  G4UnionSolid *solidGuideTubOut     = new G4UnionSolid("GuideTubOut", solidGuide1Out, solidGuide2Out, rotationGuide2, guide2PosVec - *fGuideOutPosVec);
  G4UnionSolid *solidGuideTubIn      = new G4UnionSolid("GuideTubIn", solidGuide1In, solidGuide2In, rotationGuide2, guide2PosVec - *fGuideInPosVec);
  
  fSolidGuideOut  = new G4UnionSolid("GuideOut", solidGuideTubOut, solidTiltOut, rotationTilt, tiltOutPosVec - *fGuideOutPosVec);
  fSolidGuideIn   = new G4UnionSolid("GuideIn", solidGuideTubIn, solidTiltInHole, rotationTilt, tiltOutPosVec- *fGuideInPosVec);
  fSolidGuideCase = new G4SubtractionSolid("GuideCase", fSolidGuideOut, fSolidGuideIn, 0, *fGuideInPosVec - *fGuideOutPosVec);
  fLogicGuideCase = new G4LogicalVolume(fSolidGuideCase, fMatBox, "GuideCase");
  fPhysGuideCase  = new G4PVPlacement(0,    *fGuideOutPosVec,  fLogicGuideCase,  "GuideCase",  fLogicWorld, false, 0, fCheckOverlaps);
 
  fScorerPosVec = new G4ThreeVector(tiltOutPosVec + G4ThreeVector(0.*mm, -tiltOutY, tiltOutX) + fScorerZ*G4ThreeVector(0.*mm, -cos(fAngleTilt), sin(fAngleTilt)));
}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4ACDetectorConstruction8::ConstructScorer()
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
void G4ACDetectorConstruction8::ConstructAerogel()
{

  const G4double aerogelX = 92.*mm/2;
  const G4double aerogelY = fAerogelThickness;
  const G4double aerogelZ = 92.*mm/2;
  
  fSolidAerogel1 = new G4Box("Aerogel1", aerogelX - fAerogelGap, aerogelY - fAerogelGap, aerogelZ - fAerogelGap);
  fSolidAerogel2 = new G4Box("Aerogel2", aerogelX - fAerogelGap, aerogelY - fAerogelGap, aerogelZ - fAerogelGap);
  fSolidAerogel3 = new G4Box("Aerogel3", aerogelX - fAerogelGap, aerogelY - fAerogelGap, aerogelZ - fAerogelGap);
  fSolidAerogel4 = new G4Box("Aerogel4", aerogelX - fAerogelGap, aerogelY - fAerogelGap, aerogelZ - fAerogelGap);
  fSolidAerogel5 = new G4Box("Aerogel5", aerogelX - fAerogelGap, aerogelY - fAerogelGap, aerogelZ - fAerogelGap);
  fSolidAerogel6 = new G4Box("Aerogel6", aerogelX - fAerogelGap, aerogelY - fAerogelGap, aerogelZ - fAerogelGap);

  fLogicAerogel1 = new G4LogicalVolume(fSolidAerogel1, fMatAerogel, "Aerogel1");
  fLogicAerogel2 = new G4LogicalVolume(fSolidAerogel2, fMatAerogel, "Aerogel2");
  fLogicAerogel3 = new G4LogicalVolume(fSolidAerogel3, fMatAerogel, "Aerogel3");
  fLogicAerogel4 = new G4LogicalVolume(fSolidAerogel4, fMatAerogel, "Aerogel4");
  fLogicAerogel5 = new G4LogicalVolume(fSolidAerogel5, fMatAerogel, "Aerogel5");
  fLogicAerogel6 = new G4LogicalVolume(fSolidAerogel6, fMatAerogel, "Aerogel6");
  
  fAerogelPosVec1 = new G4ThreeVector(0.*mm, fBoxThickness + aerogelY, fHolderInPosVec->z() + aerogelZ);
  fAerogelPosVec2 = new G4ThreeVector(0.*mm, fBoxThickness + aerogelY, fHolderInPosVec->z() + 3*aerogelZ);
  fAerogelPosVec3 = new G4ThreeVector(0.*mm, fBoxThickness + 3*aerogelY, fHolderInPosVec->z() + aerogelZ);
  fAerogelPosVec4 = new G4ThreeVector(0.*mm, fBoxThickness + 3*aerogelY, fHolderInPosVec->z() + 3*aerogelZ);
  fAerogelPosVec5 = new G4ThreeVector(0.*mm, fBoxThickness + 5*aerogelY, fHolderInPosVec->z() + aerogelZ);
  fAerogelPosVec6 = new G4ThreeVector(0.*mm, fBoxThickness + 5*aerogelY, fHolderInPosVec->z() + 3*aerogelZ);

  new G4PVPlacement(0, *fAerogelPosVec1, fLogicAerogel1, "Aerogel1", fLogicIn, true, 0, fCheckOverlaps);
  new G4PVPlacement(0, *fAerogelPosVec2, fLogicAerogel1, "Aerogel2", fLogicIn, true, 1, fCheckOverlaps);
  new G4PVPlacement(0, *fAerogelPosVec3, fLogicAerogel1, "Aerogel3", fLogicIn, true, 2, fCheckOverlaps);
  new G4PVPlacement(0, *fAerogelPosVec4, fLogicAerogel1, "Aerogel4", fLogicIn, true, 3, fCheckOverlaps);
  new G4PVPlacement(0, *fAerogelPosVec5, fLogicAerogel1, "Aerogel5", fLogicIn, true, 4, fCheckOverlaps);
  new G4PVPlacement(0, *fAerogelPosVec6, fLogicAerogel1, "Aerogel6", fLogicIn, true, 5, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction8::ConstructAerogelBoundary()
{
  G4OpticalSurface *opAerogel1AirSurface = new G4OpticalSurface("Aerogel1AirSurface");
  opAerogel1AirSurface->SetType(dielectric_dielectric);
  opAerogel1AirSurface->SetFinish(groundair);
  opAerogel1AirSurface->SetModel(unified);
}