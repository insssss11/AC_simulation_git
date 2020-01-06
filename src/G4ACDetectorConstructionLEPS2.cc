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
// $Id: G4ACDetectorConstructionLEPS2.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACDetectorConstructionLEPS2.cc
/// \brief Definition of the G4ACDetectorConstructionLEPS2 class

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

#include "G4ACDetectorConstructionLEPS2.hh"
#include "G4ACAerogelSD.hh"

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern const G4double TARGET_LENGTH, TARGET_CENTER_POSITION,
  AC_DETECTOR_R, TPC_EDGE_R, TPC_DOWNSTREAM_Z;

G4ACDetectorConstructionLEPS2::G4ACDetectorConstructionLEPS2()
: G4ACBaseDetectorConstruction(0.1*mm),
  fTPCradius(TPC_EDGE_R), fOffsetAngle(-24.*deg), fOffsetTPCrocation(TPC_DOWNSTREAM_Z), fNelem(30),
  fSolidGuideCaseA(0), fSolidGuideCaseB(0), fSolidGuideCaseC(0),
  fSolidGuideInA(0), fSolidGuideInB(0), fSolidGuideInC(0),
  fSolidGuideOutA(0), fSolidGuideOutB(0), fSolidGuideOutC(0)
{
  fRotationElement = new G4RotationMatrix*[fNelem];
  fPhysHolderCaseArr = new G4PVPlacement*[fNelem];
  fPhysGuideCaseArr = new G4PVPlacement*[fNelem];
  fPhysInArr = new G4PVPlacement*[fNelem];
  fPhysScorerArr = new G4PVPlacement*[fNelem];
  fPhysGlassArr = new G4PVPlacement*[fNelem];
  fPhysCathodeArr = new G4PVPlacement*[fNelem];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstructionLEPS2::~G4ACDetectorConstructionLEPS2()
{
  delete[] fPhysHolderCaseArr;
  delete[] fPhysGuideCaseArr;
  delete[] fPhysInArr;
  delete[] fPhysScorerArr;
  delete[] fPhysGlassArr;
  delete[] fPhysCathodeArr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACDetectorConstructionLEPS2::Construct()
{
  //.......oooOO0OOooo........oooOO0OOooo......
  

  /////////////////////////////////////////////////////////////////
  // world
  /////////////////////////////////////////////////////////////////
  G4double worldHalfSizeX = 1000.*mm;
  G4double worldHalfSizeY = 1000.*mm;
  G4double worldHalfSizeZ = 1000.*mm;
  fSolidWorld   = new G4Box("World", worldHalfSizeX, worldHalfSizeY, worldHalfSizeZ);
  fLogicWorld   = new G4LogicalVolume(fSolidWorld, fMatAir, "World");
  fPhysWorld    = new G4PVPlacement(0,    G4ThreeVector(),  fLogicWorld,  "World",  0, false, 0, fCheckOverlaps);

  // rotation
  for(G4int copyNum = 0;copyNum < fNelem;copyNum++)
  {
    fRotationElement[copyNum] = new G4RotationMatrix();
    fRotationElement[copyNum]->rotateZ(-(fOffsetAngle + copyNum*360.*deg/fNelem));
    fRotationElement[copyNum]->rotateX(-90.*deg);
  }

  /////////////////////////////////////////////////////////////////
  // construct Aerogel Holder and Light Guide, Scorer, Aerogel
  /////////////////////////////////////////////////////////////////
  ConstructHolders();
  /*
  ConstructGuides();
  ConstructScorers();
  */
  ConstructInnerAirs();
  ConstructAerogels();
  /*
  ConstructHolderBoundaries(); // construct boundary processes
  ConstructGuideBoundaries();
  ConstructAerogelBoundaries();
  ConstructGlassBoundaries();   
  */
  return fPhysWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstructionLEPS2::ConstructHolders()
{
  /////////////////////////////////////////////////////////////////
  // Aerogel Holder part 1 (inside airs)
  /////////////////////////////////////////////////////////////////  
  // bottom face
  G4double holderIn1X = fAerogel11X;
  G4double holderIn1Y = 60.*mm/2;
  
  // top face
  G4double holderIn2X = fCrossX;
  G4double holderIn2Y = fCrossY;
  
  G4double holderInZ = fAerogel1Z;

  std::vector<G4TwoVector> holderInPtVec(8);

  holderInPtVec[0] = G4TwoVector(-holderIn1X,  fBoxThickness);
  holderInPtVec[1] = G4TwoVector(-holderIn1X,  fBoxThickness + 2*holderIn1Y);
  holderInPtVec[2] = G4TwoVector(holderIn1X,  fBoxThickness + 2*holderIn1Y);
  holderInPtVec[3] = G4TwoVector(holderIn1X, fBoxThickness);
  holderInPtVec[4] = G4TwoVector(-holderIn2X,  fBoxThickness);
  holderInPtVec[5] = G4TwoVector(-holderIn2X,  fBoxThickness + 2*holderIn2Y);  
  holderInPtVec[6] = G4TwoVector(holderIn2X,  fBoxThickness + 2*holderIn2Y);
  holderInPtVec[7] = G4TwoVector(holderIn2X, fBoxThickness);

  fSolidHolderIn = new G4GenericTrap("HolderIn", holderInZ, holderInPtVec);
  fHolderInPosVec = new G4ThreeVector(0.*mm, 0.*mm, -holderInZ);
  
  /////////////////////////////////////////////////////////////////
  // Aerogel Holder part 2(volume including case)
  /////////////////////////////////////////////////////////////////
  G4double holderOutZ = holderInZ + fBoxThickness/2;
  
  // bottom face
  G4double segmentAngle = 6.*deg; // total 30 segments
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
  fSolidHolderOut = new G4GenericTrap("HolderOut", holderOutZ, holderOutPtVec);
  fHolderOutPosVec = new G4ThreeVector(0.*mm, 0.*mm, -holderOutZ);
  fSolidHolderCase  = new G4SubtractionSolid("HolderCase", fSolidHolderOut, fSolidHolderIn, 0, *fHolderInPosVec - *fHolderOutPosVec);
  fLogicHolderCase = new G4LogicalVolume(fSolidHolderCase, fMatBox, "HolderCase");
  
  // locating AC1 elements in ring shape
  G4ThreeVector posVecInLabs;
  for(G4int copyNum = 0;copyNum < fNelem;copyNum++)
  {
    posVecInLabs = G4ThreeVector(0.*mm, holderOutZ - fTPCradius, fOffsetTPCrocation);
    posVecInLabs.rotateZ(fOffsetAngle + copyNum*360.*deg/fNelem);
    fPhysHolderCaseArr[copyNum] = new G4PVPlacement(fRotationElement[copyNum],    posVecInLabs,  fLogicHolderCase,  "HolderCase",  fLogicWorld, true, copyNum, fCheckOverlaps);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void G4ACDetectorConstructionLEPS2::ConstructGuides()
{
  ConstructGuideA();
  ConstructGuideB();
  ConstructGuideC();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstructionLEPS2::ConstructScorers()
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
*/
void G4ACDetectorConstructionLEPS2::ConstructInnerAirs()
{
  // fSolidIn        = new G4UnionSolid("In", fSolidHolderIn, fSolidGuideIn, 0, *fGuideInPosVec - *fHolderInPosVec);
  fLogicIn        = new G4LogicalVolume(fSolidHolderIn, fMatAir,"In");
  // locating AC1 elements in ring shape
  G4ThreeVector posVecInLabs;
  for(G4int copyNum = 0;copyNum < fNelem;copyNum++)
  {
    posVecInLabs = G4ThreeVector(0.*mm, -fHolderInPosVec->z() - fTPCradius, fOffsetTPCrocation);
    posVecInLabs.rotateZ(fOffsetAngle + copyNum*360.*deg/fNelem);
    fPhysIn         = new G4PVPlacement(fRotationElement[copyNum],    posVecInLabs, fLogicIn,   "In", fLogicWorld, true, copyNum, fCheckOverlaps);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstructionLEPS2::ConstructAerogels()
{
  /////////////////////////////////////////////////////////////////
  // aerogel1
  /////////////////////////////////////////////////////////////////
  // The # of aerogels are six pieces and there are all different refractive indice
  // The aerogels numbering and position
  // -------------------------------------------|
  // |  3              6       |                |------
  // |  2              5       |  light guide   |  pmt|
  // |  1              4       |                |------
  // -------------------------------------------|
  //  ^
  //  |
  //  | beam
  // 
  // also aerogel pieces are daughter of aerogel and located with distance between each others 
  
  // whole Aerogel region 
  // aerogel region is daughter of caseIn,
  // the position coordinate of center of aerogels 
  std::vector<G4TwoVector> aerogel1Pt(8);

  aerogel1Pt[0] = G4TwoVector(-fAerogel11X + fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel1Pt[1] = G4TwoVector(-fAerogel11X + fAerogelGap, fBoxThickness + 2*fAerogelThickness);
  aerogel1Pt[2] = G4TwoVector(fAerogel11X - fAerogelGap, fBoxThickness + 2*fAerogelThickness);
  aerogel1Pt[3] = G4TwoVector(fAerogel11X - fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel1Pt[4] = G4TwoVector(-fAerogel12X + fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel1Pt[5] = G4TwoVector(-fAerogel12X + fAerogelGap, fBoxThickness + 2*fAerogelThickness);
  aerogel1Pt[6] = G4TwoVector(fAerogel12X - fAerogelGap, fBoxThickness + 2*fAerogelThickness);
  aerogel1Pt[7] = G4TwoVector(fAerogel12X - fAerogelGap, fBoxThickness + fAerogelGap);
  
  /////////////////////////////////////////////////////////////////
  // aerogel2
  /////////////////////////////////////////////////////////////////
  std::vector<G4TwoVector> aerogel2Pt(8);

  aerogel2Pt[0] = G4TwoVector(-fAerogel21X + fAerogelGap, fBoxThickness + 2*fAerogelThickness);
  aerogel2Pt[1] = G4TwoVector(-fAerogel21X + fAerogelGap, fBoxThickness + 4*fAerogelThickness);
  aerogel2Pt[2] = G4TwoVector(fAerogel21X - fAerogelGap, fBoxThickness + 4*fAerogelThickness);
  aerogel2Pt[3] = G4TwoVector(fAerogel21X - fAerogelGap, fBoxThickness + 2*fAerogelThickness);
  aerogel2Pt[4] = G4TwoVector(-fAerogel22X + fAerogelGap, fBoxThickness + 2*fAerogelThickness);
  aerogel2Pt[5] = G4TwoVector(-fAerogel22X + fAerogelGap, fBoxThickness + 4*fAerogelThickness);
  aerogel2Pt[6] = G4TwoVector(fAerogel22X + fAerogelGap, fBoxThickness + 4*fAerogelThickness);
  aerogel2Pt[7] = G4TwoVector(fAerogel22X + fAerogelGap, fBoxThickness + 2*fAerogelThickness);
  
  /////////////////////////////////////////////////////////////////
  // aerogel3
  /////////////////////////////////////////////////////////////////
  std::vector<G4TwoVector> aerogel3Pt(8);

  aerogel3Pt[0] = G4TwoVector(-fAerogel31X + fAerogelGap, fBoxThickness + 4*fAerogelThickness);
  aerogel3Pt[1] = G4TwoVector(-fAerogel31X + fAerogelGap, fBoxThickness + 6*fAerogelThickness - fAerogelGap);
  aerogel3Pt[2] = G4TwoVector(fAerogel31X - fAerogelGap, fBoxThickness + 6*fAerogelThickness - fAerogelGap);
  aerogel3Pt[3] = G4TwoVector(fAerogel31X - fAerogelGap, fBoxThickness + 4*fAerogelThickness);
  aerogel3Pt[4] = G4TwoVector(-fAerogel32X + fAerogelGap, fBoxThickness + 4*fAerogelThickness);
  aerogel3Pt[5] = G4TwoVector(-fAerogel32X + fAerogelGap, fBoxThickness + 6*fAerogelThickness - fAerogelGap);
  aerogel3Pt[6] = G4TwoVector(fAerogel32X - fAerogelGap, fBoxThickness + 6*fAerogelThickness - fAerogelGap);
  aerogel3Pt[7] = G4TwoVector(fAerogel32X - fAerogelGap, fBoxThickness + 4*fAerogelThickness);

  fSolidAerogel1 = new G4GenericTrap("Aerogel1", fAerogel1Z - fAerogelGap, aerogel1Pt);
  fSolidAerogel2 = new G4GenericTrap("Aerogel2", fAerogel2Z - fAerogelGap, aerogel2Pt);
  fSolidAerogel3 = new G4GenericTrap("Aerogel3", fAerogel3Z - fAerogelGap, aerogel3Pt);
  
  G4VSolid *solidAerogel12 = new G4UnionSolid("Aerogel12", fSolidAerogel1, fSolidAerogel2, 0, G4ThreeVector(0.*mm, 0.*mm, fAerogel1Z - fAerogel2Z));
  G4VSolid *solidAerogel = new G4UnionSolid("Aerogel", solidAerogel12, fSolidAerogel3, 0, G4ThreeVector(0.*mm, 0.*mm, fAerogel1Z - fAerogel3Z));


  fLogicAerogel1 = new G4LogicalVolume(solidAerogel, fMatAerogel, "Aerogel");
  
  // locating AC1 elements in ring shape
  /*
  for(G4int copyNum = 0;copyNum < 1;copyNum++)
  {
    fPhysAerogel1Arr[copyNum] = new G4PVPlacement(0, G4ThreeVector(), fLogicAerogel1, "Aerogel",  fLogicIn, false, copyNum, fCheckOverlaps);
  }
  */
  fPhysAerogel1Arr = new G4PVPlacement(0, G4ThreeVector(), fLogicAerogel1, "Aerogel",  fLogicIn, false, 0, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void G4ACDetectorConstructionLEPS2::ConstructGuideA()
{
  /////////////////////////////////////////////////////////////////
  // Light Guide (inside air)
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


  G4double guide1InBoxX = guide1In2X, guide1InBoxY = 2*fBoxThickness/2, guide1InBoxZ = guide1In2X*cos(fAngleTilt);

  G4VSolid *solidGuide1In1 = new G4GenericTrap("Guide1In1A", guide1InZ, guide1InPtVec);
  G4VSolid *solidGuide1In2 = new G4Box("Guide1In2A", guide1InBoxX, guide1InBoxY, guide1InBoxZ);
  
  auto guide1In1PosVec = new G4ThreeVector(0.*mm, 0.*mm, guide1InZ);
  auto guide1In2PosVec = new G4ThreeVector(0.*mm, fBoxThickness, 2*guide1InZ - guide1InBoxZ);

  G4UnionSolid *solidGuide1In = new G4UnionSolid("Guide1InA", solidGuide1In1, solidGuide1In2, 0, *guide1In2PosVec - *guide1In1PosVec);
  // G4VSolid *solidGuide1In = solidGuide1In1;
  fGuideInPosVecA = guide1In1PosVec;

  /////////////////////////////////////////////////////////////////
  // Light Guide Outside
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
  G4VSolid *solidGuide1Out = new G4GenericTrap("Guide1OutA", guide1OutZ, guide1OutPtVec);
  fGuideOutPosVecA = new G4ThreeVector(0.*mm, 0.*mm, guide1OutZ);

  /////////////////////////////////////////////////////////////////
  // Tilt board In
  /////////////////////////////////////////////////////////////////
  // sakak-kidung
  G4double tiltInX = guide1InBoxZ, tiltInY = tiltInX*tan(fAngleTilt), tiltInZ = guide1InBoxX;
  std::vector <G4TwoVector> tiltInPolygon(3);
  tiltInPolygon[0] = G4TwoVector(-2.*tiltInX - fBoxThickness, 0.*mm);
  tiltInPolygon[1] = G4TwoVector(-fBoxThickness, 0.*mm);
  tiltInPolygon[2] = G4TwoVector(-fBoxThickness, -2.*tiltInY);
  G4ExtrudedSolid *solidTiltIn = new G4ExtrudedSolid("TiltInA", tiltInPolygon, tiltInZ, 
    G4TwoVector(), 1., G4TwoVector(), 1.);

  /////////////////////////////////////////////////////////////////
  // Tilt board Out
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
  G4ExtrudedSolid *solidTiltOut = new G4ExtrudedSolid("TiltOutA", tiltOutPolygon, tiltOutZ, 
    G4TwoVector(), 1., G4TwoVector(), 1.);


  ///////////////////////////////////////////W//////////////////////
  // Tubs : PMT hole to sutract with box and (tiltInHole = In =  PMT hole)
  /////////////////////////////////////////////////////////////////
  G4double tiltHoleZ = fBoxThickness*sin(fAngleTilt);
  G4ThreeVector tiltHolePosVec(-tiltOutX + tiltHoleZ*sin(fAngleTilt), -tiltOutY + tiltHoleZ*cos(fAngleTilt), 0.*mm);
  fRotationScorerA = new G4RotationMatrix();
  fRotationScorerA->rotate(-90.*deg, G4ThreeVector(cos(180.*deg - fAngleTilt), sin(180.*deg - fAngleTilt)));

  G4Tubs *solidTiltHole = new G4Tubs("TiltHoleA", 0.*mm, fCathodeRadius, tiltHoleZ, 0.*deg, 360.*deg);
  G4UnionSolid  *solidTiltInHole = new G4UnionSolid("TiltInHoleA", solidTiltIn, solidTiltHole, fRotationScorerA, tiltHolePosVec);

  fSolidGuideOutA  = new G4UnionSolid("GuideOutA", solidGuide1Out, solidTiltOut, rotationTilt, tiltOutPosVec - *fGuideOutPosVecA);
  fSolidGuideInA   = new G4UnionSolid("GuideInA", solidGuide1In, solidTiltInHole, rotationTilt, tiltOutPosVec- *fGuideInPosVecA);
  fSolidGuideCaseA = new G4SubtractionSolid("GuideCaseA", fSolidGuideOutA, fSolidGuideInA, 0, *fGuideInPosVecA - *fGuideOutPosVecA);
  fLogicGuideCaseA = new G4LogicalVolume(fSolidGuideCaseA, fMatBox, "GuideCaseA");
  fPhysGuideCaseA  = new G4PVPlacement(0,    *fGuideOutPosVecA,  fLogicGuideCaseA,  "GuideCaseA",  fLogicWorld, false, 0, fCheckOverlaps);
 
  fScorerPosVecA = new G4ThreeVector(tiltOutPosVec + G4ThreeVector(0.*mm, -tiltOutY, tiltOutX) + fScorerZ*G4ThreeVector(0.*mm, -cos(fAngleTilt), sin(fAngleTilt)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void G4ACDetectorConstructionLEPS2::constrcutGuideB()
{
  G4double angleWithTPCedge = 360.*deg/fNelem;
  G4double distantFromCenterPMT = 150.*mm;
  /////////////////////////////////////////////////////////////////
  // Light Guide1 (inside air)
  /////////////////////////////////////////////////////////////////
  G4double guide1InZ1 = 80.*mm/2; // distance between top face and TPC edge
  G4double guide1InZ2 = 2.*fCrossX*sin(angleWithTPCedge)/2; // distance between bottom face and TPC edge
  
  // l is leftward distance, r is rightward distance
  // top face 
  G4double guide1In2Xl = 70.*mm/2;
  G4double guide1In2Xr = 70.*mm/2;
  G4double guide1In2Yl = 90.*(1.-87.150428425705/101.704246782948)*mm/2; // = 12.878947*mm/2
  G4double guide1In2Yr = 0.*mm;
  // cross-section cut by TPC edge
  G4double guide1In3Xl = distantFromCenterPMT - fCrossX;

  // bottop face
  G4double guide1In1Xl = (guide1InZ1 + guide1InZ2)*(guide1In3Xl-guide1In2Xl)/guide1InZ1 + guide1In2Xl;
  G4double guide1In1Xr = 2*fCrossX*cos(angleWithTPCedge) - guide1In3Xl;
  G4double guide1In1Yl = (guide1InZ1 + guide1InZ2)*(45.*mm-guide1In2Yl)/guide1InZ1 + guide1In2Yl;
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

  G4VSolid *solidGuide1In1 = new G4GenericTrap("Guide1In1B", guide1InZ1 + guide1InZ2, guide1InPtVec);
  G4double guide1InBoxX = guide1In2Xl, guide1InBoxY = 1.01*fBoxThickness/2, guide1InBoxZ = guide1InZ1;
  G4Box *solidGuide1In2 = new G4Box("Guide1In2B", guide1InBoxX, guide1InBoxY, guide1InBoxZ);
  G4UnionSolid *solidGuide1In = new G4UnionSolid("Guide1InB", solidGuide1In1, solidGuide1In2, 0, G4ThreeVector(0*mm, guide1InBoxY, guide1InZ2));
  fGuideInPosVecB =
    new G4ThreeVector(G4ThreeVector(0.*mm, 0.*mm, guide1InZ1 + guide1InZ2) - G4ThreeVector(guide1In3Xl - cos(angleWithTPCedge)*fCrossX, 0.*mm, sin(angleWithTPCedge)*fCrossX));
  fGuideInPosVecB->rotateY(angleWithTPCedge);
  /////////////////////////////////////////////////////////////////
  // Light Guide Outside part 1
  /////////////////////////////////////////////////////////////////  
  G4double guide1OutZ = guide1InZ1 + guide1InZ2 + fBoxThickness;

  // l is leftward distance, r is rightward distance
  // top face 
  G4double guide1Out2Xl = guide1In2Xl + fBoxThickness;
  G4double guide1Out2Xr = guide1In2Xr + fBoxThickness;
  G4double guide1Out2Yl = guide1In2Yl + fBoxThickness;
  G4double guide1Out2Yr = guide1In2Yr + fBoxThickness;

  // bottop face
  G4double guide1Out1Xl = guide1In1Xl + fBoxThickness;
  G4double guide1Out1Xr = guide1In1Xr + fBoxThickness;
  G4double guide1Out1Yl = guide1In1Yl + fBoxThickness;
  G4double guide1Out1Yr = guide1In1Yr + fBoxThickness;

  std::vector<G4TwoVector> guide1OutPtVec(8);
  guide1OutPtVec[0] = G4TwoVector(-guide1Out1Xr, 0.*mm);
  guide1OutPtVec[1] = G4TwoVector(-guide1Out1Xr, 2*guide1Out1Yr);
  guide1OutPtVec[2] = G4TwoVector(guide1Out1Xl, 2*guide1Out1Yl);
  guide1OutPtVec[3] = G4TwoVector(guide1Out1Xl, 0.*mm);
  guide1OutPtVec[4] = G4TwoVector(-guide1Out2Xr, 0.*mm);
  guide1OutPtVec[5] = G4TwoVector(-guide1Out2Xr, 2*guide1Out2Yr);
  guide1OutPtVec[6] = G4TwoVector(guide1Out2Xl, 2*guide1Out2Yl);
  guide1OutPtVec[7] = G4TwoVector(guide1Out2Xl, 0.*mm);

  G4VSolid *solidGuide1Out = new G4GenericTrap("Guide1OutB", guide1OutZ, guide1OutPtVec);
  fGuideOutPosVecB = new G4ThreeVector(*fGuideInPosVecB);
  

  // to cut the light guide volume to fit in TPC edge
  std::vector <G4TwoVector> tpcEdgeCutPolygon(3);
  G4double tpcEdgeCutXr = 2.*guide1Out1Xr, tpcEdgeCutXl = 2.*guide1Out1Xl, tpcEdgeCutY = (tpcEdgeCutXl + tpcEdgeCutXr)*tan(angleWithTPCedge), tpcEdgeCutZ = 3.*fCrossY;
  tpcEdgeCutPolygon[0] = G4TwoVector(-tpcEdgeCutXr, 0.*mm);
  tpcEdgeCutPolygon[1] = G4TwoVector(tpcEdgeCutXl,0.*mm);
  tpcEdgeCutPolygon[2] = G4TwoVector(tpcEdgeCutXl, -tpcEdgeCutY);

  G4ExtrudedSolid *solidTpcEdgeCut = new G4ExtrudedSolid("SolidTpcEdgeCutB", tpcEdgeCutPolygon, tpcEdgeCutZ,
    G4TwoVector(), 1., G4TwoVector(), 1.);
  G4RotationMatrix *rotationTpcEdgeCut = new G4RotationMatrix();
  rotationTpcEdgeCut->rotateX(90*deg);
  G4ThreeVector tpcEdgeCutPosVec(0.*mm, 0.*mm, guide1InZ2 - guide1InZ1 - tan(angleWithTPCedge)*(guide1In3Xl + tpcEdgeCutXr));

  /////////////////////////////////////////////////////////////////
  // Tilt board In
  /////////////////////////////////////////////////////////////////
  // sakak-kidung
  G4double tiltInX = guide1InBoxX, tiltInY = tiltInX*tan(fAngleTilt), tiltInZ = guide1InBoxZ;
  std::vector <G4TwoVector> tiltInPolygon(3);
  tiltInPolygon[0] = G4TwoVector(-2.*tiltInX - fBoxThickness, -2.*tiltInY);
  tiltInPolygon[1] = G4TwoVector(- fBoxThickness, 0.*mm);
  tiltInPolygon[2] = G4TwoVector(-2.*tiltInX - fBoxThickness);
  G4ExtrudedSolid *solidTiltIn = new G4ExtrudedSolid("TiltInB", tiltInPolygon, tiltInZ, 
    G4TwoVector(), 1., G4TwoVector(), 1.);

  /////////////////////////////////////////////////////////////////
  // Tilt board Out
  /////////////////////////////////////////////////////////////////
  // samkak-kidung
  std::vector <G4TwoVector> tiltOutPolygon(3);
  G4double tiltOutX = tiltInX + fBoxThickness, tiltOutY = tiltOutX*tan(fAngleTilt), tiltOutZ = tiltInZ + fBoxThickness;
  G4ThreeVector tiltOutPosVec(tiltOutX, 0*mm, guide1InZ2);
  tiltOutPolygon[0] = G4TwoVector(-2*tiltOutX, -2.*tiltOutY);
  tiltOutPolygon[1] = G4TwoVector(0.*mm, 0.*mm);
  tiltOutPolygon[2] = G4TwoVector(-2*tiltOutX, 0.*mm);

  G4ExtrudedSolid *solidTiltOut = new G4ExtrudedSolid("TiltOutB", tiltOutPolygon, tiltOutZ, 
    G4TwoVector(), 1., G4TwoVector(), 1.);


  ///////////////////////////////////////////-//////////////////////
  // Tubs : PMT hole to sutract with box and (tiltInHole = In =  PMT hole)
  /////////////////////////////////////////////////////////////////
  G4double tiltHoleZ = fBoxThickness*sin(fAngleTilt);
  G4ThreeVector tiltHolePosVec(-tiltOutX - tiltHoleZ*sin(fAngleTilt), -tiltOutY + tiltHoleZ*cos(fAngleTilt), 0.*mm);

  G4RotationMatrix *rotationTiltHole = new G4RotationMatrix();
  rotationTiltHole->rotate(-90.*deg, G4ThreeVector(cos(fAngleTilt), sin(fAngleTilt)));

  G4Tubs *solidTiltHole = new G4Tubs("TiltHoleB", 0.*mm, fCathodeRadius, tiltHoleZ, 0.*deg, 360.*deg);
  G4UnionSolid  *solidTiltInHole = new G4UnionSolid("TiltInHoleB", solidTiltIn, solidTiltHole, rotationTiltHole, tiltHolePosVec);

  G4UnionSolid *solidGuideOutBeforeTPCedgeCut = new G4UnionSolid("GuideOutBeforeTPCedgeCutB", solidGuide1Out, solidTiltOut, 0, tiltOutPosVec);
  G4UnionSolid *solidGuideInBeforeTPCedgeCut = new G4UnionSolid("GuideInBeforeTPCedgeCutB", solidGuide1In, solidTiltInHole, 0, tiltOutPosVec);
  
  fSolidGuideOutB  = new G4SubtractionSolid("GuideOutB", solidGuideOutBeforeTPCedgeCut, solidTpcEdgeCut, rotationTpcEdgeCut, tpcEdgeCutPosVec);
  fSolidGuideInB   = new G4SubtractionSolid("GuideInB", solidGuideInBeforeTPCedgeCut, solidTpcEdgeCut, rotationTpcEdgeCut, tpcEdgeCutPosVec);
  fSolidGuideCaseB = new G4SubtractionSolid("GuideCaseB", fSolidGuideOutB, fSolidGuideInB, 0, *fGuideInPosVecB - *fGuideOutPosVecB);
  fLogicGuideCaseB = new G4LogicalVolume(fSolidGuideCaseB, fMatBox, "GuideCaseB");

  fRotationGuideB = new G4RotationMatrix();
  fRotationGuideB->rotateY(-angleWithTPCedge);
  // fPhysGuideCase  = new G4PVPlacement(fRotationGuide, *fGuideOutPosVec,  fLogicGuideCase,  "GuideCaseB",  fLogicWorld, false, 0, fCheckOverlaps);
 
  fRotationScorerB = new G4RotationMatrix();
  fRotationScorerB->rotateY(-angleWithTPCedge);
  fRotationScorerB->rotate(90.*deg, G4ThreeVector(cos(fAngleTilt), sin(fAngleTilt)));

  fScorerPosVecB = 
    new G4ThreeVector(G4ThreeVector(0.*mm, -tiltOutY, guide1InZ2) + fScorerZ*G4ThreeVector(sin(fAngleTilt), -cos(fAngleTilt), 0.*mm));
  fScorerPosVecB->rotateY(angleWithTPCedge);
  *fScorerPosVecB = *fScorerPosVecB + *fGuideOutPosVecB;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstructionLEPS2::constrcutGuideC()
{
  G4double angleWithTPCedge = 2*360.*deg/fNelem;
  G4double distantFromCenterPMT = 150.*mm;
  /////////////////////////////////////////////////////////////////
  // Light Guide1 (inside air)
  /////////////////////////////////////////////////////////////////
  G4double guide1InZ1 = 80.*mm/2; // distance between top face and TPC edge
  G4double guide1InZ2 = 2.*fCrossX*sin(angleWithTPCedge)/2; // distance between bottom face and TPC edge
  G4double guide1InZ3 = 2.*fCrossX*sin(angleWithTPCedge/2)/2;

  // l is leftward distance, r is rightward distance
  // top face 
  G4double guide1In2Xl = 70.*mm/2;
  G4double guide1In2Xr = 70.*mm/2;
  G4double guide1In2Yl = 90.*(1.-126.753555907414/155.225120922720)*mm/2; // = 16.507901*mm/2
  G4double guide1In2Yr = 0.*mm;

  // bottop face
  G4double guide1In1Xl = guide1In2Xl;
  G4double guide1In1Xr = 2*fCrossX*cos(angleWithTPCedge) - guide1In2Xl;
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

  G4VSolid *solidGuide1In1 = new G4GenericTrap("Guide1In1C", guide1InZ1 + guide1InZ2 + guide1InZ3, guide1InPtVec);
  G4double guide1InBoxX = guide1In2Xl, guide1InBoxY = 1.01*fBoxThickness/2, guide1InBoxZ = guide1InZ1;
  G4Box *solidGuide1In2 = new G4Box("Guide1In2C", guide1InBoxX, guide1InBoxY, guide1InBoxZ);
  G4UnionSolid *solidGuide1In = new G4UnionSolid("Guide1InC", solidGuide1In1, solidGuide1In2, 0, G4ThreeVector(0*mm, guide1InBoxY, guide1InZ2 + guide1InZ3));
  fGuideInPosVecC =
    new G4ThreeVector(G4ThreeVector(0.*mm, 0.*mm, guide1InZ1 + guide1InZ3 - guide1InZ2) - G4ThreeVector(guide1In2Xl - cos(angleWithTPCedge)*fCrossX, 0.*mm, -sin(angleWithTPCedge)*fCrossX));
  fGuideInPosVecC->rotateY(angleWithTPCedge);
  /////////////////////////////////////////////////////////////////
  // Light Guide Outside part 1
  /////////////////////////////////////////////////////////////////  
  G4double guide1OutZ = guide1InZ1 + guide1InZ2 + guide1InZ3 + fBoxThickness;

  // l is leftward distance, r is rightward distance
  // top face 
  G4double guide1Out2Xl = guide1In2Xl + fBoxThickness;
  G4double guide1Out2Xr = guide1In2Xr + fBoxThickness;
  G4double guide1Out2Yl = guide1In2Yl + fBoxThickness;
  G4double guide1Out2Yr = guide1In2Yr + fBoxThickness;

  // bottop face
  G4double guide1Out1Xl = guide1In1Xl + fBoxThickness;
  G4double guide1Out1Xr = guide1In1Xr + fBoxThickness;
  G4double guide1Out1Yl = guide1In1Yl + fBoxThickness;
  G4double guide1Out1Yr = guide1In1Yr + fBoxThickness;

  std::vector<G4TwoVector> guide1OutPtVec(8);
  guide1OutPtVec[0] = G4TwoVector(-guide1Out1Xr, 0.*mm);
  guide1OutPtVec[1] = G4TwoVector(-guide1Out1Xr, 2*guide1Out1Yr);
  guide1OutPtVec[2] = G4TwoVector(guide1Out1Xl, 2*guide1Out1Yl);
  guide1OutPtVec[3] = G4TwoVector(guide1Out1Xl, 0.*mm);
  guide1OutPtVec[4] = G4TwoVector(-guide1Out2Xr, 0.*mm);
  guide1OutPtVec[5] = G4TwoVector(-guide1Out2Xr, 2*guide1Out2Yr);
  guide1OutPtVec[6] = G4TwoVector(guide1Out2Xl, 2*guide1Out2Yl);
  guide1OutPtVec[7] = G4TwoVector(guide1Out2Xl, 0.*mm);

  G4VSolid *solidGuide1Out = new G4GenericTrap("Guide1OutC", guide1OutZ, guide1OutPtVec);
  fGuideOutPosVecC = new G4ThreeVector(*fGuideInPosVecC);
  

  // to cut the light guide volume to fit in TPC edge
  std::vector <G4TwoVector> tpcEdgeCutPolygon(3);
  G4double tpcEdgeCutXr = 2.*guide1Out1Xr, tpcEdgeCutXl = 2.*guide1Out1Xl, tpcEdgeCutY = (tpcEdgeCutXl + tpcEdgeCutXr)*tan(angleWithTPCedge), tpcEdgeCutZ = 3.*fCrossY;
  tpcEdgeCutPolygon[0] = G4TwoVector(-tpcEdgeCutXr, 0.*mm);
  tpcEdgeCutPolygon[1] = G4TwoVector(tpcEdgeCutXl,0.*mm);
  tpcEdgeCutPolygon[2] = G4TwoVector(tpcEdgeCutXl, -tpcEdgeCutY);

  G4ExtrudedSolid *solidTpcEdgeCut = new G4ExtrudedSolid("SolidTpcEdgeCutC", tpcEdgeCutPolygon, tpcEdgeCutZ,
    G4TwoVector(), 1., G4TwoVector(), 1.);
  G4RotationMatrix *rotationTpcEdgeCut = new G4RotationMatrix();
  rotationTpcEdgeCut->rotateX(90*deg);
  G4ThreeVector tpcEdgeCutPosVec(0.*mm, 0.*mm, guide1InZ2 - guide1InZ3 - guide1InZ1 - tan(angleWithTPCedge)*(guide1In2Xl + tpcEdgeCutXr));

  /////////////////////////////////////////////////////////////////
  // Tilt board In
  /////////////////////////////////////////////////////////////////
  // sakak-kidung
  G4double tiltInX = guide1InBoxX, tiltInY = tiltInX*tan(fAngleTilt), tiltInZ = guide1InBoxZ;
  std::vector <G4TwoVector> tiltInPolygon(3);
  tiltInPolygon[0] = G4TwoVector(-2.*tiltInX - fBoxThickness, -2.*tiltInY);
  tiltInPolygon[1] = G4TwoVector(- fBoxThickness, 0.*mm);
  tiltInPolygon[2] = G4TwoVector(-2.*tiltInX - fBoxThickness);
  G4ExtrudedSolid *solidTiltIn = new G4ExtrudedSolid("TiltInC", tiltInPolygon, tiltInZ, 
    G4TwoVector(), 1., G4TwoVector(), 1.);

  /////////////////////////////////////////////////////////////////
  // Tilt board Out
  /////////////////////////////////////////////////////////////////
  // samkak-kidung
  std::vector <G4TwoVector> tiltOutPolygon(3);
  G4double tiltOutX = tiltInX + fBoxThickness, tiltOutY = tiltOutX*tan(fAngleTilt), tiltOutZ = tiltInZ + fBoxThickness;
  G4ThreeVector tiltOutPosVec(tiltOutX, 0*mm, guide1InZ2 + guide1InZ3);
  tiltOutPolygon[0] = G4TwoVector(-2*tiltOutX, -2.*tiltOutY);
  tiltOutPolygon[1] = G4TwoVector(0.*mm, 0.*mm);
  tiltOutPolygon[2] = G4TwoVector(-2*tiltOutX, 0.*mm);

  G4ExtrudedSolid *solidTiltOut = new G4ExtrudedSolid("TiltOutC", tiltOutPolygon, tiltOutZ, 
    G4TwoVector(), 1., G4TwoVector(), 1.);


  ///////////////////////////////////////////-//////////////////////
  // Tubs : PMT hole to sutract with box and (tiltInHole = In =  PMT hole)
  /////////////////////////////////////////////////////////////////
  G4double tiltHoleZ = fBoxThickness*sin(fAngleTilt);
  G4ThreeVector tiltHolePosVec(-tiltOutX - tiltHoleZ*sin(fAngleTilt), -tiltOutY + tiltHoleZ*cos(fAngleTilt), 0.*mm);

  G4RotationMatrix *rotationTiltHole = new G4RotationMatrix();
  rotationTiltHole->rotate(-90.*deg, G4ThreeVector(cos(fAngleTilt), sin(fAngleTilt)));

  G4Tubs *solidTiltHole = new G4Tubs("TiltHoleC", 0.*mm, fCathodeRadius, tiltHoleZ, 0.*deg, 360.*deg);
  G4UnionSolid  *solidTiltInHole = new G4UnionSolid("TiltInHoleC", solidTiltIn, solidTiltHole, rotationTiltHole, tiltHolePosVec);

  G4UnionSolid *solidGuideOutBeforeTPCedgeCut = new G4UnionSolid("GuideOutBeforeTPCedgeCutC", solidGuide1Out, solidTiltOut, 0, tiltOutPosVec);
  G4UnionSolid *solidGuideInBeforeTPCedgeCut = new G4UnionSolid("GuideInBeforeTPCedgeCutC", solidGuide1In, solidTiltInHole, 0, tiltOutPosVec);
  
  fSolidGuideOutC  = new G4SubtractionSolid("GuideOutC", solidGuideOutBeforeTPCedgeCut, solidTpcEdgeCut, rotationTpcEdgeCut, tpcEdgeCutPosVec);
  fSolidGuideInC   = new G4SubtractionSolid("GuideInC", solidGuideInBeforeTPCedgeCut, solidTpcEdgeCut, rotationTpcEdgeCut, tpcEdgeCutPosVec);
  fSolidGuideCaseC = new G4SubtractionSolid("GuideCaseC", fSolidGuideOutC, fSolidGuideInC, 0, *fGuideInPosVecC - *fGuideOutPosVecC);
  fLogicGuideCaseC = new G4LogicalVolume(fSolidGuideCaseC, fMatBox, "GuideCaseC");

  fRotationGuideC = new G4RotationMatrix();
  fRotationGuideC->rotateY(-angleWithTPCedge);
 
  fRotationScorerC = new G4RotationMatrix();
  fRotationScorerC->rotateY(-angleWithTPCedge);
  fRotationScorerC->rotate(90.*deg, G4ThreeVector(cos(fAngleTilt), sin(fAngleTilt)));

  fScorerPosVecC = 
    new G4ThreeVector(G4ThreeVector(0.*mm, -tiltOutY, guide1InZ2 + guide1InZ3) + fScorerZ*G4ThreeVector(sin(fAngleTilt), -cos(fAngleTilt), 0.*mm));
  fScorerPosVecC->rotateY(angleWithTPCedge);
  *fScorerPosVecC = *fScorerPosVecC + *fGuideOutPosVecC;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
*/
void G4ACDetectorConstructionLEPS2::ConstructSDandField()
{
  auto magField = new G4UniformMagField(1.*tesla, 0., 0.);
  auto fieldMan = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fieldMan->SetDetectorField(magField);
  fieldMan->CreateChordFinder(magField);

  // register HolderCase as Sensitve detector
  G4ACAerogelSD *AerogelSD= new G4ACAerogelSD("AerogelSD", fOffsetAngle);
  G4SDManager *sdMan  = G4SDManager::GetSDMpointer();
  sdMan->AddNewDetector(AerogelSD);
  fLogicAerogel1->SetSensitiveDetector(AerogelSD);
}
