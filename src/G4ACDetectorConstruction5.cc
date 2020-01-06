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
// $Id: G4ACDetectorConstruction5.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACDetectorConstruction5.cc
/// \brief Implementation of the G4ACDetectorConstruction5 class

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
#include "G4Trap.hh"
#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "G4GenericTrap.hh"

#include "globals.hh"
#include <math.h>

#include "G4ACDetectorConstruction5.hh"
#include "G4ACFineMeshSD.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction5::G4ACDetectorConstruction5()
: G4ACBaseDetectorConstruction(144.4*mm/2, 90.*mm/2)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction5::~G4ACDetectorConstruction5()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACDetectorConstruction5::Construct()
{
  //.......oooOO0OOooo........oooOO0OOooo......
  

  /////////////////////////////////////////////////////////////////
  // world
  /////////////////////////////////////////////////////////////////
  G4double worldHalfSizeX = 200.*mm;
  G4double worldHalfSizeY = 300.*mm;
  G4double worldHalfSizeZ = 400.*mm;
  fSolidWorld = new G4Box("World", worldHalfSizeX, worldHalfSizeY, worldHalfSizeZ);
  fLogicWorld   = new G4LogicalVolume(fSolidWorld, fMatAir, "World");
  fPhysWorld      = new G4PVPlacement(0,    G4ThreeVector(),  fLogicWorld,  "World",  0, false, 0, fCheckOverlaps);

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


void G4ACDetectorConstruction5::ConstructGuide()
{
  /////////////////////////////////////////////////////////////////
  // Light Guide1 (inside air)
  /////////////////////////////////////////////////////////////////
  G4double guide1InZ = 100.*mm/2;
  // bottop face
  G4double guide1In1X = fCrossX;
  G4double guide1In1Y = fCrossY;

  // top fae
  G4double guide1In2X = 65.*mm/2;
  G4double guide1In2Y = 65.*mm/2;

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

  G4double guide2OutRmin = 0, guide2OutRmax = 2.*guide1Out2X;
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
  G4double tiltOutX = guide1Out2X, tiltOutY = tiltOutX*tan(fAngleTilt), tiltOutZ = guide1Out2Y;
  G4ThreeVector tiltOutPosVec(tiltOutX, 0.*mm, 2.*guide1OutZ + tiltOutZ);

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
  
  fSolidGuideOut  = new G4UnionSolid("GuideOut", solidGuideTubOut, solidTiltOut, 0, tiltOutPosVec - *fGuideOutPosVec);
  fSolidGuideIn   = new G4UnionSolid("GuideIn", solidGuideTubIn, solidTiltInHole, 0, tiltOutPosVec- *fGuideInPosVec);
  fSolidGuideCase = new G4SubtractionSolid("GuideCase", fSolidGuideOut, fSolidGuideIn, 0, *fGuideInPosVec - *fGuideOutPosVec);
  fLogicGuideCase = new G4LogicalVolume(fSolidGuideCase, fMatBox, "GuideCase");
  fPhysGuideCase  = new G4PVPlacement(0,    *fGuideOutPosVec,  fLogicGuideCase,  "GuideCase",  fLogicWorld, false, 0, fCheckOverlaps);
 
  fScorerPosVec = new G4ThreeVector(tiltOutPosVec + G4ThreeVector(-tiltOutX, -tiltOutY, 0.*mm) - fScorerZ*G4ThreeVector(sin(fAngleTilt), cos(fAngleTilt), 0.*mm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......