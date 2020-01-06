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
// $Id: G4ACDetectorConstruction4.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACDetectorConstruction4.cc
/// \brief Implementation of the G4ACDetectorConstruction4 class

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
#include "G4GenericTrap.hh"
#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "globals.hh"
#include <math.h>

#include "G4ACDetectorConstruction4.hh"
#include "G4ACFineMeshSD.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction4::G4ACDetectorConstruction4()
: G4ACBaseDetectorConstruction(110.*mm/2, 110.*mm/2, 20.*mm/2, 0.05*mm, 30.*deg),
  fAerogelPosVec4(0), fAerogelPosVec5(0), fAerogelPosVec6(0),
  fSolidAerogel4(0), fSolidAerogel5(0), fSolidAerogel6(0),
  fLogicAerogel4(0), fLogicAerogel5(0), fLogicAerogel6(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction4::~G4ACDetectorConstruction4()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACDetectorConstruction4::Construct()
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

void G4ACDetectorConstruction4::ConstructHolder()
{
  /////////////////////////////////////////////////////////////////
  // Aerogel Holder part 1 (inside airs)
  /////////////////////////////////////////////////////////////////  
  // bottom face
  G4double holderIn1X = 110.*mm/2;
  G4double holderIn1Y = 70.*mm/2;
  
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

void G4ACDetectorConstruction4::ConstructAerogel()
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

void G4ACDetectorConstruction4::ConstructAerogelBoundary()
{
  G4OpticalSurface *opAerogel1AirSurface = new G4OpticalSurface("Aerogel1AirSurface");
  opAerogel1AirSurface->SetType(dielectric_dielectric);
  opAerogel1AirSurface->SetFinish(groundair);
  opAerogel1AirSurface->SetModel(unified);
  new G4LogicalSkinSurface("Aerogel1AirSurface", fLogicAerogel1, opAerogel1AirSurface);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
