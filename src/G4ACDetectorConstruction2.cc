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
// $Id: G4ACDetectorConstruction2.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACDetectorConstruction2.cc
/// \brief Implementation of the G4ACDetectorConstruction2 class

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
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

#include "G4Tubs.hh"
#include "G4GenericTrap.hh"

#include "globals.hh"
#include <math.h>

#include "G4ACDetectorConstruction2.hh"
#include "G4ACFineMeshSD.hh"


//....oooOO0OOooo........ooo
G4ACDetectorConstruction2::G4ACDetectorConstruction2()
: G4ACDetectorConstruction1()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction2::~G4ACDetectorConstruction2()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACDetectorConstruction2::Construct()
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
void G4ACDetectorConstruction2::ConstructAerogel()
{
  const G4double aerogelX = 92.*mm/2;
  const G4double aerogelY = fAerogelThickness;
  const G4double aerogelZ = 92.*mm/2;
  
  fSolidAerogel1 = new G4Box("Aerogel1", aerogelX - fAerogelGap, aerogelY - fAerogelGap, aerogelZ - fAerogelGap);
  fSolidAerogel2 = new G4Box("Aerogel2", aerogelX - fAerogelGap, aerogelY - fAerogelGap, aerogelZ - fAerogelGap);
  fSolidAerogel3 = new G4Box("Aerogel3", aerogelX - fAerogelGap, aerogelY - fAerogelGap, aerogelZ - fAerogelGap);
  fSolidAerogel4 = new G4Box("Aerogel4", aerogelX - fAerogelGap, aerogelY - fAerogelGap, aerogelZ - fAerogelGap);

  fLogicAerogel1 = new G4LogicalVolume(fSolidAerogel1, fMatAerogel, "Aerogel1");
  fLogicAerogel2 = new G4LogicalVolume(fSolidAerogel2, fMatAerogel, "Aerogel2");
  fLogicAerogel3 = new G4LogicalVolume(fSolidAerogel3, fMatAerogel, "Aerogel3");
  fLogicAerogel4 = new G4LogicalVolume(fSolidAerogel4, fMatAerogel, "Aerogel4");

  fAerogelPosVec1 = new G4ThreeVector(0.*mm, fBoxThickness + aerogelY, fHolderInPosVec->z() + aerogelZ);
  fAerogelPosVec2 = new G4ThreeVector(0.*mm, fBoxThickness + aerogelY, fHolderInPosVec->z() + 3*aerogelZ);
  fAerogelPosVec3 = new G4ThreeVector(0.*mm, fBoxThickness + 3*aerogelY, fHolderInPosVec->z() + aerogelZ);
  fAerogelPosVec4 = new G4ThreeVector(0.*mm, fBoxThickness + 3*aerogelY, fHolderInPosVec->z() + 3*aerogelZ);

  new G4PVPlacement(0, *fAerogelPosVec1, fLogicAerogel1, "Aerogel1", fLogicIn, true, 0, fCheckOverlaps);
  new G4PVPlacement(0, *fAerogelPosVec2, fLogicAerogel1, "Aerogel2", fLogicIn, true, 1, fCheckOverlaps);
  new G4PVPlacement(0, *fAerogelPosVec3, fLogicAerogel1, "Aerogel3", fLogicIn, true, 2, fCheckOverlaps);
  new G4PVPlacement(0, *fAerogelPosVec4, fLogicAerogel1, "Aerogel4", fLogicIn, true, 3, fCheckOverlaps);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......