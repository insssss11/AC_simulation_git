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
// $Id: G4ACDetectorConstruction7.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACDetectorConstruction7.cc
/// \brief Implementation of the G4ACDetectorConstruction7 class

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

#include "G4ACDetectorConstruction7.hh"
#include "G4ACFineMeshSD.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction7::G4ACDetectorConstruction7()
: G4ACDetectorConstruction5()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction7::~G4ACDetectorConstruction7()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACDetectorConstruction7::Construct()
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

void G4ACDetectorConstruction7::ConstructHolder()
{
  /////////////////////////////////////////////////////////////////
  // Aerogel Holder part 1 (1 in)
  /////////////////////////////////////////////////////////////////  
  // bottom face
  G4double holder1In1X = fAerogel31X;
  G4double holder1In1Y = 60.*mm/2;
  
  // top face
  G4double holder1In2X = fCrossX;
  G4double holder1In2Y = fCrossY;
  
  G4double holder1InZ = fAerogel3Z;

  std::vector<G4TwoVector> holder1InPtVec(8);

  holder1InPtVec[0] = G4TwoVector(-holder1In1X,  fBoxThickness);
  holder1InPtVec[1] = G4TwoVector(-holder1In1X,  fBoxThickness + 2*holder1In1Y);
  holder1InPtVec[2] = G4TwoVector(holder1In1X,  fBoxThickness + 2*holder1In1Y);
  holder1InPtVec[3] = G4TwoVector(holder1In1X, fBoxThickness);
  holder1InPtVec[4] = G4TwoVector(-holder1In2X,  fBoxThickness);
  holder1InPtVec[5] = G4TwoVector(-holder1In2X,  fBoxThickness + 2*holder1In2Y);  
  holder1InPtVec[6] = G4TwoVector(holder1In2X,  fBoxThickness + 2*holder1In2Y);
  holder1InPtVec[7] = G4TwoVector(holder1In2X, fBoxThickness);

  G4GenericTrap *solidHolder1In = new G4GenericTrap("Holder1In", holder1InZ, holder1InPtVec);

  /////////////////////////////////////////////////////////////////
  // Aerogel Holder part 2 (2 in)
  /////////////////////////////////////////////////////////////////
  // bottom face
  G4double holder2In1X = fAerogel11X;
  G4double holder2In1Y = fAerogelThickness;
  
  // top face
  G4double holder2In2X = holder1In1X;
  G4double holder2In2Y = holder1In1Y;
  
  G4double holder2InZ = fAerogel1Z - fAerogel3Z;

  std::vector<G4TwoVector> holder2InPtVec(8);

  holder2InPtVec[0] = G4TwoVector(-holder2In1X,  fBoxThickness);
  holder2InPtVec[1] = G4TwoVector(-holder2In1X,  fBoxThickness + 2*holder2In1Y);
  holder2InPtVec[2] = G4TwoVector(holder2In1X,  fBoxThickness + 2*holder2In1Y);
  holder2InPtVec[3] = G4TwoVector(holder2In1X, fBoxThickness);
  holder2InPtVec[4] = G4TwoVector(-holder2In2X,  fBoxThickness);
  holder2InPtVec[5] = G4TwoVector(-holder2In2X,  fBoxThickness + 2*holder2In2Y);  
  holder2InPtVec[6] = G4TwoVector(holder2In2X,  fBoxThickness + 2*holder2In2Y);
  holder2InPtVec[7] = G4TwoVector(holder2In2X, fBoxThickness);

  G4GenericTrap *solidHolder2In = new G4GenericTrap("Holder2In", holder2InZ, holder2InPtVec);
  fHolderInPosVec = new G4ThreeVector(0.*mm, 0.*mm, -holder1InZ);
  fSolidHolderIn = new G4UnionSolid("HolderIn", solidHolder1In, solidHolder2In, 0, G4ThreeVector(0.*mm, 0.*mm, -(holder1InZ + holder2InZ)));

  /////////////////////////////////////////////////////////////////
  // Aerogel Holder part 3 (1 out)
  /////////////////////////////////////////////////////////////////  
  // bottom face
  G4double holder1Out1X = fAerogel31X;
  G4double holder1Out1Y = holder1In1Y + fBoxThickness;
  
  // top face
  G4double holder1Out2X = fCrossX + fBoxThickness;
  G4double holder1Out2Y = fCrossY + fBoxThickness;
  
  G4double holder1OutZ = fAerogel3Z;

  std::vector<G4TwoVector> holder1OutPtVec(8);

  holder1OutPtVec[0] = G4TwoVector(-holder1Out1X,  0);
  holder1OutPtVec[1] = G4TwoVector(-holder1Out1X,  2*holder1Out1Y);
  holder1OutPtVec[2] = G4TwoVector(holder1Out1X,  2*holder1Out1Y);
  holder1OutPtVec[3] = G4TwoVector(holder1Out1X, 0);
  holder1OutPtVec[4] = G4TwoVector(-holder1Out2X,  0);
  holder1OutPtVec[5] = G4TwoVector(-holder1Out2X,  2*holder1Out2Y);  
  holder1OutPtVec[6] = G4TwoVector(holder1Out2X,  2*holder1Out2Y);
  holder1OutPtVec[7] = G4TwoVector(holder1Out2X, 0);
  
  fHolderOutPosVec = new G4ThreeVector(*fHolderInPosVec);
  G4GenericTrap *solidHolder1Out = new G4GenericTrap("Holder1Out", holder1OutZ, holder1OutPtVec);

  /////////////////////////////////////////////////////////////////
  // Aerogel Holder part 4 (2 out)
  /////////////////////////////////////////////////////////////////
  // bottom face
  G4double holder2Out1X = fAerogel11X + fBoxThickness*(1 - tan(5.*deg));
  G4double holder2Out1Y = fAerogelThickness + fBoxThickness*(1 - tan(4.46*deg)/2);
  
  // top face
  G4double holder2Out2X = holder1Out1X;
  G4double holder2Out2Y = holder1Out1Y;
  
  G4double holder2OutZ = fAerogel1Z - fAerogel3Z + fBoxThickness/2;

  std::vector<G4TwoVector> holder2OutPtVec(8);

  holder2OutPtVec[0] = G4TwoVector(-holder2Out1X,  0);
  holder2OutPtVec[1] = G4TwoVector(-holder2Out1X,  2*holder2Out1Y);
  holder2OutPtVec[2] = G4TwoVector(holder2Out1X,  2*holder2Out1Y);
  holder2OutPtVec[3] = G4TwoVector(holder2Out1X, 0);
  holder2OutPtVec[4] = G4TwoVector(-holder2Out2X,  0);
  holder2OutPtVec[5] = G4TwoVector(-holder2Out2X,  2*holder2Out2Y);  
  holder2OutPtVec[6] = G4TwoVector(holder2Out2X,  2*holder2Out2Y);
  holder2OutPtVec[7] = G4TwoVector(holder2Out2X, 0);

  G4GenericTrap *solidHolder2Out = new G4GenericTrap("Holder2Out", holder2OutZ, holder2OutPtVec);
  fSolidHolderOut = new G4UnionSolid("HolderOut", solidHolder1Out, solidHolder2Out, 0, G4ThreeVector(0.*mm, 0.*mm, -(holder1OutZ + holder2OutZ)));

  fSolidHolderCase  = new G4SubtractionSolid("HolderCase", fSolidHolderOut, fSolidHolderIn, 0, *fHolderInPosVec - *fHolderOutPosVec);
  fLogicHolderCase = new G4LogicalVolume(fSolidHolderCase, fMatBox, "HolderCase");
  fPhysHolderCase = new G4PVPlacement(0,    *fHolderOutPosVec,  fLogicHolderCase,  "HolderCase",  fLogicWorld, false, 0, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

