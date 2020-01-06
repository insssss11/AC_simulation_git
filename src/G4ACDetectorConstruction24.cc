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
// $Id: G4ACDetectorConstruction24.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACDetectorConstruction24.cc
/// \brief Definition of the G4ACDetectorConstruction24 class

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

#include "G4ACDetectorConstruction24.hh"
#include "G4ACFineMeshSD.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4ACDetectorConstruction24::G4ACDetectorConstruction24()
: G4ACBaseDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction24::~G4ACDetectorConstruction24()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACDetectorConstruction24::Construct()
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

void G4ACDetectorConstruction24::ConstructGuide()
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

  G4VSolid *solidGuide1In1 = new G4GenericTrap("Guide1In1", guide1InZ, guide1InPtVec);
  G4VSolid *solidGuide1In2 = new G4Box("Guide1In2", guide1InBoxX, guide1InBoxY, guide1InBoxZ);
  
  auto guide1In1PosVec = new G4ThreeVector(0.*mm, 0.*mm, guide1InZ);
  auto guide1In2PosVec = new G4ThreeVector(0.*mm, fBoxThickness, 2*guide1InZ - guide1InBoxZ);

  G4UnionSolid *solidGuide1In = new G4UnionSolid("Guide1In", solidGuide1In1, solidGuide1In2, 0, *guide1In2PosVec - *guide1In1PosVec);
  // G4VSolid *solidGuide1In = solidGuide1In1;
  fGuideInPosVec = guide1In1PosVec;

  /////////////////////////////////////////////////////////////////
  // Light Guide Outside
  /////////////////////////////////////////////////////////////////  
  G4double guide1OutZ = guide1InZ + fBoxThickness/2;

  // bottom face half length
  G4double guide1Out1X = fCrossX + fBoxThickness;
  G4double guide1Out1Y = fCrossY + fBoxThickness/2;

  // top face
  G4double guide1Out2X = guide1In2X + fBoxThickness;
  G4double guide1Out2Y = guide1In2Y + fBoxThickness/2;
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
  // Tilt board In
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
  G4ExtrudedSolid *solidTiltOut = new G4ExtrudedSolid("TiltOut", tiltOutPolygon, tiltOutZ, 
    G4TwoVector(), 1., G4TwoVector(), 1.);


  ///////////////////////////////////////////W//////////////////////
  // Tubs : PMT hole to sutract with box and (tiltInHole = In =  PMT hole)
  /////////////////////////////////////////////////////////////////
  G4double tiltHoleZ = fBoxThickness*sin(fAngleTilt);
  G4ThreeVector tiltHolePosVec(-tiltOutX + tiltHoleZ*sin(fAngleTilt), -tiltOutY + tiltHoleZ*cos(fAngleTilt), 0.*mm);
  fRotationScorer = new G4RotationMatrix();
  fRotationScorer->rotate(-90.*deg, G4ThreeVector(cos(180.*deg - fAngleTilt), sin(180.*deg - fAngleTilt)));

  G4Tubs *solidTiltHole = new G4Tubs("TiltHole", 0.*mm, fCathodeRadius, tiltHoleZ, 0.*deg, 360.*deg);
  G4UnionSolid  *solidTiltInHole = new G4UnionSolid("TiltInHole", solidTiltIn, solidTiltHole, fRotationScorer, tiltHolePosVec);

  fSolidGuideOut  = new G4UnionSolid("GuideOut", solidGuide1Out, solidTiltOut, rotationTilt, tiltOutPosVec - *fGuideOutPosVec);
  fSolidGuideIn   = new G4UnionSolid("GuideIn", solidGuide1In, solidTiltInHole, rotationTilt, tiltOutPosVec- *fGuideInPosVec);
  fSolidGuideCase = new G4SubtractionSolid("GuideCase", fSolidGuideOut, fSolidGuideIn, 0, *fGuideInPosVec - *fGuideOutPosVec);
  fLogicGuideCase = new G4LogicalVolume(fSolidGuideCase, fMatBox, "GuideCase");
  fPhysGuideCase  = new G4PVPlacement(0,    *fGuideOutPosVec,  fLogicGuideCase,  "GuideCase",  fLogicWorld, false, 0, fCheckOverlaps);
 
  fScorerPosVec = new G4ThreeVector(tiltOutPosVec + G4ThreeVector(0.*mm, -tiltOutY, tiltOutX) + fScorerZ*G4ThreeVector(0.*mm, -cos(fAngleTilt), sin(fAngleTilt)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction24::ConstructScorer()
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

void G4ACDetectorConstruction24::ConstructHolderBoundary()
{
  // photon energy definition
  /*
  const G4int nEntries = 20;
  const G4double pEmin = 1.9074*eV;
  const G4double pEmax = 6.1992*eV;
  G4double pE[nEntries] = {};

  for(G4int j = 0; j < nEntries;j++)
  {
    G4double dE = (pEmax - pEmin)/(nEntries- 1);
    pE[j] = pEmin + j*dE;
  }

  G4double reflectivity[nEntries] = {}; // reflection is decided by reflectivity only in dielectric_metal surface type!
  for(G4int j = 0;j < nEntries;j++)
  {
    reflectivity[j]= fSplineAlMylarRef->Eval(pE[j]);
  }
  */
  // G4MaterialPropertiesTable *MPTsurface = new G4MaterialPropertiesTable();
  // MPTsurface->AddProperty("REFLECTIVITY", pE, reflectivity, nEntries);
  G4OpticalSurface *opInHolderCaseSurface = new G4OpticalSurface("InHolderCaseSurface");
  opInHolderCaseSurface->SetType(dielectric_LUTDAVIS);
  opInHolderCaseSurface->SetFinish(RoughTeflon_LUT);
  opInHolderCaseSurface->SetModel(DAVIS);
  
  G4LogicalBorderSurface *logicSurfaceInHolderCase = new G4LogicalBorderSurface("InHolderCaseSurface",
    fPhysIn, fPhysHolderCase, opInHolderCaseSurface);  
  // opInHolderCaseSurface->SetMaterialPropertiesTable(MPTsurface);
  opInHolderCaseSurface = dynamic_cast<G4OpticalSurface*>
    (logicSurfaceInHolderCase->GetSurface(fPhysIn, fPhysHolderCase));
  if(opInHolderCaseSurface) opInHolderCaseSurface->DumpInfo();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
