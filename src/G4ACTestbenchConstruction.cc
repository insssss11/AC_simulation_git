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
// $Id:G4ACTestbenchConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \fileG4ACTestbenchConstruction.cc
/// \brief Definition of theG4ACTestbenchConstruction class

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

#include "G4ACTestbenchConstruction.hh"
#include "G4ACFineMeshSD.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACTestbenchConstruction::G4ACTestbenchConstruction()
:G4ACBaseDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACTestbenchConstruction::~G4ACTestbenchConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *G4ACTestbenchConstruction::Construct()
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
  ConstructGuide();
  ConstructInnerAir();

  ConstructGuideBoundary();
  ConstructScorer();
  // ConstructGlassBoundary();  

  #ifndef AC_TESTBENCH4 // in set 4 no aerogel block is used
    ConstructAerogel();
    ConstructAerogelBoundary();
  #endif

  
  return fPhysWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACTestbenchConstruction::ConstructGuide()
{
  /////////////////////////////////////////////////////////////////
  // Light Guide Part
  /////////////////////////////////////////////////////////////////    
  G4double guide1Out1X = fAerogel11X + fBoxThickness;
  G4double guide1Out2X = fAerogel12X + fBoxThickness;
  G4double guide1OutY = fAerogel1Y + fBoxThickness;
  G4double guide2Out1X = fAerogel21X + fBoxThickness;
  G4double guide2Out2X = fAerogel22X + fBoxThickness;  
  G4double guide2OutY = fAerogel2Y + fBoxThickness;
  G4double guide1OutZ = fAerogel1Z + fBoxThickness;
  G4double guide2OutZ = fAerogel2Z + fBoxThickness;
  
  // aerogel slot close to LED
  std::vector<G4TwoVector> guide1OutPtVec(8);
  guide1OutPtVec[0] = G4TwoVector(-guide1Out1X, 0.*mm);
  guide1OutPtVec[1] = G4TwoVector(-guide1Out2X, 2*guide1OutY);
  guide1OutPtVec[2] = G4TwoVector(guide1Out2X, 2*guide1OutY);
  guide1OutPtVec[3] = G4TwoVector(guide1Out1X, 0.*mm);
  guide1OutPtVec[4] = G4TwoVector(-guide1Out1X, 0.*mm);
  guide1OutPtVec[5] = G4TwoVector(-guide1Out2X, 2*guide1OutY);
  guide1OutPtVec[6] = G4TwoVector(guide1Out2X, 2*guide1OutY);
  guide1OutPtVec[7] = G4TwoVector(guide1Out1X, 0.*mm);

  // aerogel slot far to LED
  std::vector<G4TwoVector> guide2OutPtVec(8);
  guide2OutPtVec[0] = G4TwoVector(-guide2Out1X, 0.*mm);
  guide2OutPtVec[1] = G4TwoVector(-guide2Out2X, 2*guide2OutY);
  guide2OutPtVec[2] = G4TwoVector(guide2Out2X, 2*guide2OutY);
  guide2OutPtVec[3] = G4TwoVector(guide2Out1X, 0.*mm);
  guide2OutPtVec[4] = G4TwoVector(-guide2Out1X, 0.*mm);
  guide2OutPtVec[5] = G4TwoVector(-guide2Out2X, 2*guide2OutY);
  guide2OutPtVec[6] = G4TwoVector(guide2Out2X, 2*guide2OutY);
  guide2OutPtVec[7] = G4TwoVector(guide2Out1X, 0.*mm);

  G4VSolid *solidGuide1Out = new G4GenericTrap("Guide1Out", guide1OutZ, guide1OutPtVec);
  G4VSolid *solidGuide2Out = new G4GenericTrap("Guide2Out", guide2OutZ, guide2OutPtVec);

  fGuideOutPosVec = new G4ThreeVector(0.*mm, 0.*mm, guide1OutZ);
  G4UnionSolid *fSolidGuideOut = new G4UnionSolid("GuideOut", solidGuide1Out, solidGuide2Out, 0, G4ThreeVector(0.*mm, 0.*mm, 2.*guide2OutZ));
  
  /////////////////////////////////////////////////////////////////
  // Light Guide In
  /////////////////////////////////////////////////////////////////  
  G4double guide1InZ = fAerogelThickness;
  G4double guide2InZ = guide1InZ;
  G4double guide3InX = 10.*mm/2;
  G4double guide3InY = guide3InX;
  G4double guide3InZ = 2.*(fAerogelThickness + fBoxThickness);

  std::vector<G4TwoVector> guide1InPtVec(8);
  guide1InPtVec[0] = G4TwoVector(-fAerogel11X, fBoxThickness);
  guide1InPtVec[1] = G4TwoVector(-fAerogel12X, fBoxThickness + 2.*fAerogel1Y);
  guide1InPtVec[2] = G4TwoVector(fAerogel12X, fBoxThickness + 2.*fAerogel1Y);
  guide1InPtVec[3] = G4TwoVector(fAerogel11X, fBoxThickness);
  guide1InPtVec[4] = G4TwoVector(-fAerogel11X, fBoxThickness);
  guide1InPtVec[5] = G4TwoVector(-fAerogel12X, fBoxThickness + 2.*fAerogel1Y);
  guide1InPtVec[6] = G4TwoVector(fAerogel12X, fBoxThickness + 2.*fAerogel1Y);
  guide1InPtVec[7] = G4TwoVector(fAerogel11X,fBoxThickness);

  std::vector<G4TwoVector> guide2InPtVec(8);
  guide2InPtVec[0] = G4TwoVector(-fAerogel21X, fBoxThickness);
  guide2InPtVec[1] = G4TwoVector(-fAerogel22X, fBoxThickness + 2.*fAerogel2Y);
  guide2InPtVec[2] = G4TwoVector(fAerogel22X, fBoxThickness + 2.*fAerogel2Y);
  guide2InPtVec[3] = G4TwoVector(fAerogel21X, fBoxThickness);
  guide2InPtVec[4] = G4TwoVector(-fAerogel21X, fBoxThickness);
  guide2InPtVec[5] = G4TwoVector(-fAerogel22X, fBoxThickness + 2.*fAerogel2Y);
  guide2InPtVec[6] = G4TwoVector(fAerogel22X, fBoxThickness + 2.*fAerogel2Y);
  guide2InPtVec[7] = G4TwoVector(fAerogel21X, fBoxThickness);

  fGuideInPosVec = new G4ThreeVector(0.*mm, 0.*mm, guide1OutZ);
  G4VSolid *solidGuide1In = new G4GenericTrap("Guide1In", guide1InZ, guide1InPtVec);
  G4VSolid *solidGuide2In = new G4GenericTrap("Guide2In", guide2InZ, guide2InPtVec);
  G4VSolid *solidGuide3In = new G4Box("Guide3In", guide3InX, guide3InY, guide3InZ);
  G4UnionSolid *solidGuide13In = new G4UnionSolid("Guide13In", solidGuide1In, solidGuide3In, 0,
    G4ThreeVector(0.*mm, 30.*mm + fBoxThickness, guide3InZ - guide1OutZ));
  fSolidGuideIn = new G4UnionSolid("GuideIn", solidGuide13In, solidGuide2In, 0,
    G4ThreeVector(0.*mm, 0.*mm, guide1OutZ + guide2OutZ));

  /////////////////////////////////////////////////////////////////
  // Light Guide Reflector
  /////////////////////////////////////////////////////////////////  
  G4double guideReflector1Z = fAerogelThickness + fReflectorThickness;
  G4double guideReflector2Z = guideReflector1Z;
  G4double guideReflector3X = 10.*mm/2 + fReflectorThickness;
  G4double guideReflector3Y = guideReflector3X;
  G4double guideReflector3Z = 2.*(fAerogelThickness + fBoxThickness);

  std::vector<G4TwoVector> guideReflector1PtVec(8);
  guideReflector1PtVec[0] = G4TwoVector(-fAerogel11X - fReflectorThickness/cos(15.*deg), fBoxThickness - fReflectorThickness);
  guideReflector1PtVec[1] = G4TwoVector(-fAerogel12X - fReflectorThickness/cos(15.*deg), fBoxThickness + fReflectorThickness + 2.*fAerogel1Y);
  guideReflector1PtVec[2] = G4TwoVector(fAerogel12X + fReflectorThickness/cos(15.*deg), fBoxThickness + fReflectorThickness + 2.*fAerogel1Y);
  guideReflector1PtVec[3] = G4TwoVector(fAerogel11X + fReflectorThickness/cos(15.*deg), fBoxThickness - fReflectorThickness);
  guideReflector1PtVec[4] = G4TwoVector(-fAerogel11X - fReflectorThickness/cos(15.*deg), fBoxThickness - fReflectorThickness);
  guideReflector1PtVec[5] = G4TwoVector(-fAerogel12X - fReflectorThickness/cos(15.*deg), fBoxThickness + fReflectorThickness + 2.*fAerogel1Y);
  guideReflector1PtVec[6] = G4TwoVector(fAerogel12X + fReflectorThickness/cos(15.*deg), fBoxThickness + fReflectorThickness + 2.*fAerogel1Y);
  guideReflector1PtVec[7] = G4TwoVector(fAerogel11X + fReflectorThickness/cos(15.*deg), fBoxThickness - fReflectorThickness);

  std::vector<G4TwoVector> guideReflector2PtVec(8);
  guideReflector2PtVec[0] = G4TwoVector(-fAerogel21X - fReflectorThickness/cos(15.*deg), fBoxThickness - fReflectorThickness);
  guideReflector2PtVec[1] = G4TwoVector(-fAerogel22X - fReflectorThickness/cos(15.*deg), fBoxThickness + fReflectorThickness + 2.*fAerogel2Y);
  guideReflector2PtVec[2] = G4TwoVector(fAerogel22X + fReflectorThickness/cos(15.*deg), fBoxThickness + fReflectorThickness + 2.*fAerogel2Y);
  guideReflector2PtVec[3] = G4TwoVector(fAerogel21X + fReflectorThickness/cos(15.*deg), fBoxThickness - fReflectorThickness);
  guideReflector2PtVec[4] = G4TwoVector(-fAerogel21X - fReflectorThickness/cos(15.*deg), fBoxThickness - fReflectorThickness);
  guideReflector2PtVec[5] = G4TwoVector(-fAerogel22X - fReflectorThickness/cos(15.*deg), fBoxThickness + fReflectorThickness + 2.*fAerogel2Y);
  guideReflector2PtVec[6] = G4TwoVector(fAerogel22X + fReflectorThickness/cos(15.*deg), fBoxThickness + fReflectorThickness + 2.*fAerogel2Y);
  guideReflector2PtVec[7] = G4TwoVector(fAerogel21X + fReflectorThickness/cos(15.*deg), fBoxThickness - fReflectorThickness);

  // to make the centers of masses locate at Origin
  fGuideReflectorPosVec = new G4ThreeVector(0.*mm, 0.*mm, guide1OutZ);
  G4VSolid *solidGuideReflector1 = new G4GenericTrap("GuideReflector1", guideReflector1Z, guideReflector1PtVec);
  G4VSolid *solidGuideReflector2 = new G4GenericTrap("GuideReflector2", guideReflector2Z, guideReflector2PtVec);
  G4VSolid *solidGuideReflector3 = new G4Box("GuideReflector3", guideReflector3X, guideReflector3Y, guideReflector3Z);
  G4UnionSolid *solidGuideReflector13 = new G4UnionSolid("GuideReflector13Out", solidGuideReflector1, solidGuideReflector3, 0,
    G4ThreeVector(0.*mm, 30.*mm + fBoxThickness, guideReflector3Z - guide1OutZ));
  G4UnionSolid *solidGuideReflector123 = new G4UnionSolid("GuideReflector123", solidGuideReflector13, solidGuideReflector2, 0,
    G4ThreeVector(0.*mm, 0.*mm, guide1OutZ + guide2OutZ));
  
  fSolidGuideCase = new G4SubtractionSolid("GuideCase", fSolidGuideOut, solidGuideReflector123, 0, *fGuideOutPosVec - *fGuideReflectorPosVec);
  fSolidGuideReflector = new G4SubtractionSolid("GuideReflector", solidGuideReflector123, fSolidGuideIn, 0, *fGuideReflectorPosVec - *fGuideInPosVec);

  fLogicGuideCase = new G4LogicalVolume(fSolidGuideCase, fMatBox, "GuideCase");
  fLogicGuideReflector = new G4LogicalVolume(fSolidGuideReflector, fMatReflector, "GuideReflector");

  fPhysGuideCase  = new G4PVPlacement(0,    *fGuideOutPosVec,  fLogicGuideCase,  "GuideCase",  fLogicWorld, false, 0, fCheckOverlaps);
  fPhysGuideReflector  = new G4PVPlacement(0,    *fGuideReflectorPosVec,  fLogicGuideReflector,  "GuideReflector",  fLogicWorld, false, 0, fCheckOverlaps);
  
  fScorerPosVec = new G4ThreeVector(0.*mm, 30.*mm + fBoxThickness, 150.*mm/2 + 2.*(guide1OutZ + guide2OutZ));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACTestbenchConstruction::ConstructScorer()
{
  fRotationScorer = new G4RotationMatrix();
  fRotationScorer->rotateX(180.*deg);
  /////////////////////////////////////////////////////////////////
  // PMT = Scorer volume (this has Glass, cathode volume as daughter volumes)
  /////////////////////////////////////////////////////////////////
  fSolidScorer = new G4Tubs("Scorer", 0.*mm, 60.*mm/2, 150.*mm/2, 0.*deg, 360.*deg);
  // the photon go to the cathode of FMPMT through UV silicon cookie(glass)
  fSolidGlass  = new G4Tubs("Glass", 0.*mm, 53.*mm/2, 1.*mm/2, 0.*deg, 360.*deg);
  // cathode
  fSolidCathode  = new G4Tubs("Cathode", 0.*mm, 46.*mm/2, 10.*mm/2, 0.*deg, 360.*deg);
  fLogicScorer  = new G4LogicalVolume(fSolidScorer, fMatScorer, "Scorer");
  fLogicGlass   = new G4LogicalVolume(fSolidGlass, fMatGlass, "Glass");
  fLogicCathode = new G4LogicalVolume(fSolidCathode, fMatAir, "Cathode");

  fPhysScorer     = new G4PVPlacement(fRotationScorer, *fScorerPosVec, fLogicScorer, "Scorer", fLogicWorld, false, 0, fCheckOverlaps);
  fPhysGlass      = new G4PVPlacement(0,    G4ThreeVector(0.*mm, 0.*mm, 150.*mm/2 - 1.*mm/2) ,    fLogicGlass,    "Glass", fLogicScorer, false, 0, fCheckOverlaps);
  fPhysCathode    = new G4PVPlacement(0,    G4ThreeVector(0.*mm, 0.*mm, 150.*mm/2 - 1.*mm - 10.*mm/2)  ,  fLogicCathode,  "Cathode", fLogicScorer, false, 0, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACTestbenchConstruction::ConstructInnerAir()
{
  fLogicIn        = new G4LogicalVolume(fSolidGuideIn, fMatAir,"In");
  fPhysIn         = new G4PVPlacement(0,    *fGuideInPosVec, fLogicIn,   "In", fLogicWorld, false, 0, fCheckOverlaps);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void G4ACTestbenchConstruction::ConstructAerogel()
{
  #ifdef AC_TESTBENCH1
    ConstructAerogel115a();
    ConstructAerogel115b();
  #elif AC_TESTBENCH2
    ConstructAerogel115a();
  #elif AC_TESTBENCH3
    ConstructAerogel115b();
  #endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACTestbenchConstruction::ConstructAerogel115a()
{
  std::vector<G4TwoVector> aerogel1Pt(8);
  aerogel1Pt[0] = G4TwoVector(-fAerogel11X + fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel1Pt[1] = G4TwoVector(-fAerogel12X + fAerogelGap, fBoxThickness + 2*fAerogel1Y - fAerogelGap);
  aerogel1Pt[2] = G4TwoVector(fAerogel12X - fAerogelGap, fBoxThickness + 2*fAerogel1Y - fAerogelGap);
  aerogel1Pt[3] = G4TwoVector(fAerogel11X - fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel1Pt[4] = G4TwoVector(-fAerogel11X + fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel1Pt[5] = G4TwoVector(-fAerogel12X + fAerogelGap, fBoxThickness + 2*fAerogel1Y - fAerogelGap);
  aerogel1Pt[6] = G4TwoVector(fAerogel12X - fAerogelGap, fBoxThickness + 2*fAerogel1Y - fAerogelGap);
  aerogel1Pt[7] = G4TwoVector(fAerogel11X - fAerogelGap, fBoxThickness + fAerogelGap);
  fAerogelPosVec1 = new G4ThreeVector(0.*mm, 0.*mm, 0.*mm);
  fSolidAerogel1 = new G4GenericTrap("Aerogel1", fAerogel1Z - fAerogelGap, aerogel1Pt);
  fLogicAerogel1 = new G4LogicalVolume(fSolidAerogel1, fMatAerogel, "Aerogel1");
  fPhysAerogel1 = new G4PVPlacement(0, *fAerogelPosVec1, fLogicAerogel1, "Aerogel1", fLogicIn, false, 0, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACTestbenchConstruction::ConstructAerogel115b()
{
  std::vector<G4TwoVector> aerogel2Pt(8);
  aerogel2Pt[0] = G4TwoVector(-fAerogel21X + fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel2Pt[1] = G4TwoVector(-fAerogel22X + fAerogelGap, fBoxThickness + 2*fAerogel2Y - fAerogelGap);
  aerogel2Pt[2] = G4TwoVector(fAerogel22X - fAerogelGap, fBoxThickness + 2*fAerogel2Y - fAerogelGap);
  aerogel2Pt[3] = G4TwoVector(fAerogel21X - fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel2Pt[4] = G4TwoVector(-fAerogel21X + fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel2Pt[5] = G4TwoVector(-fAerogel22X + fAerogelGap, fBoxThickness + 2*fAerogel2Y - fAerogelGap);
  aerogel2Pt[6] = G4TwoVector(fAerogel22X - fAerogelGap, fBoxThickness + 2*fAerogel2Y - fAerogelGap);
  aerogel2Pt[7] = G4TwoVector(fAerogel21X - fAerogelGap, fBoxThickness + fAerogelGap);
  fAerogelPosVec2 = new G4ThreeVector(0.*mm, 0.*mm, fAerogel1Z + fAerogel2Z + 2.*fBoxThickness);
  fSolidAerogel2 = new G4GenericTrap("Aerogel2", fAerogel2Z - fAerogelGap, aerogel2Pt);
  fLogicAerogel2 = new G4LogicalVolume(fSolidAerogel2, fMatAerogel, "Aerogel2");
  fPhysAerogel2 = new G4PVPlacement(0, *fAerogelPosVec2, fLogicAerogel2, "Aerogel2", fLogicIn, false, 0, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACTestbenchConstruction::ConstructGuideBoundary()
{
  // photon energy definition
  const G4int nEntries = 20;
  const G4double pEmin = 1.9074*eV;
  const G4double pEmax = 6.1992*eV;
  G4double pE[nEntries] = {};
  // Reflectivity is included in DAVIS LUT model
  for(G4int j = 0; j < nEntries;j++)
  {
    G4double dE = (pEmax - pEmin)/(nEntries- 1);
    pE[j] = pEmin + j*dE;
  }

  G4double reflectivity[nEntries] = {};
  for(G4int j = 0;j < nEntries;j++)
  {
    reflectivity[j]= fSplineAlMylarRef->Eval(pE[j]);
  }
  G4MaterialPropertiesTable *MPTsurface = new G4MaterialPropertiesTable();
  MPTsurface->AddProperty("REFLECTIVITY", pE, reflectivity, nEntries);
  
  
  G4OpticalSurface *opInGuideCaseSurface = new G4OpticalSurface("InGuideCaseSurface");
  opInGuideCaseSurface->SetType(dielectric_metal);
  opInGuideCaseSurface->SetFinish(polished); // mechanically grounded surface, with esr film & meltmount
  opInGuideCaseSurface->SetModel(unified);
  
  G4LogicalBorderSurface *logicSurfaceInGuideCase = new G4LogicalBorderSurface("InGuideCaseSurface",
    fPhysIn, fPhysGuideReflector, opInGuideCaseSurface);  

  opInGuideCaseSurface = dynamic_cast<G4OpticalSurface*>
    (logicSurfaceInGuideCase->GetSurface(fPhysIn, fPhysGuideReflector));
  if(opInGuideCaseSurface) opInGuideCaseSurface->DumpInfo();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACTestbenchConstruction::ConstructAerogelBoundary()
{
  G4OpticalSurface *opAerogel1AirSurface = new G4OpticalSurface("Aerogel1AirSurface");
  opAerogel1AirSurface->SetType(dielectric_dielectric);
  opAerogel1AirSurface->SetFinish(polished);
  opAerogel1AirSurface->SetModel(unified);
  G4OpticalSurface *opAerogel2AirSurface = new G4OpticalSurface("Aerogel2AirSurface");
  opAerogel2AirSurface->SetType(dielectric_dielectric);
  opAerogel2AirSurface->SetFinish(polished);
  opAerogel2AirSurface->SetModel(unified);
  new G4LogicalSkinSurface("LogicalAerogel1AirSurface", fLogicAerogel1, opAerogel1AirSurface);
  new G4LogicalSkinSurface("LogicalAerogel2AirSurface", fLogicAerogel2, opAerogel2AirSurface);
}

