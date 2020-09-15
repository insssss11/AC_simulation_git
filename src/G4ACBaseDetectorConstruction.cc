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
// $Id: G4ACBaseDetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACBaseDetectorConstruction.cc
/// \brief Implementation of the G4ACBaseDetectorConstruction class



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

#include "G4ACBaseDetectorConstruction.hh"
#include "G4ACFineMeshSD.hh"


extern const G4double TARGET_LENGTH, TARGET_CENTER_POSITION,
  AC_DETECTOR_R, TPC_EDGE_R, TPC_DOWNSTREAM_Z;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ACBaseDetectorConstruction::G4ACBaseDetectorConstruction()
: G4VUserDetectorConstruction(),
  fCrossX(144.0*mm/2), fCrossY(90.*mm/2),
  
  fAerogelThickness(20.*mm/2),
  fAerogel1Z(AC_DETECTOR_R/2),
  fAerogel11X(91.0*mm/2), fAerogel11Y(fAerogelThickness),
  fAerogel12X(144.*mm/2), fAerogel12Y(fAerogelThickness),
  fAerogel2Z(238.85*mm/2),
  fAerogel21X(95.1*mm/2), fAerogel21Y(fAerogelThickness),
  fAerogel22X(144.*mm/2), fAerogel22Y(fAerogelThickness),
  fAerogel3Z(221.9*mm/2),
  fAerogel31X(98.6*mm/2), fAerogel31Y(fAerogelThickness),
  fAerogel32X(144.*mm/2), fAerogel32Y(fAerogelThickness),

  fBoxThickness(0.5*mm), fReflectorThickness(65.*um), fAngleTilt(30.*deg), fAerogelGap(0.05*mm), fPMTradius(78.*mm/2), fCathodeRadius(64.*mm/2), fGlassThickness(3.*mm), fScorerZ(123.7*mm/2), fCathodeZ(5. *mm),
  fCheckOverlaps(true),
  fSolidWorld(0),
  fSolidHolderCase(0), fSolidHolderOut(0), fSolidHolderIn(0),
  fSolidGuideCase(0), fSolidGuideOut(0), fSolidGuideIn(0),
  fSolidOut(0), fSolidIn(0),
  fSolidScorer(0), fSolidGlass(0), fSolidCathode(0),
  fSolidAerogel1(0), fSolidAerogel2(0), fSolidAerogel3(0),
  fSolidHolderReflector(0), fSolidGuideReflector(0),
  fLogicWorld(0),
  fLogicHolderCase(0),
  fLogicGuideCase(0),
  fLogicIn(0),
  fLogicScorer(0), fLogicGlass(0), fLogicCathode(0),
  fLogicAerogel1(0),  fLogicAerogel2(0),  fLogicAerogel3(0),
  fLogicHolderReflector(0), fLogicGuideReflector(0),
  fPhysWorld(0),
  fPhysHolderCase(0),
  fPhysGuideCase(0),
  fPhysIn(0),
  fPhysScorer(0), fPhysGlass(0), fPhysCathode(0),
  fPhysAerogel1(0), fPhysAerogel2(0), fPhysAerogel3(0),
  fPhysHolderReflector(0), fPhysGuideReflector(0),
  fHolderOutPosVec(0), fHolderInPosVec(0),
  fGuideOutPosVec(0), fGuideInPosVec(0),
  fScorerPosVec(0),
  fAerogelPosVec1(0), fAerogelPosVec2(0), fAerogelPosVec3(0),
  fHolderReflectorPosVec(0), fGuideReflectorPosVec(0),

  fRotationScorer(0),
  fMatAir(0), fMatBox(0), fMatGlass(0), fMatScorer(0), fMatAerogel(0), fMatReflector(0),
  fSplineAeroAbs(0), fSplineESRref(0), fSplineAlMylarRef(0), fSplineTyvekRef(0)
{
  ConstructSpline(); // data spline definitions
  ConstructMaterial(); // material definitions
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ACBaseDetectorConstruction::G4ACBaseDetectorConstruction(G4double BoxThickness)
: G4VUserDetectorConstruction(),
  fCrossX(140.0*mm/2), fCrossY(90.*mm/2),
  fAerogelThickness(20.*mm/2),
  fAerogel1Z(AC_DETECTOR_R/2),
  fAerogel11X(91.5*mm/2), fAerogel11Y(fAerogelThickness),
  fAerogel12X(144.4*mm/2), fAerogel12Y(fAerogelThickness),
  fAerogel2Z(238.85*mm/2),
  fAerogel21X(95.1*mm/2), fAerogel21Y(fAerogelThickness),
  fAerogel22X(144.4*mm/2), fAerogel22Y(fAerogelThickness),
  fAerogel3Z(221.9*mm/2),
  fAerogel31X(98.6*mm/2), fAerogel31Y(fAerogelThickness),
  fAerogel32X(144.4*mm/2), fAerogel32Y(fAerogelThickness),
  fBoxThickness(BoxThickness), fReflectorThickness(65.*um), fAngleTilt(30.*deg), fAerogelGap(0.05*mm), fPMTradius(78.*mm/2), fCathodeRadius(64.*mm/2), fGlassThickness(3.*mm), fScorerZ(123.7*mm/2), fCathodeZ(5. *mm),
  fCheckOverlaps(true),
  fSolidWorld(0),
  fSolidHolderCase(0), fSolidHolderOut(0), fSolidHolderIn(0),
  fSolidGuideCase(0), fSolidGuideOut(0), fSolidGuideIn(0),
  fSolidOut(0), fSolidIn(0),
  fSolidScorer(0), fSolidGlass(0), fSolidCathode(0),
  fSolidAerogel1(0), fSolidAerogel2(0), fSolidAerogel3(0),
  fSolidHolderReflector(0), fSolidGuideReflector(0),
  fLogicWorld(0),
  fLogicHolderCase(0),
  fLogicGuideCase(0),
  fLogicIn(0),
  fLogicScorer(0), fLogicGlass(0), fLogicCathode(0),
  fLogicAerogel1(0),  fLogicAerogel2(0),  fLogicAerogel3(0),
  fLogicHolderReflector(0), fLogicGuideReflector(0),
  fPhysWorld(0),
  fPhysHolderCase(0), fPhysGuideCase(0),
  fPhysIn(0),
  fPhysScorer(0), fPhysGlass(0), fPhysCathode(0),
  fPhysAerogel1(0), fPhysAerogel2(0), fPhysAerogel3(0),
  fPhysHolderReflector(0), fPhysGuideReflector(0),
  fHolderOutPosVec(0), fHolderInPosVec(0),
  fGuideOutPosVec(0), fGuideInPosVec(0),
  fScorerPosVec(0),
  fAerogelPosVec1(0), fAerogelPosVec2(0), fAerogelPosVec3(0),
  fHolderReflectorPosVec(0), fGuideReflectorPosVec(0),
  fRotationScorer(0),
  fMatAir(0), fMatBox(0), fMatGlass(0), fMatScorer(0), fMatAerogel(0), fMatReflector(0),
  fSplineAeroAbs(0), fSplineESRref(0), fSplineAlMylarRef(0), fSplineTyvekRef(0)
{
  ConstructSpline(); // data spline definitions
  ConstructMaterial(); // material definitions
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ACBaseDetectorConstruction::G4ACBaseDetectorConstruction(G4double crossX, G4double crossY, G4double aerogelThickness, G4double aerogelGap, G4double angleTilt)
: G4VUserDetectorConstruction(),
  fCrossX(crossX), fCrossY(crossY),
  
  fAerogelThickness(aerogelThickness),
  fAerogel1Z(256.3*mm/2),
  fAerogel11X(91.5*mm/2), fAerogel11Y(fAerogelThickness),
  fAerogel12X(144.4*mm/2), fAerogel12Y(fAerogelThickness),
  fAerogel2Z(238.85*mm/2),
  fAerogel21X(95.1*mm/2), fAerogel21Y(fAerogelThickness),
  fAerogel22X(144.4*mm/2), fAerogel22Y(fAerogelThickness),
  fAerogel3Z(221.9*mm/2),
  fAerogel31X(98.6*mm/2), fAerogel31Y(fAerogelThickness),
  fAerogel32X(144.4*mm/2), fAerogel32Y(fAerogelThickness),

  fBoxThickness(0.5*mm), fReflectorThickness(65.*um), fAngleTilt(angleTilt), fAerogelGap(aerogelGap), fPMTradius(78.*mm/2), fCathodeRadius(64.*mm/2), fGlassThickness(3.*mm), fScorerZ(123.7*mm/2), fCathodeZ(5. *mm),
  fCheckOverlaps(true),
  fSolidWorld(0),
  fSolidHolderCase(0), fSolidHolderOut(0), fSolidHolderIn(0),
  fSolidGuideCase(0), fSolidGuideOut(0), fSolidGuideIn(0),
  fSolidOut(0), fSolidIn(0),
  fSolidScorer(0), fSolidGlass(0), fSolidCathode(0),
  fSolidAerogel1(0), fSolidAerogel2(0), fSolidAerogel3(0),
  fSolidHolderReflector(0), fSolidGuideReflector(0),
  fLogicWorld(0),
  fLogicHolderCase(0),
  fLogicGuideCase(0),
  fLogicIn(0),
  fLogicScorer(0), fLogicGlass(0), fLogicCathode(0),
  fLogicAerogel1(0),  fLogicAerogel2(0),  fLogicAerogel3(0),
  fLogicHolderReflector(0), fLogicGuideReflector(0),
  fPhysWorld(0),
  fPhysHolderCase(0),
  fPhysGuideCase(0),
  fPhysIn(0),
  fPhysScorer(0), fPhysGlass(0), fPhysCathode(0),
  fPhysAerogel1(0), fPhysAerogel2(0), fPhysAerogel3(0),
  fPhysHolderReflector(0), fPhysGuideReflector(0),
  fHolderOutPosVec(0), fHolderInPosVec(0),
  fGuideOutPosVec(0), fGuideInPosVec(0),
  fScorerPosVec(0),
  fAerogelPosVec1(0), fAerogelPosVec2(0), fAerogelPosVec3(0),
  fHolderReflectorPosVec(0), fGuideReflectorPosVec(0),
  fRotationScorer(0),

  fMatAir(0), fMatBox(0), fMatGlass(0), fMatScorer(0), fMatAerogel(0), fMatReflector(0),
  fSplineAeroAbs(0), fSplineESRref(0), fSplineAlMylarRef(0), fSplineTyvekRef(0)
{
  ConstructSpline(); // data spline definitions
  ConstructMaterial(); // material definitions
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACBaseDetectorConstruction::~G4ACBaseDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACBaseDetectorConstruction::Construct()
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

void G4ACBaseDetectorConstruction::ConstructMaterial()
{
  G4NistManager *nist = G4NistManager::Instance();

  /////////////////////////////////////////////////////////////////
  // array of energy of photon to defining Material
  /////////////////////////////////////////////////////////////////

  // photon energy definition
  const G4int nEntries = 20;
  const G4double pEmin = 1.9074*eV; // 650 nm
  const G4double pEmax = 6.1992*eV; // 200 nm
  G4double pE[nEntries] = {};
  for(G4int j = 0; j < nEntries;j++)
  {
    G4double dE = (pEmax - pEmin)/(nEntries- 1);
    pE[j] = pEmin + j*dE;
  }
  /////////////////////////////////////////////////////////////////
  // Aerogel Material
  ///////////////////////////////////////////////////////////////// 
  G4double rIndexAerogel = 1.05;
  G4double densityAerogel = (rIndexAerogel - 1.)/0.19;

  // define wavelength dependent properties(or maybe later.....)
  G4double scatAerogelArray[nEntries] = {};
  G4double absAerogelArray[nEntries] = {};
  G4double rIndexAerogelArray[nEntries] = {};

  G4Material *elC = nist->FindOrBuildMaterial("G4_C");
  G4Material *matSiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material *matH2O  = nist->FindOrBuildMaterial("G4_WATER");
  
  // finally define aerogels
	for(G4int j = 0;j < nEntries;j++)
	{
		absAerogelArray[j] = GetAerogelAbs(pE[j]);
		scatAerogelArray[j] = GetAerogelScat(pE[j]);
		rIndexAerogelArray[j] = rIndexAerogel;
  }
  // matarial property table of aerogel
  G4MaterialPropertiesTable *MPTaerogel = new G4MaterialPropertiesTable();
  MPTaerogel->AddProperty("RINDEX", pE, rIndexAerogelArray, nEntries);
  MPTaerogel->AddProperty("ABSLENGTH", pE, absAerogelArray, nEntries);
  MPTaerogel->AddProperty("RAYLEIGH", pE, scatAerogelArray, nEntries);

  fMatAerogel = new G4Material("Aerogel", densityAerogel, 3);
  fMatAerogel->AddMaterial(elC, 0.1*perCent);
  fMatAerogel->AddMaterial(matSiO2, 62.5*perCent);
  fMatAerogel->AddMaterial(matH2O, 37.4*perCent);
  fMatAerogel->SetMaterialPropertiesTable(MPTaerogel);
  /////////////////////////////////////////////////////////////////
  // Air Material
  /////////////////////////////////////////////////////////////////   
  G4MaterialPropertiesTable *MPTair = new G4MaterialPropertiesTable();

  fMatAir  = nist->FindOrBuildMaterial("G4_AIR");
  
  // refractive index
  G4double rIndexAir[nEntries] = {};
  for(G4int j = 0;j < nEntries;j++)
  {
    rIndexAir[j] = 1.0;
  }

  MPTair = new G4MaterialPropertiesTable();
  MPTair->AddProperty("RINDEX", pE, rIndexAir, nEntries);
  fMatAir->SetMaterialPropertiesTable(MPTair);

  /////////////////////////////////////////////////////////////////
  // Detector Box Material (polypropylene)
  /////////////////////////////////////////////////////////////////  
  G4double densityBox = 0.9*g/cm3;
  //refractive index
  G4double rIndexBox[nEntries] = {};
  G4double absBoxArray[nEntries] = {};
  for(G4int j = 0;j < nEntries;j++)
  {
    rIndexBox[j] = 1.4900;
    absBoxArray[j] = 1.*nm;
  }

  fMatBox = new G4Material("Box", densityBox, 1);;
  G4MaterialPropertiesTable *MPTBox = new G4MaterialPropertiesTable();

  G4Material *matPP = nist->FindOrBuildMaterial("G4_POLYPROPYLENE");
  MPTBox->AddProperty("RINDEX", pE, rIndexBox, nEntries);
  MPTBox->AddProperty("ABSLENGTH", pE, absBoxArray, nEntries);
  fMatBox->SetMaterialPropertiesTable(MPTBox);
  fMatBox->AddMaterial(matPP, 100.*perCent);
  /////////////////////////////////////////////////////////////////
  // Reflector Material(ESR)
  ///////////////////////////////////////////////////////////////// 
  G4double pEreflector[] = {
    2.069010169*eV, 2.108284284*eV, 2.144795492*eV, 2.176688782*eV, 2.22324978*eV,
    2.267069056*eV, 2.306523886*eV, 2.353199393*eV, 2.401784385*eV, 2.450033554*eV,
    2.514933163*eV, 2.56011235*eV,  2.604277714*eV, 2.6121898*eV,   2.658169296*eV,
    2.700068882*eV, 2.740369754*eV, 2.79102241*eV,  2.846801716*eV, 2.894980375*eV,
    2.916012667*eV, 2.941890097*eV, 2.986894366*eV, 3.029679815*eV, 3.066373008*eV,
    3.088824826*eV, 3.119350312*eV, 3.162184455*eV, 3.206136733*eV, 3.272127761*eV,
    3.325292422*eV, 3.363605838*eV, 3.417199331*eV, 3.479994302*eV, 3.5451302*eV,
    3.620875049*eV, 3.693237327*eV, 3.797418437*eV
  };

  G4double scintReflector[] = {
    0.002606491, 0.003237212, 0.003218877, 0.003529104, 0.004810347, 
    0.005116174, 0.006582967, 0.008616675, 0.011471787, 0.015559739, 
    0.02292744,  0.029071102, 0.035626933, 0.038088211, 0.045054011, 
    0.052432713, 0.059812881, 0.066778681, 0.072101672, 0.073345513, 
    0.071710039, 0.070432462, 0.062207421, 0.054393081, 0.045760272, 
    0.040416012, 0.031373234, 0.022328255, 0.014926085, 0.007930216, 
    0.004792012, 0.002502349, 0.00151373,  0.001175635, 0.0011639, 
    0.001149966, 0.002511149, 0.003725654
  };
  assert(sizeof(pEreflector)==sizeof(pEreflector));
  const G4int nEntriesReflector = sizeof(pEreflector)/sizeof(G4double);
  G4double densityReflector = 1.29*g/cm3;
  G4double rIndexReflector[nEntriesReflector] = {};
  G4double absReflector[nEntriesReflector] = {};
  for(G4int j = 0;j < nEntriesReflector;j++)
  {
    rIndexReflector[j] = 1.4900;
    absReflector[j] = 35.*cm;
  }

  fMatReflector = new G4Material("Reflector", densityReflector, 1);;
  G4MaterialPropertiesTable *MPTreflector = new G4MaterialPropertiesTable();
  MPTreflector->AddProperty("RINDEX", pEreflector, rIndexReflector, nEntriesReflector);
  MPTreflector->AddProperty("ABSLENGTH", pEreflector, absReflector, nEntriesReflector);
  MPTreflector->AddProperty("FASTCOMPONENT", pEreflector, scintReflector, nEntriesReflector);
  MPTreflector->AddConstProperty("SCINTILLATIONYIELD", 5000./MeV);
  MPTreflector->AddConstProperty("YIELDRATIO", 1.);
  MPTreflector->AddConstProperty("SLOWTIMECONSTANT", 10.*ns);
  MPTreflector->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  MPTreflector->AddConstProperty("RESOLUTIONSCALE", 1.0);

  fMatReflector->SetMaterialPropertiesTable(MPTreflector);
  G4Material *matReflector = nist->FindOrBuildMaterial("G4_POLYPROPYLENE");
  fMatReflector->AddMaterial(matReflector, 100.*perCent);
  
  /////////////////////////////////////////////////////////////////
  // Glass Material(Sylgard 527)
  /////////////////////////////////////////////////////////////////
  G4double densityGlass = 0.95*g/cm3;

  G4double rIndexGlass[nEntries] = {};
  for(G4int j = 0;j < nEntries;j++)
  {
    rIndexGlass[j] = 1.4074;
  }
  
  G4MaterialPropertiesTable *MPTGlass = new G4MaterialPropertiesTable();
  
  fMatGlass = new G4Material("Glass", densityGlass, 3);
  fMatGlass->AddMaterial(elC, 0.1*perCent);
  fMatGlass->AddMaterial(matSiO2, 62.5*perCent);
  fMatGlass->AddMaterial(matH2O, 37.4*perCent);

  MPTGlass->AddProperty("RINDEX", pE, rIndexGlass, nEntries);

  fMatGlass->SetMaterialPropertiesTable(MPTGlass);

  /////////////////////////////////////////////////////////////////
  // FMPMT Box Material(case of PMT)
  /////////////////////////////////////////////////////////////////
  fMatScorer  = nist->FindOrBuildMaterial("G4_AIR");
  // the refractive index is not set to kill photon
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACBaseDetectorConstruction::ConstructHolder()
{
  /////////////////////////////////////////////////////////////////
  // Aerogel Holder Part 1 In
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
  // Aerogel Holder Part 1 Out
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

  /////////////////////////////////////////////////////////////////
  // Aerogel Holder Part 1 Reflector
  /////////////////////////////////////////////////////////////////
  std::vector<G4TwoVector> holderReflectorOutPtVec(8);
  G4double holderReflectorOutZ = fAerogel1Z + fReflectorThickness/2;
  holderReflectorOutPtVec[0] = holderInPtVec[0] + G4TwoVector(-fReflectorThickness,  -fReflectorThickness);
  holderReflectorOutPtVec[1] = holderInPtVec[1] + G4TwoVector(-fReflectorThickness,  fReflectorThickness);
  holderReflectorOutPtVec[2] = holderInPtVec[2] + G4TwoVector(fReflectorThickness,  fReflectorThickness);
  holderReflectorOutPtVec[3] = holderInPtVec[3] + G4TwoVector(fReflectorThickness,  -fReflectorThickness);
  holderReflectorOutPtVec[4] = holderInPtVec[4] + G4TwoVector(-fReflectorThickness,  -fReflectorThickness);
  holderReflectorOutPtVec[5] = holderInPtVec[5] + G4TwoVector(-fReflectorThickness,  fReflectorThickness);
  holderReflectorOutPtVec[6] = holderInPtVec[6] + G4TwoVector(fReflectorThickness,  fReflectorThickness);
  holderReflectorOutPtVec[7] = holderInPtVec[7] + G4TwoVector(fReflectorThickness,  -fReflectorThickness);
  fHolderReflectorPosVec = new G4ThreeVector(0.*mm, 0.*mm, -holderReflectorOutZ);

  // position vector of the centers of masses(center of mass)
  fSolidHolderOut = new G4GenericTrap("HolderOut", holderOutZ, holderOutPtVec);
  fHolderOutPosVec = new G4ThreeVector(0.*mm, 0.*mm, -holderOutZ);
  G4VSolid *solidHolderReflectorOut = new G4GenericTrap("HolderReflectorOut", holderReflectorOutZ, holderReflectorOutPtVec);
  
  fSolidHolderCase  = new G4SubtractionSolid("HolderCase", fSolidHolderOut, solidHolderReflectorOut, 0, *fHolderReflectorPosVec - *fHolderOutPosVec);
  fSolidHolderReflector = new G4SubtractionSolid("HolderReflector", solidHolderReflectorOut, fSolidHolderIn, 0, *fHolderInPosVec - *fHolderReflectorPosVec);  
  
  fLogicHolderCase = new G4LogicalVolume(fSolidHolderCase, fMatBox, "HolderCase");
  fLogicHolderReflector = new G4LogicalVolume(fSolidHolderReflector, fMatReflector, "HolderReflector");

  fPhysHolderCase = new G4PVPlacement(0,    *fHolderOutPosVec,  fLogicHolderCase,  "HolderCase",  fLogicWorld, false, 0, fCheckOverlaps);
  fPhysHolderReflector = new G4PVPlacement(0, *fHolderReflectorPosVec, fLogicHolderReflector, "HolderReflector", fLogicWorld, false, 0, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACBaseDetectorConstruction::ConstructGuide()
{
  /////////////////////////////////////////////////////////////////
  // Light Guide Part 1 In
  /////////////////////////////////////////////////////////////////
  G4double guide1InZ = 200.*mm/2;
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
  guide1InPtVec[4] = G4TwoVector(-guide1In2X, fBoxThickness + 2*(guide1In1Y - guide1In2Y));
  guide1InPtVec[5] = G4TwoVector(-guide1In2X, fBoxThickness + 2*guide1In1Y);
  guide1InPtVec[6] = G4TwoVector(guide1In2X, fBoxThickness + 2*guide1In1Y);
  guide1InPtVec[7] = G4TwoVector(guide1In2X, fBoxThickness + 2*(guide1In1Y - guide1In2Y));


  // to make the centers of masses locate at Origin
  G4VSolid *solidGuide1In = new G4GenericTrap("Guide1In", guide1InZ, guide1InPtVec);
  fGuideInPosVec = new G4ThreeVector(0.*mm, 0.*mm, guide1InZ);

  /////////////////////////////////////////////////////////////////
  // Light Guide Part 1 Out
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
  guide1OutPtVec[4] = G4TwoVector(-guide1Out2X, 2*(guide1Out1Y - guide1Out2Y));
  guide1OutPtVec[5] = G4TwoVector(-guide1Out2X, 2*guide1Out1Y);
  guide1OutPtVec[6] = G4TwoVector(guide1Out2X, 2*guide1Out1Y);
  guide1OutPtVec[7] = G4TwoVector(guide1Out2X, 2*(guide1Out1Y - guide1Out2Y));

  // to make the centers of masses locate at Origin
  G4VSolid *solidGuide1Out = new G4GenericTrap("Guide1Out", guide1OutZ, guide1OutPtVec);
  fGuideOutPosVec = fGuideInPosVec;
  
  /////////////////////////////////////////////////////////////////
  // Light Guide Part 1 Reflector
  /////////////////////////////////////////////////////////////////  
  G4double guideReflector1OutZ = guide1InZ;
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
  G4VSolid *solidGuideReflector1Out = new G4GenericTrap("GuideReflector1Out", guideReflector1OutZ, guideReflector1OutPtVec);
  fGuideReflectorPosVec = fGuideInPosVec;

  /////////////////////////////////////////////////////////////////
  // Light Guide Part2 Out(tub)
  /////////////////////////////////////////////////////////////////
  G4RotationMatrix *rotationGuide2 = new G4RotationMatrix();
  rotationGuide2->rotateY(90.*deg);
  G4ThreeVector guide2PosVec(0.*mm, 2*(guide1Out1Y - guide1Out2Y), 2.*guide1OutZ);

  G4double guide2OutRmin = 0, guide2OutRmax = 2.*guide1Out2X;
  G4double guide2OutZ = guide1Out2X;
  G4double guide2OutSphi = 0.*deg, guide2OutDphi = 90.*deg;

  G4Tubs *solidGuide2Out = new G4Tubs("Guide2Out", guide2OutRmin, guide2OutRmax, guide2OutZ, guide2OutSphi, guide2OutDphi);

  /////////////////////////////////////////////////////////////////
  // Light Guide Part2 In(tub)
  /////////////////////////////////////////////////////////////////
  G4double guide2InRmin = fBoxThickness, guide2InRmax = guide2OutRmax - fBoxThickness;
  G4double guide2InZ = guide2OutZ - fBoxThickness;
  G4double guide2InSphi = 0.*deg, guide2InDphi = 90.*deg;

  G4Tubs *solidGuide2In = new G4Tubs("Guide2In", guide2InRmin, guide2InRmax, guide2InZ, guide2InSphi, guide2InDphi);

  /////////////////////////////////////////////////////////////////
  // Light Guide Part2 Reflector
  /////////////////////////////////////////////////////////////////
  G4Tubs *solidGuideReflector2Out = new G4Tubs("GuideReflector2Out", guide2InRmin, guide2InRmax + fReflectorThickness, guide2InZ + fReflectorThickness, guide2InSphi, guide2InDphi);

  /////////////////////////////////////////////////////////////////
  // Light Guide Part3(tilt) Out
  /////////////////////////////////////////////////////////////////
  // samkak-kidung

  std::vector <G4TwoVector> tiltOutPolygon(3);
  G4double tiltOutX = guide1Out2X, tiltOutY = tiltOutX*tan(fAngleTilt), tiltOutZ = guide1Out2Y;
  G4ThreeVector tiltOutPosVec(tiltOutX, 2*(guide1Out1Y - guide1Out2Y), 2.*guide1OutZ + tiltOutZ);

  // define bottom face before rotation 90 deg with X axis
  tiltOutPolygon[0] = G4TwoVector(-2*tiltOutX, 0.*mm);
  tiltOutPolygon[1] = G4TwoVector(0.*mm,0.*mm);
  tiltOutPolygon[2] = G4TwoVector(0.*mm, -2.*tiltOutY);

  // before lotation with x axis
  G4ExtrudedSolid *solidTiltOut = new G4ExtrudedSolid("TiltOut", tiltOutPolygon, tiltOutZ, 
    G4TwoVector(), 1., G4TwoVector(), 1.);

  /////////////////////////////////////////////////////////////////
  // Light Guide Part3(tilt) In
  /////////////////////////////////////////////////////////////////
  // sakak-kidung
  std::vector <G4TwoVector> tiltInPolygon(3);
  tiltInPolygon[0] = G4TwoVector(-2.*tiltOutX + fBoxThickness, 0.*mm);
  tiltInPolygon[1] = G4TwoVector(-fBoxThickness, 0.*mm);
  tiltInPolygon[2] = G4TwoVector(-fBoxThickness, 2.*(-tiltOutX + fBoxThickness)*tan(fAngleTilt));
  G4ExtrudedSolid *solidTiltIn = new G4ExtrudedSolid("TiltIn", tiltInPolygon, tiltOutZ - fBoxThickness, 
    G4TwoVector(), 1., G4TwoVector(), 1.);

  /////////////////////////////////////////////////////////////////
  // Light Guide Part3(tilt) Reflector
  /////////////////////////////////////////////////////////////////
  // sakak-kidung
  std::vector <G4TwoVector> guideReflectorOutPolygon(3);
  guideReflectorOutPolygon[0] = tiltInPolygon[0] + G4TwoVector(-fReflectorThickness, 0.*mm);
  guideReflectorOutPolygon[1] = tiltInPolygon[1] + G4TwoVector(fReflectorThickness, 0.*mm);
  guideReflectorOutPolygon[2] = tiltInPolygon[2] + G4TwoVector(fReflectorThickness, -2.*fReflectorThickness*tan(fAngleTilt));
  G4ExtrudedSolid *solidGuideReflector3Out = new G4ExtrudedSolid("GuideReflector3Out", guideReflectorOutPolygon, tiltOutZ - fBoxThickness + fReflectorThickness, 
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

  G4UnionSolid *solidGuideTubOut     = new G4UnionSolid("GuideTubOut", solidGuide1Out, solidGuide2Out, rotationGuide2, guide2PosVec - *fGuideOutPosVec);
  G4UnionSolid *solidGuideTubIn      = new G4UnionSolid("GuideTubIn", solidGuide1In, solidGuide2In, rotationGuide2, guide2PosVec - *fGuideInPosVec);
  G4UnionSolid *solidGuideReflector12Out = new G4UnionSolid("GuideReflector12Out", solidGuideReflector1Out, solidGuideReflector2Out, rotationGuide2, guide2PosVec - *fGuideReflectorPosVec);

  fSolidGuideOut  = new G4UnionSolid("GuideOut", solidGuideTubOut, solidTiltOut, 0, tiltOutPosVec - *fGuideOutPosVec);
  fSolidGuideIn   = new G4UnionSolid("GuideIn", solidGuideTubIn, solidTiltInHole, 0, tiltOutPosVec- *fGuideInPosVec);
  G4UnionSolid *solidGuideReflectorOutWithoutHole = new G4UnionSolid("GuideReflectorOutWithoutHole", solidGuideReflector12Out, solidGuideReflector3Out,0,  tiltOutPosVec - *fGuideReflectorPosVec);
  G4SubtractionSolid *solidGuideCaseWithoutHole = new G4SubtractionSolid("GuideCaseWithoutHole", fSolidGuideOut, solidGuideReflectorOutWithoutHole, 0, *fGuideReflectorPosVec - *fGuideOutPosVec);
  fSolidGuideCase = new G4SubtractionSolid("GuideCase", solidGuideCaseWithoutHole, fSolidGuideIn, 0, *fGuideInPosVec - *fGuideOutPosVec);
  fSolidGuideReflector = new G4SubtractionSolid("GuideReflector", solidGuideReflectorOutWithoutHole, fSolidGuideIn, 0, *fGuideReflectorPosVec - *fGuideOutPosVec);

  fLogicGuideCase = new G4LogicalVolume(fSolidGuideCase, fMatBox, "GuideCase");
  fLogicGuideReflector = new G4LogicalVolume(fSolidGuideReflector, fMatReflector, "GuideReflector");

  fPhysGuideCase  = new G4PVPlacement(0,    *fGuideOutPosVec,  fLogicGuideCase,  "GuideCase",  fLogicWorld, false, 0, fCheckOverlaps);
  fPhysGuideReflector  = new G4PVPlacement(0,    *fGuideReflectorPosVec,  fLogicGuideReflector,  "GuideReflector",  fLogicWorld, false, 0, fCheckOverlaps);

  fScorerPosVec = new G4ThreeVector(tiltOutPosVec + G4ThreeVector(-tiltOutX, -tiltOutY, 0.*mm) - fScorerZ*G4ThreeVector(sin(fAngleTilt), cos(fAngleTilt), 0.*mm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACBaseDetectorConstruction::ConstructScorer()
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

void G4ACBaseDetectorConstruction::ConstructInnerAir()
{
  fSolidIn        = new G4UnionSolid("In", fSolidHolderIn, fSolidGuideIn, 0, *fGuideInPosVec - *fHolderInPosVec);
  fLogicIn        = new G4LogicalVolume(fSolidIn, fMatAir,"In");
  fPhysIn         = new G4PVPlacement(0,    *fHolderInPosVec, fLogicIn,   "In", fLogicWorld, false, 0, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACBaseDetectorConstruction::ConstructAerogel()
{
  /////////////////////////////////////////////////////////////////
  // aerogel1
  /////////////////////////////////////////////////////////////////
  // The # of aerogels are six pieces and there are all different refractive indice
  // The aerogels numbering and position
  // -------------------------------------------|
  // |  6              3       |                |------
  // |  5              2       |  light guide   |  pmt|
  // |  4              1       |                |------
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

  G4double aerogel31X = (fAerogel11X + fCrossX)/2; 
  aerogel1Pt[0] = G4TwoVector(-fAerogel11X + fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel1Pt[1] = G4TwoVector(-fAerogel11X + fAerogelGap, fBoxThickness + 2*fAerogelThickness - fAerogelGap);
  aerogel1Pt[2] = G4TwoVector(fAerogel11X - fAerogelGap, fBoxThickness + 2*fAerogelThickness - fAerogelGap);
  aerogel1Pt[3] = G4TwoVector(fAerogel11X - fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel1Pt[4] = G4TwoVector(-aerogel31X + fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel1Pt[5] = G4TwoVector(-aerogel31X + fAerogelGap, fBoxThickness + 2*fAerogelThickness - fAerogelGap);
  aerogel1Pt[6] = G4TwoVector(aerogel31X - fAerogelGap, fBoxThickness + 2*fAerogelThickness - fAerogelGap);
  aerogel1Pt[7] = G4TwoVector(aerogel31X - fAerogelGap, fBoxThickness + fAerogelGap);
  
  /////////////////////////////////////////////////////////////////
  // aerogel2
  /////////////////////////////////////////////////////////////////
  std::vector<G4TwoVector> aerogel2Pt(8);

  aerogel2Pt[0] = G4TwoVector(-aerogel31X + fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel2Pt[1] = G4TwoVector(-aerogel31X + fAerogelGap, fBoxThickness + 2*fAerogelThickness - fAerogelGap);
  aerogel2Pt[2] = G4TwoVector(aerogel31X - fAerogelGap, fBoxThickness + 2*fAerogelThickness - fAerogelGap);
  aerogel2Pt[3] = G4TwoVector(aerogel31X - fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel2Pt[4] = G4TwoVector(-fAerogel22X + fAerogelGap, fBoxThickness + fAerogelGap);
  aerogel2Pt[5] = G4TwoVector(-fAerogel22X + fAerogelGap, fBoxThickness + 2*fAerogelThickness - fAerogelGap);
  aerogel2Pt[6] = G4TwoVector(fAerogel22X - fAerogelGap, fBoxThickness + 2*fAerogelThickness - fAerogelGap);
  aerogel2Pt[7] = G4TwoVector(fAerogel22X - fAerogelGap, fBoxThickness + fAerogelGap);

  fSolidAerogel1 = new G4GenericTrap("Aerogel1", fAerogel1Z/2 - fAerogelGap, aerogel1Pt);
  fSolidAerogel2 = new G4GenericTrap("Aerogel2", fAerogel1Z/2 - fAerogelGap, aerogel2Pt);
  
  fLogicAerogel1 = new G4LogicalVolume(fSolidAerogel1, fMatAerogel, "Aerogel1");
  fLogicAerogel2 = new G4LogicalVolume(fSolidAerogel2, fMatAerogel, "Aerogel2");

  new G4PVPlacement(0, G4ThreeVector(0.*mm, 0.*mm, -fAerogel1Z/2), fLogicAerogel1, "Aerogel1", fLogicIn, true, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.*mm, 2.*fAerogelThickness, -fAerogel1Z/2), fLogicAerogel1, "Aerogel2", fLogicIn, true, 1, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.*mm, 4.*fAerogelThickness, -fAerogel1Z/2), fLogicAerogel1, "Aerogel3", fLogicIn, true, 2, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.*mm, 0.*mm, fAerogel1Z/2), fLogicAerogel2, "Aerogel4", fLogicIn, true, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.*mm, 2.*fAerogelThickness, fAerogel1Z/2), fLogicAerogel2, "Aerogel5", fLogicIn, true, 1, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.*mm, 4.*fAerogelThickness, fAerogel1Z/2), fLogicAerogel2, "Aerogel6", fLogicIn, true, 2, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACBaseDetectorConstruction::ConstructHolderBoundary()
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
    reflectivity[j]= fSplineESRref->Eval(pE[j]);
  }
  */
  // G4MaterialPropertiesTable *MPTsurface = new G4MaterialPropertiesTable();
  // MPTsurface->AddProperty("REFLECTIVITY", pE, reflectivity, nEntries);
  G4OpticalSurface *opInHolderCaseSurface = new G4OpticalSurface("InHolderCaseSurface");
  opInHolderCaseSurface->SetType(dielectric_LUTDAVIS);
  opInHolderCaseSurface->SetFinish(PolishedESR_LUT); // mechanically grounded surface, with esr film & meltmount
  opInHolderCaseSurface->SetModel(DAVIS);
  
  G4LogicalBorderSurface *logicSurfaceInHolderCase = new G4LogicalBorderSurface("InHolderCaseSurface",
    fPhysIn, fPhysHolderReflector, opInHolderCaseSurface);  
  // opInHolderCaseSurface->SetMaterialPropertiesTable(MPTsurface);
  opInHolderCaseSurface = dynamic_cast<G4OpticalSurface*>
    (logicSurfaceInHolderCase->GetSurface(fPhysIn, fPhysHolderReflector));
  if(opInHolderCaseSurface) opInHolderCaseSurface->DumpInfo();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACBaseDetectorConstruction::ConstructGuideBoundary()
{
  // photon energy definition
  const G4int nEntries = 20;
  const G4double pEmin = 1.9074*eV;
  const G4double pEmax = 6.1992*eV;
  G4double pE[nEntries] = {};
  G4double reflectivity[nEntries] = {};
  // Reflectivity is included in DAVIS LUT model
  for(G4int j = 0; j < nEntries;j++)
  {
    G4double dE = (pEmax - pEmin)/(nEntries- 1);
    pE[j] = pEmin + j*dE;
  }
  G4MaterialPropertiesTable *MPTsurface = new G4MaterialPropertiesTable();
  G4OpticalSurface *opInGuideCaseSurface = new G4OpticalSurface("InGuideCaseSurface");
  #if defined(AC_REFLECTOR_DAVISLUT_ESR)
  opInGuideCaseSurface->SetModel(DAVIS);
  opInGuideCaseSurface->SetType(dielectric_LUTDAVIS);
  opInGuideCaseSurface->SetFinish(PolishedESR_LUT);
  #elif defined(AC_REFLECTOR_LUT_TYVEK)
  opInGuideCaseSurface->SetModel(LUT);
  opInGuideCaseSurface->SetType(dielectric_LUT);
  opInGuideCaseSurface->SetFinish(groundtyvekair);
  #elif (defined(AC_REFLECTOR_UNIFIED_ESR) || defined(AC_REFLECTOR_UNIFIED_ALMYLAR))
  opInGuideCaseSurface->SetModel(unified);
  opInGuideCaseSurface->SetType(dielectric_metal);
  opInGuideCaseSurface->SetFinish(polished);
  for(G4int j = 0;j < nEntries;j++)
  {
    #if defined(AC_REFLECTOR_UNIFIED_ESR)
    reflectivity[j] = fSplineESRref->Eval(pE[j]);
    // reflectivity[j] = 1.;
    #elif defined(AC_REFLECTOR_UNIFIED_ALMYLAR)
    reflectivity[j] = fSplineAlMylarRef->Eval(pE[j]);
    #endif
  }
  MPTsurface->AddProperty("REFLECTIVITY", pE, reflectivity, nEntries);
  opInGuideCaseSurface->SetMaterialPropertiesTable(MPTsurface);
  #endif

  G4LogicalBorderSurface *logicSurfaceInGuideCase = new G4LogicalBorderSurface("InGuideCaseSurface",
    fPhysIn, fPhysGuideReflector, opInGuideCaseSurface);
  opInGuideCaseSurface = dynamic_cast<G4OpticalSurface*>
    (logicSurfaceInGuideCase->GetSurface(fPhysIn, fPhysGuideReflector));
  if(opInGuideCaseSurface) opInGuideCaseSurface->DumpInfo();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACBaseDetectorConstruction::ConstructAerogelBoundary()
{
  G4OpticalSurface *opAerogel1AirSurface = new G4OpticalSurface("Aerogel1AirSurface");
  opAerogel1AirSurface->SetType(dielectric_dielectric);
  opAerogel1AirSurface->SetFinish(polished);
  opAerogel1AirSurface->SetModel(unified);
  G4OpticalSurface *opAerogel2AirSurface = new G4OpticalSurface("Aerogel2AirSurface");
  opAerogel2AirSurface->SetType(dielectric_dielectric);
  opAerogel2AirSurface->SetFinish(polished);
  opAerogel2AirSurface->SetModel(unified);
  G4OpticalSurface *opAerogel3AirSurface = new G4OpticalSurface("Aerogel3AirSurface");
  opAerogel3AirSurface->SetType(dielectric_dielectric);
  opAerogel3AirSurface->SetFinish(polished);
  opAerogel3AirSurface->SetModel(unified);
  new G4LogicalSkinSurface("LogicalAerogel1AirSurface", fLogicAerogel1, opAerogel1AirSurface);
  new G4LogicalSkinSurface("LogicalAerogel2AirSurface", fLogicAerogel2, opAerogel2AirSurface);
  new G4LogicalSkinSurface("LogicalAerogel3AirSurface", fLogicAerogel3, opAerogel3AirSurface);  
}


void G4ACBaseDetectorConstruction::ConstructGlassBoundary()
{
  // photon energy definition
  const G4int nEntries = 20;
  const G4double pEmin = 1.9074*eV;
  const G4double pEmax = 6.1992*eV;
  G4double pE[nEntries] = {};
  G4double rIndexGlass[nEntries] = {};
  for(G4int j = 0; j < nEntries;j++)
  {
    G4double dE = (pEmax - pEmin)/(nEntries- 1);
    pE[j] = pEmin + j*dE;
  }
  for(G4int j = 0;j < nEntries;j++)
  {
    rIndexGlass[j] = 1.4074;
  }
  G4MaterialPropertiesTable *MPTglassSurface = new G4MaterialPropertiesTable();
  MPTglassSurface->AddProperty("RINDEX", pE, rIndexGlass, nEntries);

  G4OpticalSurface *opInGlassSurface = new G4OpticalSurface("InGlassSurface");
  opInGlassSurface->SetType(dielectric_dielectric);
  opInGlassSurface->SetFinish(polished); 
  opInGlassSurface->SetModel(glisur);
  
  G4LogicalBorderSurface *logicSurfaceInGlass = new G4LogicalBorderSurface("InGlassSurface",
    fPhysIn, fPhysGlass, opInGlassSurface);  
  opInGlassSurface->SetMaterialPropertiesTable(MPTglassSurface);
  opInGlassSurface = dynamic_cast<G4OpticalSurface*>
    (logicSurfaceInGlass->GetSurface(fPhysIn, fPhysGlass));
  if(opInGlassSurface) opInGlassSurface->DumpInfo();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACBaseDetectorConstruction::ConstructSDandField()
{
  // register Cathode as Sensitve detector
  G4ACFineMeshSD *fineMeshSD= new G4ACFineMeshSD("FMPMT");
  G4SDManager *sdMan  = G4SDManager::GetSDMpointer();
  sdMan->AddNewDetector(fineMeshSD);
  fLogicCathode->SetSensitiveDetector(fineMeshSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACBaseDetectorConstruction::ConstructSpline()
{
  
  G4double aeroPeV[] ={
    1.413758474*eV, 1.459908551*eV, 1.509173306*eV, 1.564166667*eV, 1.620855685*eV, 1.681831141*eV,
    1.718745231*eV, 1.784528954*eV, 1.865805169*eV, 1.929806524*eV, 2.020682793*eV, 2.11630251*eV,
    2.221461105*eV, 2.337572647*eV, 2.455296799*eV, 2.597923875*eV, 2.758204569*eV, 2.939562917*eV,
    3.146529171*eV, 3.374152254*eV, 3.686964286*eV, 3.886737991*eV, 4.062637327*eV, 4.156416709*eV,
    4.27135124*eV,  4.374982342*eV, 4.507258505*eV, 4.82820173*eV,  5.006344716*eV, 5.148449838*eV,
    5.352661597*eV, 5.692060283*eV
    };

  G4double aeroAbs[] ={
    18.22930033*cm,  18.01357865*cm,  18.01357865*cm,  17.58935847*cm,  17.58935847*cm,
    17.58935847*cm,  17.55294532*cm,  17.38801475*cm,  17.06200245*cm,  17.38080952*cm,
    17.58935847*cm,  17.58935847*cm,  17.79999991*cm,  17.79999991*cm,  18.22930033*cm,
    17.79999991*cm,  18.22930033*cm,  19.3481747*cm,   20.78212726*cm,  21.28335064*cm,
    19.3481747*cm,   15.42872364*cm,  11.87134701*cm,  9.579878477*cm,  7.197472265*cm,
    5.407543229*cm,  4.13066524*cm,   3.169786417*cm,  2.598842737*cm,  2.186552956*cm,
    1.662455389*cm,  1.294493875*cm
    };
/*
  G4double aeroPeV[] ={
    2.349253992*eV,  2.462128734*eV,  2.586392131*eV,  2.723871333*eV,  2.753134972*eV,
    2.909447609*eV,  3.047884083*eV,  3.208160392*eV,  3.333378144*eV,  3.506737416*eV,
    3.688464845*eV,  3.843375772*eV,  3.950030427*eV,  4.088699518*eV,  4.237459028*eV,
    4.397467501*eV,  4.570017325*eV,  4.756661359*eV,  4.902274785*eV,  5.158890327*eV,
    5.375360402*eV,  5.635454545*eV,  5.949507407*eV
  };

  G4double aeroAbs[] ={
    1021.409752*cm,  1021.409752*cm,  1021.409752*cm,  1021.409752*cm,  1021.409752*cm, 
    1021.409752*cm,  759.2241255*cm,  564.3259223*cm,  456.5727631*cm,  339.3673053*cm, 
    236.7173008*cm,  151.6980506*cm,  110.3926096*cm,  69.26119813*cm,  42.54318497*cm, 
    26.6919314*cm,  16.39532387*cm,  10.28655612*cm,  7.328751566*cm,  4.696561352*cm, 
    3.799706047*cm,  2.884828501*cm,  2.190231396*cm

  };*/
  // Reflectivity Spectra for Commonly Used Reflectors, IEEE
  G4double esrPeV[] = {
    1.549775*eV,    1.585849102*eV, 1.642278383*eV, 1.680865717*eV, 1.738036995*eV,
    1.799236955*eV, 1.850430661*eV, 1.972906784*eV, 2.069954588*eV, 2.193731034*eV,
    2.314383743*eV, 2.453277084*eV, 2.58632055*eV,  2.64855985*eV,  2.724202561*eV,
    2.804299344*eV, 2.954879118*eV, 3.062395487*eV, 3.178039634*eV, 3.220964247*eV,
    3.250230695*eV, 3.265064271*eV, 3.272536267*eV, 3.280033863*eV, 3.287574611*eV,
    3.295150111*eV, 3.310397492*eV, 3.325786636*eV, 3.349135987*eV, 3.372824688*eV,
    3.429417053*eV, 3.496459321*eV, 3.714315673*eV, 3.939300546*eV, 4.14471155*eV,
    4.413175955*eV, 4.597528099*eV, 4.7501772*eV,   4.913291591*eV, 6.1991*eV
  };

  G4double esrRef[] = {
    0.99, 0.99, 0.98, 0.98, 0.98,
    0.98, 0.97,   0.97,   0.97, 0.98,
    0.98, 0.99, 0.98, 0.99, 0.97,
    0.98, 0.99, 0.98, 0.96, 0.88,
    0.74, 0.62, 0.53, 0.45, 0.39,
    0.34, 0.27, 0.17, 0.13, 0.13,
    0.13, 0.12, 0.12, 0.14, 0.16,
    0.14, 0.16, 0.21, 0.26, 0.26
  };

  G4double alMylarPeV[] = {
    1.9074*eV, 1.572137663*eV, 1.616538756*eV, 1.661588103*eV, 1.703118905*eV,
    1.753197819*eV, 1.806311018*eV, 1.865172122*eV, 1.925403054*eV, 1.992425569*eV,
    2.061307268*eV, 2.135122244*eV, 2.217854083*eV, 2.303539649*eV, 2.400133169*eV,
    2.50518209*eV, 2.619847457*eV, 2.740252766*eV, 2.872258614*eV, 3.017626345*eV,
    3.18557244*eV, 3.357485568*eV, 3.557830527*eV, 3.783602356*eV, 4.039982797*eV,
    4.320527383*eV, 4.613026592*eV, 4.96516377*eV
  };

  G4double alMylarRef[] ={
    0.70000, 0.75414,  0.76123,  0.76832,  0.77541,
    0.77778, 0.78251,  0.78487,  0.78723,  0.7896,
    0.7896, 0.7896, 0.79196,  0.7896, 0.7896,
    0.78723, 0.78487,  0.78487,  0.78251,  0.77778,
    0.77305, 0.76596,  0.7565, 0.74704,  0.73995,
    0.73759, 0.70922,  0.68558
  };

  G4double tyvekPeV[] = {
    1.560135027*eV, 1.602034066*eV, 1.648142782*eV, 1.696981972*eV, 1.746673221*eV,
    1.783667076*eV, 1.824586599*eV, 1.879685001*eV, 1.943479235*eV, 2.006115416*eV,
    2.078942824*eV, 2.154015589*eV, 2.203750997*eV, 2.273752462*eV, 2.407584668*eV,
    2.513245208*eV, 2.604697277*eV, 2.71845171*eV, 2.749781015*eV, 2.876869965*eV,
    3.022650377*eV, 3.241660567*eV, 3.419765139*eV, 3.618567923*eV, 3.841911537*eV,
    4.106381781*eV, 4.343074897*eV, 4.608723124*eV, 4.942877871*eV
  };

  G4double tyvekRef[] = {
    0.97872, 0.97163, 0.97636, 0.974, 0.97636, 
    0.974, 0.97163, 0.97636, 0.974, 0.974, 
    0.97872, 0.97636, 0.97636, 0.974, 0.974, 
    0.974, 0.974, 0.974, 0.974, 0.97163,
    0.97163, 0.96454, 0.95508, 0.94326, 0.92199, 
    0.90307, 0.87234, 0.83452, 0.80615
  };

  assert(sizeof(aeroPeV) == sizeof(aeroAbs));
  assert(sizeof(esrPeV) == sizeof(esrRef));
  assert(sizeof(alMylarPeV) == sizeof(alMylarRef));
  assert(sizeof(tyvekPeV) == sizeof(tyvekRef));

  fSplineAeroAbs  = new TSpline3("Aerogel Absortion Length", aeroPeV, aeroAbs, sizeof(aeroPeV)/sizeof(G4double));
  fSplineESRref   = new TSpline3("ESR Reflectivity", esrPeV, esrRef, sizeof(esrPeV)/sizeof(G4double));
  fSplineAlMylarRef   = new TSpline3("Aluminized Mylar Reflectivity", alMylarPeV, alMylarRef, sizeof(alMylarPeV)/sizeof(G4double));
  fSplineTyvekRef = new TSpline3("Tyvek Reflectivity", tyvekPeV, tyvekRef, sizeof(tyvekPeV)/sizeof(G4double));
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// A threshold Cherenkokv detector for K+/pi+ separation using silica aerogel, R. Siudak
G4double  G4ACBaseDetectorConstruction::GetAerogelAbs(G4double peV)
{
  if(fSplineAeroAbs == nullptr)
  {
    G4cerr << "fSplineAeroAbs in G4ACBaseDetectorConstruction class is null pointer!!" << G4endl
    << "It returns 0." << G4endl;
    return 0;
  }

  return fSplineAeroAbs->Eval(peV) < 0? 0.1 * cm : fSplineAeroAbs->Eval(peV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double  G4ACBaseDetectorConstruction::GetAerogelScat(G4double peV)
{
  peV = peV/eV;
  return 313.39*cm/(peV*peV*peV*peV);
  // return 1.*cm;
}
