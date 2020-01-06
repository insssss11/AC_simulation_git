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
// $Id: G4ACDetectorConstruction29.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACDetectorConstruction29.cc
/// \brief Definition of the G4ACDetectorConstruction29 class

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

#include "G4ACDetectorConstruction29.hh"
#include "G4ACFineMeshSD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction29::G4ACDetectorConstruction29()
: G4ACDetectorConstruction21()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4ACDetectorConstruction29::~G4ACDetectorConstruction29()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACDetectorConstruction29::Construct()
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

void G4ACDetectorConstruction29::ConstructHolderBoundary()
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
    fPhysIn, fPhysHolderReflector, opInHolderCaseSurface);  
  // opInHolderCaseSurface->SetMaterialPropertiesTable(MPTsurface);
  opInHolderCaseSurface = dynamic_cast<G4OpticalSurface*>
    (logicSurfaceInHolderCase->GetSurface(fPhysIn, fPhysHolderReflector));
  if(opInHolderCaseSurface) opInHolderCaseSurface->DumpInfo();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction29::ConstructGuideBoundary()
{ 
  G4OpticalSurface *opInGuideCaseSurface = new G4OpticalSurface("InGuideCaseSurface");
  opInGuideCaseSurface->SetType(dielectric_LUT);
  opInGuideCaseSurface->SetFinish(groundtyvekair);
  opInGuideCaseSurface->SetModel(LUT);
  
  G4LogicalBorderSurface *logicSurfaceInGuideCase = new G4LogicalBorderSurface("InGuideCaseSurface",
    fPhysIn, fPhysGuideReflector, opInGuideCaseSurface);  
  
  // opInGuideCaseSurface->SetMaterialPropertiesTable(MPTsurface);
  opInGuideCaseSurface = dynamic_cast<G4OpticalSurface*>
    (logicSurfaceInGuideCase->GetSurface(fPhysIn, fPhysGuideReflector));
  if(opInGuideCaseSurface) opInGuideCaseSurface->DumpInfo();
}

void G4ACDetectorConstruction29::ConstructAerogel()
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

void G4ACDetectorConstruction29::ConstructMaterial()
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
  MPTreflector->AddConstProperty("SCINTILLATIONYIELD", 0.801*5000./MeV);
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