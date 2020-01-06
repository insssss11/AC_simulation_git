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
// $Id: G4ACDetectorConstruction39.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACDetectorConstruction39.cc
/// \brief Definition of the G4ACDetectorConstruction39 class

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

#include "G4ACDetectorConstruction39.hh"
#include "G4ACFineMeshSD.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction39::G4ACDetectorConstruction39()
: G4ACDetectorConstruction22()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction39::~G4ACDetectorConstruction39()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACDetectorConstruction39::Construct()
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

void G4ACDetectorConstruction39::ConstructHolderBoundary()
{
  // photon energy definition
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
    reflectivity[j]= fSplineTyvekRef->Eval(pE[j]);
  }
  G4MaterialPropertiesTable *MPTsurface = new G4MaterialPropertiesTable();
  MPTsurface->AddProperty("REFLECTIVITY", pE, reflectivity, nEntries);
  G4OpticalSurface *opInHolderCaseSurface = new G4OpticalSurface("InHolderCaseSurface");
  opInHolderCaseSurface->SetType(dielectric_metal);
  opInHolderCaseSurface->SetFinish(ground); // mechanically grounded surface, with esr film & meltmount
  opInHolderCaseSurface->SetModel(unified);
  
  G4LogicalBorderSurface *logicSurfaceInHolderCase = new G4LogicalBorderSurface("InHolderCaseSurface",
    fPhysIn, fPhysHolderReflector, opInHolderCaseSurface);  
  opInHolderCaseSurface->SetMaterialPropertiesTable(MPTsurface);
  opInHolderCaseSurface = dynamic_cast<G4OpticalSurface*>
    (logicSurfaceInHolderCase->GetSurface(fPhysIn, fPhysHolderReflector));
  if(opInHolderCaseSurface) opInHolderCaseSurface->DumpInfo();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction39::ConstructGuideBoundary()
{
  // photon energy definition
  const G4int nEntries = 20;
  const G4double pEmin = 1.9074*eV;
  const G4double pEmax = 6.1992*eV;
  G4double pE[nEntries] = {};
  for(G4int j = 0; j < nEntries;j++)
  {
    G4double dE = (pEmax - pEmin)/(nEntries- 1);
    pE[j] = pEmin + j*dE;
  }

  G4double reflectivity[nEntries] = {};
  for(G4int j = 0;j < nEntries;j++)
  {
    reflectivity[j]= fSplineESRref->Eval(pE[j]);
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
    (logicSurfaceInGuideCase->GetSurface(fPhysIn, fPhysGuideCase));
  if(opInGuideCaseSurface) opInGuideCaseSurface->DumpInfo();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
