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
// $Id: G4ACDetectorConstruction40.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file G4ACDetectorConstruction40.cc
/// \brief Definition of the G4ACDetectorConstruction40 class


#include "G4ACDetectorConstruction40.hh"
#include "G4ACFineMeshSD.hh"
#include "G4Box.hh"
#include "G4ExtrudedSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4RotationMatrix.hh"
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction40::G4ACDetectorConstruction40()
: G4ACDetectorConstruction20()
{
  fNxMppc = 3;  fXmppc = 10.1*mm/2; fXsens = 6.*mm/2; fXrIn = fCrossX; fXrOut = fXmppc*fNxMppc;
  fNyMppc = 2;  fYmppc = 8.9*mm/2;  fYsens = 6.*mm/2; fYrIn = fCrossY;  fYrOut = fYmppc*fNyMppc;
  fL = 200.*mm/2;fNptr = 150;
  fGlassZ = 0.4*mm/2;fMppcZ  = 2.0*mm/2;fSensZ = (fMppcZ - fGlassZ)/2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction40::G4ACDetectorConstruction40(G4String name)
: G4ACDetectorConstruction20(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACDetectorConstruction40::~G4ACDetectorConstruction40()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G4ACDetectorConstruction40::Construct()
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

void G4ACDetectorConstruction40::ConstructGuide()
{
  // define intersaction solid of two solids which have Winston cross-sections
  std::vector<G4TwoVector> guideOutXptr, guideReflectorXptr, guideInXptr;
  std::vector<G4TwoVector> guideOutYptr, guideReflectorYptr, guideInYptr;
  // paraboloid approximation with 100 points
  FillWinston(guideInXptr, fNptr, fXrIn, fXrOut, fL);
  FillWinstonDistance(guideReflectorXptr, guideInXptr, fReflectorThickness, fNptr, fXrIn, fXrOut, fL);
  FillWinstonDistance(guideOutXptr, guideInXptr, fBoxThickness, fNptr, fXrIn, fXrOut, fL);
  FillWinston(guideInYptr, fNptr, fYrIn, fYrOut, fL);
  FillWinstonDistance(guideReflectorYptr, guideInYptr, fReflectorThickness, fNptr, fYrIn, fYrOut, fL);
  FillWinstonDistance(guideOutYptr, guideInYptr, fBoxThickness, fNptr, fYrIn, fYrOut, fL);
    
  G4ExtrudedSolid *solidGuideInX = new G4ExtrudedSolid("GuideInX", guideInXptr, fYrIn, G4TwoVector(), 1., G4TwoVector(), 1.);
  G4ExtrudedSolid *solidGuideInY = new G4ExtrudedSolid("GuideInY", guideInYptr, fXrIn, G4TwoVector(), 1., G4TwoVector(), 1.);
  G4ExtrudedSolid *solidGuideOutX = new G4ExtrudedSolid("GuideOutX", guideOutXptr, fYrIn, G4TwoVector(), 1., G4TwoVector(), 1.);
  G4ExtrudedSolid *solidGuideOutY = new G4ExtrudedSolid("GuideOutY", guideOutYptr, fXrIn, G4TwoVector(), 1., G4TwoVector(), 1.);
  G4ExtrudedSolid *solidGuideReflectorX = new G4ExtrudedSolid("GuideReflectorX", guideReflectorXptr, fYrIn, G4TwoVector(), 1., G4TwoVector(), 1.);
  G4ExtrudedSolid *solidGuideReflectorY = new G4ExtrudedSolid("GuideReflectorX", guideReflectorYptr, fXrIn, G4TwoVector(), 1., G4TwoVector(), 1.);
  G4RotationMatrix *rotationGuideY = new G4RotationMatrix();rotationGuideY->rotateY(90.*deg);
  G4RotationMatrix *rotationGuide = new G4RotationMatrix();rotationGuide->rotateX(-90.*deg);
  fSolidGuideIn = new G4IntersectionSolid("GuideIn", solidGuideInX, solidGuideInY, rotationGuideY, G4ThreeVector());
  fSolidGuideOut = new G4IntersectionSolid("GuideOut", solidGuideOutX, solidGuideOutY, rotationGuideY, G4ThreeVector());
  G4VSolid *solidGuideReflectorOut = new G4IntersectionSolid("GuideReflectorOut", solidGuideReflectorX, solidGuideReflectorY, rotationGuideY, G4ThreeVector());
  fSolidGuideReflector = new G4SubtractionSolid("GuideReflector", solidGuideReflectorOut, fSolidGuideIn);
  fSolidGuideCase = new G4SubtractionSolid("GuideCase", fSolidGuideOut, solidGuideReflectorOut);
  /*
  G4RotationMatrix *rotationGuideX = new G4RotationMatrix();rotationGuideX->rotate(90.*deg, G4ThreeVector(1., 0., 0.));
  G4RotationMatrix *rotationGuideY = new G4RotationMatrix();rotationGuideY->rotate(90.*deg, G4ThreeVector(0., 0., 1.));rotationGuideY->rotate(-90.*deg, G4ThreeVector(1., 0., 0.));  
  auto boxxx = new G4Box("boxxx", fYrIn, 200.*mm, fXrIn);
  fSolidGuideIn = new G4UnionSolid("GuideIn", solidGuideInY, boxxx, nullptr, G4ThreeVector(0.*mm, -199.*mm, 0));
  fSolidGuideIn = solidGuideInY;
  fSolidGuideOut = solidGuideOutY;
  fSolidGuideReflector = new G4SubtractionSolid("GuideCase", solidGuideReflectorY, fSolidGuideIn); 
  fSolidGuideCase = new G4SubtractionSolid("GuideCase", fSolidGuideOut, fSolidGuideReflector); 
  */
  fLogicGuideCase = new G4LogicalVolume(fSolidGuideCase, fMatBox, "GuideCase");
  fLogicGuideReflector = new G4LogicalVolume(fSolidGuideReflector, fMatReflector, "GuideReflector");
  
  fGuideOutPosVec = new G4ThreeVector(0.*mm, fYrIn + fBoxThickness, 0.*mm);
  fGuideInPosVec = new G4ThreeVector(*fGuideOutPosVec);
  fGuideReflectorPosVec = new G4ThreeVector(*fGuideOutPosVec);
  fPhysGuideCase  = new G4PVPlacement(rotationGuide,    *fGuideOutPosVec,  fLogicGuideCase,  "GuideCase",  fLogicWorld, false, 0, fCheckOverlaps);
  fPhysGuideReflector  = new G4PVPlacement(rotationGuide,    *fGuideReflectorPosVec,  fLogicGuideReflector,  "GuideReflector",  fLogicWorld, false, 0, fCheckOverlaps);

  fRotationScorer = new G4RotationMatrix();
  fScorerPosVec = new G4ThreeVector(0., fYrIn + fBoxThickness, 2*fL + fMppcZ);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction40::ConstructScorer()
{
  /////////////////////////////////////////////////////////////////
  // MPPC sensitive area = Scorer volume (this has Glass, cathode volume as daughter volumes)
  /////////////////////////////////////////////////////////////////
  fSolidScorer = new G4Box("Scorer", fNxMppc*fXmppc, fNyMppc*fYmppc, fMppcZ);
  fSolidGlass  = new G4Box("Glass", fNxMppc*fXmppc, fNyMppc*fYmppc, fGlassZ);
  // cathode
  fSolidCathode  = new G4Box("Cathode", fXsens, fYsens, fSensZ);
  fLogicScorer  = new G4LogicalVolume(fSolidScorer, fMatScorer, "Scorer");
  fLogicGlass   = new G4LogicalVolume(fSolidGlass, fMatGlass, "Glass");
  fLogicCathode = new G4LogicalVolume(fSolidCathode, fMatAir, "Cathode");
  fPhysScorer     = new G4PVPlacement(fRotationScorer, *fScorerPosVec, fLogicScorer, "Scorer", fLogicWorld, false, 0, fCheckOverlaps);
  fPhysGlass      = new G4PVPlacement(0,    G4ThreeVector(0.*mm, 0.*mm, fGlassZ - fMppcZ) , fLogicGlass, "Glass", fLogicScorer, false, 0, fCheckOverlaps);
  for(Int_t x = 0;x < fNxMppc; x++)
    for(Int_t y = 0;y < fNyMppc; y++)
      new G4PVPlacement(0, G4ThreeVector((2*x - fNxMppc + 1)*fXmppc, (2*y - fNyMppc + 1)*fYmppc, 2*fGlassZ + fSensZ - fMppcZ), 
      fLogicCathode, "Cathode", fLogicScorer, false, y + x*fNyMppc, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction40::ConstructInnerAir()
{
  
  G4RotationMatrix *rotationGuide = new G4RotationMatrix();rotationGuide->rotateX(-90.*deg);
  //G4RotationMatrix *rotationGuideY = new G4RotationMatrix();rotationGuideY->rotate(90.*deg, G4ThreeVector(0., 0., 1.));rotationGuideY->rotate(-90.*deg, G4ThreeVector(1., 0., 0.));  
  fSolidIn        = new G4UnionSolid("In", fSolidHolderIn, fSolidGuideIn, rotationGuide, *fGuideInPosVec - *fHolderInPosVec);
  fLogicIn        = new G4LogicalVolume(fSolidIn, fMatAir,"In");
  fPhysIn         = new G4PVPlacement(0,    *fHolderInPosVec, fLogicIn,   "In", fLogicWorld, false, 0, fCheckOverlaps);
  // fLogicIn        = new G4LogicalVolume(fSolidGuideIn, fMatAir,"In");
  // fPhysIn         = new G4PVPlacement(rotationGuide, *fHolderInPosVec, fLogicIn, "In", fLogicWorld, false, 0, fCheckOverlaps);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction40::FillWinston(std::vector<G4TwoVector> &vec, G4int N, G4double r1, G4double r2, G4double l)
{
  vec.resize(2*(N + 1));
  vec.clear();
  G4double x, y;
  vec.emplace_back(-r1, 0.);
  for(G4int i = 1;i < N + 1;i++)
  {
    WinstonA1L(x, y, r1, l, i, N);
    vec.emplace_back(x, y);
  }
  for(G4int i = 0;i < N + 1;i++)
  {
    vec.emplace_back(-vec[N - i].x(), vec[N - i].y());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACDetectorConstruction40::FillWinstonDistance(std::vector<G4TwoVector> &dst, const std::vector<G4TwoVector> &src, G4double dist, G4int N, G4double r1, G4double r2, G4double l)
{
  G4double dy, x2, y2; // derivative of paraboloid
  const G4double k1 = (r2*r2 - r1*r1)/(2.*l);
  dst.resize(2*(N + 1));
  dst.clear();
  dst.emplace_back(src[0].x() - dist, src[0].y());
  for(G4int i = 1;i < N + 1;i++)
  {
    dst.emplace_back(src[i].x() - dist, src[i].y());
  }
  auto size = dst.size();
  for(G4int i = 0;i < size;i++)
  {
    dst.emplace_back(-dst[size - i - 1].x(), dst[size - i - 1].y());
  }
  /*
  std::cout << "----------------------------------------" << std::endl;
  std::cout << dst.size() << "\t" << dst.capacity() << std::endl;
  for(G4int i =0;i < 2*(N + 1);i++)
    std::cout << dst[i].x() << "\t" << dst[i].y() << std::endl;
  std::cout << "----------------------------------------" << std::endl;*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ACDetectorConstruction40::WinstonA1A2(G4double &x, G4double &y, G4double a1, G4double a2, G4int n, G4int N)
{
  // tan, sin, cos
  const G4double pi = 3.1415926535897932384;  
  G4double theta = asin(a2/a1);
  G4double alpha = (pi/2 - theta)*(N - n)/N;
  G4double l = (a1 + a2)*sqrt(a1*a1/(a2*a2) - 1);
  G4double f = a2*(1. + a2/a1);
  y = l - 2*f*sin(alpha)/(1. - sin(alpha - theta));
  x = a2 - 2*f*cos(alpha)/(1. - sin(alpha - theta));
  return l/2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ACDetectorConstruction40::WinstonA1L(G4double &x, G4double &y, G4double a1, G4double l, G4int n, G4int N)
{
  // tan, sin, cos
  const G4double pi = 3.1415926535897932384;  
  G4double coeff[] = {1., 2., 4*l*l/(a1*a1), -2., -1.};
  G4double a2 = a1*GetRootNewton(this->Poly4, this->Poly4D, coeff, 0.5);
  G4double f = a2*(1. + a2/a1);
  G4double theta = asin(a2/a1);
  G4double alpha = (pi/2 - theta)*(N - n)/N;
  y = 2*(l - f*sin(alpha)/(1. - sin(alpha - theta)));
  x = a2 - 2*f*cos(alpha)/(1. - sin(alpha - theta));
  return a2;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ACDetectorConstruction40::GetRootNewton(G4double (*f)(G4double*, G4double*), G4double (*df)(G4double*, G4double*), G4double *p, G4double x_init, G4double tol)
{
  G4double x_old = x_init, x, err;
  do{
    x = x_old - f(&x_old, p)/df(&x_old, p);
    err = x_old > x?(x_old - x)/x:(x - x_old)/x;
    x_old = x;
  }while(err > tol);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ACDetectorConstruction40::Poly4(G4double *x, G4double *p)
{
  G4double xin = x[0];
  return p[0]*xin*xin*xin*xin + p[1]*xin*xin*xin  + p[2]*xin*xin + p[3]*xin + p[4];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ACDetectorConstruction40::Poly4D(G4double *x, G4double *p)
{
  G4double xin = x[0];
  return 4.*p[0]*xin*xin*xin + 3.*p[1]*xin*xin  + 2.*p[2]*xin + p[3];  
}