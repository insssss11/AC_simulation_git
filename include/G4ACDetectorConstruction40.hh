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
// $Id: G4ACDetectorConstruction40.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file G4ACDetectorConstruction40.hh
/// \brief Definition of the G4ACDetectorConstruction40 class
// Design TYPE39  : (TYPEM, ESR, davis)

#ifndef G4ACDetectorConstruction40_h
#define G4ACDetectorConstruction40_h 1

#include <vector>
#include "G4PVPlacement.hh"
#include "G4ACDetectorConstruction20.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4TwoVector.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

/// Detector construction class to define materials and geometry.

class G4ACDetectorConstruction40 : public G4ACDetectorConstruction20
{
  public:
    G4ACDetectorConstruction40();
    G4ACDetectorConstruction40(G4String file);
    virtual ~G4ACDetectorConstruction40();
    
    virtual G4VPhysicalVolume *Construct();
  protected:
    virtual void ConstructGuide();
    virtual void ConstructScorer();
    virtual void ConstructInnerAir();
  private:
    G4int fNxMppc, fNyMppc;
    G4double fXmppc, fYmppc; // xy size of a single MPPC module
    G4double fXsens, fYsens; // xy size of sensitive area of a single MPPC module
    G4double fGlassZ, fSensZ, fMppcZ;
    G4double fXrIn, fXrOut; // parabolid parameters(crosssection is parallel to x axis)
    G4double fYrIn, fYrOut; // parabolid parameters(crosssection is parallel to x axis)
    G4double fL;
    G4int fNptr; // how many points to approximate with?
    void FillWinston(std::vector<G4TwoVector> &vec, G4int N, G4double r1, G4double r2, G4double l);
    void ReadConfigureFile(G4String file);
    void FillWinstonDistance(std::vector<G4TwoVector> &dst, const std::vector<G4TwoVector> &src, G4double dist, G4int N, G4double r1, G4double r2, G4double l);// to make a shell shape solid with Winston cross-section
    // Get n-th point on winston curve with two independent variables(input radius a1, output radius a2, length l)
    // it returns a value of variable fixed by two independent variables
    G4double WinstonA1A2(G4double &x, G4double &y, G4double a1, G4double a2, G4int n, G4int N);
    G4double WinstonA1L(G4double &x, G4double &y, G4double a1, G4double l, G4int n, G4int N);
    // root finding algorithm used in WinstonA1L
    G4double GetRootNewton(G4double (*f)(G4double*, G4double*), G4double (*df)(G4double*, G4double*), G4double *p, G4double x_init, G4double tol = 1e-8);
    static G4double Poly4(G4double*, G4double*);
    static G4double Poly4D(G4double*, G4double*);
};


namespace ac_detector_type40
{
  typedef G4ACDetectorConstruction40 G4ACDetectorConstruction;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

