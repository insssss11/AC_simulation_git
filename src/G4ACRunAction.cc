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
// $Id: G4ACRunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file G4ACRunAction.cc
/// \brief Implementation of the G4ACRunAction class

#include "G4ACRunAction.hh"
#include "G4ACPrimaryGeneratorAction.hh"
#include "G4ACAnalysis.hh"

// #include "G4ACRun.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACRunAction::G4ACRunAction()
: G4UserRunAction(),
  fTimer(0), fMessenger(0), fFileName(0)
{
  fTimer = new G4Timer();
  fFileName =  new G4String("simulateAC");
  fMessenger = new G4ACRunActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ACRunAction::~G4ACRunAction()
{
  delete fTimer;
  delete fMessenger;
  delete fFileName;
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACRunAction::BeginOfRunAction(const G4Run* aRun)
{
  auto man = G4AnalysisManager::Instance();
  man->SetNtupleMerging(true);
  man->SetNtupleRowWise(false);
  man->OpenFile(fFileName->data());
  man->SetVerboseLevel(0);
  #ifndef AC_DETECTOR_DESIGN_LEPS2  
    // tree t0 and t1 will be filled if gunmode is single and t2 will be filled otherwise 
    man->CreateNtuple("t0", "Total Npe(per an event)");
    man->CreateNtupleIColumn(0, "evtID");
    man->CreateNtupleIColumn(0, "Ncheren");
    man->CreateNtupleIColumn(0, "Nscint");
    man->CreateNtupleIColumn(0, "Ncol");
    man->CreateNtupleIColumn(0, "Npe");
    man->FinishNtuple();

    man->CreateNtuple("t1", "Photon statistics(per an detected photon)");
    man->CreateNtupleIColumn(1, "photonID");
    man->CreateNtupleIColumn(1, "evtID");
    man->CreateNtupleFColumn(1, "wavelen");
    man->CreateNtupleFColumn(1, "t_signal");
    man->CreateNtupleIColumn(1, "creator"); // 0 : cherenkov, 1 : scintillation
    man->FinishNtuple();
  #endif
  #ifdef AC_DETECTOR_DESIGN_LEPS2
    man->CreateNtuple("t2", "Hit data of Holder Case(Particle ID order : pi+, pi-, K+, K-. Position is in mm and momentum is in MeV/c)");
    man->CreateNtupleFColumn(2, "x");
    man->CreateNtupleFColumn(2, "y");
    man->CreateNtupleFColumn(2, "z");
    man->CreateNtupleFColumn(2, "px");
    man->CreateNtupleFColumn(2, "py");
    man->CreateNtupleFColumn(2, "pz");
    man->CreateNtupleIColumn(2, "PID");
    man->CreateNtupleIColumn(2, "ch"); // decay channel
    man->CreateNtupleIColumn(2, "multi"); // multiplicity (if one particle hits multiple detectors)
    man->CreateNtupleIColumn(2, "volumeNum");
    man->FinishNtuple();

    man->CreateNtuple("t3", "Hit data of Holder Case(collected once in an event)");
    man->CreateNtupleIColumn(3, "ch"); // decay channel
    man->CreateNtupleIColumn(3, "eventID");
    man->CreateNtupleIColumn(3, "NofParticlesDet");
    man->FinishNtuple();
  #endif
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  fTimer->Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACRunAction::EndOfRunAction(const G4Run* aRun)
{
  fTimer->Stop();
  auto man = G4AnalysisManager::Instance();
    G4cout << "number of event = " << aRun->GetNumberOfEvent()
          << " " << *fTimer << G4endl;  
  man->Write();
  G4cout << "File " << man->GetFileName()+".root has been saved." << G4endl;
  man->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ACRunAction::SetFileName(G4String fileName)
{
  (*fFileName) = fileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const char* G4ACRunAction::GetFileName() const
{
  return fFileName->data();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
