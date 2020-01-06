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
// $Id: exampleG4AC.cc 86065 2014-11-07 08:51:15Z gcosmo $
//
/// \file exampleG4AC.cc
/// \brief Main program of the G4AC example




#include "G4ACDetectorConstruction.hh"
#include "G4ACActionInitialization.hh"
#include "G4ACPhysicsList.hh"
#include "G4ACAnalysis.hh"
#include "time.h"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "QBBC.hh"
#include "G4String.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool ACverbose = false; // print output?

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrintUsage()
{
    G4cerr << " Usage: " << G4endl;
    G4cerr << " simulateAC [-m macro(without .mac) ] [-r seed] [-t n threads]"
           << G4endl;
}


int main(int argc,char** argv)
{
  if ( argc > 9 )
  {
    PrintUsage();
    return 1;
  }
  
  G4String file;
  
  // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(time(NULL));
  // default number of threads
  G4int nThreads = 4;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv); // without option
  }
  else{
    for(G4int i = 1;i < argc;i = i+2) // with option
    {
        if(G4String(argv[i]) == "-m") 
        {
          file = argv[i+1];
        }
        else if(G4String(argv[i]) == "-r")
        {
          G4Random::setTheSeed(atoi(argv[i+1]));
        }
        else if(G4String(argv[i]) == "-t")
        {
          nThreads = atoi(argv[i+1]);
        }
        else  {
          PrintUsage();
          return 1;
        }
    }
  }

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new G4ACDetectorConstruction());

  // Physics list
  runManager->SetUserInitialization(new G4ACPhysicsList());
    
  // User action initialization
  runManager->SetUserInitialization(new G4ACActionInitialization());
  
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/ACverbose guidance.
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + file + ".mac");
  }
  else { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
    
  delete visManager;
  delete runManager; // the AnalysisManager is deleted here!
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
