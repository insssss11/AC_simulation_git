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
// $Id: G4ACBaseDetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// This header file select the design of AC detector

#ifndef G4ACDetectorConstruction_h
#define G4ACDetectorConstruction_h 1
#ifdef AC_DETECTOR_DESIGN_DEFAULT
    #include "G4ACBaseDetectorConstruction.hh"
    using namespace ac_detector_base;
#elif defined AC_DETECTOR_DESIGN_TYPE1
    #include "G4ACDetectorConstruction1.hh"
    using namespace ac_detector_type1;
#elif defined AC_DETECTOR_DESIGN_TYPE2
    #include "G4ACDetectorConstruction2.hh"
    using namespace ac_detector_type2;
#elif defined AC_DETECTOR_DESIGN_TYPE3
    #include "G4ACDetectorConstruction3.hh"
    using namespace ac_detector_type3;
#elif defined AC_DETECTOR_DESIGN_TYPE4
    #include "G4ACDetectorConstruction4.hh"
    using namespace ac_detector_type4;
#elif defined AC_DETECTOR_DESIGN_TYPE5
    #include "G4ACDetectorConstruction5.hh"
    using namespace ac_detector_type5;
#elif defined AC_DETECTOR_DESIGN_TYPE6
    #include "G4ACDetectorConstruction6.hh"
    using namespace ac_detector_type6;
#elif defined AC_DETECTOR_DESIGN_TYPE7
    #include "G4ACDetectorConstruction7.hh"
    using namespace ac_detector_type7;    
#elif defined AC_DETECTOR_DESIGN_TYPE8
    #include "G4ACDetectorConstruction8.hh"
    using namespace ac_detector_type8;    
#elif defined AC_DETECTOR_DESIGN_TYPE20
    #include "G4ACDetectorConstruction20.hh"
    using namespace ac_detector_type20;    
#elif defined AC_DETECTOR_DESIGN_TYPE21
    #include "G4ACDetectorConstruction21.hh"
    using namespace ac_detector_type21;    
#elif defined AC_DETECTOR_DESIGN_TYPE22
    #include "G4ACDetectorConstruction22.hh"
    using namespace ac_detector_type22;
#elif defined AC_DETECTOR_DESIGN_TYPE23
    #include "G4ACDetectorConstruction23.hh"
    using namespace ac_detector_type23;
#elif defined AC_DETECTOR_DESIGN_TYPE24
    #include "G4ACDetectorConstruction24.hh"
    using namespace ac_detector_type24;
#elif defined AC_DETECTOR_DESIGN_TYPE25
    #include "G4ACDetectorConstruction25.hh"
    using namespace ac_detector_type25;
#elif defined AC_DETECTOR_DESIGN_TYPE26
    #include "G4ACDetectorConstruction26.hh"
    using namespace ac_detector_type26;
#elif defined AC_DETECTOR_DESIGN_TYPE27
    #include "G4ACDetectorConstruction27.hh"
    using namespace ac_detector_type27;  
#elif defined AC_DETECTOR_DESIGN_TYPE28
    #include "G4ACDetectorConstruction28.hh"
    using namespace ac_detector_type28;
#elif defined AC_DETECTOR_DESIGN_TYPE29
    #include "G4ACDetectorConstruction29.hh"
    using namespace ac_detector_type29;      
#elif defined AC_DETECTOR_DESIGN_TYPE30
    #include "G4ACDetectorConstruction30.hh"
    using namespace ac_detector_type30;
#elif defined AC_DETECTOR_DESIGN_TYPE31
    #include "G4ACDetectorConstruction31.hh"
    using namespace ac_detector_type31;
#elif defined AC_DETECTOR_DESIGN_TYPE32
    #include "G4ACDetectorConstruction32.hh"
    using namespace ac_detector_type32;
#elif defined AC_DETECTOR_DESIGN_TYPE33
    #include "G4ACDetectorConstruction33.hh"
    using namespace ac_detector_type33;
#elif defined AC_DETECTOR_DESIGN_TYPE34
    #include "G4ACDetectorConstruction34.hh"
    using namespace ac_detector_type34;
#elif defined AC_DETECTOR_DESIGN_TYPE35
    #include "G4ACDetectorConstruction35.hh"
    using namespace ac_detector_type35;
#elif defined AC_DETECTOR_DESIGN_TYPE36
    #include "G4ACDetectorConstruction36.hh"
    using namespace ac_detector_type36;   
#elif defined AC_DETECTOR_DESIGN_TYPE37
    #include "G4ACDetectorConstruction37.hh"
    using namespace ac_detector_type37;
#elif defined AC_DETECTOR_DESIGN_TYPE38
    #include "G4ACDetectorConstruction38.hh"
    using namespace ac_detector_type38;
#elif defined AC_DETECTOR_DESIGN_TYPE39
    #include "G4ACDetectorConstruction39.hh"
    using namespace ac_detector_type39;
#elif defined AC_DETECTOR_DESIGN_TYPE40
    #include "G4ACDetectorConstruction40.hh"
    using namespace ac_detector_type40;
#elif defined AC_DETECTOR_DESIGN_LEPS2
    #include "G4ACDetectorConstructionLEPS2.hh"
    using namespace ac_detector_type_leps2;
#elif defined AC_TESTBENCH1
    #include "G4ACTestbenchConstruction.hh"
    using namespace ac_testbench;
#elif defined AC_TESTBENCH2
    #include "G4ACTestbenchConstruction.hh"
    using namespace ac_testbench;
#elif defined AC_TESTBENCH3
    #include "G4ACTestbenchConstruction.hh"
    using namespace ac_testbench;
#elif defined AC_TESTBENCH4
    #include "G4ACTestbenchConstruction.hh"
    using namespace ac_testbench;                                      
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif