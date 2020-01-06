// written by Hyunmin Yang in HANUL
// supervisor : Prof Ahn
// This class inherits G4VSensitveDetector and play a role as Fine mesh pmt

#ifndef G4ACFINEMESHSD_HH
#define G4ACFINEMESHSD_HH 1

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "Randomize.hh"
#include "TSpline.h"

class G4VSensitiveDetector;
class G4Step;

class G4ACFineMeshSD : public G4VSensitiveDetector{
    public:
    // constructor
    G4ACFineMeshSD(G4String name);
    // destructor
    virtual ~G4ACFineMeshSD();

    private:
    // G4double fTravelLength;
    G4int fEvtID;
    G4int fNphotonCheren;
    G4int fNphotonScint;
    G4int fNphotonCollected;
    G4int fNphotonDetected;
    TSpline3 *fQEspline;
    TSpline3 *GetQuantumEfficiency(); // interpolate wavelength dependent quantum efficiency of Hamamatsu R5543

    // TT + TTS : time taken between emitting photoelectron at cathode and pulse peak
    // rise time : time to rise from 10 % to 90 %  of pulse peak amplitude        
    const G4double fTimeTrans; // Electron Transition Time
    const G4double fTimeTransSpread; // Electron Transition Time Spread (sigma of gaussian)
    const G4double fTimeRise; // rise time of signal
    G4RandGauss *fRandGause;

    G4bool IsDetected(G4double ePhoton = 0); // determine PMT detects incoming photons or not
    G4double GetTransitTime(G4Step *aStep) const;
    public:
    G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);

    void Initialize(G4HCofThisEvent* HCE);
    void EndOfEvent(G4HCofThisEvent *HCE);
};


#endif