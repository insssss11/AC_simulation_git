// written by Hyunmin Yang in HANUL
// supervisor : Prof Ahn
// This class inherits G4VSensitveDetector to extract data about incident position of pi+, pi-, K+ 

#ifndef G4ACAerogelSD_HH
#define G4ACAerogelSD_HH 1

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"

class G4VSensitiveDetector;
class G4Step;

class G4ACAerogelSD : public G4VSensitiveDetector{
    public:
    // constructor
    G4ACAerogelSD(G4String name, G4double foffsetAngle);
    // destructor
    virtual ~G4ACAerogelSD();

    private:
    G4double fPosition[3];
    G4double fMomentum[3];
    G4double fOffsetAngle;
    
    // to take information of multiplicity
    G4int fMultiplicity;
    G4int fPreviousVolumeID;
    G4int fNofParticlesDetected;
    G4String fCurrentParticleName;
    void CalculateIncidentInfo(const G4Step *aStep, G4double *position, G4double *momentum); // calculate incident position and momentum from the origin of Holder volume

    public:
    G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);

    void Initialize(G4HCofThisEvent* HCE);
    void EndOfEvent(G4HCofThisEvent *HCE);
};


#endif
