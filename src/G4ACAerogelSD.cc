#include "G4ACAerogelSD.hh"
#include "G4ACAnalysis.hh"

#include "G4StepPoint.hh" 
#include "Randomize.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4ACPrimaryGeneratorAction.hh"
#include "G4ACBaseDetectorConstruction.hh"

G4ACAerogelSD::G4ACAerogelSD(G4String name, G4double offsetAngle)
    :G4VSensitiveDetector(name), fOffsetAngle(offsetAngle), fMultiplicity(0), fPreviousVolumeID(-1), fNofParticlesDetected(0), fCurrentParticleName("")
{
    G4cout << "Creating SD with name " + name << G4endl;
}

G4ACAerogelSD::~G4ACAerogelSD()
{
}

G4bool G4ACAerogelSD::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    // G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    // G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
    G4Track *aTrack = aStep->GetTrack();
    
    G4String particleName = aTrack->GetParticleDefinition()->GetParticleName();
    G4int volumeID = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);

    auto man = G4AnalysisManager::Instance();
    auto runMan = G4RunManager::GetRunManager();
    auto primaryAction = (G4ACPrimaryGeneratorAction *)runMan->GetUserPrimaryGeneratorAction();
    
    // this if-else statement is for excluding multiple counting of particle hits in the same detector volume(I want only data about first incident step into detector volmue)
    if(particleName == fCurrentParticleName)
    {
        if(volumeID == fPreviousVolumeID)
            return true;
        else
        {
            fMultiplicity++;
        }
    }
    else
    {
        fNofParticlesDetected++;
        fMultiplicity = 0;
    }

    // fill tuple
    if(particleName == "pi+")
        man->FillNtupleIColumn(2, 6, 0);
    else if(particleName == "pi-")
        man->FillNtupleIColumn(2, 6, 1);
    else if(particleName == "kaon+")
        man->FillNtupleIColumn(2, 6, 2);
    else if(particleName == "kaon-")
        man->FillNtupleIColumn(2, 6, 3);
    else{
        fNofParticlesDetected--;
        aTrack->SetTrackStatus(fStopAndKill);
        return true;
    }

    fPreviousVolumeID = volumeID;
    fCurrentParticleName = particleName;

    CalculateIncidentInfo(aStep, fPosition, fMomentum);
    man->FillNtupleFColumn(2, 0, fPosition[0]);
    man->FillNtupleFColumn(2, 1, fPosition[1]);
    man->FillNtupleFColumn(2, 2, fPosition[2]);
    man->FillNtupleFColumn(2, 3, fMomentum[0]);
    man->FillNtupleFColumn(2, 4, fMomentum[1]);
    man->FillNtupleFColumn(2, 5, fMomentum[2]);
    man->FillNtupleIColumn(2, 7, primaryAction->GetEventGenerator()->GetDecayChannel());
    man->FillNtupleIColumn(2, 8, fMultiplicity);
    man->FillNtupleIColumn(2, 9, volumeID);
    man->AddNtupleRow(2);
    return true;
}

void G4ACAerogelSD::CalculateIncidentInfo(const G4Step *aStep, G4double *position, G4double *momentum)
{
    
    // G4VPhysicalVolume *vol = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    G4int copyNum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
    
    G4ThreeVector incidentVector = aStep->GetPreStepPoint()->GetPosition();
    incidentVector.rotateZ(-(fOffsetAngle + copyNum*12.*deg + 180.*deg));
    incidentVector[1] = incidentVector[1] - TPC_EDGE_R + AC_DETECTOR_R;
    incidentVector[2] = incidentVector[2] - TPC_DOWNSTREAM_Z;
    G4ThreeVector incidentMomentum = aStep->GetPreStepPoint()->GetMomentum();
    incidentMomentum.rotateZ(-(fOffsetAngle + copyNum*12.*deg + 180.*deg));

    position[0] = incidentVector[0]/mm;
    position[1] = incidentVector[1]/mm;
    position[2] = incidentVector[2]/mm;
    momentum[0] = incidentMomentum[0]/MeV;
    momentum[1] = incidentMomentum[1]/MeV;
    momentum[2] = incidentMomentum[2]/MeV;
   
}

void G4ACAerogelSD::Initialize(G4HCofThisEvent*)
{
    fPosition[0] = 0., fPosition[1] = 0., fPosition[2] = 0.;
    fMomentum[0] = 0., fMomentum[1] = 0., fMomentum[2] = 0.;
    fPreviousVolumeID = -1;
    fCurrentParticleName = "";
    fMultiplicity = 0;
    fNofParticlesDetected = 0;
}

void G4ACAerogelSD::EndOfEvent(G4HCofThisEvent *)
{
    auto runMan = G4RunManager::GetRunManager();
    auto primaryAction = (G4ACPrimaryGeneratorAction *)runMan->GetUserPrimaryGeneratorAction();
    auto man = G4AnalysisManager::Instance();
    man->FillNtupleIColumn(3, 0, primaryAction->GetEventGenerator()->GetDecayChannel());
    man->FillNtupleIColumn(3, 1, runMan->GetCurrentEvent()->GetEventID());
    man->FillNtupleIColumn(3, 2, fNofParticlesDetected);
    man->AddNtupleRow(3);
}
