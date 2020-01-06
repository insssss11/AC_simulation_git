#include "G4ACFineMeshSD.hh"
#include "G4ACAnalysis.hh"

#include "G4StepPoint.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"

G4ACFineMeshSD::G4ACFineMeshSD(G4String name)
    : G4VSensitiveDetector(name),
      fEvtID(0),
      fNphotonCheren(0), fNphotonScint(0), fNphotonCollected(0), fNphotonDetected(0), fQEspline(0),
      fTimeTrans(12.3 * ns), fTimeTransSpread(0.5 * ns), fTimeRise(2.9 * ns)
{
  G4cout << "Creating SD with name " + name << G4endl;
  fQEspline = GetQuantumEfficiency();
  fRandGause = new G4RandGauss(G4Random::getTheEngine(), 0., 2.35 * fTimeTransSpread);
}

G4ACFineMeshSD::~G4ACFineMeshSD()
{
  delete fQEspline;
  delete fRandGause;
}

G4bool G4ACFineMeshSD::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
  G4double eP = aStep->GetTrack()->GetTotalEnergy();
  G4String procName;
  if(aStep->GetTrack()->GetCreatorProcess() != nullptr)
    procName = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
  else{
    procName = "It a primary particle :(";
  }
  auto man = G4AnalysisManager::Instance();

  fNphotonCollected++;
  if (IsDetected(eP))
  {
    man->FillNtupleIColumn(1, 0, fNphotonDetected);
    man->FillNtupleIColumn(1, 1, fEvtID);
    // man->FillNtupleFColumn(1, 1, 1.23984193e-3/eP);
    man->FillNtupleFColumn(1, 3, GetTransitTime(aStep));
    if (procName == "Scintillation")
    {
      fNphotonScint++;
      man->FillNtupleIColumn(1, 4, 1);
    }
    else if (procName == "Cerenkov")
    {
      fNphotonCheren++;
      man->FillNtupleIColumn(1, 4, 0);
    }
    else
      man->FillNtupleIColumn(1, 4, 2);
    man->AddNtupleRow(1);
    fNphotonDetected++;
  }
  aStep->GetTrack()->SetTrackStatus(fStopAndKill);
  return true;
}

void G4ACFineMeshSD::Initialize(G4HCofThisEvent *)
{
  // prepare for new event
  fNphotonScint = 0;
  fNphotonCheren = 0;
  fNphotonCollected = 0;
  fNphotonDetected = 0;
  fEvtID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
}

void G4ACFineMeshSD::EndOfEvent(G4HCofThisEvent *)
{
  G4cout << "Number of photons reached to Cathode : "
         << fNphotonCollected << G4endl;
  G4cout << "Number of photons reached to Cathode and made signal : "
         << fNphotonDetected << G4endl;
  auto man = G4AnalysisManager::Instance();
  man->FillNtupleIColumn(0, 0, fEvtID);
  man->FillNtupleIColumn(0, 1, fNphotonCheren);
  man->FillNtupleIColumn(0, 2, fNphotonScint);
  man->FillNtupleIColumn(0, 3, fNphotonCollected);
  man->FillNtupleIColumn(0, 4, fNphotonDetected);
}

G4bool G4ACFineMeshSD::IsDetected(G4double ePhoton)
{
  return G4UniformRand() < fQEspline->Eval(ePhoton) / 100;
}

TSpline3 *G4ACFineMeshSD::GetQuantumEfficiency()
{
  const int nBins = 63;
  double pE[nBins] =
  {
    1.742949282 * eV, 1.750743809 * eV, 1.75860587 * eV, 1.770532256 * eV, 1.774543747 * eV,
    1.778573457 * eV, 1.79487959 * eV, 1.807304786 * eV, 1.815684264 * eV, 1.854376744 * eV,
    1.867645198 * eV, 1.876594933 * eV, 1.890181527 * eV, 1.899351992 * eV, 1.917956556 * eV,
    1.922664785 * eV, 1.971053569 * eV, 1.986050997 * eV, 2.006401895 * eV, 2.02194089 * eV,
    2.043041229 * eV, 2.064583167 * eV, 2.086587736 * eV, 2.114757454 * eV, 2.137850528 * eV,
    2.149585371 * eV, 2.17344583 * eV, 2.204030697 * eV, 2.235488634 * eV, 2.274439159 * eV,
    2.328539179 * eV, 2.363680038 * eV, 2.414695709 * eV, 2.437237939 * eV, 2.467957203 * eV,
    2.531784117 * eV, 2.607653188 * eV, 2.688209811 * eV, 2.706796464 * eV, 2.783761392 * eV,
    2.875758173 * eV, 2.985386573 * eV, 3.10369688 * eV, 3.166435538 * eV, 3.23177137 * eV,
    3.356425271 * eV, 3.506711694 * eV, 3.671087345 * eV, 3.832782557 * eV, 3.889889375 * eV,
    4.009377981 * eV, 4.203032676 * eV, 4.416344006 * eV, 4.624991607 * eV, 4.824430332 * eV,
    4.85433168 * eV, 4.977737808 * eV, 5.107603072 * eV, 5.244403837 * eV, 5.502329916 * eV,
    5.661447417 * eV, 6.1992 * eV
  };

  double QE[nBins] =
  {
    0.02030673, 0.022284865, 0.026834901, 0.031720498, 0.03749212,
    0.038197066, 0.057481752, 0.06794852, 0.081820099, 0.151072355,
    0.249454303, 0.294858736, 0.355077375, 0.427573716, 0.505448602,
    0.47808487, 1.166310539, 1.404504063, 1.660309298, 1.962605448,
    2.320058843, 2.794016716, 3.302896768, 3.904649, 4.530895022,
    4.880721197, 5.356634326, 6.332555431, 7.348538943, 8.217051332,
    10.27211801, 12.1436712, 14.09321249, 14.62884371, 15.75977453,
    17.62422384, 19.3472837, 20.46491758, 20.46680256, 21.24956285,
    22.06327611, 22.90920403, 22.92028429, 23.35554468, 23.79907075,
    23.81003315, 23.8215491, 23.83307063, 22.97417897, 22.55485295,
    21.7390219, 19.82028979, 18.07090904, 16.17186101, 14.47171496,
    14.47271467, 12.24757197, 10.36453924, 8.77119972, 5.619012523,
    4.755147578, 2.828870065
  };
  return new TSpline3("Quantum Efficiency", pE, QE, nBins);
}

G4double G4ACFineMeshSD::GetTransitTime(G4Step *aStep) const
{
  return aStep->GetTrack()->GetGlobalTime() + (fTimeTrans + fRandGause->fire()) / ns;
}