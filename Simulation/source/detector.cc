// include detector header file
#include "detector.hh"

// action initialisation constructor
SensitiveDetector::SensitiveDetector(G4String name) : G4VSensitiveDetector(name) {}

// action initialisation destructor
SensitiveDetector::~SensitiveDetector() {}

// method to process particle hits in detector
G4bool SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist) {

  // get track of particle in sensitive detector
  G4Track *track = aStep->GetTrack();

  // get step points when particle enters and leaves sensitive detector
  G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint *postStepPoint = aStep->GetPostStepPoint();

  // get kinetic energy when particle enters sensitive detector
  G4double kenParticle = preStepPoint->GetKineticEnergy();

  // get kinetic energy when particle enters sensitive detector
  G4ThreeVector momParticle = preStepPoint->GetMomentum();

  // get particle position vector when enters sensitive detector
  G4ThreeVector posParticle = preStepPoint->GetPosition();

  // call instance of analysis manager
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  // get current event number from run manager
  G4int Nevnt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  // write event number and particle energy in keV
  man->FillNtupleIColumn(0, Nevnt);
  man->FillNtupleDColumn(1, kenParticle / keV);
  man->FillNtupleDColumn(2, momParticle[0] / keV);
  man->FillNtupleDColumn(3, momParticle[1] / keV);
  man->FillNtupleDColumn(4, momParticle[2] / keV);

  // add new row to root file
  man->AddNtupleRow(0);

return true;
}
