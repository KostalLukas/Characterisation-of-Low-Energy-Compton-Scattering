// include run header file
#include "run.hh"

// run action initialisation constructor
RunAction::RunAction() {

  // define messenger command directory
  fMessenger = new G4GenericMessenger(this, "/output/", "prepare root output parameters");

  // define messenger commands
  fMessenger->DeclareProperty("simNum", simNum, "set output simulation number");
  fMessenger->DeclareProperty("simAng", simAng, "set scattering angle");

  // specify generic values
  simNum = 0;
  simAng = 20;

}

// run action initialisation destructor
RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run*) {

  // convert simulation number to string
  strSimNum = std::to_string(simNum);

  // round scattering angle to integer and convert to string
  intSimAng = std::round(simAng);
  strSimAng = std::to_string(intSimAng);

  while (strSimAng.length() < 3) {
    strSimAng = "0" + strSimAng;
  }

  // initialise analysis manager singeltonian
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  // set thread merging for output
  man->SetNtupleMerging(true);

  // open root output file
  man->OpenFile("../output/Sim" + strSimNum + "_" + strSimAng + ".root");

  // define output columns
  man->CreateNtuple("events", "events");

  // define output rows
  man->CreateNtupleIColumn("Nevnt");
  man->CreateNtupleDColumn("Etot");
  man->CreateNtupleDColumn("Xmom");
  man->CreateNtupleDColumn("Ymom");
  man->CreateNtupleDColumn("Zmom");

  // finish writing rows to 0th column
  man->FinishNtuple(0);
}

void RunAction::EndOfRunAction(const G4Run*) {

  // initialise analysis manager singelton
  G4AnalysisManager *man = G4AnalysisManager::Instance();

  // write out root output file and close it
  man->Write();
  man->Reset();
  //man->CloseFile();
}
