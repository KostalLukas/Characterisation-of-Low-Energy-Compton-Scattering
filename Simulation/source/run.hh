#ifndef RUN_HH
#define RUN_HH

// include all used header files
#include <string>
#include "G4UserRunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4GenericMessenger.hh"

// run action class
class RunAction : public G4UserRunAction {
public:
  // constructor
  RunAction();
  // destructor
  ~RunAction();

  // functions to start and end run action from G4UserRunAction class
  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

  // declare variables for generic messenger
  G4int simNum;
  G4double simAng;

  // declate string to store simulation number and angle
  int intSimAng;
  std::string strSimNum, strSimAng;

private:
  // initialise generic messenger
  G4GenericMessenger *fMessenger;

};

#endif
