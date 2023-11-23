#ifndef ACTION_HH
#define ACTION_HH

// include all used header files
#include "G4VUserActionInitialization.hh"

// include particle generator header file
#include "generator.hh"

// include run header file
#include "run.hh"

// action initialization class
class ActionInitialization : public G4VUserActionInitialization {
public:
  // constructor
  ActionInitialization();
  // destructor
  ~ActionInitialization();

  // method to build and run simulation for master thread
  virtual void BuildForMaster() const;

  // method to build and run simulation for worker thread
  virtual void Build() const;
};

#endif
