// include physics header file
#include "physics.hh"

// physics list constructor
PhysicsList::PhysicsList() {

  // load physics list to implement EM interactions
  RegisterPhysics (new G4EmStandardPhysics());

  // load physics list to implement optical photon interactions
  //RegisterPhysics (new G4OpticalPhysics());

}

// physics list destructor
PhysicsList::~PhysicsList() {}
