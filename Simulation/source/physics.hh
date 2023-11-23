#ifndef PHYSICS_HH
#define PHYSICS_HH

// include all used header files
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"

// physics list class
class PhysicsList : public G4VModularPhysicsList {
public:
  // constructor
  PhysicsList();
  // destructor
  ~PhysicsList();
};

#endif
