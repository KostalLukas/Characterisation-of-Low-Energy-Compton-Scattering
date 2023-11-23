#ifndef GENERATOR_HH
#define GENERATOR_HH

// include all used header files
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4GenericMessenger.hh"
#include "Randomize.hh"

// primary generator action class
class PrimaryGenerator : public G4VUserPrimaryGeneratorAction {
public:
  // constructor
  PrimaryGenerator();
  // destructor
  ~PrimaryGenerator();

  // method to generate primaries
  virtual void GeneratePrimaries(G4Event*);

private:
  // initialise particle gun
  G4ParticleGun *fParticleGun;

  // initialise generic messenger
  G4GenericMessenger *fMessenger;

  // declare variables for generic messenger
  G4bool momRand, momRest;
  G4double partE, scatAng;
  G4ThreeVector pos, mom;

};

#endif
