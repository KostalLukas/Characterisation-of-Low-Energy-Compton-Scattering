// include generator header file
#include "generator.hh"

// primary generator constructor
PrimaryGenerator::PrimaryGenerator() {

  // define messenger command directory
  fMessenger = new G4GenericMessenger(this, "/setup/", "prepare experimental setup");

  // define messenger commands
  fMessenger->DeclareProperty("scatAng", scatAng, "set scattering angle");
  fMessenger->DeclareProperty("momRand", momRand, "set random momentum direction vector");
  fMessenger->DeclareProperty("momRest", momRest, "restrict momentum direction vector to front hemisphere");
  fMessenger->DeclareProperty("partE", partE, "set particle energy in keV");

  fMessenger->DeclareProperty("momX", mom[0], "set momentum direction vector x component");
  fMessenger->DeclareProperty("momY", mom[1], "set momentum direction vector y component");
  fMessenger->DeclareProperty("momZ", mom[2], "set momentum direction vector z component");

  // specify generic values
  scatAng = 20;
  momRand = false;
  momRest = false;
  partE = 1000;

  // initial particle position and momentum direction vector
  pos = G4ThreeVector(0, 0, -180*mm);
  mom = G4ThreeVector(0, 0, 1);

  // define particle gun with 1 particle per 1 event
  fParticleGun = new G4ParticleGun();
}

// primary generator destructor
PrimaryGenerator::~PrimaryGenerator() {

  // delete particle gun
  delete fParticleGun;
}

// define method to generate primaries
void PrimaryGenerator::GeneratePrimaries(G4Event *anEvent) {

  // initialise particle table
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

  // define particle and get properties from particle table
  G4String particleName = "gamma";
  G4ParticleDefinition *particle = particleTable->FindParticle(particleName);

  // change angle direction and to radians
  double radAng = -scatAng / 180 * M_PI;

  // rotate source position coordinates
  G4double posRotX = pos[0] * std::cos(radAng) - pos[2] * std::sin(radAng);
  G4double posRotY = pos[1];
  G4double posRotZ = pos[0] * std::sin(radAng) + pos[2] * std::cos(radAng);
  G4ThreeVector posRot(posRotX*mm, posRotY*mm, posRotZ*mm);

  // check if momentum direction vector specified to random
  if (momRand == true && momRest == false ) {
    // generate components for random momentum direction vector
    G4double momX = 2 * G4UniformRand() - 1;
    G4double momY = 2 * G4UniformRand() - 1;
    G4double momZ = 2 * G4UniformRand() - 1;
    mom = G4ThreeVector(momX, momY, momZ);
  }
  else if (momRand == true && momRest == true) {
    // generate components restricted to front hemisphere
    G4double momX = 2 * G4UniformRand() - 1;
    G4double momY = 2 * G4UniformRand() - 1;
    G4double momZ = G4UniformRand();
    mom = G4ThreeVector(momX, momY, momZ);
  }

  // rotate momentum direction vector
  G4double momRotX = mom[0] * std::cos(radAng) - mom[2] * std::sin(radAng);
  G4double momRotY = mom[1];
  G4double momRotZ = mom[0] * std::sin(radAng) + mom[2] * std::cos(radAng);
  G4ThreeVector momRot(momRotX, momRotY, momRotZ);

  // assign particle properties to particle gun
  fParticleGun->SetParticlePosition(posRot);
  fParticleGun->SetParticleMomentumDirection(momRot);
  fParticleGun->SetParticleEnergy(partE*keV);
  fParticleGun->SetParticleDefinition(particle);

  // generate primary vertex and pass it to event
  fParticleGun->GeneratePrimaryVertex(anEvent);

}
