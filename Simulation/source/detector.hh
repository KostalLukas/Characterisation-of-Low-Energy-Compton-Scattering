#ifndef DETECTOR_HH
#define DETECTOR_HH

// include all used header files
#include <fstream>
#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

// sensitive detector class
class SensitiveDetector : public G4VSensitiveDetector {
public:
  // constructor
  SensitiveDetector(G4String);
  // destructor
  ~SensitiveDetector();

private:
  // method to process particle hits in detector
  virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);

};

#endif
