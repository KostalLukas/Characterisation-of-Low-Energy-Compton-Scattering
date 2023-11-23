#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

// include all used header files
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4GenericMessenger.hh"
#include "G4VisAttributes.hh"

// include sensitive detetor header file
#include "detector.hh"

// detector construction class
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  // constructor
  DetectorConstruction();
  // destructor
  ~DetectorConstruction();

  // method to construct a physical volume
  virtual G4VPhysicalVolume *Construct();

  // initialise generic messenger
  G4GenericMessenger *fMessenger;

  // declare variables for generic messenger
  G4double scatAng, shieldAng, shieldDist;
  G4bool block, targ, shield;

  // declare varaible for rotating vector since redefined
  G4double posRotX, posRotY, posRotZ;

  // declare source and shield position vectors for rotating
  G4ThreeVector posHous, posCap, posShield, posRotHous, posRotCap, posRotShield;

  // declare variables for angles converted to radians
  double radScatAng, radShieldAng;

  // declare variables for angles for rotating physical volumes
  G4double rotShieldAng, rotScatAng;

  // declare rotation matrices for rotating physical volumes
  G4RotationMatrix *targRot, *shieldRot, *sourceRot;

  // declare pointer to nist database manager isntance
  G4NistManager *nist;

  // declare material pointers
  G4Material *worldMat, *targMat, *blockMat, *scintMat;

  // declare box solid volume pointers
  G4Box *solidWorld, *solidBlock, *solidShield, *solidScintBlock;

  // declare tubs solid volume pointers
  G4Tubs *solidSourceCavity, *solidSourceAperture, *solidSourceTub,
         *solidSourceCap1, *solidSourceCap2, *solidScintTub, *solidScintCap1,
         *solidScintCap2, *solidTarg, *solidScint;

  // declare variable solid volume pointers
  G4VSolid *solidSourceSub, *solidSourceHousUni1, *solidSourceHousUni2,
           *solidSourceHous, *solidSourceCapUni, *solidSourceCap,
           *solidScintCavity1, *solidScintCavity2, *solidScintAperture,
           *solidScintSub1, *solidScintSub2, *solidScintHous, *solidScintUni1,
           *solidScintCap;

  // declare logical volume pointers
  G4LogicalVolume *logicWorld, *logicTarg, *logicSourceHous, *logicSourceCap,
                  *logicScintHous, *logicScintCap, *logicShield, *logicScint;

  // declare physical volume pointers
  G4VPhysicalVolume *physWorld, *physTarg, *physSourceHous, *physSourceCap,
                    *physScintHous, *physScintCap, *physScint, *physShield;

  // visual attributes
  G4VisAttributes *rVis, *gVis, *yVis, *bVis;

  // method to define materials
  void DefineMaterials();

  // method to construct sensitive detector or fields
  virtual void ConstructSDandField();
};

#endif
