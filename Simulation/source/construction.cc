// include construction header file
#include "construction.hh"

// detector constructor
DetectorConstruction::DetectorConstruction() {

    // define messenger command directory
    fMessenger = new G4GenericMessenger(this, "/setup/", "prepare experimental setup");

    // define messenger commands
    fMessenger->DeclareProperty("scatAng", scatAng, "set scattering angle");
    fMessenger->DeclareProperty("block", block, "block source output");
    fMessenger->DeclareProperty("targ", targ, "construct scattering target");
    fMessenger->DeclareProperty("shield", shield, "construct shielding block");
    fMessenger->DeclareProperty("shieldAng", shieldAng, "set shielding block angle");
    fMessenger->DeclareProperty("shieldDist", shieldDist, "set shielding block edge distance");

    // specify messenger variable generic values
    scatAng = 30;
    block = false;
    targ = true;
    shield = true;
    shieldAng = 105; //120;
    shieldDist = 18;

    // call function to define materials
    DefineMaterials();
}

// detector destructor
DetectorConstruction::~DetectorConstruction() {}

// method for defining materials
void DetectorConstruction::DefineMaterials() {

  // get instance of nist material database manager
  nist = G4NistManager::Instance();

  // get materials from nist manager
  worldMat = nist->FindOrBuildMaterial("G4_AIR");
  targMat = nist->FindOrBuildMaterial("G4_Al");
  blockMat = nist->FindOrBuildMaterial("G4_Pb");

  // create scintillator NaI(Ti) material
  scintMat = new G4Material("scintMat", 3670*kg/m3, 2);
  scintMat->AddElement(nist->FindOrBuildElement("Na"), 1);
  scintMat->AddElement(nist->FindOrBuildElement("I"), 1);
}

// method for constucting detector volume
G4VPhysicalVolume *DetectorConstruction::Construct() {

  // print current geometry setup
  G4cout << G4endl;
  G4cout << "#---------------------------------------------------------------------#" << G4endl;
  G4cout << "scattering angle   = " << scatAng << " deg" << G4endl;
  G4cout << "shielding angle    = " << shieldAng << " deg" << G4endl;
  G4cout << "shielding distance = " << shieldDist << " mm" << G4endl;
  G4cout << "block source       = " << block << G4endl;
  G4cout << "#---------------------------------------------------------------------#" << G4endl;
  G4cout << G4endl;

  // specify source and shield position vectors for rotating
  posHous = G4ThreeVector(0, 0, -175*mm);
  posCap = G4ThreeVector(0, 0, -189*mm);
  posShield = G4ThreeVector(0, 0, -(shieldDist + 100)*mm);

  // change angle direction and to radians
  radScatAng = -scatAng / 180 * M_PI;
  radShieldAng = -shieldAng / 180 * M_PI;

  // rotate source housing position coordinates
  posRotX = posHous[0] * std::cos(radScatAng) - posHous[2] * std::sin(radScatAng);
  posRotY = posHous[1];
  posRotZ = posHous[0] * std::sin(radScatAng) + posHous[2] * std::cos(radScatAng);
  posRotHous = G4ThreeVector(posRotX*mm, posRotY*mm, posRotZ*mm);

  // rotate source housing position coordinates
  posRotX = posCap[0] * std::cos(radScatAng) - posCap[2] * std::sin(radScatAng);
  posRotY = posCap[1];
  posRotZ = posCap[0] * std::sin(radScatAng) + posCap[2] * std::cos(radScatAng);
  posRotCap = G4ThreeVector(posRotX*mm, posRotY*mm, posRotZ*mm);

  // rotate shield block position coordinates
  posRotX = posShield[0] * std::cos(radShieldAng) - posShield[2] * std::sin(radShieldAng);
  posRotY = posShield[1];
  posRotZ = posShield[0] * std::sin(radShieldAng) + posShield[2] * std::cos(radShieldAng);
  posRotShield = G4ThreeVector(posRotX*mm, posRotY*mm, posRotZ*mm);

  // create world solid volume realisstically 200*mm, 50*mm, 300*mm
  solidWorld = new G4Box("solidWorld", 400*mm, 50*mm, 500*mm);

  // create universal block solid volume
  solidBlock = new G4Box("solidBlock", 48.5*mm, 50*mm, 25*mm);

  // create shielding block solid volume
  solidShield = new G4Box("solidBlock", 25*mm, 50*mm, 100*mm);

  // create source housing solid volume
  solidSourceCavity = new G4Tubs("solidSourceCavity", 0, 6*mm, 16*mm, 0, 2*M_PI*rad);
  solidSourceAperture = new G4Tubs("solidSourceAperture", 0, 4*mm, 11*mm, 0, 2*M_PI*rad);
  solidSourceSub = new G4UnionSolid("solidSourceSub", solidSourceCavity, solidSourceAperture, 0, G4ThreeVector(0,0,27*mm));
  solidSourceHousUni1 = new G4SubtractionSolid("solidSourceHousUni1", solidBlock, solidSourceSub, 0, G4ThreeVector(0,0,-10*mm));
  solidSourceHousUni2 = new G4UnionSolid("solidSourceHousUni2", solidSourceHousUni1, solidBlock, 0, G4ThreeVector(0,0,-54*mm));

  // if specified also add volume to block source
  if (block == true) {
    solidSourceHous = new G4UnionSolid("solidSource", solidSourceHousUni2, solidBlock, 0, G4ThreeVector(0,0,52*mm));
  }
  else {
    solidSourceHous = solidSourceHousUni2;
  }

  // create source cap solid volume
  solidSourceTub = new G4Tubs("solidSourceTub", 1*mm, 3.8*mm, 11*mm, 0, 2*M_PI*rad);
  solidSourceCap1 = new G4Tubs("solidSourceCap1", 0, 3.8*mm, 1*mm, 0, 2*M_PI*rad);
  solidSourceCap2 = new G4Tubs("solidSourceCap2", 0, 10*mm, 2*mm, 0, 2*M_PI*rad);
  solidSourceCapUni = new G4UnionSolid("solidSourceCapUni", solidSourceTub, solidSourceCap1, 0, G4ThreeVector(0,0,12*mm));
  solidSourceCap = new G4UnionSolid("solidSourceCap", solidSourceCapUni, solidSourceCap2, 0, G4ThreeVector(0,0,-13*mm));

  // create scintillator housing solid volume
  solidScintBlock = new G4Box("solidBlock", 48.5*mm, 50*mm, 50*mm);
  solidScintCavity1 = new G4Tubs("solidScintCavity1", 0, 22.3*mm, 22.5*mm, 0, 2*M_PI*rad);
  solidScintCavity2 = new G4Tubs("solidScintCavity2", 0, 31*mm, 21*mm, 0, 2*M_PI*rad);
  solidScintAperture = new G4Tubs("solidScintAperture", 0, 10*mm, 8*mm, 0, 2*M_PI*rad);
  solidScintSub1 = new G4UnionSolid("solidScintSub1", solidScintCavity1, solidScintAperture, 0, G4ThreeVector(0,0,-30.5*mm));
  solidScintSub2 = new G4UnionSolid("solidScintSub2", solidScintSub1, solidScintCavity2, 0, G4ThreeVector(0,0,42.5*mm));
  solidScintHous = new G4SubtractionSolid("solidScintBlock", solidScintBlock, solidScintSub2, 0, G4ThreeVector(0,0,-12.5));

  // create scintillator cap solid volume
  solidScintTub = new G4Tubs("solidScintTub", 19.5*mm, 21.5*mm, 25.5*mm, 0, 2*M_PI*rad);
  solidScintCap1 = new G4Tubs("solidScintCap1", 0, 21.5*mm, 1*mm, 0, 2*M_PI*rad);
  solidScintCap2 = new G4Tubs("solidScintCap2", 0, 21.5*mm, 0.2*mm, 0, 2*M_PI*rad);
  solidScintUni1 = new G4UnionSolid("solidScintUni1", solidScintTub, solidScintCap1, 0, G4ThreeVector(0,0,26.5*mm));
  solidScintCap = new G4UnionSolid("solidScintCap", solidScintUni1, solidScintCap2, 0, G4ThreeVector(0,0,-25.7*mm));

  // create cylindrical target volume
  solidTarg = new G4Tubs("solidTarg", 0, 10*mm, 50*mm, 0, 2*M_PI*rad);

  // crete scintillator crystal solid volume
  solidScint = new G4Tubs("solidScint", 0, 19.05*mm, 25.4*mm, 0, 2*M_PI*rad);

  // create logical volumes
  logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
  logicTarg = new G4LogicalVolume(solidTarg, targMat, "logicTarg");
  logicScint = new G4LogicalVolume(solidScint, scintMat, "logicScint");
  logicSourceHous = new G4LogicalVolume(solidSourceHous, blockMat, "logicSource");
  logicSourceCap = new G4LogicalVolume(solidSourceCap, targMat, "logicSourceCap");
  logicScintHous = new G4LogicalVolume(solidScintHous, blockMat, "logicScintHous");
  logicScintCap = new G4LogicalVolume(solidScintCap, targMat, "logicScintCap");
  logicShield = new G4LogicalVolume(solidShield, blockMat, "logicShield");

  // create visualisation attributes
  rVis = new G4VisAttributes(G4Colour(1, 0.2, 0.2));
  gVis = new G4VisAttributes(G4Colour(0.2, 1, 0.2));
  bVis = new G4VisAttributes(G4Colour(0.2, 0.2, 1));
  yVis = new G4VisAttributes(G4Colour(1, 1, 0.2));

  // set visualisation attributes
  logicSourceCap->SetVisAttributes(rVis);
  logicScint->SetVisAttributes(bVis);
  logicScintCap->SetVisAttributes(yVis);
  logicTarg->SetVisAttributes(gVis);

  // create world physical volume
  physWorld = new G4PVPlacement(0, G4ThreeVector(0,0,0), logicWorld, "physWorld", 0, false, 0, true);

  // if specified create scattering target physical volume
  if (targ == true) {
    // create rotation matrxi for rotating target volume
    targRot = new G4RotationMatrix();
    targRot->rotateX(90*deg);
    targRot->rotateY(0*deg);
    targRot->rotateZ(0*deg);

    // create target physical volume
    physTarg = new G4PVPlacement(targRot, G4ThreeVector(0,0,0), logicTarg, "physTarg", logicWorld, false, 0, true);
  }

  // create rotation matrxi for rotating source housing volume
  sourceRot = new G4RotationMatrix();
  rotScatAng = radScatAng;
  sourceRot->rotateX(0*deg);
  sourceRot->rotateY(rotScatAng*rad);
  sourceRot->rotateZ(0*deg);

  // create source housing physical volume
  physSourceHous = new G4PVPlacement(sourceRot, posRotHous, logicSourceHous, "physSourceHous", logicWorld, false, 0, true);

  // create source cap physical volume
  physSourceCap = new G4PVPlacement(sourceRot, posRotCap, logicSourceCap, "physSourceCap", logicWorld, false, 0, true);

  // create scintillator housing physical volume
  physScintHous = new G4PVPlacement(0, G4ThreeVector(0,0,200*mm), logicScintHous, "physScintHous", logicWorld, false, 0, true);

  // create scintillator cap physical volume
  physScintCap = new G4PVPlacement(0, G4ThreeVector(0,0,193.4*mm), logicScintCap, "physScintHous", logicWorld, false, 0, true);

  // create scintillator physical volume
  physScint = new G4PVPlacement(0, G4ThreeVector(0,0,193.4*mm), logicScint, "physScint", logicWorld, false, 0, true);

  // if specified create shielding block physical volume
  if (shield == true) {
    // create rotation matrxi for rotating shielding block volume
    shieldRot = new G4RotationMatrix();
    rotShieldAng = radShieldAng;
    shieldRot->rotateX(0*deg);
    shieldRot->rotateY(rotShieldAng*rad);
    shieldRot->rotateZ(0*deg);

    // create shielding block physical volume
    physShield = new G4PVPlacement(shieldRot, posRotShield, logicShield, "physShield", logicWorld, false, 0, true);
  }

  return physWorld;
}

// method for constructing sensitive detector
void DetectorConstruction::ConstructSDandField() {

  // initialise instance of sensitive detector
  SensitiveDetector *sensDet = new SensitiveDetector("SensitiveDetector");

  // set logic detector volume to sensitive detector
  logicScint->SetSensitiveDetector(sensDet);
}
