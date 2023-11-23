#include <iostream>

// include used header files
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"

// include simulation header files
#include "construction.hh"
#include "physics.hh"
#include "action.hh"

int main(int argc, char** argv) {

	#ifdef G4MULTITHREADED
		G4MTRunManager *runManager = new G4MTRunManager();
	#else
		G4RunManager *runManager = new G4RunManager();
	#endif

	// initialsie new detector construction
	runManager->SetUserInitialization(new DetectorConstruction());

	// initialsie physics list
	runManager->SetUserInitialization(new PhysicsList());

	// initialsie action initialisation
	runManager->SetUserInitialization(new ActionInitialization());

	// initialise user interface
	G4UIExecutive *ui = 0;

	// if simulation is executed without arguments run with user interface
	if(argc == 1) {
		ui = new G4UIExecutive(argc, argv);
	}

	// initialise user interface manager
	G4UImanager *uiManager = G4UImanager::GetUIpointer();

	// if user interface initialised set it up
	if(ui) {
		// initialise run manager
		runManager->Initialize();

		// initialise visualisation manager
		G4VisManager *visManager = new G4VisExecutive("quite");
		visManager->Initialize();

		// create OGL user interface
		uiManager->ApplyCommand("/vis/open OGL");

		// set user interface to refresh at each event
		uiManager->ApplyCommand("/vis/set/autoRefresh true");

		// set initial view vector
		uiManager->ApplyCommand("/vis/viewer/set/viewpointVector 0 0 0");

		// show defined volumes
		uiManager->ApplyCommand("/vis/drawVolume");
		uiManager->ApplyCommand("/vis/geometry/set/forceCloud all -1 true 100000");

		// show particle trajectories
		uiManager->ApplyCommand("/vis/scene/add/trajectories smooth");

		// draw axis with m units and scale
		uiManager->ApplyCommand("/vis/scene/add/axes 0 0 0 100 mm");
		uiManager->ApplyCommand("/vis/scene/add/scale 100 mm");

		// show event id
		uiManager->ApplyCommand("/vis/scene/add/eventID");

		// accumulate at the end of action
		uiManager->ApplyCommand("vis/scene/endOfEventAction accumulate");

		// start user interface session
		ui->SessionStart();
	}

	// otherwise run batch mode with macros as argument
	else {
		G4String macName = argv[1];
		uiManager->ApplyCommand("/control/execute " + macName);
	}

	return 0;
}
