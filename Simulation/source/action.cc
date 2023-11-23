// include action header file
#include "action.hh"

// action initialisation constructor
ActionInitialization::ActionInitialization() {}

// action initialisation destructor
ActionInitialization::~ActionInitialization() {}

// define method to build simulation for master thread
void ActionInitialization::BuildForMaster() const {

  // initialise and set run action from run.hh
  RunAction *runAction = new RunAction();
  SetUserAction(runAction);
}

// define method to build simulation for worker thread
void ActionInitialization::Build() const {

  // create primary generator action
  PrimaryGenerator *generator = new PrimaryGenerator();
  SetUserAction(generator);

  // initialise and set run action from run.hh
  RunAction *runAction = new RunAction();
  SetUserAction(runAction);
}
