#include "CLYCDetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "G4OpticalPhysics.hh"
#include "G4PhysListFactory.hh"
#include "PhysicsList.hh"

#include "G4UImanager.hh"
#include "QBBC.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

#include "G4ParticleHPManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// const G4double det_dist = 5 * CLHEP::m;

int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Optionally: choose a different Random engine...
  G4Random::setTheEngine(new CLHEP::MTwistEngine);
  
  // Construct the default run manager
  //

  CLHEP::HepRandom::setTheSeed((unsigned)clock()); 

  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);


  // Set mandatory initialization classes
  //

    // runManager->SetUserInitialization(new CLYCDetectorConstruction);


  // Physics list QGSP_BERT_HP
    G4PhysListFactory *physListFactory = new G4PhysListFactory(); 
    G4VModularPhysicsList *physicsList =  physListFactory->GetReferencePhysList("QGSP_BERT_HP");   
    runManager->SetUserInitialization(physicsList);
    physicsList->SetVerboseLevel(2);


    // // const std::vector<G4String> v = physListFactory->AvailablePhysLists();
    // // for(auto i = 0; i< v.size(); i++)
    // // {
    // //   G4cout << v.at(i) << G4endl;
    // // }
    // physicsList->SetVerboseLevel(2);
    // physicsList->RemovePhysics(physicsList->GetPhysics("G4RadioactiveDecayPhysics"));
    // runManager->SetUserInitialization(physicsList);// to set parameters in code, if wanted

  // G4VUserPhysicsList *physicsList =  physListFactory->GetReferencePhysList("FTFP_BERT_HP"); 


////////////////////////////////////////////////////////////////////////////////////

  // Physics list QGSP_BERT_HP_MOD - removed G4RadioactiveDecayPhysics

////////////////////////////////////////////////////////////////////////////////////
  // PhysicsList* phys = new PhysicsList;
  // runManager->SetUserInitialization(phys);
  // phys->SetVerboseLevel(2);
////////////////////////////////////////////////////////////////////////////////////

  // Replaced HP environmental variables with C++ calls
  // G4ParticleHPManager::GetInstance()->SetSkipMissingIsotopes( false );
  // G4ParticleHPManager::GetInstance()->SetDoNotAdjustFinalState( false );
  // G4ParticleHPManager::GetInstance()->SetUseOnlyPhotoEvaporation( false );
  // G4ParticleHPManager::GetInstance()->SetNeglectDoppler( false );
  // G4ParticleHPManager::GetInstance()->SetProduceFissionFragments( false );
  // G4ParticleHPManager::GetInstance()->SetUseWendtFissionModel( false );
  // G4ParticleHPManager::GetInstance()->SetUseNRESP71Model( false );
  //G4ParticleHPManager::GetInstance()->SetSkipMissingIsotopes( true );
  //G4ParticleHPManager::GetInstance()->SetDoNotAdjustFinalState( true );
  //G4ParticleHPManager::GetInstance()->SetUseOnlyPhotoEvaporation( true );
  //G4ParticleHPManager::GetInstance()->SetNeglectDoppler( true );
  //G4ParticleHPManager::GetInstance()->SetProduceFissionFragments( true );
  //G4ParticleHPManager::GetInstance()->SetUseWendtFissionModel( true );
  //G4ParticleHPManager::GetInstance()->SetUseNRESP71Model( true );

  //Old Physics List
  // G4VModularPhysicsList* physicsList = new QBBC;
  // physicsList->SetVerboseLevel(1);
  // G4OpticalPhysics*opticalPhysics =new G4OpticalPhysics();
  // physicsList->RegisterPhysics(opticalPhysics);
  // auto opticalParams = G4OpticalParameters::Instance();
  // opticalParams->SetWLSTimeProfile("delta");
  // runManager->SetUserInitialization(physicsList);// to set parameters in code, if wanted

  // User action initialization
  // Detector construction

  auto detConstruction = new CLYCDetectorConstruction();

  runManager->SetUserInitialization(detConstruction);
  runManager->SetUserInitialization(new ActionInitialization(detConstruction));
  
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  
  delete visManager;
  delete runManager;
}

/*
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"


#include "CLYCDetectorConstruction.hh"
// #include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

int main(int argc, char** argv)
{
  // construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  runManager->SetUserInitialization(new CLYCDetectorConstruction);
  // runManager->SetUserInitialization(new PhysicsList);

  // set mandatory user action class
  runManager->SetUserAction(new PrimaryGeneratorAction);

  // initialize G4 kernel
  runManager->Initialize();

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( argc == 1 ) 
  {
    // interactive mode : define UI session
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/execute init.mac");
    ui->SessionStart();
    delete ui;
  }
  else 
  {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  // job termination
  delete visManager;
  delete runManager;
  return 0;
}

*/