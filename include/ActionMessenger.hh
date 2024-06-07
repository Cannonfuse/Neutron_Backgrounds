#ifndef ActionMessenger_h
#define ActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class ActionInitialization;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Messenger class that defines commands for B2bDetectorConstruction.
///
/// It implements commands:
/// - /B2/det/setTargetMaterial name
/// - /B2/det/setChamberMaterial name
/// - /B2/det/stepMax value unit

class ActionMessenger: public G4UImessenger
{
  public:
    ActionMessenger(ActionInitialization* );
    virtual ~ActionMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    ActionInitialization*  fActInit;

    G4UIdirectory*           fActionDirectory;

    G4UIcmdWithAString* fEnergyAngleDist;
    G4UIcmdWithAString* fEnergyZDist;
    G4UIcmdWithAString* fEnergyAngleZBins;
    G4UIcmdWithAString* fNeutrons;



    G4UIcmdWithABool* fUseDists;
    G4UIcmdWithABool* fUseNeutrons;


};

#endif
