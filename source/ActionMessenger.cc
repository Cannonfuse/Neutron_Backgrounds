#include "ActionMessenger.hh"
#include "ActionInitialization.hh"
#include "G4String.hh"
// #include "CLYCDetectorConstruction.hh"

#include "G4UIdirectory.hh"
// #include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"


ActionMessenger::ActionMessenger(ActionInitialization* ActInit)
 : G4UImessenger(),
   fActInit(ActInit)
{
  fActionDirectory = new G4UIdirectory("/PartDist/");
  fActionDirectory->SetGuidance("Files for user particle distributions");

  fEnergyAngleDist = new G4UIcmdWithAString("/PartDist/EnergyAngleDist",this);
  fEnergyAngleDist->SetGuidance("Defines the 2D Energy Angle distribution file name");
  fEnergyAngleDist->SetParameterName("EnergyAngleDist",true);
  fEnergyAngleDist->SetDefaultValue("none");
  fEnergyAngleDist->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fEnergyZDist = new G4UIcmdWithAString("/PartDist/EnergyZDist",this);
  fEnergyZDist->SetGuidance("Defines the 2D Energy Z-position distribution file name");
  fEnergyZDist->SetParameterName("EnergyZDist",true);
  fEnergyZDist->SetDefaultValue("none");
  fEnergyZDist->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fEnergyAngleZBins = new G4UIcmdWithAString("/PartDist/EnergyAngleZBins",this);
  fEnergyAngleZBins->SetGuidance("Defines the Energy, Angle, and Z-position bins file name");
  fEnergyAngleZBins->SetParameterName("EnergyAngleZBins",true);
  fEnergyAngleZBins->SetDefaultValue("none");
  fEnergyAngleZBins->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fNeutrons = new G4UIcmdWithAString("/PartDist/NeutronsFile",this);
  fNeutrons->SetGuidance("Defines the Energy, Angle, and Z-position neutrons file name");
  fNeutrons->SetParameterName("NeutronsFile",true);
  fNeutrons->SetDefaultValue("none");
  fNeutrons->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseDists = new G4UIcmdWithABool("/PartDist/UseDists", this);
  fUseDists->SetGuidance("Enables the use of the 2D distributions for particle generation");
  fUseDists->SetParameterName("UseDists",true);
  fUseDists->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseNeutrons = new G4UIcmdWithABool("/PartDist/UseNeutrons", this);
  fUseNeutrons->SetGuidance("Enables the use of known neutrons for particle generation");
  fUseNeutrons->SetParameterName("UseNeutrons",true);
  fUseNeutrons->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

}

ActionMessenger::~ActionMessenger()
{
  delete fEnergyAngleDist;
  delete fEnergyZDist;
  delete fUseDists;
  delete fEnergyAngleZBins;
  delete fNeutrons;
  delete fUseNeutrons;

}

void ActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

  if( command == fEnergyAngleDist ) {
    fActInit->SetEnergyAngleDist(newValue);
  }
  if( command == fEnergyZDist ) {
    fActInit->SetEnergyZDist(newValue);
  }
  if( command == fEnergyAngleZBins ) {
    fActInit->SetEnergyAngleZBins(newValue);
  }
  if( command == fUseDists) {
    fActInit->SetUseDists(fUseDists->GetNewBoolValue(newValue));
  }
  if( command == fNeutrons ) {
    fActInit->SetNeutronsData(newValue);
  }
  if( command == fUseNeutrons) {
    fActInit->SetUseNeutrons(fUseNeutrons->GetNewBoolValue(newValue));
  }
}
