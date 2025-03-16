//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// 
/// \file B2bDetectorMessenger.cc
/// \brief Implementation of the B2bDetectorMessenger class

#include "DetectorMessenger.hh"
#include "CLYCDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(CLYCDetectorConstruction* Det)
 : G4UImessenger(),
   fDetectorConstruction(Det)
{
  fDetDirectory = new G4UIdirectory("/Detectors/");
  fDetDirectory->SetGuidance("Detector construction control");

  fC7LYCDistance = new G4UIcmdWithADoubleAndUnit("/Detectors/C7LYCDistance",this);
  fC7LYCDistance->SetGuidance("Define the distance of the C7LYC detector");
  fC7LYCDistance->SetParameterName("C7LYCDistance",false);
  fC7LYCDistance->SetUnitCategory("Length");
  fC7LYCDistance->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fC7LYC_TopcapPosition = new G4UIcmdWithADoubleAndUnit("/Detectors/C7LYCTopcapPosition",this);
  fC7LYC_TopcapPosition->SetGuidance("Define the distance of the C7LYC detector topcap");
  fC7LYC_TopcapPosition->SetParameterName("C7LYCTopcapPosition",true);
  fC7LYC_TopcapPosition->SetUnitCategory("Length");
  fC7LYC_TopcapPosition->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fC6LYCDistance = new G4UIcmdWithADoubleAndUnit("/Detectors/C6LYCDistance",this);
  fC6LYCDistance->SetGuidance("Define the distance of the C6LYC detector");
  fC6LYCDistance->SetParameterName("C6LYCDistance",false);
  fC6LYCDistance->SetUnitCategory("Length");
  fC6LYCDistance->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fC7LYC_X = new G4UIcmdWithADoubleAndUnit("/Detectors/C7LYC_X",this);
  fC7LYC_X->SetGuidance("Define the X position of the C7LYC detector relative to the beam axis");
  fC7LYC_X->SetParameterName("C7LYC_X",false);
  fC7LYC_X->SetUnitCategory("Length");
  fC7LYC_X->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fC7LYC_Y = new G4UIcmdWithADoubleAndUnit("/Detectors/C7LYC_Y",this);
  fC7LYC_Y->SetGuidance("Define the Y position of the C7LYC detector relative to the beam axis");
  fC7LYC_Y->SetParameterName("C7LYC_Y",false);
  fC7LYC_Y->SetUnitCategory("Length");
  fC7LYC_Y->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fC6LYC_X = new G4UIcmdWithADoubleAndUnit("/Detectors/C6LYC_X",this);
  fC6LYC_X->SetGuidance("Define the X position of the C6LYC detector relative to the beam axis");
  fC6LYC_X->SetParameterName("C6LYC_X",false);
  fC6LYC_X->SetUnitCategory("Length");
  fC6LYC_X->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fC6LYC_Y = new G4UIcmdWithADoubleAndUnit("/Detectors/C6LYC_Y",this);
  fC6LYC_Y->SetGuidance("Define the Y position of the C6LYC detector relative to the beam axis");
  fC6LYC_Y->SetParameterName("C6LYC_Y",false);
  fC6LYC_Y->SetUnitCategory("Length");
  fC6LYC_Y->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseC6LYC = new G4UIcmdWithABool("/Detectors/UseC6LYC", this);
  fUseC6LYC->SetGuidance("Enable the Li6-enriched CLYC detector (1 in x 1 in)");
  fUseC6LYC->SetParameterName("UseC6LYC",false);
  fUseC6LYC->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseC6LYC_Case = new G4UIcmdWithABool("/Detectors/UseC6LYC_Case", this);
  fUseC6LYC_Case->SetGuidance("Enable the Li6-enriched CLYC detector casing");
  fUseC6LYC_Case->SetParameterName("UseC6LYC_Case",true);
  fUseC6LYC_Case->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseC7LYC = new G4UIcmdWithABool("/Detectors/UseC7LYC", this);
  fUseC7LYC->SetGuidance("Enable the Li6-depleted CLYC detector (75 mm x 10 mm)");
  fUseC7LYC->SetParameterName("UseC7LYC",false);
  fUseC7LYC->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseC7LYC_Case = new G4UIcmdWithABool("/Detectors/UseC7LYC_Case", this);
  fUseC7LYC_Case->SetGuidance("Enable the Li6-depleted CLYC detector casing");
  fUseC7LYC_Case->SetParameterName("UseC7LYC_Case",true);
  fUseC7LYC_Case->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fC6LYCSlices = new G4UIcmdWithAnInteger("/Detectors/C6LYC_Slices", this);
  fC6LYCSlices->SetGuidance("Slices the C6LYC detector into n miniature detectors along the beam axis");
  fC6LYCSlices->SetParameterName("C6LYC_Slices",true);
  fC6LYCSlices->SetDefaultValue(1);
  fC6LYCSlices->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fC7LYCSlices = new G4UIcmdWithAnInteger("/Detectors/C7LYC_Slices", this);
  fC7LYCSlices->SetGuidance("Slices the C7LYC detector into n miniature detectors along the beam axis");
  fC7LYCSlices->SetParameterName("C7LYC_Slices",true);
  fC7LYCSlices->SetDefaultValue(1);
  fC7LYCSlices->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseStructure = new G4UIcmdWithABool("/Detectors/UseStructure", this);
  fUseStructure->SetGuidance("Enables the concrete structure of the TOF facility");
  fUseStructure->SetParameterName("UseStructure",true);
  fUseStructure->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);
  
  fUseDummy = new G4UIcmdWithABool("/Detectors/UseDummy", this);
  fUseDummy->SetGuidance("Enables the dummy detector for checking the subtended solid angle inside the tunnel");
  fUseDummy->SetParameterName("UseDummy",true);
  fUseDummy->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseBe9Target = new G4UIcmdWithABool("/Detectors/UseBe9Target", this);
  fUseBe9Target->SetGuidance("Enables Be9 target inside of the target chamber");
  fUseBe9Target->SetParameterName("UseBe9Target",true);
  fUseBe9Target->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseLargeChamber = new G4UIcmdWithABool("/Detectors/UseLargeChamber", this);
  fUseLargeChamber->SetGuidance("Enables Large Target Chamber");
  fUseLargeChamber->SetParameterName("UseLargeChamber",true);
  fUseLargeChamber->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseGasCell = new G4UIcmdWithABool("/Detectors/UseGasCell", this);
  fUseGasCell->SetGuidance("Enables Gas Cell");
  fUseGasCell->SetParameterName("UseGasCell",true);
  fUseGasCell->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseLTC = new G4UIcmdWithABool("/Detectors/UseLTC", this);
  fUseLTC->SetGuidance("Enable the large tunnel collimator");
  fUseLTC->SetParameterName("UseBe9Target",true);
  fUseLTC->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseMTC = new G4UIcmdWithABool("/Detectors/UseMTC", this);
  fUseMTC->SetGuidance("Enables the medium tunnel collimator");
  fUseMTC->SetParameterName("UseMTC",true);
  fUseMTC->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fUseFTC = new G4UIcmdWithABool("/Detectors/UseFTC", this);
  fUseFTC->SetGuidance("Enables the fine tunnel collimator");
  fUseFTC->SetParameterName("UseFTC",true);
  fUseFTC->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fGasCellDirectory = new G4UIdirectory("/GasCell/");
  fGasCellDirectory->SetGuidance("Gas Cell Control");

  fGasCellPressure = new G4UIcmdWithADoubleAndUnit("/GasCell/GasCellPressure",this);
  fGasCellPressure->SetGuidance("Define a gas cell pressure in pascal");
  fGasCellPressure->SetParameterName("GasCellPressure",false);
  fGasCellPressure->SetUnitCategory("Pressure");
  fGasCellPressure->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fGasCellPosition = new G4UIcmdWithADoubleAndUnit("/GasCell/GasCellPosition",this);
  fGasCellPosition->SetGuidance("Define the center of the gas cell volume");
  fGasCellPosition->SetParameterName("GasCellPosition",false);
  fGasCellPosition->SetUnitCategory("Length");
  fGasCellPosition->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fGasCellLength = new G4UIcmdWithADoubleAndUnit("/GasCell/GasCellLength",this);
  fGasCellLength->SetGuidance("Define the length of the gas volume");
  fGasCellLength->SetParameterName("GasCellLength",false);
  fGasCellLength->SetUnitCategory("Length");
  fGasCellLength->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fGasCellDiameter = new G4UIcmdWithADoubleAndUnit("/GasCell/GasCellDiameter",this);
  fGasCellDiameter->SetGuidance("Define the internal diameter of the gas cell");
  fGasCellDiameter->SetParameterName("GasCellDiameter",false);
  fGasCellDiameter->SetUnitCategory("Length");
  fGasCellDiameter->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);

  fWorldMaterial = new G4UIcmdWithAString("/Detectors/WorldMaterial", this);
  fWorldMaterial->SetGuidance("Sets the WorldMaterial the world is filled with");
  fWorldMaterial->SetDefaultValue("G4_Galactic");
  fWorldMaterial->SetParameterName("WorldMaterial",true);
  fWorldMaterial->AvailableForStates(G4State_PreInit);//,G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fDetDirectory;
  delete fGasCellDirectory;
  delete fC7LYCDistance;
  delete fC7LYC_TopcapPosition;
  delete fC6LYCDistance;
  delete fC7LYC_X;
  delete fC7LYC_Y;
  delete fC6LYC_X;
  delete fC6LYC_Y;
  delete fGasCellPressure;
  delete fGasCellPosition;
  delete fGasCellLength;
  delete fGasCellDiameter;
  delete fUseC6LYC;
  delete fUseC7LYC;
  delete fUseStructure;
  delete fUseDummy;
  delete fUseBe9Target;
  delete fUseLargeChamber;
  delete fUseGasCell;
  delete fUseLTC;
  delete fUseMTC;
  delete fUseFTC;
  delete fC6LYCSlices;
  delete fC7LYCSlices;
  delete fWorldMaterial;
  delete fUseC6LYC_Case;
  delete fUseC7LYC_Case;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

  if( command == fC7LYCDistance ) {
    fDetectorConstruction->SetC7LYCDistance(fC7LYCDistance->GetNewDoubleValue(newValue));
  }
  if( command == fC7LYC_TopcapPosition ) {
    fDetectorConstruction->SetC7LYC_TopcapPosition(fC7LYC_TopcapPosition->GetNewDoubleValue(newValue));
  }
  if( command == fC6LYCDistance ) {
    fDetectorConstruction->SetC6LYCDistance(fC6LYCDistance->GetNewDoubleValue(newValue));
  }
  if( command == fC7LYC_X) {
    fDetectorConstruction->SetC7LYC_X(fC7LYC_X->GetNewDoubleValue(newValue));
  }
  if( command == fC7LYC_Y) {
    fDetectorConstruction->SetC7LYC_Y(fC7LYC_Y->GetNewDoubleValue(newValue));
  }
    if( command == fC6LYC_X) {
    fDetectorConstruction->SetC6LYC_X(fC6LYC_X->GetNewDoubleValue(newValue));
  }
  if( command == fC6LYC_Y) {
    fDetectorConstruction->SetC6LYC_Y(fC6LYC_Y->GetNewDoubleValue(newValue));
  }
  if( command == fUseC6LYC) {
    fDetectorConstruction->SetUseC6LYC(fUseC6LYC->GetNewBoolValue(newValue));
  } 
  if( command == fUseC7LYC) {
    fDetectorConstruction->SetUseC7LYC(fUseC7LYC->GetNewBoolValue(newValue));
  } 
  if( command == fUseStructure) {
    fDetectorConstruction->SetUseStructure(fUseStructure->GetNewBoolValue(newValue));
  } 
  if( command == fUseDummy) {
    fDetectorConstruction->SetUseDummy(fUseDummy->GetNewBoolValue(newValue));
  } 
  if( command == fUseBe9Target) {
    fDetectorConstruction->SetUseBe9target(fUseBe9Target->GetNewBoolValue(newValue));
  }
  if( command == fUseLargeChamber) {
    fDetectorConstruction->SetUseLargeTarget(fUseLargeChamber->GetNewBoolValue(newValue));
  }
  if( command == fUseGasCell) {
    fDetectorConstruction->SetUseGasCell(fUseGasCell->GetNewBoolValue(newValue));
  }
  if( command == fUseLTC) {
    fDetectorConstruction->SetUseLTC(fUseLTC->GetNewBoolValue(newValue));
  } 
  if( command == fUseMTC) {
    fDetectorConstruction->SetUseMTC(fUseMTC->GetNewBoolValue(newValue));
  } 
  if( command == fUseFTC) {
    fDetectorConstruction->SetUseFTC(fUseFTC->GetNewBoolValue(newValue));
  }  
  if( command == fGasCellPressure ) {
    fDetectorConstruction->setGasCellPressure(fGasCellPressure->GetNewDoubleValue(newValue));
  }  
  if( command ==  fGasCellPosition) {
    fDetectorConstruction->setGasCellPosition(fGasCellPosition->GetNewDoubleValue(newValue));
  }  
  if( command ==  fGasCellLength) {
    fDetectorConstruction->setGasCellLength(fGasCellLength->GetNewDoubleValue(newValue));
  }  
  if( command ==  fGasCellDiameter) {
    fDetectorConstruction->setGasCellDiameter(fGasCellDiameter->GetNewDoubleValue(newValue));
  }  
  if( command ==  fC6LYCSlices) {
    fDetectorConstruction->setC6LYC_Slices(fC6LYCSlices->GetNewIntValue(newValue));
  }  
  if( command ==  fC7LYCSlices) {
    fDetectorConstruction->setC7LYC_Slices(fC7LYCSlices->GetNewIntValue(newValue));
  }  
  if( command ==  fWorldMaterial) {
    fDetectorConstruction->setWorldMaterial(newValue);
  }
  if( command == fUseC6LYC_Case) {
    fDetectorConstruction->SetUseC6LYC_Case(fUseC6LYC_Case->GetNewBoolValue(newValue));
  } 
  if( command == fUseC7LYC_Case) {
    fDetectorConstruction->SetUseC7LYC_Case(fUseC7LYC_Case->GetNewBoolValue(newValue));
  } 
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
