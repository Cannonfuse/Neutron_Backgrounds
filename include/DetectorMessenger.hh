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
/// \file B2bDetectorMessenger.hh
/// \brief Definition of the B2bDetectorMessenger class

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class CLYCDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Messenger class that defines commands for B2bDetectorConstruction.
///
/// It implements commands:
/// - /B2/det/setTargetMaterial name
/// - /B2/det/setChamberMaterial name
/// - /B2/det/stepMax value unit

class DetectorMessenger : public G4UImessenger
{
  public:
    DetectorMessenger(CLYCDetectorConstruction* );
    virtual ~DetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    CLYCDetectorConstruction*  fDetectorConstruction;

    G4UIdirectory*           fDetDirectory;
    G4UIdirectory*           fGasCellDirectory;

    G4UIcmdWithADoubleAndUnit* fC7LYCDistance;
    G4UIcmdWithADoubleAndUnit* fC7LYC_TopcapPosition;

    G4UIcmdWithADoubleAndUnit* fC6LYCDistance;
    G4UIcmdWithADoubleAndUnit* fC7LYC_X;
    G4UIcmdWithADoubleAndUnit* fC7LYC_Y;
    G4UIcmdWithADoubleAndUnit* fC6LYC_X;
    G4UIcmdWithADoubleAndUnit* fC6LYC_Y;
    G4UIcmdWithADoubleAndUnit* fGasCellPressure;
    G4UIcmdWithADoubleAndUnit* fGasCellPosition;
    G4UIcmdWithADoubleAndUnit* fGasCellLength;
    G4UIcmdWithADoubleAndUnit* fGasCellDiameter;
    G4UIcmdWithADoubleAndUnit* fBe9TgtThickness;
    G4UIcmdWithABool* fUseC6LYC;
    G4UIcmdWithABool* fUseC7LYC;
    G4UIcmdWithABool* fUseC6LYC_Case;
    G4UIcmdWithABool* fUseC7LYC_Case;
    G4UIcmdWithABool* fUseStructure;
    G4UIcmdWithABool* fUseDummy;
    G4UIcmdWithABool* fUseBe9Target;
    G4UIcmdWithABool* fUseLargeChamber;
    G4UIcmdWithABool* fUseGasCell;
    G4UIcmdWithABool* fUseLTC;
    G4UIcmdWithABool* fUseMTC;
    G4UIcmdWithABool* fUseFTC;
    G4UIcmdWithAnInteger *fC6LYCSlices;
    G4UIcmdWithAnInteger *fC7LYCSlices;
    G4UIcmdWithAString *fWorldMaterial;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
