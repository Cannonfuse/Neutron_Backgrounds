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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "CLYCDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(const CLYCDetectorConstruction* detectorConstruction, EventAction* eventAction)
: G4UserSteppingAction(),
  fDetConstruction(detectorConstruction),
  fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

  // if (!fScoringVolume) { 
  //   const CLYCDetectorConstruction* detectorConstruction
  //     = static_cast<const CLYCDetectorConstruction*>
  //       (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  //   fScoringVolume = detectorConstruction->GetScoringVolume();   
  // }

  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  auto sl = step->GetStepLength();
  // auto deltae = step->GetDeltaEnergy();
  auto deltat = step->GetDeltaTime();
      
  // check if we are in scoring volume
  if (volume == fDetConstruction->GetCLYCVolume())
  {
    fEventAction->AddEdep(edep);
    fEventAction->AddLstep(sl);
    // fEventAction->AdddeltaE(deltae);
    fEventAction->AdddeltaT(deltat);
  }

  // collect energy deposited in this step
  // G4double edepStep = step->GetTotalEnergyDeposit();
  // G4double lStep = step->GetStepLength();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
