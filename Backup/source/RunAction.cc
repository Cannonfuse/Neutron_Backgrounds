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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "CLYCDetectorConstruction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "g4root.hh"
// #include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
: G4UserRunAction()
{ 
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  auto analysisManager = G4AnalysisManager::Instance();

  analysisManager->SetVerboseLevel(1);
  
  analysisManager->CreateH2("Edep_Lstep","Energy Deposited vs. Particle distance",4000,0.,0.8*MeV,250,0,25*mm);
  analysisManager->CreateH1("Lstep","Particle distance",10000,0,1000*mm);
  analysisManager->CreateH1("Edep","Energy Deposited",80000,0.,0.8*MeV);
  // analysisManager->CreateH1("Deltae","Particle Energy Change",10000,0.,10.*MeV);
  analysisManager->CreateH1("Deltat","Particle Time Change",10000,0.,1e6*microsecond);

  analysisManager->CreateNtuple("Data", "Data");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("Lstep");
  // analysisManager->CreateNtupleDColumn("Deltae");
  analysisManager->CreateNtupleDColumn("Deltat");
  analysisManager->FinishNtuple();

  analysisManager->SetNtupleMerging(true);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{ 
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile(); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  auto analysisManager = G4AnalysisManager::Instance();
  // analysisManager->GetH2(analysisManager->GetFirstH2Id())->
  // G4cout << "WRITING FILE"
  analysisManager->Write();
  analysisManager->CloseFile(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void RunAction::AddEdep(G4double edep)
// {
//   fEdep  += edep;
// }

// void RunAction::AddLstep(G4double lstep)
// {
//   fLstep  += lstep;
// }



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

