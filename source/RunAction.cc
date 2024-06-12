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
#include "EventAction.hh"
// #include "PrimaryGeneratorAction.hh"
// #include "CLYCDetectorConstruction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4RNGHelper.hh"
#include "G4Run.hh"
// #include "G4AnalysisManager.hh"
#include "g4root.hh"
// #include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <ctime>
// #include "RunMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction* eventAction)
: G4UserRunAction(),  fEventAction(eventAction)
{ 

  // fMessenger = new RunMessenger(this);

  // printf("RunAction it is %u\n",GetSaveVectors());

  G4RunManager::GetRunManager()->SetPrintProgress(1);

  auto analysisManager = G4AnalysisManager::Instance();

  analysisManager->SetVerboseLevel(1);
  
  // analysisManager->CreateH2("Edep_Lstep","Energy Deposited vs. Particle distance",150,0.,15*MeV,100,0,250*mm); //0
  // analysisManager->CreateH2("Edep_DeltaT","TOF vs. E",200,0,2000*nanosecond,200,0.,20*MeV); //1
  // analysisManager->CreateH2("GunEnergy_Edep","Neutron Energy vs. Energy Deposited",200,0.,20000*keV,200,0.,20000*keV); //2
  // analysisManager->CreateH2("GunEnergy_Lstep","Neutron Energy vs. Particle distance",200,0.,20*MeV,250,0,250*mm); //3
  // analysisManager->CreateH2("A_Z","A vs Z",200,0,200,400,0,400);//4
  // analysisManager->CreateH2("GunEnergy_FKE_C6LYC","Neutron Energy vs. Final Kinetic Energy (C6LYC)",2000,0.,20*MeV,2000,0.,20*MeV);//5
  // analysisManager->CreateH2("GunEnergy_FKE_C7LYC","Neutron Energy vs. Final Kinetic Energy (C7LYC)",2000,0.,20*MeV,2000,0.,20*MeV);//6
  // analysisManager->CreateH2("Step_FKE","Steps vs. Final Kinetic Energy",100,0,1000,200,0.,20*MeV);//7
  // analysisManager->CreateH2("GunEnergy_DeltaT","TOF vs. E",1000,0,2000*nanosecond,2000,0.,20*MeV);//8
  // analysisManager->CreateH2("InDetDeltaT_GunEnergy","TOF vs. E",100,0.,20*MeV,100,0,100*nanosecond);//9
  // analysisManager->CreateH2("InDetDeltaD_GunEnergy","TOF vs. E",100,0.,20*MeV,100,0,100*mm);//10


  // analysisManager->CreateH2("Edep_Lstep","Energy Deposited vs. Particle distance",150,0.,15*MeV,100,0,250*mm); //0
  // analysisManager->CreateH2("Edep_DeltaT","TOF vs. E",200,0,2000*nanosecond,200,0.,20*MeV); //1
  // analysisManager->CreateH2("GunEnergy_Edep","Neutron Energy vs. Energy Deposited",200,0.,20000*keV,200,0.,20000*keV); //2
  // analysisManager->CreateH2("GunEnergy_Lstep","Neutron Energy vs. Particle distance",200,0.,20*MeV,250,0,250*mm); //3
  // analysisManager->CreateH2("A_Z","A vs Z",200,0,200,400,0,400);//4
  analysisManager->CreateH2("GunEnergy_FKE_C6LYC","Neutron Energy vs. Final Kinetic Energy (C6LYC)",2000,0.,20*MeV,2000,0.,20*MeV);//0
  analysisManager->CreateH2("GunEnergy_FKE_C7LYC","Neutron Energy vs. Final Kinetic Energy (C7LYC)",2000,0.,20*MeV,2000,0.,20*MeV);//1
  analysisManager->CreateH2("C6LYC_Final_Position","XY Neutron Final Position",400,-200*mm,200*mm,400,-200*mm,200*mm);//2
  analysisManager->CreateH2("C7LYC_Final_Position","XY Neutron Final Position",400,-200*mm,200*mm,400,-200*mm,200*mm);//3
  analysisManager->CreateH2("C7LYC_Energy_Angle","Energy vs. Angle",100,1*MeV,5*MeV,100,0,M_PI);//4
  analysisManager->CreateH2("C6LYC_Slice_Distance","DetectorSlice vs. Distance Traveled",100,0,100,1000,0,100*cm);//5
  analysisManager->CreateH2("C7LYC_Slice_Distance","DetectorSlice vs. Distance Traveled",100,0,100,1000,0,100*cm);//6


  // analysisManager->CreateH2("Step_FKE","Steps vs. Final Kinetic Energy",100,0,1000,200,0.,20*MeV);//7
  // analysisManager->CreateH2("GunEnergy_DeltaT","TOF vs. E",1000,0,2000*nanosecond,2000,0.,20*MeV);//8
  // analysisManager->CreateH2("InDetDeltaT_GunEnergy","TOF vs. E",100,0.,20*MeV,100,0,100*nanosecond);//9
  // analysisManager->CreateH2("InDetDeltaD_GunEnergy","TOF vs. E",100,0.,20*MeV,100,0,100*mm);//10

  // analysisManager->CreateH3("Reaction_X_Y_Z","Reaction_Position",75,-3.75*cm,3.75*cm,75,-3.75*cm,3.75*cm,75,100*cm,101*cm);

  analysisManager->CreateH3("Start_X_Y_Z","Start_Position",75,-3.75*cm,3.75*cm,75,-3.75*cm,3.75*cm,75,-22.15*cm,-19.12*cm);





  analysisManager->CreateH1("Lstep","Particle distance",1000,0,250*mm);
  analysisManager->CreateH1("Edep","Energy Deposited",1500,0.,15*MeV);
  // analysisManager->CreateH1("Deltae","Particle Energy Change",10000,0.,10.*MeV);
  analysisManager->CreateH1("Deltat","Particle Time Change",10000,0., 10000*nanosecond);
  analysisManager->CreateH1("GunEnergy","Neutron Energy",1500,0.,15*MeV);

  analysisManager->CreateH1("C6LYC_InDetDeltaT","Particle Time Change",10000,0., 100*nanosecond);
  analysisManager->CreateH1("C6LYC_InDetDeltaD","Particle Distance",10000,0., 100*mm);
  analysisManager->CreateH1("C7LYC_InDetDeltaT","Particle Time Change",10000,0., 100*nanosecond);
  analysisManager->CreateH1("C7LYC_InDetDeltaD","Particle Distance",10000,0., 100*mm);
  analysisManager->CreateH1("C6LYC_FKE_NEQ_GE","Neutron Energy",2000,0., 20*MeV);
  analysisManager->CreateH1("C7LYC_FKE_NEQ_GE","Neutron Energy",2000,0., 20*MeV);
  analysisManager->CreateH1("C7LYC_TOF","Time of flight",1000,0., 1000*ns);


  analysisManager->CreateNtuple("Data", "Data");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("Lstep");
  analysisManager->CreateNtupleDColumn("Deltat");
  analysisManager->CreateNtupleDColumn("GunEnergy");
  analysisManager->CreateNtupleIColumn("Z");
  analysisManager->CreateNtupleIColumn("A");
  analysisManager->CreateNtupleDColumn("fKE");
  analysisManager->CreateNtupleIColumn("Steps");
  analysisManager->CreateNtupleDColumn("PosX");
  analysisManager->CreateNtupleDColumn("PosY");
  analysisManager->CreateNtupleDColumn("PosZ");
  analysisManager->CreateNtupleIColumn("Detector");
  analysisManager->CreateNtupleDColumn("PreDetKE");
  // analysisManager->CreateNtupleDColumn("VecX",eventAction->Get_pos_x_vector());
  // analysisManager->CreateNtupleDColumn("VecY",eventAction->Get_pos_y_vector());
  // analysisManager->CreateNtupleDColumn("VecZ",eventAction->Get_pos_z_vector());
  // analysisManager->CreateNtupleDColumn("VecLStep",eventAction->GetLstepVector());
  // analysisManager->CreateNtupleDColumn("VecDeltaT",eventAction->GetDeltaTVector());
  // analysisManager->CreateNtupleDColumn("VecEndKE",eventAction->GetEndKEVector());
  // analysisManager->CreateNtupleDColumn("VecEdep",eventAction->GetEdepVector());
  // analysisManager->CreateNtupleIColumn("VecParticleA",eventAction->GetAVector());
  // analysisManager->CreateNtupleIColumn("VecParticleZ",eventAction->GetZVector());

  analysisManager->FinishNtuple(0);

  analysisManager->CreateNtuple("Reacs_C6", "Reacs_C6");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("Lstep");
  analysisManager->CreateNtupleDColumn("Deltat");
  analysisManager->CreateNtupleDColumn("GunEnergy");
  analysisManager->CreateNtupleIColumn("Z");
  analysisManager->CreateNtupleIColumn("A");
  analysisManager->CreateNtupleDColumn("fKE");
  analysisManager->CreateNtupleIColumn("Steps");
  analysisManager->CreateNtupleDColumn("PosX");
  analysisManager->CreateNtupleDColumn("PosY");
  analysisManager->CreateNtupleDColumn("PosZ");
  analysisManager->CreateNtupleDColumn("PreDetKE");
  analysisManager->CreateNtupleIColumn("Detector");
  analysisManager->CreateNtupleDColumn("GlobalTime");
  analysisManager->CreateNtupleIColumn("Slice");

#ifdef SAVEVECTORS
  analysisManager->CreateNtupleDColumn("VecX",eventAction->Get_pos_x_vector());
  analysisManager->CreateNtupleDColumn("VecY",eventAction->Get_pos_y_vector());
  analysisManager->CreateNtupleDColumn("VecZ",eventAction->Get_pos_z_vector());
  analysisManager->CreateNtupleDColumn("VecLStep",eventAction->GetLstepVector());
  analysisManager->CreateNtupleDColumn("VecDeltaT",eventAction->GetDeltaTVector());
  analysisManager->CreateNtupleDColumn("VecEndKE",eventAction->GetEndKEVector());
  analysisManager->CreateNtupleDColumn("VecEdep",eventAction->GetEdepVector());
  analysisManager->CreateNtupleIColumn("VecParticleA",eventAction->GetAVector());
  analysisManager->CreateNtupleIColumn("VecParticleZ",eventAction->GetZVector());
  analysisManager->CreateNtupleIColumn("DetectorSlice",eventAction->GetSliceVector());
#endif


  analysisManager->FinishNtuple(1);

  analysisManager->CreateNtuple("Reacs_C7", "Reacs_C7");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("Lstep");
  analysisManager->CreateNtupleDColumn("Deltat");
  analysisManager->CreateNtupleDColumn("GunEnergy");
  analysisManager->CreateNtupleIColumn("Z");
  analysisManager->CreateNtupleIColumn("A");
  analysisManager->CreateNtupleDColumn("fKE");
  analysisManager->CreateNtupleIColumn("Steps");
  analysisManager->CreateNtupleDColumn("PosX");
  analysisManager->CreateNtupleDColumn("PosY");
  analysisManager->CreateNtupleDColumn("PosZ");
  analysisManager->CreateNtupleDColumn("PreDetKE");
  analysisManager->CreateNtupleIColumn("Detector");
  analysisManager->CreateNtupleDColumn("GlobalTime");
  analysisManager->CreateNtupleIColumn("Slice");

 #ifdef SAVEVECTORS
  analysisManager->CreateNtupleDColumn("VecX",eventAction->Get_pos_x_vector());
  analysisManager->CreateNtupleDColumn("VecY",eventAction->Get_pos_y_vector());
  analysisManager->CreateNtupleDColumn("VecZ",eventAction->Get_pos_z_vector());
  analysisManager->CreateNtupleDColumn("VecLStep",eventAction->GetLstepVector());
  analysisManager->CreateNtupleDColumn("VecDeltaT",eventAction->GetDeltaTVector());
  analysisManager->CreateNtupleDColumn("VecEndKE",eventAction->GetEndKEVector());
  analysisManager->CreateNtupleDColumn("VecEdep",eventAction->GetEdepVector());
  analysisManager->CreateNtupleIColumn("VecParticleA",eventAction->GetAVector());
  analysisManager->CreateNtupleIColumn("VecParticleZ",eventAction->GetZVector());
  analysisManager->CreateNtupleIColumn("DetectorSlice",eventAction->GetSliceVector());
#endif
  analysisManager->FinishNtuple(2);

  analysisManager->CreateNtuple("DummyDetector", "DummyDetector");
  analysisManager->CreateNtupleDColumn("PosX");
  analysisManager->CreateNtupleDColumn("PosY");
  analysisManager->CreateNtupleDColumn("PosZ");

  analysisManager->FinishNtuple(3);


  // analysisManager->CreateNtuple("ReacsParticles", "ReacsParticles");
  // analysisManager->CreateNtupleDColumn("VecParticleA",eventAction->GetAVector());
  // analysisManager->CreateNtupleDColumn("VecParticleZ",eventAction->GetZVector());

  // analysisManager->FinishNtuple(2);

  analysisManager->SetNtupleMerging(true);
  


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();
  // printf("Gozer RunAction it is %u\n",GetSaveVectors());

  // delete fMessenger;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{ 
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  // printf("BeingOfRunAction it is %u\n",GetSaveVectors());
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


void RunAction::BuildAnalysis()
{
  auto analysisManager = G4AnalysisManager::Instance();

  analysisManager->SetVerboseLevel(1);
  
  // analysisManager->CreateH2("Edep_Lstep","Energy Deposited vs. Particle distance",150,0.,15*MeV,100,0,250*mm); //0
  // analysisManager->CreateH2("Edep_DeltaT","TOF vs. E",200,0,2000*nanosecond,200,0.,20*MeV); //1
  // analysisManager->CreateH2("GunEnergy_Edep","Neutron Energy vs. Energy Deposited",200,0.,20000*keV,200,0.,20000*keV); //2
  // analysisManager->CreateH2("GunEnergy_Lstep","Neutron Energy vs. Particle distance",200,0.,20*MeV,250,0,250*mm); //3
  // analysisManager->CreateH2("A_Z","A vs Z",200,0,200,400,0,400);//4
  // analysisManager->CreateH2("GunEnergy_FKE_C6LYC","Neutron Energy vs. Final Kinetic Energy (C6LYC)",2000,0.,20*MeV,2000,0.,20*MeV);//5
  // analysisManager->CreateH2("GunEnergy_FKE_C7LYC","Neutron Energy vs. Final Kinetic Energy (C7LYC)",2000,0.,20*MeV,2000,0.,20*MeV);//6
  // analysisManager->CreateH2("Step_FKE","Steps vs. Final Kinetic Energy",100,0,1000,200,0.,20*MeV);//7
  // analysisManager->CreateH2("GunEnergy_DeltaT","TOF vs. E",1000,0,2000*nanosecond,2000,0.,20*MeV);//8
  // analysisManager->CreateH2("InDetDeltaT_GunEnergy","TOF vs. E",100,0.,20*MeV,100,0,100*nanosecond);//9
  // analysisManager->CreateH2("InDetDeltaD_GunEnergy","TOF vs. E",100,0.,20*MeV,100,0,100*mm);//10


  // analysisManager->CreateH2("Edep_Lstep","Energy Deposited vs. Particle distance",150,0.,15*MeV,100,0,250*mm); //0
  // analysisManager->CreateH2("Edep_DeltaT","TOF vs. E",200,0,2000*nanosecond,200,0.,20*MeV); //1
  // analysisManager->CreateH2("GunEnergy_Edep","Neutron Energy vs. Energy Deposited",200,0.,20000*keV,200,0.,20000*keV); //2
  // analysisManager->CreateH2("GunEnergy_Lstep","Neutron Energy vs. Particle distance",200,0.,20*MeV,250,0,250*mm); //3
  // analysisManager->CreateH2("A_Z","A vs Z",200,0,200,400,0,400);//4
  analysisManager->CreateH2("GunEnergy_FKE_C6LYC","Neutron Energy vs. Final Kinetic Energy (C6LYC)",2000,0.,20*MeV,2000,0.,20*MeV);//0
  analysisManager->CreateH2("GunEnergy_FKE_C7LYC","Neutron Energy vs. Final Kinetic Energy (C7LYC)",2000,0.,20*MeV,2000,0.,20*MeV);//1
  analysisManager->CreateH2("C6LYC_Final_Position","XY Neutron Final Position",400,-200*mm,200*mm,400,-200*mm,200*mm);//2
  analysisManager->CreateH2("C7LYC_Final_Position","XY Neutron Final Position",400,-200*mm,200*mm,400,-200*mm,200*mm);//3
  analysisManager->CreateH2("C7LYC_Energy_Angle","Energy vs. Angle",100,1*MeV,5*MeV,100,0,M_PI);//4
  analysisManager->CreateH2("C6LYC_Slice_Distance","DetectorSlice vs. Distance Traveled",100,0,100,1000,0,100*cm);//5
  analysisManager->CreateH2("C7LYC_Slice_Distance","DetectorSlice vs. Distance Traveled",100,0,100,1000,0,100*cm);//6


  // analysisManager->CreateH2("Step_FKE","Steps vs. Final Kinetic Energy",100,0,1000,200,0.,20*MeV);//7
  // analysisManager->CreateH2("GunEnergy_DeltaT","TOF vs. E",1000,0,2000*nanosecond,2000,0.,20*MeV);//8
  // analysisManager->CreateH2("InDetDeltaT_GunEnergy","TOF vs. E",100,0.,20*MeV,100,0,100*nanosecond);//9
  // analysisManager->CreateH2("InDetDeltaD_GunEnergy","TOF vs. E",100,0.,20*MeV,100,0,100*mm);//10

  // analysisManager->CreateH3("Reaction_X_Y_Z","Reaction_Position",75,-3.75*cm,3.75*cm,75,-3.75*cm,3.75*cm,75,100*cm,101*cm);

  analysisManager->CreateH3("Start_X_Y_Z","Start_Position",75,-3.75*cm,3.75*cm,75,-3.75*cm,3.75*cm,75,-22.15*cm,-19.12*cm);





  analysisManager->CreateH1("Lstep","Particle distance",1000,0,250*mm);
  analysisManager->CreateH1("Edep","Energy Deposited",1500,0.,15*MeV);
  // analysisManager->CreateH1("Deltae","Particle Energy Change",10000,0.,10.*MeV);
  analysisManager->CreateH1("Deltat","Particle Time Change",10000,0., 10000*nanosecond);
  analysisManager->CreateH1("GunEnergy","Neutron Energy",1500,0.,15*MeV);

  analysisManager->CreateH1("C6LYC_InDetDeltaT","Particle Time Change",10000,0., 100*nanosecond);
  analysisManager->CreateH1("C6LYC_InDetDeltaD","Particle Distance",10000,0., 100*mm);
  analysisManager->CreateH1("C7LYC_InDetDeltaT","Particle Time Change",10000,0., 100*nanosecond);
  analysisManager->CreateH1("C7LYC_InDetDeltaD","Particle Distance",10000,0., 100*mm);
  analysisManager->CreateH1("C6LYC_FKE_NEQ_GE","Neutron Energy",2000,0., 20*MeV);
  analysisManager->CreateH1("C7LYC_FKE_NEQ_GE","Neutron Energy",2000,0., 20*MeV);
  analysisManager->CreateH1("C7LYC_TOF","Time of flight",1000,0., 1000*ns);


  analysisManager->CreateNtuple("Data", "Data");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("Lstep");
  analysisManager->CreateNtupleDColumn("Deltat");
  analysisManager->CreateNtupleDColumn("GunEnergy");
  analysisManager->CreateNtupleIColumn("Z");
  analysisManager->CreateNtupleIColumn("A");
  analysisManager->CreateNtupleDColumn("fKE");
  analysisManager->CreateNtupleIColumn("Steps");
  analysisManager->CreateNtupleDColumn("PosX");
  analysisManager->CreateNtupleDColumn("PosY");
  analysisManager->CreateNtupleDColumn("PosZ");
  analysisManager->CreateNtupleIColumn("Detector");
  analysisManager->CreateNtupleDColumn("PreDetKE");
  // analysisManager->CreateNtupleDColumn("VecX",eventAction->Get_pos_x_vector());
  // analysisManager->CreateNtupleDColumn("VecY",eventAction->Get_pos_y_vector());
  // analysisManager->CreateNtupleDColumn("VecZ",eventAction->Get_pos_z_vector());
  // analysisManager->CreateNtupleDColumn("VecLStep",eventAction->GetLstepVector());
  // analysisManager->CreateNtupleDColumn("VecDeltaT",eventAction->GetDeltaTVector());
  // analysisManager->CreateNtupleDColumn("VecEndKE",eventAction->GetEndKEVector());
  // analysisManager->CreateNtupleDColumn("VecEdep",eventAction->GetEdepVector());
  // analysisManager->CreateNtupleIColumn("VecParticleA",eventAction->GetAVector());
  // analysisManager->CreateNtupleIColumn("VecParticleZ",eventAction->GetZVector());

  analysisManager->FinishNtuple(0);

  analysisManager->CreateNtuple("Reacs_C6", "Reacs_C6");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("Lstep");
  analysisManager->CreateNtupleDColumn("Deltat");
  analysisManager->CreateNtupleDColumn("GunEnergy");
  analysisManager->CreateNtupleIColumn("Z");
  analysisManager->CreateNtupleIColumn("A");
  analysisManager->CreateNtupleDColumn("fKE");
  analysisManager->CreateNtupleIColumn("Steps");
  analysisManager->CreateNtupleDColumn("PosX");
  analysisManager->CreateNtupleDColumn("PosY");
  analysisManager->CreateNtupleDColumn("PosZ");
  analysisManager->CreateNtupleDColumn("PreDetKE");
  analysisManager->CreateNtupleIColumn("Detector");
  analysisManager->CreateNtupleDColumn("GlobalTime");
  if(GetSaveVectors())
  {
    analysisManager->CreateNtupleDColumn("VecX",fEventAction->Get_pos_x_vector());
    analysisManager->CreateNtupleDColumn("VecY",fEventAction->Get_pos_y_vector());
    analysisManager->CreateNtupleDColumn("VecZ",fEventAction->Get_pos_z_vector());
    analysisManager->CreateNtupleDColumn("VecLStep",fEventAction->GetLstepVector());
    analysisManager->CreateNtupleDColumn("VecDeltaT",fEventAction->GetDeltaTVector());
    analysisManager->CreateNtupleDColumn("VecEndKE",fEventAction->GetEndKEVector());
    analysisManager->CreateNtupleDColumn("VecEdep",fEventAction->GetEdepVector());
    analysisManager->CreateNtupleIColumn("VecParticleA",fEventAction->GetAVector());
    analysisManager->CreateNtupleIColumn("VecParticleZ",fEventAction->GetZVector());
    analysisManager->CreateNtupleIColumn("DetectorSlice",fEventAction->GetSliceVector());
  }


  analysisManager->FinishNtuple(1);

  analysisManager->CreateNtuple("Reacs_C7", "Reacs_C7");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("Lstep");
  analysisManager->CreateNtupleDColumn("Deltat");
  analysisManager->CreateNtupleDColumn("GunEnergy");
  analysisManager->CreateNtupleIColumn("Z");
  analysisManager->CreateNtupleIColumn("A");
  analysisManager->CreateNtupleDColumn("fKE");
  analysisManager->CreateNtupleIColumn("Steps");
  analysisManager->CreateNtupleDColumn("PosX");
  analysisManager->CreateNtupleDColumn("PosY");
  analysisManager->CreateNtupleDColumn("PosZ");
  analysisManager->CreateNtupleDColumn("PreDetKE");
  analysisManager->CreateNtupleIColumn("Detector");
  analysisManager->CreateNtupleDColumn("GlobalTime");
  if(GetSaveVectors())
  {
    // printf("Saving vectors ho ho ho\n");
    analysisManager->CreateNtupleDColumn("VecX",fEventAction->Get_pos_x_vector());
    analysisManager->CreateNtupleDColumn("VecY",fEventAction->Get_pos_y_vector());
    analysisManager->CreateNtupleDColumn("VecZ",fEventAction->Get_pos_z_vector());
    analysisManager->CreateNtupleDColumn("VecLStep",fEventAction->GetLstepVector());
    analysisManager->CreateNtupleDColumn("VecDeltaT",fEventAction->GetDeltaTVector());
    analysisManager->CreateNtupleDColumn("VecEndKE",fEventAction->GetEndKEVector());
    analysisManager->CreateNtupleDColumn("VecEdep",fEventAction->GetEdepVector());
    analysisManager->CreateNtupleIColumn("VecParticleA",fEventAction->GetAVector());
    analysisManager->CreateNtupleIColumn("VecParticleZ",fEventAction->GetZVector());
    analysisManager->CreateNtupleIColumn("DetectorSlice",fEventAction->GetSliceVector());
  }


  analysisManager->FinishNtuple(2);

  analysisManager->CreateNtuple("DummyDetector", "DummyDetector");
  analysisManager->CreateNtupleDColumn("PosX");
  analysisManager->CreateNtupleDColumn("PosY");
  analysisManager->CreateNtupleDColumn("PosZ");

  analysisManager->FinishNtuple(3);


  // analysisManager->CreateNtuple("ReacsParticles", "ReacsParticles");
  // analysisManager->CreateNtupleDColumn("VecParticleA",eventAction->GetAVector());
  // analysisManager->CreateNtupleDColumn("VecParticleZ",eventAction->GetZVector());

  // analysisManager->FinishNtuple(2);

  analysisManager->SetNtupleMerging(true);
  
  return;
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

