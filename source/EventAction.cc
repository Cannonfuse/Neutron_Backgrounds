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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
// #include "G4AnalysisManager.hh"
#include "g4root.hh"


#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
// #include "PhysicalConstants.hh"
#include <random>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double energyRes(G4double edep)
{
  std::random_device DEVICE;
  std::mt19937_64 GENERATOR(DEVICE());
  // std::default_random_engine GENERATOR;
  double EDEP = static_cast<double> (edep);
  EDEP *= 1000.;
  // G4cout << edep * MeV << G4endl;
  // G4cout << EDEP << G4endl;
   
  std::normal_distribution<double> resA(-1.15446028e-04, 1.85495156e-05);
  std::normal_distribution<double> resB(2.14550165e-01, 1.69963576e-02);
  double stddev = (  resA(GENERATOR) * EDEP + resB(GENERATOR)) / (2 * std::sqrt(2. * std::log(2)));
  stddev *= EDEP;
  // G4cout << stddev << G4endl;

  // std::normal_distribution<double> edepRes(0, stddev);
  // double edepAdd = edepRes(GENERATOR);
  // G4double modedep = edep + edepAdd * keV;

  std::normal_distribution<double> edepRes(static_cast<double> (edep) * 1000., stddev);
  return edepRes(GENERATOR) * keV;
}

G4double detEnergyResponse(G4double edep)
{
  std::random_device DEVICE;
  std::mt19937_64 GENERATOR(DEVICE());
  // std::default_random_engine GENERATOR;
  double EDEP = static_cast<double> (edep);
  EDEP *= 1000.;
  std::normal_distribution<double> resA(0, 1);
  if(EDEP <= 300.)
  {
    if(resA(GENERATOR) >= 0)
    {
      return EDEP * keV;
    }
    return 0 * keV;
  }
  return EDEP * keV;
}

G4double calcTime(G4double start_KE, G4double end_KE, G4double DIST)
{
  // G4cout << DIST << G4endl;
  // lambda E,M: constants.c * np.sqrt(1 - np.power(1+E/M,-2))
  G4double vstart = CLHEP::c_light * sqrt(1 - pow(1 + start_KE/(939.550*MeV),-2));
  // G4double vend = CLHEP::c_light * sqrt(1 - pow(1 + end_KE/(939.550*MeV),-2));
  // G4double vavg = (vstart - vend)/2;

  // if(DIST > 0)
  // {
  //   return ((DIST/mm)/v);
  // }

  // G4cout << (DIST/mm)/vstart << ", " << (DIST/mm)/vavg  << G4endl;

  if(DIST > 0)
  {
    return ((DIST/mm)/vstart);
  }
    return 0;
}


EventAction::EventAction()
: G4UserEventAction(),
  fEdep(0.),
  fLstep(0.),
  fGunEnergy(0.)

{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{    
  fEdep = 0.;
  fPreDetectorEnergy = 0.;
  fLstep = 0.;
  fSteps = 0;
  fDeltat = 0.;
  inDetDeltaT = 0.;
  inDetDeltaD = 0.;
  A = 0;
  Z = 0;
  fdummyXPosition = 0;
  fdummyYPosition = 0;
  fdummyZPosition = 0;
  fDummy=false;
  hasCl35=false;
  hasLi6=false;
  goodPreDet=false;
  Detector=0;
  LstepVector.clear();
  DeltaTVector.clear();
  endKEVector.clear();
  EdepVector.clear();
  particleAVector.clear();
  particleZVector.clear();
  fGunEnergy = event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();
  pos_x.clear();
  pos_y.clear();
  pos_z.clear();
  auto startposition = event->GetPrimaryVertex()->GetPosition();
  AddTo_pos_x(startposition.getX());
  AddTo_pos_y(startposition.getY());
  AddTo_pos_z(startposition.getZ());
  AddToLstepVector(0);
  AddToEndKEVector(fGunEnergy);
  AddToDeltaTVector(0);
  EdepVector.push_back(0);
  int startA = event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetAtomicMass();
  int startZ = event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetAtomicNumber();
  AddToAVector(startA);
  AddToZVector(startZ);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{   
  // This is the analysis manager
  auto analysisManager = G4AnalysisManager::Instance(); 

  auto  dataNtuple = analysisManager->GetFirstNtupleId();
  auto reacs6Ntuple = analysisManager->GetFirstNtupleId() + 1;
  auto reacs7Ntuple = analysisManager->GetFirstNtupleId() + 2;
  auto dummyNtuple = analysisManager->GetFirstNtupleId() + 3;

  // auto reacsParticlesNtuple = analysisManager->GetFirstNtupleId() + 2;

  
  // G4double new_Deltat{0};
  // G4int ntupleID = analysisManager->GetFirstNtupleId();

  // if(LstepVector.size() > 0)
  // {
  //   for(auto i = 0; i < LstepVector.size(); ++i)
  //   {
  //     new_Deltat += calcTime(startKEVector.at(i)/MeV,endKEVector.at(i)/MeV,LstepVector.at(i)/mm);
  //   }
  //   new_Deltat/=ns;
  // }

  // G4cout << fDeltat << ", " << new_Deltat << G4endl;

  // analysisManager->FillNtupleDColumn(dataNtuple, 0, fEdep);
  // analysisManager->FillNtupleDColumn(dataNtuple, 1, fLstep);
  // analysisManager->FillNtupleDColumn(dataNtuple, 2, fDeltat);
  // analysisManager->FillNtupleDColumn(dataNtuple, 3, fGunEnergy);
  // analysisManager->FillNtupleDColumn(dataNtuple, 4, fPreDetectorEnergy);
  // analysisManager->FillNtupleIColumn(dataNtuple, 5, Detector);

  // analysisManager->FillNtupleDColumn(dataNtuple, 0, fEdep);
  // analysisManager->FillNtupleDColumn(dataNtuple, 1, fLstep);
  // analysisManager->FillNtupleDColumn(dataNtuple, 2, fDeltat);
  // analysisManager->FillNtupleDColumn(dataNtuple, 3, fGunEnergy);
  // analysisManager->FillNtupleIColumn(dataNtuple, 4, Z);
  // analysisManager->FillNtupleIColumn(dataNtuple, 5, A);
  // analysisManager->FillNtupleDColumn(dataNtuple, 6, fKineticEnergy);
  // analysisManager->FillNtupleIColumn(dataNtuple, 7, fSteps);
  // analysisManager->FillNtupleDColumn(dataNtuple, 8, fXPosition);
  // analysisManager->FillNtupleDColumn(dataNtuple, 9, fYPosition);
  // analysisManager->FillNtupleDColumn(dataNtuple, 10, fZPosition);
  // analysisManager->FillNtupleDColumn(dataNtuple, 11, fPreDetectorEnergy);
  // analysisManager->FillNtupleIColumn(dataNtuple, 12, Detector);

  // analysisManager->AddNtupleRow(dataNtuple);

  if(fDummy)
  {  
    analysisManager->FillNtupleDColumn(dummyNtuple, 0, fdummyXPosition);
    analysisManager->FillNtupleDColumn(dummyNtuple, 1, fdummyYPosition);
    analysisManager->FillNtupleDColumn(dummyNtuple, 2, fdummyZPosition);

    analysisManager->AddNtupleRow(dummyNtuple);
  }  

  // G4cout << analysisManager->GetFirstNtupleId() << " , " << analysisManager->GetFirstNtupleId()+1 << G4endl;


  // G4cout << G4BestUnit(fLstep, "Length")<< ", " << G4BestUnit(fEdep, "Energy") << G4endl;

  // G4cout << analysisManager->GetNofNtuples() << G4endl;

  analysisManager->FillH1(0,fLstep);
  analysisManager->FillH1(1,fEdep);
  analysisManager->FillH1(2,fDeltat);
  analysisManager->FillH1(3,fGunEnergy);


  // analysisManager->FillH2(analysisManager->GetFirstH2Id(),fEdep,fLstep);
  // analysisManager->FillH2(analysisManager->GetFirstH2Id()+1,fDeltat,fEdep);
  // analysisManager->FillH2(analysisManager->GetFirstH2Id()+2,fGunEnergy,fEdep);
  // analysisManager->FillH2(analysisManager->GetFirstH2Id()+3,fGunEnergy,fLstep);
  // analysisManager->FillH2(analysisManager->GetFirstH2Id()+4,Z,A);
  // analysisManager->FillH2(analysisManager->GetFirstH2Id()+7,fSteps,fKineticEnergy);


  if(Detector == 6)
  {

    analysisManager->FillNtupleDColumn(reacs6Ntuple, 0, fEdep);
    analysisManager->FillNtupleDColumn(reacs6Ntuple, 1, fLstep);
    analysisManager->FillNtupleDColumn(reacs6Ntuple, 2, fDeltat);
    analysisManager->FillNtupleDColumn(reacs6Ntuple, 3, fGunEnergy);
    analysisManager->FillNtupleIColumn(reacs6Ntuple, 4, Z);
    analysisManager->FillNtupleIColumn(reacs6Ntuple, 5, A);
    analysisManager->FillNtupleDColumn(reacs6Ntuple, 6, fKineticEnergy);
    analysisManager->FillNtupleIColumn(reacs6Ntuple, 7, fSteps);
    analysisManager->FillNtupleDColumn(reacs6Ntuple, 8, fXPosition);
    analysisManager->FillNtupleDColumn(reacs6Ntuple, 9, fYPosition);
    analysisManager->FillNtupleDColumn(reacs6Ntuple, 10, fZPosition);
    analysisManager->FillNtupleDColumn(reacs6Ntuple, 11, fPreDetectorEnergy);
    analysisManager->FillNtupleIColumn(reacs6Ntuple, 12, Detector);

    analysisManager->AddNtupleRow(reacs6Ntuple);

    if(issethasCl35() or issethasLi6())
    {
      analysisManager->FillH1(analysisManager->GetFirstH1Id()+4,inDetDeltaT);
      analysisManager->FillH1(analysisManager->GetFirstH1Id()+5,inDetDeltaD);
      if(issetGoodPreDetEn())
      {
        analysisManager->FillH1(analysisManager->GetFirstH1Id()+8,fPreDetectorEnergy);
      }


      analysisManager->FillH2(analysisManager->GetFirstH2Id(),fGunEnergy,fPreDetectorEnergy);
      analysisManager->FillH2(analysisManager->GetFirstH2Id()+2,fXPosition,fYPosition);

      // analysisManager->FillH2(analysisManager->GetFirstH2Id()+8,fDeltat,fGunEnergy);
      // analysisManager->FillH2(analysisManager->GetFirstH2Id()+9,fGunEnergy,inDetDeltaT);
      // analysisManager->FillH2(analysisManager->GetFirstH2Id()+10,fGunEnergy,inDetDeltaD);
      // analysisManager->FillH3(analysisManager->GetFirstH3Id(),particleX,particleY,particleZ);
    }
  }

  if(Detector == 7)
    {

      analysisManager->FillNtupleDColumn(reacs7Ntuple, 0, fEdep);
      analysisManager->FillNtupleDColumn(reacs7Ntuple, 1, fLstep);
      analysisManager->FillNtupleDColumn(reacs7Ntuple, 2, fDeltat);
      analysisManager->FillNtupleDColumn(reacs7Ntuple, 3, fGunEnergy);
      analysisManager->FillNtupleIColumn(reacs7Ntuple, 4, Z);
      analysisManager->FillNtupleIColumn(reacs7Ntuple, 5, A);
      analysisManager->FillNtupleDColumn(reacs7Ntuple, 6, fKineticEnergy);
      analysisManager->FillNtupleIColumn(reacs7Ntuple, 7, fSteps);
      analysisManager->FillNtupleDColumn(reacs7Ntuple, 8, fXPosition);
      analysisManager->FillNtupleDColumn(reacs7Ntuple, 9, fYPosition);
      analysisManager->FillNtupleDColumn(reacs7Ntuple, 10, fZPosition);
      analysisManager->FillNtupleDColumn(reacs7Ntuple, 11, fPreDetectorEnergy);
      analysisManager->FillNtupleIColumn(reacs7Ntuple, 12, Detector);

      analysisManager->AddNtupleRow(reacs7Ntuple);


      if(issethasCl35())
      {
        analysisManager->FillH1(analysisManager->GetFirstH1Id()+6,inDetDeltaT);
        analysisManager->FillH1(analysisManager->GetFirstH1Id()+7,inDetDeltaD);
        if(issetGoodPreDetEn())
        {
          analysisManager->FillH1(analysisManager->GetFirstH1Id()+9,fPreDetectorEnergy);
        }


        analysisManager->FillH2(analysisManager->GetFirstH2Id()+1,fGunEnergy,fPreDetectorEnergy);
        analysisManager->FillH2(analysisManager->GetFirstH2Id()+3,fXPosition,fYPosition);

        // analysisManager->FillH2(analysisManager->GetFirstH2Id()+8,fDeltat,fGunEnergy);
        // analysisManager->FillH2(analysisManager->GetFirstH2Id()+9,fGunEnergy,inDetDeltaT);
        // analysisManager->FillH2(analysisManager->GetFirstH2Id()+10,fGunEnergy,inDetDeltaD);
        // analysisManager->FillH3(analysisManager->GetFirstH3Id(),particleX,particleY,particleZ);
      }
    }

  // else
  // {
  //     analysisManager->SetNtupleActivation(analysisManager->GetFirstNtupleId()+1,false);
  //     analysisManager->SetNtupleActivation(analysisManager->GetFirstNtupleId(),true);

  //     analysisManager->FillNtupleDColumn(0, fEdep);
  //     analysisManager->FillNtupleDColumn(1, fLstep);
  //     analysisManager->FillNtupleDColumn(2, new_Deltat);
  //     analysisManager->FillNtupleDColumn(3, fGunEnergy);
  //     analysisManager->FillNtupleDColumn(4, Z);
  //     analysisManager->FillNtupleDColumn(5, A);
  //     analysisManager->FillNtupleDColumn(6, fKineticEnergy);
  //     analysisManager->FillNtupleIColumn(7, fSteps);
  // }
  // accumulate statistics in run action

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
