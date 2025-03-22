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
#include "G4ParticleDefinition.hh"
#include "G4HadronicProcess.hh"
#include "G4SteppingManager.hh"



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

  auto particle = step->GetTrack()->GetParticleDefinition();
  auto track = step->GetTrack();
  G4SteppingManager*  steppingManager = fpSteppingManager;


  // if(abs(particle->GetPDGEncoding()) < 100)
  // {
  //   step->GetTrack()->SetTrackStatus(fStopAndKill);
  // }
    // auto detC7LYCvolume_position = fDetConstruction->GetC7LYCVolume()->GetLogicalVolume()->GetMasterSolid()->GetExtent();
    // detC7LYCvolume_position.GetZmax

  auto zPos = step->GetPreStepPoint()->GetPosition().getZ();
  
  auto C6LYCDist = fDetConstruction->GetC6LYCDistance();
  auto C7LYCDist = fDetConstruction->GetC7LYCDistance();

  auto UseC6LYC = fDetConstruction->GetUseC6LYC();
  auto UseC7LYC = fDetConstruction->GetUseC7LYC();

  // if((zPos >=  C6LYCDist+1*m && UseC6LYC) || (zPos >=  C7LYCDist+1*m && UseC7LYC))
  //   {step->GetTrack()->SetTrackStatus(fStopAndKill);}


  // if (!fScoringVolume) { 
  //   const CLYCDetectorConstruction* detectorConstruction
  //     = static_cast<const CLYCDetectorConstruction*>
  //       (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  //   fScoringVolume = detectorConstruction->GetScoringVolume();   
  // }

  // get volume of the current step
  auto ivolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  auto fvolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();

  // Get the detector volume
  auto detC6LYCvolume = fDetConstruction->GetC6LYCVolume();
  auto detC7LYCvolume = fDetConstruction->GetC7LYCVolume();
  auto detdummyvolume = fDetConstruction->GetdummydetectorVolume();
  auto Be9targetvol = fDetConstruction->GetBe9TargetVolume();


  if(ivolume == detdummyvolume && !fEventAction->GetDummy())
  {
    auto dummyposition = step->GetPreStepPoint()->GetPosition();
    fEventAction->SetDummy();
    fEventAction->SetDummyPosX(dummyposition.getX());
    fEventAction->SetDummyPosY(dummyposition.getY());
    fEventAction->SetDummyPosZ(dummyposition.getZ());
  }

  // G4double estart = step->GetPreStepPoint()->GetTr;
  // G4double estop =  step->GetPostStepPoint()->GetTotalEnergy();


  // G4cout << estart << ", " << estop << G4endl;

  // energy deposit
  auto ske = step->GetPostStepPoint()->GetKineticEnergy();
  auto fke = step->GetPreStepPoint()->GetKineticEnergy();
  auto edep = step->GetTotalEnergyDeposit();
  auto non_ion_edep = step->GetNonIonizingEnergyDeposit();
  auto sl = step->GetStepLength();
  auto deltat = step->GetDeltaTime();
  G4int particleA = particle->GetAtomicMass();
  G4int particleZ = particle->GetAtomicNumber();
  // bool isProton = particleA == 1 && particleZ == 1;
  // bool isAlpha = particleA == 4 && particleZ == 2;
  // bool isTriton = particleA == 3 && particleZ == 1;
  // bool isSulfur35 = particleA == 35 && particleZ == 16;
  // bool isPhosphorus32 = particleA == 32 && particleZ == 15;
  // bool isElectronGamma = particleA == 0 && particleZ == 0;
  // bool isChlorine35 = particleA == 35 && particleZ == 17;
  // bool isLithium6 = particleA == 6 && particleZ == 3;
  // bool isChlorine36 = particleA == 36 && particleZ == 17;
  bool isProton{false}, isTriton{false}, isAlpha{false}, isPhosphorus32{false}, isSulfur35{false}, isDeuteron{false};


  // if(isProton)
  // {
  //   fEventAction->sethasH1();
  // }
  // if(isAlpha)
  // {
  //   fEventAction->sethasHe4();
  // }
  // if(isPhosphorus32)
  // {
  //   fEventAction->sethasP32();
  // }
  // if(isSulfur35)
  // {
  //   fEventAction->sethasS35();
  // }

  

  // Not for doing backgrounds!
  // if you only track reactions, then you aren't tracking all of the particles that interacted with Cl35 or Li6
  // additionally, saving everything increases file sizes. But you also dont care about neutrons that don't interact with
  // Cl35 or Li6
  // if(isSulfur35 or isPhosphorus32)
  // {
  //   fEventAction->sethasCl35();
  // }

  // if(isChlorine35)
  // {
  //   fEventAction->sethasCl35();
  //   // if(particle->GetPDGEncoding() == 2112)
  //   // {fEventAction->clearsethasLi6();}

  // }
  // else if(isLithium6)
  // {
  //   fEventAction->sethasLi6();
  //   // if(particle->GetPDGEncoding() == 2112)
  //   // {fEventAction->clearsethasCl35();}
  // }

  // if(particleA == 35 and particleZ == 17)
  // {
  //   fEventAction->sethasCl35();
  // }

  // if(particleA == 6 and particleZ == 3)
  // {
  //   fEventAction->sethasLi6();
  // }

  auto position = step->GetPostStepPoint()->GetPosition();
  fEventAction->AddTo_pos_x(position.getX());
  fEventAction->AddTo_pos_y(position.getY());
  fEventAction->AddTo_pos_z(position.getZ());
  fEventAction->AddToLstepVector(sl);
  fEventAction->AddToDeltaTVector(deltat/ns);
  fEventAction->AddToEndKEVector(ske);
  fEventAction->AddToEdepVector(edep);
  fEventAction->AddToAVector(particleA);
  fEventAction->AddToZVector(particleZ);
  fEventAction->AddLstep(sl/mm);
  fEventAction->AdddeltaT(deltat/ns);
  fEventAction->SetGlobalTime(step->GetPostStepPoint()->GetGlobalTime()/ns);
  auto det=0;

  if(particle->GetPDGEncoding() == 2112)
  {
    // Get the final position of the particle when it stops. Used for calculating
    // corrections due to scintillation light travel time in the crystal.
    auto initialposition = step->GetPreStepPoint()->GetPosition();
    auto finalposition = step->GetPostStepPoint()->GetPosition();

    fEventAction->setFKE(ske);
    fEventAction->AddInDetDeltaT(deltat/ns);
    fEventAction->AddInDetDeltaD(sl/mm);
    fEventAction->SetPosX(finalposition.getX());
    fEventAction->SetPosY(finalposition.getY());
    fEventAction->SetPosZ(finalposition.getZ());
    if(fDetConstruction->GetUseC6LYC() && fvolume ==  detC6LYCvolume)
    {
      fEventAction->AddToSliceVector(detC6LYCvolume->GetCopyNo());
      fEventAction->SetSlice(detC6LYCvolume->GetCopyNo());
    }
    else if(fDetConstruction->GetUseC7LYC() && fvolume == detC7LYCvolume)
    {
      fEventAction->AddToSliceVector(detC7LYCvolume->GetCopyNo());
      fEventAction->SetSlice(detC7LYCvolume->GetCopyNo());
    }
    else if(fDetConstruction->GetUseBe9target() &&  ivolume == Be9targetvol && track->GetTrackStatus() != fAlive)
    {
      fSecondary = steppingManager->GetfSecondary();
      if(fSecondary->size() > 0)
      {
        G4int secondaryPDG;
        G4double secondaryEn;
        G4bool G4bool{false};
        for(auto i = 0; i < fSecondary->size(); ++i)
        {
          // secondaryPDG = step->GetSecondary()->at(i)->GetDynamicParticle()->GetPDGcode();
          // secondaryPDG_En = step->GetSecondary()->at(i)->GetDynamicParticle()->GetKineticEnergy();
          secondaryPDG = fSecondary->at(i)->GetDynamicParticle()->GetPDGcode();
          secondaryEn = fSecondary->at(i)->GetKineticEnergy();
          if(secondaryPDG == 2112)// && (fvolume == detC6LYCvolume || fvolume == detC7LYCvolume))          
          {
            fEventAction->AddTo_SecondaryNeutrons(secondaryEn);
          }
        }
      }
    }
  }
  else
  {
    fEventAction->AddToSliceVector(-1);
    fEventAction->SetSlice(-1);

  }
  // step->GetNonIonizingEnergyDeposit();

  // G4cout << deltat/nanosecond << G4endl;

  // if(particle->GetParticleType() == )



  // if(particle->GetParticleType()  == "lepton")
  // {
    // G4cout << particle->GetParticleType() << ", " << particle->GetParticleSubType() << ", " << particle->GetPDGEncoding() << G4endl;
  // }
  
  // if(particle->GetPDGEncoding() == 2112) // particle code for a neutron
  // {
  //   if(ivolume ==  detC6LYCvolume && fvolume == detC6LYCvolume)
  //   {
  //     fEventAction->SetSlice(detC6LYCvolume->GetCopyNo());
  //   }
  //   else if(ivolume == detC7LYCvolume && fvolume == detC7LYCvolume)
  //   {
  //     fEventAction->SetSlice(detC7LYCvolume->GetCopyNo());
  //   }
  //   else
  //   {
  //     fEventAction->SetSlice(-1);
  //   }
  // }
  // else
  // {
  //   fEventAction->AddToSliceVector(-1);
  //   fEventAction->SetSlice(-1);

  // }

  if( (fvolume == detC6LYCvolume && ivolume != detC6LYCvolume) || (fvolume == detC7LYCvolume && ivolume != detC7LYCvolume))
  {
    if(particle->GetPDGEncoding() == 2112 && fEventAction->GetDetector() < 1)
    {
      // auto KE = step->GetTrack()->GetDynamicParticle()->GetKineticEnergy();
      // if(KE != fEventAction->retGunEnergy())
      // {printf("KE = %f, fke = %f, GE = %f, PDG_Encoding = %i\n",KE,fke,fEventAction->retGunEnergy(),particle->GetPDGEncoding());}
      fEventAction->SetPreDetectorEnergy(fke);
      // fEventAction->SetPreDetectorEnergy(fke);
      G4double GE{fEventAction->retGunEnergy()}, PDE{fEventAction->retPreDetEnergy()};
      if(PDE >= (GE*0.99))
      {
        fEventAction->setGoodPreDetEn();
      }
    }
  }

// D2Gas->GetTotNbOfAtomsPerVolume()

  if(ivolume == detC6LYCvolume )
  {
    fEventAction->setC6LYC();
    if (fvolume == detC6LYCvolume && ivolume == detC6LYCvolume)
    {
      // fEventAction->setC6LYC();
      fEventAction->SetDetector(6);
      det=6;
      fEventAction->AddEdep(edep);
      if(step->GetSecondary()->size() > 0)
      {
        G4int secondaryZ,secondaryA;
        for(auto i = 0; i < step->GetSecondary()->size(); ++i)
        {
          secondaryA = step->GetSecondary()->at(i)->GetDynamicParticle()->GetDefinition()->GetAtomicMass();
          secondaryZ = step->GetSecondary()->at(i)->GetDynamicParticle()->GetDefinition()->GetAtomicNumber();
          isSulfur35 = (secondaryZ == 16 && secondaryA == 35) ? true : isSulfur35;
          isPhosphorus32 = (secondaryZ == 15 && secondaryA == 32) ? true : isPhosphorus32;
          isProton = (secondaryZ == 1 && secondaryA == 1) ? true : isProton;
          isAlpha = (secondaryZ == 2 && secondaryA == 4) ? true : isAlpha;
          isTriton = (secondaryZ == 1 && secondaryA == 3) ? true : isTriton;
          isTriton = (secondaryZ == 1 && secondaryA == 3) ? true : isTriton;
          isDeuteron = (secondaryZ == 1 && secondaryA == 2) ? true : isDeuteron;
  
        }
      }
      if(isProton && isSulfur35)
      {
        // printf("is Cl35(n,p)!\n");
        fEventAction->sethasH1();
        fEventAction->sethasS35();
        fEventAction->sethasCl35();
      }
      else if(isAlpha && isPhosphorus32)
      {
        // printf("is Cl35(n,alpha)\n");
        fEventAction->sethasHe4();
        fEventAction->sethasP32();
        fEventAction->sethasCl35();
      }
      else if((isAlpha && isTriton) || ((isAlpha && isDeuteron)))
      {
        fEventAction->sethasLi6();
      }
  
      
  
      // Get the final position of the particle when it stops. Used for calculating
      // corrections due to scintillation light travel time in the crystal.
      // auto initialposition = step->GetPreStepPoint()->GetPosition();
      // auto finalposition = step->GetPostStepPoint()->GetPosition();
      fEventAction->SetA(particleA);
      fEventAction->SetZ(particleZ);
    }
  }
  else if(ivolume == detC7LYCvolume)
  {
    fEventAction->setC7LYC();
    if (fvolume == detC7LYCvolume && ivolume == detC7LYCvolume)
    {
      // bool isProton{false}, isTriton{false}, isAlpha{false}, isPhosphorus32{false}, isSulfur35{false};
      // fEventAction->setC7LYC();
      fEventAction->SetDetector(7);
      det=7;
      fEventAction->AddEdep(edep);
      if(step->GetSecondary()->size() > 0)
      {
        G4int secondaryZ,secondaryA;
        for(auto i = 0; i < step->GetSecondary()->size(); ++i)
        {
          secondaryA = step->GetSecondary()->at(i)->GetDynamicParticle()->GetDefinition()->GetAtomicMass();
          secondaryZ = step->GetSecondary()->at(i)->GetDynamicParticle()->GetDefinition()->GetAtomicNumber();
          isSulfur35 = (secondaryZ == 16 && secondaryA == 35) ? true : isSulfur35;
          isPhosphorus32 = (secondaryZ == 15 && secondaryA == 32) ? true : isPhosphorus32;
          isProton = (secondaryZ == 1 && secondaryA == 1) ? true : isProton;
          isAlpha = (secondaryZ == 2 && secondaryA == 4) ? true : isAlpha;
          isTriton = (secondaryZ == 1 && secondaryA == 3) ? true : isTriton;
          isTriton = (secondaryZ == 1 && secondaryA == 3) ? true : isTriton;
          isDeuteron = (secondaryZ == 1 && secondaryA == 2) ? true : isDeuteron;
  
        }
      }
      if(isProton && isSulfur35)
      {
        // printf("is Cl35(n,p)!\n");
        fEventAction->sethasH1();
        fEventAction->sethasS35();
        fEventAction->sethasCl35();
        fEventAction->incNPcounter();
      }
      else if(isAlpha && isPhosphorus32)
      {
        // printf("is Cl35(n,alpha)\n");
        fEventAction->sethasHe4();
        fEventAction->sethasP32();
        fEventAction->sethasCl35();
        fEventAction->incNALPHAcounter();
      }
      else if((isAlpha && isTriton) || ((isAlpha && isDeuteron)))
      {
        fEventAction->sethasLi6();
      }
  
      // Get the final position of the particle when it stops. Used for calculating
      // corrections due to scintillation light travel time in the crystal.
      // auto initialposition = step->GetPreStepPoint()->GetPosition();
      // auto finalposition = step->GetPostStepPoint()->GetPosition();
      fEventAction->SetA(particleA);
      fEventAction->SetZ(particleZ);
    }
  }

  // if(fvolume == detC6LYCvolume)
  // {

  // }
  // else if(fvolume == detC7LYCvolume)
  // {

  // }
  // if (fvolume == detC6LYCvolume && ivolume == detC6LYCvolume)
  // {
  //   // fEventAction->setC6LYC();
  //   fEventAction->SetDetector(6);
  //   det=6;
  //   fEventAction->AddEdep(edep);
  //   if(step->GetSecondary()->size() > 0)
  //   {
  //     G4int secondaryZ,secondaryA;
  //     for(auto i = 0; i < step->GetSecondary()->size(); ++i)
  //     {
  //       secondaryA = step->GetSecondary()->at(i)->GetDynamicParticle()->GetDefinition()->GetAtomicMass();
  //       secondaryZ = step->GetSecondary()->at(i)->GetDynamicParticle()->GetDefinition()->GetAtomicNumber();
  //       isSulfur35 = (secondaryZ == 16 && secondaryA == 35) ? true : isSulfur35;
  //       isPhosphorus32 = (secondaryZ == 15 && secondaryA == 32) ? true : isPhosphorus32;
  //       isProton = (secondaryZ == 1 && secondaryA == 1) ? true : isProton;
  //       isAlpha = (secondaryZ == 2 && secondaryA == 4) ? true : isAlpha;
  //       isTriton = (secondaryZ == 1 && secondaryA == 3) ? true : isTriton;
  //       isTriton = (secondaryZ == 1 && secondaryA == 3) ? true : isTriton;
  //       isDeuteron = (secondaryZ == 1 && secondaryA == 2) ? true : isDeuteron;

  //     }
  //   }
  //   if(isProton && isSulfur35)
  //   {
  //     // printf("is Cl35(n,p)!\n");
  //     fEventAction->sethasH1();
  //     fEventAction->sethasS35();
  //     fEventAction->sethasCl35();
  //   }
  //   else if(isAlpha && isPhosphorus32)
  //   {
  //     // printf("is Cl35(n,alpha)\n");
  //     fEventAction->sethasHe4();
  //     fEventAction->sethasP32();
  //     fEventAction->sethasCl35();
  //   }
  //   else if((isAlpha && isTriton) || ((isAlpha && isDeuteron)))
  //   {
  //     fEventAction->sethasLi6();
  //   }

    

  //   // Get the final position of the particle when it stops. Used for calculating
  //   // corrections due to scintillation light travel time in the crystal.
  //   // auto initialposition = step->GetPreStepPoint()->GetPosition();
  //   // auto finalposition = step->GetPostStepPoint()->GetPosition();
  //   fEventAction->SetA(particleA);
  //   fEventAction->SetZ(particleZ);
  // }
  // (((fvolume == detC6LYCvolume && ivolume == detC6LYCvolume) || (fvolume == detC7LYCvolume && ivolume == detC7LYCvolume)))
  // check if we are in scoring volume
  // else if (fvolume == detC7LYCvolume && ivolume == detC7LYCvolume)
  // {
  //   // bool isProton{false}, isTriton{false}, isAlpha{false}, isPhosphorus32{false}, isSulfur35{false};
  //   // fEventAction->setC7LYC();
  //   fEventAction->SetDetector(7);
  //   det=7;
  //   fEventAction->AddEdep(edep);
  //   if(step->GetSecondary()->size() > 0)
  //   {
  //     G4int secondaryZ,secondaryA;
  //     for(auto i = 0; i < step->GetSecondary()->size(); ++i)
  //     {
  //       secondaryA = step->GetSecondary()->at(i)->GetDynamicParticle()->GetDefinition()->GetAtomicMass();
  //       secondaryZ = step->GetSecondary()->at(i)->GetDynamicParticle()->GetDefinition()->GetAtomicNumber();
  //       isSulfur35 = (secondaryZ == 16 && secondaryA == 35) ? true : isSulfur35;
  //       isPhosphorus32 = (secondaryZ == 15 && secondaryA == 32) ? true : isPhosphorus32;
  //       isProton = (secondaryZ == 1 && secondaryA == 1) ? true : isProton;
  //       isAlpha = (secondaryZ == 2 && secondaryA == 4) ? true : isAlpha;
  //       isTriton = (secondaryZ == 1 && secondaryA == 3) ? true : isTriton;
  //       isTriton = (secondaryZ == 1 && secondaryA == 3) ? true : isTriton;
  //       isDeuteron = (secondaryZ == 1 && secondaryA == 2) ? true : isDeuteron;

  //     }
  //   }
  //   if(isProton && isSulfur35)
  //   {
  //     // printf("is Cl35(n,p)!\n");
  //     fEventAction->sethasH1();
  //     fEventAction->sethasS35();
  //     fEventAction->sethasCl35();
  //   }
  //   else if(isAlpha && isPhosphorus32)
  //   {
  //     // printf("is Cl35(n,alpha)\n");
  //     fEventAction->sethasHe4();
  //     fEventAction->sethasP32();
  //     fEventAction->sethasCl35();
  //   }
  //   else if((isAlpha && isTriton) || ((isAlpha && isDeuteron)))
  //   {
  //     fEventAction->sethasLi6();
  //   }

  //   // Get the final position of the particle when it stops. Used for calculating
  //   // corrections due to scintillation light travel time in the crystal.
  //   // auto initialposition = step->GetPreStepPoint()->GetPosition();
  //   // auto finalposition = step->GetPostStepPoint()->GetPosition();
  //   fEventAction->SetA(particleA);
  //   fEventAction->SetZ(particleZ);
  // }

    // if(particle->GetPDGEncoding() == 2112)// and not(fEventAction->issethasCl35() or fEventAction->issethasLi6()))
    // {
    //   fEventAction->setFKE(KE);
    //   fEventAction->AddInDetDeltaT(deltat/ns);
    //   fEventAction->AddInDetDeltaD(sl/mm);
    //   // printf("step length = %f\n",(double)sl/mm);
    //   fEventAction->SetPosX(finalposition.getX());
    //   fEventAction->SetPosY(finalposition.getY());
    //   fEventAction->SetPosZ(finalposition.getZ());
      // fEventAction->SetParticleX(finalposition.getX());
      // fEventAction->SetParticleY(finalposition.getY());
      // fEventAction->SetParticleZ(finalposition.getZ());

      // if(fDetConstruction->GetUseC6LYC() || fDetConstruction->GetUseC7LYC())
      // {
      //   if(ivolume ==  detC6LYCvolume && fvolume ==  detC6LYCvolume)
      //   {
      //     fEventAction->SetSlice(detC6LYCvolume->GetCopyNo());
      //   }
      //   else if(ivolume == detC7LYCvolume && fvolume == detC7LYCvolume)
      //   {
      //     fEventAction->SetSlice(detC7LYCvolume->GetCopyNo());
      //   }
      //   else
      //   {
      //     fEventAction->SetSlice(-1);
      //   }
      // }
      // else if(fDetConstruction->GetUseC7LYC())
      // {
      //   if(fvolume == detC7LYCvolume)
      //   {
      //     fEventAction->SetSlice(detC7LYCvolume->GetCopyNo());
      //   } 
      //   else
      //   {
      //   }
      // }
    

    // }



    // Adds a step to the number of steps the particle took in the crystal
    fEventAction->AddStep();

    // Used for calculating total time by hand. Not in use right now
    // fEventAction->AddToLstepVector(sl/mm); 
    // fEventAction->AddToStartKEVector(fke/MeV);
    // fEventAction->AddToEndKEVector(ske/MeV);


    // G4cout << endPoint->GetKineticEnergy()/MeV - fEventAction->retGE()/MeV << G4endl;

    // if(step->GetSecondary()->size() > 0 && edep > 0)
    // {
    //   for(auto i = 0; i < step->GetSecondary()->size(); i++)
    //   {
    //     if(particle->GetAtomicMass() == 1 && particle->GetAtomicNumber() == 1)
    //     {}
    //     // {G4cout << step->GetSecondary()->size() << ": " << "i = " << i << " - "  << step->GetSecondary()->at(i)->GetDefinition()->GetParticleName() << G4endl;}
    //   }
    // }
    // if(has35)
    // {
    //   fEventAction->AddEdep(edep);
    //   fEventAction->AddLstep(sl/mm);
    //   // fEventAction->AdddeltaE(deltae);
    //   fEventAction->AdddeltaT(deltat/ns);
    //   G4cout << step->GetSecondary()->size() << G4endl;
    //   for(auto i = 0; i < step->GetSecondary()->size(); i++)
    //   {
    //     step->G
    //     G4cout << step->GetSecondary()->at(i)->GetStepLength() << ", " << ;
    //     G4cout << step->GetSecondary()->at(i)->GetParticleDefinition()->GetAtomicNumber() << ", " << step->GetSecondary()->at(i)->GetParticleDefinition()->GetBaryonNumber() << G4endl;
    //   }

    // }
  // }



  // if(step->GetPreStepPoint()->GetPosition().getZ() >= 1010*mm)
  //   {step->GetTrack()->SetTrackStatus(fStopAndKill);}

  // if ((fvolume != detC7LYCvolume && ivolume == detC7LYCvolume))
  // {
  //   step->GetTrack()->SetTrackStatus(fStopAndKill);
  // }


  // G4cout << edep << G4endl;
  // collect energy deposited in this step
  // G4double edepStep = step->GetTotalEnergyDeposit();
  // G4double lStep = step->GetStepLength();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

