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
  bool isProton = particleA == 1 && particleZ == 1;
  bool isAlpha = particleA == 4 && particleZ == 2;
  bool isTriton = particleA == 3 && particleZ == 1;
  bool isSulfur35 = particleA == 35 && particleZ == 16;
  bool isPhosphorus32 = particleA == 32 && particleZ == 15;
  bool isElectronGamma = particleA == 0 && particleZ == 0;
  bool isChlorine35 = particleA == 35 && particleZ == 17;
  bool isLithium6 = particleA == 6 && particleZ == 3;
  bool isChlorine36 = particleA == 36 && particleZ == 17;


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

  if(isChlorine35)
  {
    fEventAction->sethasCl35();
    // if(particle->GetPDGEncoding() == 2112)
    // {fEventAction->clearsethasLi6();}

  }
  else if(isLithium6)
  {
    fEventAction->sethasLi6();
    // if(particle->GetPDGEncoding() == 2112)
    // {fEventAction->clearsethasCl35();}
  }

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


  if(fDetConstruction->GetUseC6LYC() && fvolume ==  detC6LYCvolume)
  {
    fEventAction->AddToSliceVector(detC6LYCvolume->GetCopyNo());
    if(particle->GetPDGEncoding() == 2112)
    {
    fEventAction->SetSlice(detC6LYCvolume->GetCopyNo());
    }
  }
  else if(fDetConstruction->GetUseC7LYC() && fvolume == detC7LYCvolume)
  {
    fEventAction->AddToSliceVector(detC7LYCvolume->GetCopyNo());
    if(particle->GetPDGEncoding() == 2112)
    {
      fEventAction->SetSlice(detC7LYCvolume->GetCopyNo());
      // printf("Event: %i -> Energy = %f, slice = %i\n",fEventAction->getEventID(),ske,detC7LYCvolume->GetCopyNo());
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

  if((ivolume == detC6LYCvolume || fvolume == detC6LYCvolume))
  {
    fEventAction->setC6LYC();
  }
  if((ivolume == detC7LYCvolume || fvolume == detC7LYCvolume))
  {
    fEventAction->setC7LYC();
  }

  if(fvolume == detC6LYCvolume)
  {
    fEventAction->SetDetector(6);
    det=6;
  }
  else if(fvolume == detC7LYCvolume)
  {
    fEventAction->SetDetector(7);
    det=7;
  }

      
  // check if we are in scoring volume
  if (((fvolume == detC6LYCvolume && ivolume == detC6LYCvolume) || (fvolume == detC7LYCvolume && ivolume == detC7LYCvolume)))
  {
    // bool has35 = false;
    // G4cout << edep/keV << ", " << step->GetSecondaryInCurrentStep()->at(0)->/keV <<  G4endl;
    // if(particle->GetPDGEncoding() != 2112)
    // if(issethasCl35() || issethasLi6())
    fEventAction->AddEdep(edep);


    // if(isChlorine35)
    // {
    //     fEventAction->sethasCl35();
    //     // fEventAction->clearsethasLi6();
  
    // }
    // else if(isLithium6)
    // {
    //   fEventAction->sethasLi6();
    //   // fEventAction->clearsethasCl35();
    // }

    // if(fEventAction->issethasCl35())
    // {
      if(isProton)
      {
        fEventAction->sethasH1();
      }
      if(isAlpha)
      {
        fEventAction->sethasHe4();
      }
      if(isPhosphorus32)
      {
        fEventAction->sethasP32();
      }
      if(isSulfur35)
      {
        fEventAction->sethasS35();
      }
    // }


    //   if(isProton || isAlpha || isPhosphorus32 || isSulfur35)
    //   {
    //     fEventAction->AddEdep(edep);
    //     fEventAction->setCl35Reaction();
    //   }
    //   if(not(isProton || isAlpha || isPhosphorus32 || isSulfur35 || isTriton))
    //   {
    //     fEventAction->AddNonIonEdep(edep);
    //   }
    // }
    // else if(fEventAction->issethasLi6())
    // {
    //   if(isTriton || isAlpha)// || isPhosphorus32 || isSulfur35)
    //   {
    //     fEventAction->AddEdep(edep);
    //   }
    //   if(not(isProton || isAlpha || isPhosphorus32 || isSulfur35 || isTriton))
    //   {
    //     fEventAction->AddNonIonEdep(edep);
      // }
    // if(issethasLi6())
    // {
    //   if(isTriton || isAlpha)
    //   {
    //     fEventAction->AddEdep(edep);
    //     fEventAction->AddNonIonEdep(non_ion_edep);
    //   }
    // }

    // fEventAction->SetGlobalTime(step->GetPostStepPoint()->GetGlobalTime()/ns);
    // fEventAction->AddInDetDeltaT(deltat/ns);

    // const G4StepPoint* endPoint = step->GetPostStepPoint();
    // G4VProcess* process   = 
    //                 const_cast<G4VProcess*>(endPoint->GetProcessDefinedStep());

    // G4String partName = particle->GetParticleName();
    // G4String nuclearChannel = partName;
    // G4HadronicProcess* hproc = dynamic_cast<G4HadronicProcess*>(process);
    // const G4Isotope* target = NULL;
    // if (hproc) target = hproc->GetTargetIsotope();

    // Get the final position of the particle when it stops. Used for calculating
    // corrections due to scintillation light travel time in the crystal.
    auto initialposition = step->GetPreStepPoint()->GetPosition();

    auto finalposition = step->GetPostStepPoint()->GetPosition();






    // Get the final particle produced when it stops in the crystal
    // G4ParticleDefinition* particle = step->GetTrack()->GetDefinition();
    // fEventAction->SetA(particle->GetAtomicMass());
    // fEventAction->SetZ(particle->GetAtomicNumber());
    fEventAction->SetA(particleA);
    fEventAction->SetZ(particleZ);
    // if(fvolume == detC6LYCvolume)
    // {
    //   fEventAction->SetDetector(6);
    // }
    // else if(fvolume == detC7LYCvolume)
    // {
    //   fEventAction->SetDetector(7);
    // }
    // else
    // {
    //   fEventAction->SetDetector(-1);
    // }
    // Set the kinetic energy of the beam particle prior to stopping in the final step
    // if(particle->GetParticleName() == "neutron")
          // auto KE = step->GetTrack()->GetDynamicParticle()->GetKineticEnergy();
          // auto particlemass = step->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetAtomicMass

        // G4cout << "GE:" << fEventAction->retGunEnergy() <<  ", KE: " << KE << ", Edep:" << edep << ", Det: " << det << ", Z_A: " << particleZ << "_" << particleA << G4endl;


    if(particle->GetPDGEncoding() == 2112)// and not(fEventAction->issethasCl35() or fEventAction->issethasLi6()))
    {
      auto KE = step->GetTrack()->GetDynamicParticle()->GetKineticEnergy();
      // if(edep > 0) {G4cout << "KE: " << KE << ", Edep:" << edep << ", Det: " << det << G4endl;}
    // G4cout << "GE:" << fEventAction->retGunEnergy() <<  ", KE: " << KE << ", Edep:" << edep << ", Det: " << det << G4endl;
      // if(KE < fEventAction->retGunEnergy())
      // {
      fEventAction->setFKE(KE);

                // G4cout << "GE:" << fEventAction->retGunEnergy() <<  ", KE: " << KE << ", Edep:" << edep << ", Det: " << det << ", Z_A: " << particleZ << "_" << particleA << G4endl;
      // }
      fEventAction->AddInDetDeltaT(deltat/ns);
      fEventAction->AddInDetDeltaD(sl/mm);
      // printf("step length = %f\n",(double)sl/mm);
      fEventAction->SetPosX(finalposition.getX());
      fEventAction->SetPosY(finalposition.getY());
      fEventAction->SetPosZ(finalposition.getZ());
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
    

    }



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
  }



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

