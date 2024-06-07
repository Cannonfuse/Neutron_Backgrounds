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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
// #include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
#include <random>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double uniRand(double start, double end)
{
  std::random_device DEVICE;
  std::mt19937_64 GENERATOR(DEVICE());
  std::uniform_real_distribution<double> unirandnum(start, end);
  return unirandnum(GENERATOR);
}

G4double retBeamEn(const PrimaryGeneratorAction *PGA)
{
  return PGA->fGPS->GetParticleEnergy();
}

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fGPS(nullptr)//, 
  //fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fGPS  = new G4GeneralParticleSource();

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");

  // fParticleGun->SetParticleDefinition(particle);
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  // fParticleGun->SetParticleEnergy(4*MeV);
}

// PrimaryGeneratorAction::PrimaryGeneratorAction(const G4bool usedists, 
//                        const std::vector<std::vector<G4bool>> energyangledist,
//                        const std::vector<std::vector<G4bool>> energyzdist,
//                        const std::vector<std::vector<double>> energyangzbins)
//                       //  const std::default_random_engine *generator)   
// : G4VUserPrimaryGeneratorAction(),
//   fGPS(nullptr), generator(nullptr)//, 
//   //fEnvelopeBox(0)
// {
//   SetEnergyAngle_dist(energyangledist);
//   SetEnergyZ_dist(energyzdist);
//   SetEnergyAngleZ_bins(energyangzbins);
//   SetUseDists(usedists);

//   // SetGenerator(generator);
//   G4int n_particle = 1;
//   fGPS = new G4GeneralParticleSource();
//   generator = new std::mt19937_64();

//   // default particle kinematic
//   G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
//   G4String particleName;
//   G4ParticleDefinition* particle
//     = particleTable->FindParticle(particleName="gamma");
    
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGPS;
  // delete generator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{


  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  
  // G4double envSizeX = 0;
  // G4double envSizeY = 0;
  // G4double envSizeZ = 0;

  
  // if (!fEnvelopeBox)
  // {
  //   G4LogicalVolume* envLV
  //     = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
  //   if ( envLV ) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  // }

  // if ( fEnvelopeBox ) {
  //   envSizeX = fEnvelopeBox->GetXHalfLength()*2.;
  //   envSizeY = fEnvelopeBox->GetYHalfLength()*2.;
  //   envSizeZ = fEnvelopeBox->GetZHalfLength()*2.;
  // }  
  // else  {
  //   G4ExceptionDescription msg;
  //   msg << "Envelope volume of box shape not found.\n"; 
  //   msg << "Perhaps you have changed geometry.\n";
  //   msg << "The gun will be place at the center.";
  //   G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
  //    "MyCode0002",JustWarning,msg);
  // }
  


  G4double worldZHalfLength = 0;
  G4LogicalVolume* worldLV
    = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox = NULL;
  if ( worldLV ) worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  if ( worldBox ) worldZHalfLength = worldBox->GetZHalfLength();
  else  {
    G4cerr << "World volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }

  // auto phi = ((2. * G4UniformRand()) - 1.) * static_cast<G4double> (std::atan(1.25/50.));
  // auto theta = ((2. * G4UniformRand()) - 1.) * static_cast<G4double> (std::atan(1.25/50.));



  // auto phi = G4UniformRand() * 2 * M_PI;
  // auto theta = G4UniformRand() * M_PI;

  // auto phi = uniRand(0., 2. * M_PI);
  // auto theta = uniRand(-M_PI/3, M_PI/3);
  // auto theta = uniRand(-0.018257970943255358, 0.018257970943255358);

  

  

  // auto x = std::cos(phi) * std::sin(theta);
  // auto y = std::sin(phi) * std::sin(theta);
  // auto z = std::cos(theta);
  
  // G4cout << "x = " << x << ", y = " << y << ", z = " << z << G4endl;


  // auto phi = std::uniform_real_distribution

  // Note that this particular case of starting a primary particle on the world boundary
  // requires shooting in a direction towards inside the world.
  // fParticleGun->SetParticlePosition(G4ThreeVector(0. * cm, 0. * cm, 0 * cm));
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));

  
  // G4double size = 0.8; 
  // G4double x0 = size * 10*cm * (G4UniformRand()-0.5);
  // G4double y0 = size * 10*cm * (G4UniformRand()-0.5);
  // G4double z0 = -0.5 * 10*cm;
  // std::vector<std::vector<G4bool>>

  // auto eangle = GetEnergyAngle_dist();
  // auto ez = GetEnergyZ_dist();
  // auto eanglez = GetEnergyAngleZ_bins();
  
  fGPS->GeneratePrimaryVertex(anEvent);
  
  // if(GetUseDists())
  // {
  //   std::uniform_int_distribution<int> En(0,eangle[0].size()-1);
  //   std::uniform_int_distribution<int> Ang(0,eangle.size()-1);
  //   std::uniform_int_distribution<int> Zh(0,ez.size()-1);

    

  //   int en{En(*generator)}, ang{Ang(*generator)}, z{Zh(*generator)};

  //   if(eangle[ang][en] == true && ez[z][en] == true)
  //   {
  //     std::uniform_real_distribution<double> Ebin(eanglez[0][en],eanglez[0][en+1]);
  //     std::uniform_real_distribution<double> AngleBin(eanglez[1][ang],eanglez[1][ang+1]);
  //     std::uniform_real_distribution<double> ZBin(eanglez[2][z],eanglez[2][z+1]);
  //     double E{Ebin(*generator)},Angle{AngleBin(*generator)},Zed{ZBin(*generator)};

  //     // fGPS->SetParticleEnergy((G4double)E * MeV);
  //     // auto direction = new G4ParticleMomentum();
  //     // direction->setTheta(Angle);
  //     G4cout << fGPS->GetParticleDefinition() << G4endl;
  //     G4cout << fGPS->GetParticleMomentumDirection() << G4endl;



  //     printf("Jackpot! E = %i, Ang = %i, Z = %i\n",en,ang,z);
  //     printf("Max Size: E = %i, Ang = %i, Z = %i\n",eangle[0].size(),eangle.size(),ez.size());
  //   } 
  //   else
  //   {
  //     // printf("Well, shit. E = %i, Ang = %i, Z = %i\n",en,ang,z);
  //   }


  //   // printf("E = %i, Ang = %i, Z = %i\n",En(*generator),Ang(*generator),Zh(*generator));


  //   // fGPS->GeneratePrimaryVertex(anEvent);
  //   // fGPS->
  // }

  // printf("An event!\n");
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

