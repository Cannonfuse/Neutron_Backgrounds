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
#include "G4ParticleGun.hh"
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
  return PGA->fParticleGun->GetParticleEnergy();
}

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(nullptr)//, 
  //fEnvelopeBox(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4GeneralParticleSource();

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="gamma");
  // fParticleGun->SetParticleDefinition(particle);
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  // fParticleGun->SetParticleEnergy(4*MeV);

  // Hard-coded Be9(d,n) run

  /*
  double Estart[103] = {0.    ,  0.0594,  0.0999,  0.1505,  0.2001,  0.2497,  0.3002,
        0.3502,  0.3997,  0.4495,  0.4996,  0.5473,  0.5981,  0.6468,
        0.6964,  0.77  ,  0.8487,  0.9239,  1.0005,  1.0973,  1.197 ,
        1.2976,  1.3963,  1.5407,  1.6887,  1.8365,  1.9793,  2.1393,
        2.2881,  2.4185,  2.6001,  2.7223,  2.8149,  2.9124,  3.015 ,
        3.168 ,  3.3087,  3.459 ,  3.6197,  3.7625,  3.9138,  4.1079,
        4.2808,  4.4272,  4.5812,  4.7021,  4.8279,  4.9588,  5.049 ,
        5.1418,  5.2371,  5.3351,  5.4359,  5.5396,  5.6463,  5.7008,
        5.7561,  5.8122,  5.8692,  5.9856,  6.0452,  6.1056,  6.167 ,
        6.2292,  6.2925,  6.3567,  6.4219,  6.4881,  6.6236,  6.7635,
        6.9078,  7.0568,  7.2108,  7.3698,  7.5342,  7.7042,  7.88  ,
        8.062 ,  8.2504,  8.4454,  8.5456,  8.6476,  8.7514,  8.8571,
        8.9647,  9.0743,  9.186 ,  9.2997,  9.4156,  9.5337,  9.654 ,
        9.7766,  9.9016, 10.029 , 10.159 , 10.291 , 10.426 , 10.564 ,
       10.705 , 10.848 , 10.995 , 11.144 , 11.296};
  double Ewidth[103] = {0.0594, 0.0405, 0.0506, 0.0496, 0.0496, 0.0505, 0.05  , 0.0495,
       0.0498, 0.0501, 0.0477, 0.0508, 0.0487, 0.0496, 0.0736, 0.0787,
       0.0752, 0.0766, 0.0968, 0.0997, 0.1006, 0.0987, 0.1444, 0.148 ,
       0.1478, 0.1428, 0.16  , 0.1488, 0.1304, 0.1816, 0.1222, 0.0926,
       0.0975, 0.1026, 0.153 , 0.1407, 0.1503, 0.1607, 0.1428, 0.1513,
       0.1941, 0.1729, 0.1464, 0.154 , 0.1209, 0.1258, 0.1309, 0.0902,
       0.0928, 0.0953, 0.098 , 0.1008, 0.1037, 0.1067, 0.0545, 0.0553,
       0.0561, 0.057 , 0.1164, 0.0596, 0.0604, 0.0614, 0.0622, 0.0633,
       0.0642, 0.0652, 0.0662, 0.1355, 0.1399, 0.1443, 0.149 , 0.154 ,
       0.159 , 0.1644, 0.17  , 0.1758, 0.182 , 0.1884, 0.195 , 0.1002,
       0.102 , 0.1038, 0.1057, 0.1076, 0.1096, 0.1117, 0.1137, 0.1159,
       0.1181, 0.1203, 0.1226, 0.125 , 0.1274, 0.13  , 0.132 , 0.135 ,
       0.138 , 0.141 , 0.143 , 0.147 , 0.149 , 0.152 , 0.156 };
  double Yield[103] = {     0,    841,   1793,   1899,   2232,   2680,   2671,   2757,
         2749,   3072,   2509,   3133,   2807,   3019,  11049,  17513,
        14275,  14066,  28657,  32454,  33839,  32445, 102969, 115736,
       114196, 105487, 148595, 118890,  79821, 208893,  62947,  27098,
        32280,  36898, 118103,  88295, 105612, 122415,  83301,  95871,
       192578, 135124,  76765,  86991,  40176,  45227,  50758,  15786,
        16758,  16546,  16217,  16358,  17967,  19275,   2721,   2600,
         2820,   2906,  22984,   2508,   2117,   1549,   1190,    740,
          675,    837,    735,   6260,   7033,   7350,   8786,   9573,
        10255,  11491,  12448,  11393,  13695,  15775,  16585,   2338,
         2522,   2306,   2682,   2721,   2921,   2654,   2645,   2557,
         2475,   2525,   2511,   2696,   2457,   2457,   2773,   2536,
         2479,   2498,   2197,   1776,   1329,    981,    729};

  auto phi = uniRand(0., 2. * M_PI);
  auto theta = uniRand(0., M_PI);

  auto x = std::cos(phi) * std::sin(theta);
  auto y = std::sin(phi) * std::sin(theta);
  auto z = std::cos(theta);

  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0 * cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x,y,z));
  */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
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
  
  // fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

