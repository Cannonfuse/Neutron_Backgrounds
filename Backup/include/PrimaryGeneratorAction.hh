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
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <random>

class G4ParticleGun;
class G4Event;
class G4Box;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued 
/// in front of the phantom across 80% of the (X,Y) phantom size.

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();    
    virtual ~PrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);

    G4double retBeamEn() const;
    // G4ThreeVector particleDirection() const;         
  
    // method to access particle gun
    // const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
  
  private:
    G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
    // G4GeneralParticleSource* fParticleSource;
    // std::mt19937_64* GENERATOR;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double PrimaryGeneratorAction::retBeamEn() const
{
  auto retValue = fParticleGun->GetParticleEnergy();
  return retValue;
}

// inline G4ThreeVector PrimaryGeneratorAction::particleDirection() const
// {

//   std::uniform_real_distribution<G4double> uDist(0,1);

//   auto phi = uDist(GENERATOR) * static_cast<G4double> (std::atan(1.25/50.));
//   auto theta = uDist(GENERATOR) * static_cast<G4double> (std::atan(1.25/50.));

//   auto x = std::cos(theta) * std::sin(phi);
//   auto y = std::sin(theta) * std::sin(phi);
//   auto z = std::sqrt(1 - std::pow(x,2) - std::pow(y,2));
//   auto retValue = fParticleGun->GetParticleEnergy();
//   return retValue;
// }

#endif
