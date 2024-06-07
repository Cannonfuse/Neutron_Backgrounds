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
#include <vector>

// class G4ParticleGun;
class G4GeneralParticleSource;
// class std::mt19937_64;
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
    // PrimaryGeneratorAction(const G4bool usedists, 
    //                        const std::vector<std::vector<G4bool>> energyangledist,
    //                        const std::vector<std::vector<G4bool>> energyzdist,
    //                        const std::vector<std::vector<double>> energyangzbins);//,const std::default_random_engine *generator);    
    virtual ~PrimaryGeneratorAction();

    // void SetEnergyAngle_dist(std::vector<std::vector<G4bool>> dist) {EnergyAngle_dist = dist;};
    // void SetEnergyZ_dist(std::vector<std::vector<G4bool>> dist) {EnergyZ_dist = dist;};
    // void SetEnergyAngleZ_bins(std::vector<std::vector<double>> dist) {EnergyAngleZ_bins = dist;};
    // void SetUseDists(G4bool usedists) {UseDists = usedists;};
    // // void SetGenerator(std::default_random_engine *gen) {generator = gen;};


    // std::vector<std::vector<G4bool>>& GetEnergyAngle_dist() {return EnergyAngle_dist;};
    // std::vector<std::vector<G4bool>>& GetEnergyZ_dist() {return EnergyZ_dist;};
    // std::vector<std::vector<double>>& GetEnergyAngleZ_bins() {return EnergyAngleZ_bins;};
    // G4bool GetUseDists() {return UseDists;};
    // std::default_random_engine GetGenerator() {return generator;};


    // method from the base class
    virtual void GeneratePrimaries(G4Event*);

    friend G4double retBeamEn(const PrimaryGeneratorAction *PGA);

    friend double uniRand(double start, double end);

    // G4ThreeVector particleDirection() const;         
  
    // method to access particle gun
    // const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
  
  private:
    G4GeneralParticleSource*  fGPS; // pointer a to G4 gun class
    G4double beamEn;
    
    // std::vector<std::vector<G4bool>> EnergyAngle_dist;
    // std::vector<std::vector<G4bool>> EnergyZ_dist;
    // std::vector<std::vector<double>> EnergyAngleZ_bins;

    // G4bool UseDists{false};
    // std::mt19937_64 *generator;

    // G4GeneralParticleSource* fParticleSource;
    // std::mt19937_64* GENERATOR;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
