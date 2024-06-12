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
/// \file B1ActionInitialization.hh
/// \brief Definition of the B1ActionInitialization class

#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

class ActionMessenger;
class CLYCDetectorConstruction;



/// Action initialization class.

class ActionInitialization : public G4VUserActionInitialization
{
  public:
    ActionInitialization(CLYCDetectorConstruction* detConstruction);
    virtual ~ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

    



    void SetEnergyAngleDist(G4String value) {EnergyAngleDist = value;};
    void SetEnergyZDist(G4String value) {EnergyZDist = value;};
    void SetEnergyAngleZBins(G4String value) {EnergyAngleZBins = value;};
    void SetNeutronsData(G4String value) {NeutronsData = value;};

    void SetUseDists(G4bool value) {UseDists = value;}; 
    void SetUseNeutrons(G4bool value) {UseNeutrons = value;}; 
    void SetSaveAnalysisVectors(G4bool savvecs) {fUseAnalysisVectors = savvecs;};

    
    G4String GetEnergyAngleDist() const {return EnergyAngleDist;};
    G4String GetEnergyZDist() const {return EnergyZDist;};
    G4String GetEnergyAngleZBins() const {return EnergyAngleZBins;};
    G4String GetNeutronsData() const {return NeutronsData;};
    

    G4bool GetUseDists() const {return UseDists;}; 
    G4bool GetUseNeutrons() const {return UseNeutrons;}; 
    G4bool GetSaveAnalysisVectors() const {return fUseAnalysisVectors;};

  private:
    CLYCDetectorConstruction* fDetConstruction;

    ActionMessenger*  fMessenger;


    G4String EnergyAngleDist{""};
    G4String EnergyZDist{""};
    G4String EnergyAngleZBins{""};
    G4String NeutronsData{""};

    G4bool UseDists{false};
    G4bool UseNeutrons{false};
    G4bool fUseAnalysisVectors{false};

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// inline void ActionInitialization::SetEnergyAngleDist(G4String value)
// {
//   EnergyAngleDist = value;
// }

// inline void ActionInitialization::SetEnergyZDist(G4String value)
// {
//   EnergyZDist = value;;
// }

// inline void ActionInitialization::SetUseDists(G4bool value)
// {
//   UseDists = value;
// }

#endif

    
