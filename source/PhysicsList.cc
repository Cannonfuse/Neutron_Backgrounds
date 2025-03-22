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
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <iomanip>   

#include <CLHEP/Units/SystemOfUnits.h>
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4IonPhysicsPHP.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronElasticPhysicsPHP.hh"

#include "PhysicsList.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronicParameters.hh"
#include "G4HadronicProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
{
  G4cout << "<<< Geant4 Physics List simulation engine: QGSP_BERT_HP_MOD - removed G4RadioactiveDecayPhysics"<<G4endl;
  G4cout <<G4endl<<G4endl;

  G4int ver = 1;

  defaultCutValue = 0.0*CLHEP::mm;  
  SetVerboseLevel(ver);

  // EM Physics
  RegisterPhysics( new G4EmStandardPhysics(ver) );

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays
  RegisterPhysics( new G4DecayPhysics(ver) );
  RegisterPhysics( new G4RadioactiveDecayPhysics(ver) );

   // Hadron Elastic scattering
   RegisterPhysics( new G4HadronElasticPhysicsPHP(ver) );

  // Hadron Physics
  // First enable the scaling of cross section (by default disabled):
  G4HadronicParameters::Instance()->SetApplyFactorXS( false );
  // // Scaling up the nucleon inelastic cross sections by 10%

  G4cout << "Cross Section Scale Factor: " << G4HadronicParameters::Instance()->XSFactorNucleonInelastic() << ", Enabled: "<< G4HadronicParameters::Instance()->ApplyFactorXS() << G4endl;

 
  G4HadronPhysicsQGSP_BERT_HP *QGSP_BERT_HP = new G4HadronPhysicsQGSP_BERT_HP(ver);
  RegisterPhysics( QGSP_BERT_HP);
  G4HadronicParameters::Instance()->SetXSFactorNucleonInelastic( 1000. );

  G4cout << "Cross Section Scale Factor: " << G4HadronicParameters::Instance()->XSFactorNucleonInelastic() << ", Enabled: "<< G4HadronicParameters::Instance()->ApplyFactorXS() << G4endl;
  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(ver));

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver));

  // RegisterPhysics( new G4HadronicProcess(ver));
  // RegisterPhysics( new G4IonPhysicsPHP(ver));

  // auto HadParam = new G4HadronicParameters::Instance();
  // HadParam->SetApplyFactorXS( true );
  // HadParam->SetXSFactorNucleonInelastic( 1000. );
  // First enable the scaling of cross section (by default disabled):
  // G4HadronicParameters::Instance()->SetApplyFactorXS( true );

  // // Scaling up the nucleon inelastic cross sections by 10%
  // G4HadronicParameters::Instance()->SetXSFactorNucleonInelastic( 100000. );
  // // G4HadronicParameters::Instance()->SetXSFactorNucleonElastic( 10. );

  // G4HadronicParameters::Instance()->SetVerboseLevel(2);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
  SetCutValue(0*CLHEP::mm, "neutron");
  // SetCutValue(0.1*CLHEP::radian,"neutron");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
