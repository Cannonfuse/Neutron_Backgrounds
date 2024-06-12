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
/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "ActionInitialization.hh"
#include "ActionMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "csv.hpp"
#include "random"

using namespace csv;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization(CLYCDetectorConstruction* detConstruction)
 : G4VUserActionInitialization(),
  fDetConstruction(detConstruction)
{

    fMessenger = new ActionMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
  auto eventAction = new EventAction(fDetConstruction);
  auto usevectors = GetSaveAnalysisVectors();

  SetUserAction(new RunAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
  auto PGA = new PrimaryGeneratorAction();

  SetUserAction(PGA);

  auto usedists = GetUseDists();
  auto useneutrons = GetUseNeutrons();
  auto usevectors = GetSaveAnalysisVectors();

  auto EA = new EventAction(fDetConstruction);
  auto RA = new RunAction(EA);

  if(usedists)
  {
    auto eangdist = GetEnergyAngleDist();
    auto ezdist = GetEnergyZDist();
    auto eangzbins = GetEnergyAngleZBins();

    std::vector<std::vector<bool>> eang, ez;
    std::vector<std::vector<double>> bins;

    CSVFormat format;
    format.no_header();  // Parse CSVs without a header row

    // format.delimiter(',').quote('~').no_header();  // Parse CSVs without a header row
    CSVReader eangreader(eangdist,format);
    CSVReader ezreader(ezdist,format);
    CSVReader eangzreader(eangzbins,format);

    CSVRow row;

    for (CSVRow& row: eangreader) // Input iterator
    { 
      std::vector<bool> rowdat;
      for (CSVField& field: row) 
      {
          // By default, get<>() produces a std::string.
          // A more efficient get<string_view>() is also available, where the resulting
          // string_view is valid as long as the parent CSVRow is alive
          rowdat.push_back(field.get<bool>());
          // std::cout << field.get<bool>() << ",";
      }
      eang.push_back(rowdat);
    }
    for (CSVRow& row: ezreader) // Input iterator
    { 
      std::vector<bool> rowdat;
      for (CSVField& field: row) 
      {
          // By default, get<>() produces a std::string.
          // A more efficient get<string_view>() is also available, where the resulting
          // string_view is valid as long as the parent CSVRow is alive
          rowdat.push_back(field.get<bool>());
          // std::cout << field.get<bool>() << ",";
      }
      ez.push_back(rowdat);
    }
    for (CSVRow& row: eangzreader) // Input iterator
    { 
      std::vector<double> rowdat;
      for (CSVField& field: row) 
      {
          // By default, get<>() produces a std::string.
          // A more efficient get<string_view>() is also available, where the resulting
          // string_view is valid as long as the parent CSVRow is alive
          rowdat.push_back(field.get<double>());
          // std::cout << field.get<bool>() << ",";
      }
      bins.push_back(rowdat);
    }

    auto DistsEA = new EventAction(fDetConstruction, usedists, eang, ez, bins);//,generator);
    EA = DistsEA;
  } 
  if(useneutrons)
  {
    auto neutronsdata = GetNeutronsData();
    std::vector<std::vector<double>> neutrons;

    CSVFormat format;
    format.no_header();  // Parse CSVs without a header row

    // format.delimiter(',').quote('~').no_header();  // Parse CSVs without a header row
    CSVReader neutronsreader(neutronsdata,format);

    for (CSVRow& row: neutronsreader) // Input iterator
    { 
      std::vector<double> rowdat;
      for (CSVField& field: row) 
      {
          // By default, get<>() produces a std::string.
          // A more efficient get<string_view>() is also available, where the resulting
          // string_view is valid as long as the parent CSVRow is alive
          rowdat.push_back(field.get<double>());
          // std::cout << field.get<bool>() << ",";
      }
      neutrons.push_back(rowdat);
    }
    auto NeutronsDataEA = new EventAction(fDetConstruction, useneutrons, neutrons);//,generator);
    EA = NeutronsDataEA;
  }

  SetUserAction(EA);
  SetUserAction(RA);
  SetUserAction(new SteppingAction(fDetConstruction,EA));
}  

// void ActionInitialization::SetEnergyAngleDist(G4String value)
// {
//   EnergyAngleDist = value;
// }

// void ActionInitialization::SetEnergyZDist(G4String value)
// {
//   EnergyZDist = value;;
// }

// void ActionInitialization::SetUseDists(G4bool value)
// {
//   UseDists = value;
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
