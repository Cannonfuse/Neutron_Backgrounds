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
/// \file CLYCDetectorConstruction.hh
/// \brief Definition of the CLYCDetectorConstruction class

#ifndef CLYCDetectorConstruction_h
#define CLYCDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class DetectorMessenger;


/// Detector construction class to define materials and geometry.

class CLYCDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    CLYCDetectorConstruction();
    virtual ~CLYCDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    
    const G4VPhysicalVolume* GetC6LYCVolume() const;
    const G4VPhysicalVolume* GetC7LYCVolume() const;
    const G4VPhysicalVolume* GetdummydetectorVolume() const;



    void SetDetDistance(G4double dist);
    void SetC7LYCDistance(G4double dist);
    void SetC6LYCDistance(G4double dist);
    void SetC7LYC_X(G4double dist);
    void SetC7LYC_Y(G4double dist);
    void SetC6LYC_X(G4double dist);
    void SetC6LYC_Y(G4double dist);
    void SetUseC6LYC(G4bool value) {UseC6LYC = value;};
    void SetUseC7LYC(G4bool value) {UseC7LYC = value;};
    void SetUseStructure(G4bool value) {UseStructure = value;};
    void SetUseDummy(G4bool value) {UseDummy = value;};
    void SetUseBe9target(G4bool value) {UseBe9target = value;};
    void SetUseFTC(G4bool value) {UseFTC = value;};
    void SetUseMTC(G4bool value) {UseMTC = value;};
    void SetUseLTC(G4bool value) {UseLTC = value;};


    G4double GetDetDistance();
    G4double GetC7LYCDistance() {return C7LYCDetDistance;};
    G4double GetC6LYCDistance() {return C6LYCDetDistance;};
    G4double GetC7LYC_X() {return C7LYC_X;};
    G4double GetC7LYC_Y() {return C7LYC_Y;};
    G4double GetC6LYC_X() {return C6LYC_X;};
    G4double GetC6LYC_Y() {return C6LYC_Y;};
    G4bool GetUseC6LYC() {return UseC6LYC;};
    G4bool GetUseC7LYC() {return UseC7LYC;};
    G4bool GetUseStructure() {return UseStructure;};
    G4bool GetUseDummy() {return UseDummy;};
    G4bool GetUseBe9target() {return UseBe9target;};
    G4bool GetUseFTC() {return UseFTC;};
    G4bool GetUseMTC() {return UseMTC;};
    G4bool GetUseLTC() {return UseLTC;};


  private:
    // void DefineCommands();

    DetectorMessenger*  fMessenger;

    G4bool UseFTC{false};
    G4bool UseMTC{false};
    G4bool UseLTC{false};
    G4double DetDistance;
    G4double C7LYCDetDistance{0};
    G4double C6LYCDetDistance{0};
    G4double C7LYC_X{0};
    G4double C7LYC_Y{0};
    G4double C6LYC_X{0};
    G4double C6LYC_Y{0};
    G4bool UseC6LYC{false};
    G4bool UseC7LYC{false};
    G4bool UseStructure{false};
    G4bool UseDummy{false};
    G4bool UseBe9target{false};



    G4VPhysicalVolume*  fC7LYCPV;
    G4VPhysicalVolume*  fC6LYCPV;
    G4VPhysicalVolume*  fdummydetectorPV;

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline const G4VPhysicalVolume* CLYCDetectorConstruction::GetC6LYCVolume() const
{
  return fC6LYCPV;
}

inline const G4VPhysicalVolume* CLYCDetectorConstruction::GetC7LYCVolume() const
{
  return fC7LYCPV;
}

inline const G4VPhysicalVolume* CLYCDetectorConstruction::GetdummydetectorVolume() const
{
  return fdummydetectorPV;
}

inline void CLYCDetectorConstruction::SetDetDistance(G4double dist)
{
  DetDistance = dist;
}

inline void CLYCDetectorConstruction::SetC7LYCDistance(G4double dist)
{
  C7LYCDetDistance = dist;
}

inline void CLYCDetectorConstruction::SetC6LYCDistance(G4double dist)
{
  C6LYCDetDistance = dist;
}

inline void CLYCDetectorConstruction::SetC7LYC_X(G4double dist)
{
  C7LYC_X = dist;
}

inline void CLYCDetectorConstruction::SetC7LYC_Y(G4double dist)
{
  C7LYC_Y = dist;
}

inline void CLYCDetectorConstruction::SetC6LYC_X(G4double dist)
{
  C6LYC_X = dist;
}

inline void CLYCDetectorConstruction::SetC6LYC_Y(G4double dist)
{
  C6LYC_Y = dist;
}


inline G4double CLYCDetectorConstruction::GetDetDistance()
{
  return DetDistance;
}

// inline G4double CLYCDetectorConstruction::GetC7LYCDistance()
// {
//   return C7LYCDetDistance;
// }

// inline G4double CLYCDetectorConstruction::C6LYCDetDistance()
// {
//   return C6LYCDetDistance;
// }


#endif