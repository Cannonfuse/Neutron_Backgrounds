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
// #define USE_CADMESH_TETGEN 1


#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
// #include "tetgen.h"

class G4VPhysicalVolume;
class G4LogicalVolume;
class DetectorMessenger;
class G4VisAttributes;
class G4Material;
class OUG4Materials;


/// Detector construction class to define materials and geometry.

class CLYCDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    CLYCDetectorConstruction();
    virtual ~CLYCDetectorConstruction();

    // void DefineMaterials();

    virtual G4VPhysicalVolume* Construct();
    
    const G4VPhysicalVolume* GetC6LYCVolume() const;
    const G4VPhysicalVolume* GetC7LYCVolume() const;
    const G4VPhysicalVolume* GetdummydetectorVolume() const;
    const G4VPhysicalVolume* GetGasCellVolume() const;
    const G4VPhysicalVolume* GetGasCellReplicaVolume() const;
    const G4VPhysicalVolume* GetBe9TargetVolume() const;

    void SetDetDistance(G4double dist);
    void SetC7LYCDistance(G4double dist);
    void SetC6LYCDistance(G4double dist);
    void SetC7LYC_X(G4double dist);
    void SetC7LYC_Y(G4double dist);
    void SetC6LYC_X(G4double dist);
    void SetC6LYC_Y(G4double dist);
    void SetUseC6LYC(G4bool value) {UseC6LYC = value;};
    void SetUseC7LYC(G4bool value) {UseC7LYC = value;};
    void SetUseC6LYC_Case(G4bool value) {UseC6LYC_Case = value;};
    void SetUseC7LYC_Case(G4bool value) {UseC7LYC_Case = value;};
    void SetUseStructure(G4bool value) {UseStructure = value;};
    void SetUseDummy(G4bool value) {UseDummy = value;};
    void SetUseBe9target(G4bool value) {UseBe9target = value;};
    void SetUseLargeTarget(G4bool value) {UseLargeTarget = value;};
    void SetUseGasCell(G4bool value) {UseGasCell = value;};
    void SetUseFTC(G4bool value) {UseFTC = value;};
    void SetUseMTC(G4bool value) {UseMTC = value;};
    void SetUseLTC(G4bool value) {UseLTC = value;};
    void SetOverlaps(G4bool value) {fCheckOverlaps = value;}
    void setGasCellPressure(G4double pressure) {GasCellPressure = pressure;};
    void setGasCellPosition(G4double position) {GasCellPosition = position;};
    void setGasCellLength(G4double length) {GasCellLength = length;};
    void setGasCellDiameter(G4double diameter) {GasCellDiameter = diameter;};
    void setC6LYC_Slices(G4int slices) {C6LYC_Slices = slices;};
    void setC7LYC_Slices(G4int slices) {C7LYC_Slices = slices;};
    void setWorldMaterial(G4String mat) {WorldMat = mat;};
    void SetBe9Thickness(G4double thickness) {Be9TgtThickness = thickness;};


    G4double GetDetDistance();
    G4double GetC7LYCDistance() const {return C7LYCDetDistance;};
    G4double GetC6LYCDistance() const {return C6LYCDetDistance;};
    G4double GetC7LYC_X() const {return C7LYC_X;};
    G4double GetC7LYC_Y() const {return C7LYC_Y;};
    G4double GetC6LYC_X() const {return C6LYC_X;};
    G4double GetC6LYC_Y() const {return C6LYC_Y;};
    G4bool GetUseC6LYC() const {return UseC6LYC;};
    G4bool GetUseC7LYC() const {return UseC7LYC;};
    G4bool GetUseC6LYC_Case() const {return UseC6LYC_Case;};
    G4bool GetUseC7LYC_Case() const {return UseC7LYC_Case;};
    G4bool GetUseStructure() const {return UseStructure;};
    G4bool GetUseDummy() const {return UseDummy;};
    G4bool GetUseBe9target() const {return UseBe9target;};
    G4bool GetUseLargeTarget() const {return UseLargeTarget;};
    G4bool GetUseGasCell() const {return UseGasCell;};
    G4bool GetUseFTC() const {return UseFTC;};
    G4bool GetUseMTC() const {return UseMTC;};
    G4bool GetUseLTC() const {return UseLTC;};
    G4bool GetOverlaps() const {return fCheckOverlaps;};
    G4double GetGasCellPressure() const {return GasCellPressure;};
    G4double GetGasCellPosition() const {return GasCellPosition;};
    G4double GetGasCellLength() const {return GasCellLength;};
    G4double GetGasCellDiameter() const {return GasCellDiameter;};
    G4int GetC6LYC_Slices() const {return C6LYC_Slices;};
    G4int GetC7LYC_Slices() const {return C7LYC_Slices;};
    G4String GetWorldMaterial() const {return WorldMat;};
    G4double GetBe9Thickness() {return Be9TgtThickness;};

    //C7LYC parts
    G4double GetC7LYC_TopcapPosition() {return C7LYC_TopcapPosition;};
    void SetC7LYC_TopcapPosition(G4double value) {C7LYC_TopcapPosition = value;};

    // Material Set/Get functions
    // void SetC6LYCMaterial(G4Material* mat) {C6LYC = mat;};
    // void SetC7LYCMaterial(G4Material* mat) {C7LYC = mat;};
    // void SetHAVARMaterial(G4Material* mat) {HAVAR = mat;};
    // void SetMuMetalMaterial(G4Material* mat) {MuMetal = mat;};
    // void SetQuartzMaterial(G4Material* mat) {Quartz = mat;};


    // G4Material* GetC6LYCMaterial() {return C6LYC;};
    // G4Material* GetC7LYCMaterial() {return C7LYC;};
    // G4Material* GetHAVARMaterial() {return HAVAR;};
    // G4Material* GetMuMetalMaterial() {return MuMetal;};
    // G4Material* GetQuartzMaterial() {return Quartz;};


    OUG4Materials* GetOUGeant4Mats() {return OUMats;};
    void SetOUGeant4Mats(OUG4Materials* themats) {OUMats =  themats;};



    void buildBe9Target(G4LogicalVolume* theWorld, G4Material* theMaterial, 
                        G4VisAttributes* theColor, G4bool overlaps);
    void buildMagnetStructure(G4LogicalVolume* theWorld, G4Material* theMaterial, 
                              G4VisAttributes* theColor, G4bool overlaps);
    void buildSphereDetector(G4LogicalVolume* theWorld, G4Material* theMaterial, 
                              G4VisAttributes* theColor, G4bool overlaps);





  private:
    // void DefineCommands();

    DetectorMessenger*  fMessenger;

    OUG4Materials* OUMats;



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
    G4double GasCellPressure{0};
    G4double GasCellPosition{0};
    G4double GasCellLength{0};
    G4double GasCellDiameter{0};
    G4double Be9TgtThickness{0};
    G4bool UseC6LYC{false};
    G4bool UseC7LYC{false};
    G4bool UseStructure{false};
    G4bool UseDummy{false};
    G4bool UseBe9target{false};
    G4bool UseLargeTarget{false};
    G4bool UseGasCell{false};
    G4bool UseC6LYC_Case{false};
    G4bool UseC7LYC_Case{false};

    G4int C6LYC_Slices{1};
    G4int C7LYC_Slices{1};

    G4String WorldMat{"G4_Galactic"};



    // C7LYC Parts
    G4double C7LYC_TopcapPosition{0};


    G4VPhysicalVolume* fGasCellPV;
    G4VPhysicalVolume* fGasCellReplicaPV;
    G4VPhysicalVolume* fHAVARFoilPV;
    G4VPhysicalVolume* fGasCellStainlessCasePV;
    G4VPhysicalVolume*  fC7LYCPV;
    G4VPhysicalVolume*  fC6LYCPV;
    G4VPhysicalVolume*  fdummydetectorPV;
    G4VPhysicalVolume* fSphereDetector;
    G4VPhysicalVolume* fBe9Target;



    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

    // G4Material* MuMetal;
    // G4Material* HAVAR;
    // G4Material* C7LYC;
    // G4Material* C6LYC;
    // G4Material* Quartz;

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

inline const G4VPhysicalVolume* CLYCDetectorConstruction::GetGasCellVolume() const
{
  return fGasCellPV;
}

inline const G4VPhysicalVolume* CLYCDetectorConstruction::GetGasCellReplicaVolume() const
{
  return fGasCellReplicaPV;
}

inline const G4VPhysicalVolume* CLYCDetectorConstruction::GetBe9TargetVolume() const
{
  return fBe9Target;
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