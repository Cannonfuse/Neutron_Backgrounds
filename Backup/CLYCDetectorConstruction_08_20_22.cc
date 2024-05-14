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
/// \file CLYCDetectorConstruction.cc
/// \brief Implementation of the CLYCDetectorConstruction class

#include <string>

#include "CLYCDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"

#include "G4Scintillation.hh"
#include "G4MaterialPropertiesTable.hh"

// for reading JSON files
#include "nlohmann/json.hpp"

// #define USE_CADMESH_TETGEN

#include "CADMesh.hh"

using json = nlohmann::json;


// const G4double det_dist = 5 * m;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CLYCDetectorConstruction::CLYCDetectorConstruction()
: G4VUserDetectorConstruction(),
  fCLYCPV(nullptr)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CLYCDetectorConstruction::~CLYCDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* CLYCDetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  // Detector distance; original distance was face of crystal @ 0.5 m, need to account for that with second line.
  G4double det_dist = 5 * m - 1.25/2 * cm;
  // det_dist -= 0.5 * m;

  // PE shield around detector, true = use it, false = don't use it

  bool useShield = false;


  // Specs for fine tunnel collimator (vault -> tunnel)

  G4double ftc_ID1[15] = {2.49 * cm,2.35 * cm,2.01 * cm,1.445 * cm,1.25 * cm,0.335 * cm,0.375 * cm,0.52 * cm,0.61 * cm,0.65 * cm,0.785 * cm,0.9 * cm,1.04 * cm,1.165 * cm,1.27 * cm};
  G4double ftc_ID2[15] = {2.335 * cm,2.005 * cm,1.53 * cm,1.3 * cm,1.06 * cm,0.41 * cm,0.52 * cm,0.645 * cm,0.645 * cm,0.75 * cm,0.945 * cm,0.965 * cm,1.17 * cm,1.27 * cm,1.36 * cm};
  G4double ftc_OD1[15] = {3.18 * cm,3.28 * cm,3.285 * cm,3.29 * cm,3.41 * cm,3.425 * cm,3.49 * cm,3.525 * cm,3.6 * cm,3.62 * cm,3.59 * cm,3.825 * cm,3.94 * cm,3.96 * cm,4.125 * cm};
  G4double ftc_OD2[15] = {3.32 * cm,3.225 * cm,3.32 * cm,3.37 * cm,3.45 * cm,3.45 * cm,3.52 * cm,3.605 * cm,3.62 * cm,3.715 * cm,3.825 * cm,3.94 * cm,3.98 * cm,4.105 * cm,4.235 * cm};
  G4double ftc_LEN[15] = {2 * cm,3.67 * cm,4.58 * cm,2.53 * cm,3.855 * cm,3.895 * cm,5.01 * cm,5.115 * cm,5.305 * cm,7.52 * cm,7.545 * cm,7.555 * cm,7.58 * cm,7.5 * cm,7.55 * cm};
  G4double ftc_POS[15] = {2 * cm,7.67 * cm,15.92 * cm,23.03 * cm,29.415 * cm,37.165 * cm,46.07 * cm,56.195 * cm,66.615 * cm,79.44 * cm,94.505 * cm,109.605 * cm,124.74 * cm,139.82 * cm,154.87 * cm};

  // Specs for medium tunnel collimator (vault -> tunnel)

  G4double mtc_ID1[6] = {3.125 * cm,3.3 * cm,3.4 * cm,3.65 * cm,3.85 * cm,4.05 * cm};
  G4double mtc_ID2[6] = {3.25 * cm,3.425 * cm,3.65 * cm,3.85 * cm,4.05 * cm,4.25 * cm};
  G4double mtc_OD1[6] = {6.35 * cm,6.8 * cm,7.3 * cm,8 * cm,8.7 * cm,9.5 * cm};
  G4double mtc_OD2[6] = {6.8 * cm,7.3 * cm,8.05 * cm,8.8 * cm,9.5 * cm,10.25 * cm};
  G4double mtc_LEN[6] = {10.1 * cm,10.15 * cm,15.2 * cm,15.25 * cm,15.25 * cm,15.25 * cm};
  G4double mtc_POS[6] = {10.1 * cm,30.35 * cm,55.7 * cm,86.15 * cm,116.65 * cm,147.15 * cm};

  // Specs for large tunnel collimator (vault -> tunnel)


  G4double ltc_ID1[3] = {8.1 * cm,9.55 * cm,10.45 * cm};
  G4double ltc_ID2[3] = {9.5 * cm,10.25 * cm,11 * cm};
  G4double ltc_OD1[3] = {14.9 * cm,14.9 * cm,14.9 * cm};
  G4double ltc_OD2[3] = {14.9 * cm,14.9 * cm,14.9 * cm};
  G4double ltc_LEN[3] = {30.5 * cm,15.25 * cm,14.65 * cm};
  G4double ltc_POS[3] = {30.5 * cm,76.25 * cm,106.15 * cm};

  // y height

  G4double y_height = 0 * m;

  //     
  // World
  //
  G4double base_world_unit = 1*m;
  G4double world_sizeX = 10 * base_world_unit;
  G4double world_sizeY = 10 * base_world_unit;
  G4double world_sizeZ  = 60 * base_world_unit;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       world_sizeX, world_sizeY, world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
  

  // CLYC detector material

  G4MaterialPropertiesTable* CLYC_MPT = new G4MaterialPropertiesTable();

  G4int nEntries = 8;
  G4int nEntries2 = 3;



  std::vector<G4double> CLYC_photon_energies = {2.8834 * eV, 3.0613 * eV, 3.2627 * eV, 3.3509 * eV, 3.4440 * eV, 3.8149 * eV, 4.1328 * eV, 4.5085 * eV};
  std::vector<G4double> CLYC_C1_percent = {0., 0., 0., 0., 0. , 0.1, 1., 0.1};
  std::vector<G4double> CLYC_C2_percent = {0., 0., 0., 0.1, 1., 0.1, 0., 0.};
  std::vector<G4double> CLYC_C3_percent = {0.1, 1., 0.1, 0., 0., 0., 0., 0.};
  std::vector<G4double> CLYC_refraction = {1.8, 1.81, 1.81, 1.8, 1.81, 1.81, 1.81, 1.81};
  std::vector<G4double> CLYC_absorption_length = {3.42 * cm, 3.42 * cm, 3.42 * cm, 3.42 * cm, 3.42 * cm, 3.42 * cm, 3.42 * cm, 3.42 * cm};

  G4double CLYC_time_constant[nEntries2] = {1 * ns, 50 * ns, 1000 * ns};
  G4double CLYC_scint_yield[nEntries2] = {0.15, 0.50, 0.35};


  G4double CLYC_scint_energies[nEntries] = {2.7552 * eV, 4.5085 *eV};

  CLYC_MPT->AddProperty("RINDEX",CLYC_photon_energies,CLYC_refraction);
  CLYC_MPT->AddProperty("ABSLENGTH", CLYC_photon_energies,CLYC_absorption_length);
  CLYC_MPT->AddProperty("SCINTILLATIONCOMPONENT1",CLYC_photon_energies,CLYC_C1_percent);
  CLYC_MPT->AddProperty("SCINTILLATIONCOMPONENT2",CLYC_photon_energies,CLYC_C2_percent);
  CLYC_MPT->AddProperty("SCINTILLATIONCOMPONENT3",CLYC_photon_energies,CLYC_C3_percent);


  CLYC_MPT->AddConstProperty("SCINTILLATIONYIELD",20000./MeV);
  CLYC_MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  CLYC_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1",CLYC_time_constant[0]);
  CLYC_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT2",CLYC_time_constant[1]);
  CLYC_MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT3",CLYC_time_constant[2]);
  CLYC_MPT->AddConstProperty("SCINTILLATIONYIELD1",CLYC_scint_yield[0]);
  CLYC_MPT->AddConstProperty("SCINTILLATIONYIELD2",CLYC_scint_yield[1]);
  CLYC_MPT->AddConstProperty("SCINTILLATIONYIELD3",CLYC_scint_yield[2]);



  

  G4double a, density;
  G4int z, n, ncomponents, natoms;
  G4String name, symbol;
  G4int nLi_i{2}, nCl_i{2};

  // Define Li6 enhanced 
  //G4double Li6_a,Li7_a;//,Li6enh_a;
  //Li6_a = 6.01512288742 * g/mole;
  //Li7_a = 7.01600343666 * g/mole;
  // Li6enh_a = 0.05 * Li7_a + 0.95 * Li6_a;

  // G4Isotope* Li6_i = new G4Isotope(name="Li6",z=3,n=6,a=Li6_a);
  // G4Isotope* Li7_i = new G4Isotope(name="Li7",z=3,n=7,a=Li7_a);
  G4Element* Li6Enh_el = new G4Element(name="Li6Enhanced",symbol="Li6Enh" ,nLi_i);
  G4Element* Li7Enh_el = new G4Element(name="Li7Enhanced",symbol="Li7Enh" ,1);

  // G4Element* Li6Enh_el = new G4Element(name="Li6Enhanced",symbol="Li6Enh" , z= 3., Li6enh_a);
  // Li6Enh_el->AddIsotope(Li6_i,95*perCent);
  // Li6Enh_el->AddIsotope(Li7_i,5*perCent);


  // Get the natural abundance litium from the geant4 database
  G4Element* Li_i = nist->FindOrBuildElement("Li",true);
  
  // Get the Li6 and Li7 isotopes from the geant4 database; must use const_cast<G4Isotope *> () to successfully use in construction of enriched element.
  // Get Li6 and print its properties
  G4cout << *Li_i->GetIsotopeVector()->at(0) << G4endl;
  G4Isotope* Li6_i = const_cast<G4Isotope *> (Li_i->GetIsotope(0));
  
  // Get Li7 and print its properties
  G4cout << *Li_i->GetIsotopeVector()->at(1) << G4endl;
  G4Isotope* Li7_i = const_cast<G4Isotope *> (Li_i->GetIsotope(1));

  // Add Li6 and Li7 to the enriched element
  Li6Enh_el->AddIsotope(Li6_i, 95 * perCent);
  Li6Enh_el->AddIsotope(Li7_i, 5 * perCent);
  // Li6Enh_el->AddIsotope(Li7_i, 100 * perCent);

  Li7Enh_el->AddIsotope(Li7_i, 100 * perCent);


  G4Element* Cl_el = new G4Element(name="ClNatural",symbol="Clnat" ,nCl_i);
  // G4Element* Li6Enh_el = new G4Element(name="Li6Enhanced",symbol="Li6Enh" , z= 3., Li6enh_a);
  // Li6Enh_el->AddIsotope(Li6_i,95*perCent);
  // Li6Enh_el->AddIsotope(Li7_i,5*perCent);


  // Get the natural abundance litium from the geant4 database
  G4Element* Cl_i = nist->FindOrBuildElement("Cl",true);
  
  // Get the Li6 and Li7 isotopes from the geant4 database; must use const_cast<G4Isotope *> () to successfully use in construction of enriched element.
  // Get Li6 and print its properties
  G4cout << *Cl_i->GetIsotopeVector()->at(0) << G4endl;
  G4Isotope* Cl35_i = const_cast<G4Isotope *> (Cl_i->GetIsotope(0));
  
  // Get Li7 and print its properties
  G4cout << *Cl_i->GetIsotopeVector()->at(1) << G4endl;
  G4Isotope* Cl37_i = const_cast<G4Isotope *> (Cl_i->GetIsotope(1));

  // Add Li6 and Li7 to the enriched element
  Cl_el->AddIsotope(Cl35_i, 76.76 * perCent);
  Cl_el->AddIsotope(Cl37_i, 24.24 * perCent);

  // Get the other elements from constructing the detector from the geant4 database. All other elements used in construction of the detector are of natural abundance.
  G4Element* elCl = nist->FindOrBuildElement(17);
  G4Element* elY = nist->FindOrBuildElement(39);
  G4Element* elCs = nist->FindOrBuildElement(55);
  G4Element* elCe = nist->FindOrBuildElement(58);

  // Define the CLYC molecule
  // density = 3.31 * g/cm3;
  density = 3.31 * g/cm3;
  G4Material* CLYC_molecule = new G4Material(name="CLYC",density,ncomponents=4);
  CLYC_molecule->AddElement(elCs, 20*perCent);
  CLYC_molecule->AddElement(Li6Enh_el, 10*perCent);
  CLYC_molecule->AddElement(elY, 10*perCent);
  CLYC_molecule->AddElement(elCl, 60*perCent);

  CLYC_molecule->SetMaterialPropertiesTable(CLYC_MPT);
  CLYC_molecule->GetIonisation()->SetBirksConstant(6.95e-4 * cm / MeV);

  G4Material* C7LYC_molecule = new G4Material(name="C7LYC",density,ncomponents=4);
  C7LYC_molecule->AddElement(elCs, 20*perCent);
  C7LYC_molecule->AddElement(Li7Enh_el, 10*perCent);
  C7LYC_molecule->AddElement(elY, 10*perCent);
  C7LYC_molecule->AddElement(elCl, 60*perCent);

  C7LYC_molecule->SetMaterialPropertiesTable(CLYC_MPT);
  C7LYC_molecule->GetIonisation()->SetBirksConstant(6.95e-4 * cm / MeV);


  // density = 3.31 * g/cm3;
  // G4Material* CLYC = new G4Material(name="CLYC",density,ncomponents=2);
  // CLYC->AddElementByMassFraction(elCe, 0.5 * perCent);
  // CLYC->AddMaterial(CLYC_molecule, 99.5 * perCent);





  G4Material* WinH_mat = nist->FindOrBuildMaterial("G4_Al");
  G4Material* DC_mat = nist->FindOrBuildMaterial("G4_BRASS");
  G4Material* SC_mat = nist->FindOrBuildMaterial("G4_Al");

  

  // Source casing

  /*
  G4double innerRadiusSC1= 1.087*cm;
  G4double outerRadiusSC1 = 1.25*cm;
  G4double hzSC1 = 0.25*cm;
  G4double innerRadiusSC2 = 0*cm;
  G4double outerRadiusSC2= 1.25*cm;
  G4double hzSC2 = 0.0815 * cm;
  G4double startAngleSC = 0.*deg;
  G4double spanningAngleSC = 360.*deg;
  G4ThreeVector SCpos1 = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, .485* base_world_unit + det_dist);
  G4ThreeVector SCpos2 = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, .488315* base_world_unit + det_dist);
  G4ThreeVector SCpos3 = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, .483315* base_world_unit + det_dist);


  
   G4Tubs* SC1
     = new G4Tubs("SourceCasing1",
                  innerRadiusSC1,
                  outerRadiusSC1,
                  hzSC1,
                  startAngleSC,
                  spanningAngleSC);
    
    G4Tubs* SC2
     = new G4Tubs("SourceCasing2",
                  innerRadiusSC2,
                  outerRadiusSC2,
                  hzSC2,
                  startAngleSC,
                  spanningAngleSC);

    G4Tubs* SC3
     = new G4Tubs("SourceCasing3",
                  innerRadiusSC2,
                  outerRadiusSC2,
                  hzSC2,
                  startAngleSC,
                  spanningAngleSC);


  G4LogicalVolume* logicSC1 =                         
    new G4LogicalVolume(SC1,         //its solid
                        SC_mat,          //its material
                        "SourceCasing1");           //its name

  new G4PVPlacement(0,                       //no rotation
    SCpos1,                    //at position
    logicSC1,               //its logical volume
    "SourceCasing1",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking

  G4LogicalVolume* logicSC2 =                         
    new G4LogicalVolume(SC2,         //its solid
                        SC_mat,          //its material
                        "SourceCasing2");           //its name

  new G4PVPlacement(0,                       //no rotation
    SCpos2,                    //at position
    logicSC2,               //its logical volume
    "SourceCasing2",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking

  G4LogicalVolume* logicSC3 =                         
    new G4LogicalVolume(SC3,         //its solid
                        SC_mat,          //its material
                        "SourceCasing3");           //its name

  new G4PVPlacement(0,                       //no rotation
    SCpos3,                    //at position
    logicSC3,               //its logical volume
    "SourceCasing3",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking
  */

 #ifdef USE_CASING

  // CLYC Detector Casing
  
  G4double innerRadiusDC1= 2.55*cm;
  G4double outerRadiusDC1 = 3*cm;
  G4double hzDC1 = 6.75*cm;
  G4double innerRadiusDC2 = 0*cm;
  G4double outerRadiusDC2= 3*cm;
  G4double hzDC2 = .1 * cm;
  G4double startAngleDC = 0.*deg;
  G4double spanningAngleDC = 360.*deg;
  G4ThreeVector CLYCDCpos1 = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, .5596578* base_world_unit + det_dist);
  G4ThreeVector CLYCDCpos2 = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, .4911578* base_world_unit + det_dist);

   G4Tubs* clycDC1
     = new G4Tubs("CLYCDetectorCasing1",
                  innerRadiusDC1,
                  outerRadiusDC1,
                  hzDC1,
                  startAngleDC,
                  spanningAngleDC);
    
    G4Tubs* clycDC2
     = new G4Tubs("CLYCDetectorCasing2",
                  innerRadiusDC2,
                  outerRadiusDC2,
                  hzDC2,
                  startAngleDC,
                  spanningAngleDC);

  // G4UnionSolid* solidShape3 
  //   = new G4UnionSolid("CLYCWindowHolder",;

  G4LogicalVolume* logicclycDC1 =                         
    new G4LogicalVolume(clycDC1,         //its solid
                        DC_mat,          //its material
                        "CLYCDetectorCasing1");           //its name

  new G4PVPlacement(0,                       //no rotation
    CLYCDCpos1,                    //at position
    logicclycDC1,               //its logical volume
    "CLYCDetectorCasing1",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking

  G4LogicalVolume* logicclycDC2 =                         
    new G4LogicalVolume(clycDC2,         //its solid
                        DC_mat,          //its material
                        "CLYCDetectorCasing2");           //its name

  new G4PVPlacement(0,                       //no rotation
    CLYCDCpos2,                    //at position
    logicclycDC2,               //its logical volume
    "CLYCDetectorCasing2",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking
  
  
  
  // CLYC Detector Housing 
  
  G4double innerRadiusDH1= 1.50622*cm;
  G4double outerRadiusDH1 = 1.5875*cm;
  G4double hzDH1 = 1.50622*cm;
  G4double innerRadiusDH2 = 0*cm;
  G4double outerRadiusDH2= 1.5875*cm;
  G4double hzDH2 = 0.04064*cm;
  G4double startAngleDH = 0.*deg;
  G4double spanningAngleDH = 360.*deg;
  G4ThreeVector CLYCDHpos1 = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, .5080328* base_world_unit + det_dist);
  G4ThreeVector CLYCDHpos2 = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, .4925642* base_world_unit + det_dist);


   G4Tubs* clycDH1
     = new G4Tubs("CLYCDetectorHousing1",
                  innerRadiusDH1,
                  outerRadiusDH1,
                  hzDH1,
                  startAngleDH,
                  spanningAngleDH);
    
    G4Tubs* clycDH2
     = new G4Tubs("CLYCDetectorHousing2",
                  innerRadiusDH2,
                  outerRadiusDH2,
                  hzDH2,
                  startAngleDH,
                  spanningAngleDH);

  // G4UnionSolid* solidShape3 
  //   = new G4UnionSolid("CLYCWindowHolder",;

  G4LogicalVolume* logicclycDH1 =                         
    new G4LogicalVolume(clycDH1,         //its solid
                        WinH_mat,          //its material
                        "CLYCDetectorHousing1");           //its name

  new G4PVPlacement(0,                       //no rotation
    CLYCDHpos1,                    //at position
    logicclycDH1,               //its logical volume
    "CLYCDetectorHousing1",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking

  G4LogicalVolume* logicclycDH2 =                         
    new G4LogicalVolume(clycDH2,         //its solid
                        WinH_mat,          //its material
                        "CLYCDetectorHousing2");           //its name

  new G4PVPlacement(0,                       //no rotation
    CLYCDHpos2,                    //at position
    logicclycDH2,               //its logical volume
    "CLYCDetectorHousing2",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking

  // CLYC Window Holder 

  G4double innerRadiusWH= 1.3081*cm;
  G4double outerRadiusWH = 1.47955*cm;
  G4double hzWH = 0.1905*cm;
  G4double innerRadiusWH2 = 1.50495*cm;
  G4double outerRadiusWH2= 1.6129*cm;
  G4double hzWH2 = 0.127*cm;
  G4double startAngleWH = 0.*deg;
  G4double spanningAngleWH = 360.*deg;
  G4ThreeVector CLYCWinHpos = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, .523095* base_world_unit + det_dist);
  G4ThreeVector CLYCWinHpos2 = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, 0.52627* base_world_unit + det_dist);



   G4Tubs* clycWinH
     = new G4Tubs("CLYCWindowHolder1",
                  innerRadiusWH,
                  outerRadiusWH,
                  hzWH,
                  startAngleWH,
                  spanningAngleWH);
    
    G4Tubs* clycWinH2
     = new G4Tubs("CLYCWindowHolder2",
                  innerRadiusWH2,
                  outerRadiusWH2,
                  hzWH2,
                  startAngleWH,
                  spanningAngleWH);

  // G4UnionSolid* solidShape3 
  //   = new G4UnionSolid("CLYCWindowHolder",;

  G4LogicalVolume* logicclycWinH =                         
    new G4LogicalVolume(clycWinH,         //its solid
                        WinH_mat,          //its material
                        "CLYCWindowHolder1");           //its name

  new G4PVPlacement(0,                       //no rotation
    CLYCWinHpos,                    //at position
    logicclycWinH,               //its logical volume
    "CLYCWindowHolder1",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking

  G4LogicalVolume* logicclycWinH2 =                         
    new G4LogicalVolume(clycWinH2,         //its solid
                        WinH_mat,          //its material
                        "CLYCWindowHolder2");           //its name

  new G4PVPlacement(0,                       //no rotation
    CLYCWinHpos2,                    //at position
    logicclycWinH2,               //its logical volume
    "CLYCWindowHolder2",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking

#endif

#ifdef ONEINCHCLYC
  // CLYC Detector physical volume

  G4double innerRadius = 0.*cm;
  G4double outerRadius = 1.25*cm;
  G4double hz = (2.5/2)*cm;
  G4double startAngle = 0.*deg;
  G4double spanningAngle = 360.*deg;
  // G4ThreeVector CLYCpos = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, 0.5125 * base_world_unit + det_dist);
  G4ThreeVector CLYCpos = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, det_dist + hz);



   G4Tubs* clycTube
     = new G4Tubs("CLYCcrystal",
                  innerRadius,
                  outerRadius,
                  hz,
                  startAngle,
                  spanningAngle);

  G4LogicalVolume* logicCLYC =                         
    new G4LogicalVolume(clycTube,         //its solid
                        CLYC_molecule,          //its material
                        "CLYCcrystal");           //its name

  fCLYCPV = new G4PVPlacement(0,                       //no rotation
    CLYCpos,                    //at position
    logicCLYC,               //its logical volume
    "CLYCcrystal",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking

#endif

#ifdef THREEINCHCLYC
  // CLYC Detector physical volume

  G4double innerRadius = 0.*cm;
  G4double outerRadius = 3.75*cm;
  G4double hz = (2.5/4)*cm;
  G4double startAngle = 0.*deg;
  G4double spanningAngle = 360.*deg;
  // G4ThreeVector CLYCpos = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, 0.5125 * base_world_unit + det_dist);
  G4ThreeVector CLYCpos = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit + y_height, det_dist + hz);



   G4Tubs* clycTube
     = new G4Tubs("CLYCcrystal",
                  innerRadius,
                  outerRadius,
                  hz,
                  startAngle,
                  spanningAngle);

  G4LogicalVolume* logicCLYC =                         
    new G4LogicalVolume(clycTube,         //its solid
                        C7LYC_molecule,          //its material
                        "CLYCcrystal");           //its name

  fCLYCPV = new G4PVPlacement(0,                       //no rotation
    CLYCpos,                    //at position
    logicCLYC,               //its logical volume
    "CLYCcrystal",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking

#endif

  // PE Collimator


  G4Material* collimatormat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");

  G4double col_phimin = 0.*deg, col_phimax = 360.*deg;

  G4double col_rmina, col_rmaxa, col_rminb, col_rmaxb, col_hz, col_rel_pos;

  // G4Cons* finecollimator = (G4Cons *) calloc(15, sizeof(G4Cons));;

  G4RotationMatrix rot1 = G4RotationMatrix(0,0,0);
  G4ThreeVector col_POS;
  G4Transform3D col_TR;

#ifdef USE_FTC
    G4MultiUnion* FTC = new G4MultiUnion("Fine_Collimator");
    G4ThreeVector finecollimatorpos = G4ThreeVector(0 * base_world_unit, 0* base_world_unit, 2 * base_world_unit - ftc_LEN[0]);




    // std::cout << "SIZE: " << 15*sizeof(G4Cons) <<std::endl;

    for(size_t i = 0; i < 15; i+=1)
    {
      col_rmina = ftc_ID1[i];
      col_rminb = ftc_ID2[i];
      // col_rmaxa = ftc_OD1[i];
      // col_rmaxb = ftc_OD2[i];
      col_rmaxa = ftc_OD1[i];
      col_rmaxb = ftc_OD2[i];
      col_hz = ftc_LEN[i];

      std::string newName("FTC_");
      newName+=std::to_string(15-i);

      G4Cons *newCons = new G4Cons(newName, 
      col_rmina, col_rmaxa, col_rminb, col_rmaxb, col_hz,
      col_phimin, col_phimax);


      col_POS = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit, ftc_POS[i]);
      col_TR = G4Transform3D(rot1,col_POS);

      FTC->AddNode(*newCons,col_TR);
      // std::cout << std::endl << finecollimator << ", " << (finecollimator + i)->GetRmin1() << ", " <<(finecollimator + i)->GetRmin2() << std::endl;
    }

    
    FTC->Voxelize();
    
    // G4cout << FTC->GetSurfaceArea() << G4endl;

    G4LogicalVolume* logicFTC=                         
      new G4LogicalVolume(FTC,         //its solid
                          collimatormat,          //its material
                          "Fine_Collimator");           //its name
                
    new G4PVPlacement(0,                       //no rotation
                      finecollimatorpos,                    //at position
                      logicFTC,             //its logical volume
                      "Fine_Collimator",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
#endif

#ifdef USE_MTC

    G4MultiUnion* MTC = new G4MultiUnion("Medium_Collimator");
    G4ThreeVector mediumcollimatorpos = G4ThreeVector(0 * base_world_unit, 0* base_world_unit, 2 * base_world_unit - mtc_LEN[0]);

    // std::cout << "SIZE: " << 15*sizeof(G4Cons) <<std::endl;

    for(size_t i = 0; i < 6; i+=1)
    {
      col_rmina = mtc_ID1[i];
      col_rminb = mtc_ID2[i];
      // col_rmaxa = mtc_OD1[i];
      // col_rmaxb = mtc_OD2[i];
      col_rmaxa = mtc_OD1[i];
      col_rmaxb = mtc_OD2[i];
      col_hz = mtc_LEN[i];

      std::string newName("MTC_");
      newName+=std::to_string(6-i);

      G4Cons *newCons = new G4Cons(newName, 
      col_rmina, col_rmaxa, col_rminb, col_rmaxb, col_hz,
      col_phimin, col_phimax);


      col_POS = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit, mtc_POS[i]);
      col_TR = G4Transform3D(rot1,col_POS);

      MTC->AddNode(*newCons,col_TR);
      // std::cout << std::endl << finecollimator << ", " << (finecollimator + i)->GetRmin1() << ", " <<(finecollimator + i)->GetRmin2() << std::endl;
    }

    
    MTC->Voxelize();
    
    // G4cout << FTC->GetSurfaceArea() << G4endl;

    G4LogicalVolume* logicMTC=                         
      new G4LogicalVolume(MTC,         //its solid
                          collimatormat,          //its material
                          "Medium_Collimator");           //its name
                
    new G4PVPlacement(0,                       //no rotation
                      mediumcollimatorpos,                    //at position
                      logicMTC,             //its logical volume
                      "Medium_Collimator",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

      // float *col = (float *) calloc(numberOfRows(), sizeof(float));
#endif

#ifdef USE_LTC
    G4MultiUnion* LTC = new G4MultiUnion("Large_Collimator");
    G4ThreeVector largecollimatorpos = G4ThreeVector(0 * base_world_unit, 0* base_world_unit, 2 * base_world_unit - ltc_LEN[0] + (mtc_LEN[0] + mtc_LEN[1] + mtc_LEN[2]));
    G4Material* largecollimatormat = nist->FindOrBuildMaterial("G4_NYLON-6-6");


    // std::cout << "SIZE: " << 15*sizeof(G4Cons) <<std::endl;

    for(size_t i = 0; i < 3; i+=1)
    {
      col_rmina = ltc_ID1[i];
      col_rminb = ltc_ID2[i];
      // col_rmaxa = ltc_OD1[i];
      // col_rmaxb = ltc_OD2[i];
      col_rmaxa = ltc_OD1[i];
      col_rmaxb = ltc_OD2[i];
      col_hz = ltc_LEN[i];

      std::string newName("LTC_");
      newName+=std::to_string(3-i);

      G4Cons *newCons = new G4Cons(newName, 
      col_rmina, col_rmaxa, col_rminb, col_rmaxb, col_hz,
      col_phimin, col_phimax);


      col_POS = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit, ltc_POS[i]);
      col_TR = G4Transform3D(rot1,col_POS);

      LTC->AddNode(*newCons,col_TR);
      // std::cout << std::endl << finecollimator << ", " << (finecollimator + i)->GetRmin1() << ", " <<(finecollimator + i)->GetRmin2() << std::endl;
    }

    
    LTC->Voxelize();
    
    // G4cout << FTC->GetSurfaceArea() << G4endl;

    G4LogicalVolume* logicLTC=                         
      new G4LogicalVolume(LTC,         //its solid
                          largecollimatormat,          //its material
                          "Large_Collimator");           //its name
                
    new G4PVPlacement(0,                       //no rotation
                      largecollimatorpos,                    //at position
                      logicLTC,             //its logical volume
                      "Large_Collimator",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

      // float *col = (float *) calloc(numberOfRows(), sizeof(float));

#endif
  /*
  G4Material* shape1_collimator_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4ThreeVector pos1_collimator = G4ThreeVector(0 * base_world_unit, 0* base_world_unit, 2.1664 * base_world_unit);




  // Conical section shape - front PE collimator section       
  G4double shape1_fc_rmina =  2.49*cm, shape1_fc_rmaxa = 9.4*cm;
  G4double shape1_fc_rminb =  1.06*cm, shape1_fc_rmaxb = 9.4*cm;
  G4double shape1_fc_hz = 16.64*cm;
  G4double shape1_fc_phimin = 0.*deg, shape1_fc_phimax = 360.*deg;
  G4Cons* solidShape1_front_collimator =    
    new G4Cons("Shape1_Front_Collimator", 
    shape1_fc_rmina, shape1_fc_rmaxa, shape1_fc_rminb, shape1_fc_rmaxb, shape1_fc_hz,
    shape1_fc_phimin, shape1_fc_phimax);

  // Conical section shape - rear PE collimator section       
  G4double shape1_rc_rmina =  .32*cm, shape1_rc_rmaxa = 9.4*cm;
  G4double shape1_rc_rminb =  1.25*cm, shape1_rc_rmaxb = 9.4*cm;
  G4double shape1_rc_hz = 64.58*cm;
  G4double shape1_rc_phimin = 0.*deg, shape1_rc_phimax = 360.*deg;
  G4Cons* solidShape1_rear_collimator =    
    new G4Cons("Shape1_Rear_Collimator", 
    shape1_rc_rmina, shape1_rc_rmaxa, shape1_rc_rminb, shape1_rc_rmaxb, shape1_rc_hz,
    shape1_rc_phimin, shape1_rc_phimax);

  G4MultiUnion* solidShape1_collimator
  = new G4MultiUnion("Shape1_Collimator");

  G4RotationMatrix rot1 = G4RotationMatrix(0,0,0);
  G4ThreeVector pos1_fc = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit, 0 * base_world_unit);
  G4ThreeVector pos1_rc = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit, .8121 * base_world_unit);

  G4Transform3D p1_fc_tr = G4Transform3D(rot1,pos1_fc);
  G4Transform3D p1_rc_tr = G4Transform3D(rot1,pos1_rc);
  
  solidShape1_collimator->AddNode(*solidShape1_front_collimator,p1_fc_tr);
  solidShape1_collimator->AddNode(*solidShape1_rear_collimator,p1_rc_tr);

  solidShape1_collimator->Voxelize();

  G4LogicalVolume* logicShape1_col =                         
    new G4LogicalVolume(solidShape1_collimator,         //its solid
                        shape1_collimator_mat,          //its material
                        "Shape1_Collimator");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos1_collimator,                    //at position
                    logicShape1_col,             //its logical volume
                    "Shape1_Collimator",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  */
  // G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  // G4ThreeVector pos2 = G4ThreeVector(0 * base_world_unit, 0* base_world_unit + y_height, 2.9786 * m);
        
  // // Conical section shape       

                      
  // G4LogicalVolume* logicShape2 =                         
  //   new G4LogicalVolume(solidShape2,         //its solid
  //                       shape2_mat,          //its material
  //                       "Shape2");           //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   pos2,                    //at position
  //                   logicShape2,             //its logical volume
  //                   "Shape2",                //its name
  //                   logicWorld,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

#ifdef USE_TUNNEL
    //     
  // Shape 1 -  Collimator and Housing
  //  

  // Aluminum collimator housing

  G4Material* shape1_alh_mat = nist->FindOrBuildMaterial("G4_Al");
  G4ThreeVector pos1_alh = G4ThreeVector(0 * base_world_unit, 0* base_world_unit + y_height, 2.8121* base_world_unit);


  G4double shape1_alh_rmin =  0*cm, shape1_alh_rmax = 16.17*cm;
  G4double shape1_alhs_rmin =  0*cm, shape1_alhs_rmax = 14.9*cm;

  G4Tubs* solidShape1_alhi    
    = new G4Tubs("Shape1_Al_Housing_inner",
                shape1_alhs_rmin,
                shape1_alhs_rmax,
                120*cm,
                startAngle,
                spanningAngle);

  G4Tubs* solidShape1_alho    
    = new G4Tubs("Shape1_Al_Housing_outer",
                shape1_alh_rmin,
                shape1_alh_rmax,
                100*cm,
                startAngle,
                spanningAngle);

  G4SubtractionSolid* solidShape1_alh 
    = new G4SubtractionSolid("Shape4_tub",
    solidShape1_alho,
    solidShape1_alhi);

        
  G4LogicalVolume* logicShape1_alh =                         
    new G4LogicalVolume(solidShape1_alh,         //its solid
                        shape1_alh_mat,          //its material
                        "Shape1_ALH");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos1_alh,                    //at position
                    logicShape1_alh,             //its logical volume
                    "Shape1_ALH",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  
  G4Material* shape3_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
        
  // Conical section shape       
  G4Box* solidShape3_box =    
    new G4Box("Shape3_box",
    2 * base_world_unit, 2 * base_world_unit, .8121 * base_world_unit);

  std::cout << solidShape3_box->GetName() << std::endl;

  G4Tubs* solidShape3_cylinder
    = new G4Tubs("Shape3_cylinder",
                0*cm,
                16.17*cm,
                90*cm,
                startAngle,
                spanningAngle);

  G4SubtractionSolid* solidShape3_collimator 
    = new G4SubtractionSolid("Shape3_collimator",
    solidShape3_box,
    solidShape3_cylinder);

  // G4LogicalVolume* logicShape3 =                         
  //   new G4LogicalVolume(solidShape3,         //its solid
  //                       shape3_mat,          //its material
  //                       "Shape3");           //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   pos3,                    //at position
  //                   logicShape3,             //its logical volume
  //                   "Shape3",                //its name
  //                   logicWorld,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

  G4Material* shape3a_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
  G4ThreeVector pos3a = G4ThreeVector(0 * base_world_unit, 0* base_world_unit, 30.8121* base_world_unit);
        
  // Conical section shape       
  G4Box* solidShape3_rearwall =    
    new G4Box("Shape3_rearwall",
    2 * base_world_unit, 2 * base_world_unit, .8121 * base_world_unit);

  // Conical section shape       
  G4Box* solidShape3_floor =    
    new G4Box("Shape3_floor",
    2 * base_world_unit, 0.25 * base_world_unit, 14.8121 * base_world_unit);

  G4MultiUnion* solidShape3
  = new G4MultiUnion("Shape3");

  G4RotationMatrix rot3 = G4RotationMatrix(0,0,0);
  G4ThreeVector pos3_floor = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit, 0 * base_world_unit);
  G4ThreeVector pos3_collimator = G4ThreeVector(0 * base_world_unit, 1.75* base_world_unit, -14 * base_world_unit);
  G4ThreeVector pos3_rearwall = G4ThreeVector(0 * base_world_unit, 1.75* base_world_unit, 14 * base_world_unit);

  G4Transform3D p3_floortr = G4Transform3D(rot3,pos3_floor);
  G4Transform3D p3_collimatortr = G4Transform3D(rot3,pos3_collimator);
  G4Transform3D p3_rearwalltr = G4Transform3D(rot3,pos3_rearwall);


  solidShape3->AddNode(*solidShape3_floor,p3_floortr);
  solidShape3->AddNode(*solidShape3_collimator,p3_collimatortr);
  solidShape3->AddNode(*solidShape3_rearwall,p3_rearwalltr);


  solidShape3->Voxelize();

  G4ThreeVector pos3 = G4ThreeVector(0 * base_world_unit, -1.75* base_world_unit, 16.8121* base_world_unit);


  G4LogicalVolume* logicShape3 =                         
    new G4LogicalVolume(solidShape3,         //its solid
                        shape3_mat,          //its material
                        "Shape3");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos3,                    //at position
                    logicShape3,             //its logical volume
                    "Shape3",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  // G4LogicalVolume* logicShape3a =                         
  //   new G4LogicalVolume(solidShape3a,         //its solid
  //                       shape3a_mat,          //its material
  //                       "Shape3a");           //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   pos3a,                    //at position
  //                   logicShape3a,             //its logical volume
  //                   "Shape3a",                //its name
  //                   logicWorld,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

  G4Material* shape4_mat = nist->FindOrBuildMaterial("G4_CONCRETE");

  G4ThreeVector pos4 = G4ThreeVector(0 * base_world_unit, 0* base_world_unit, 16.8121 * base_world_unit);


  G4Box* solidShape4_box =    
    new G4Box("Shape4_box",
    1.5 * base_world_unit, 1.5 * base_world_unit, 13.1879 * base_world_unit);

  G4Tubs* solidShape4_cylinder
    = new G4Tubs("Shape4_cylinder",
                0*cm,
                1.2*m,
                13.1879* m,
                startAngle,
                spanningAngle);

  G4SubtractionSolid* solidShape4_tub 
    = new G4SubtractionSolid("Shape4_tub",
    solidShape4_box,
    solidShape4_cylinder);

  // G4LogicalVolume* logicShape4_tub =                         
  //   new G4LogicalVolume(solidShape4_tub,         //its solid
  //                       shape4_mat,          //its material
  //                       "Shape4_tub");           //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   pos4,                    //at position
  //                   logicShape4_tub,             //its logical volume
  //                   "Shape4_tub",                //its name
  //                   logicWorld,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

  G4Box* solidShape4_floor =    
  new G4Box("Shape4_floor",
  1.5 * base_world_unit, 0.25 * base_world_unit, 13.1879 * base_world_unit);

  // G4LogicalVolume* logicShape4_floor =                         
  //   new G4LogicalVolume(solidShape4_floor,         //its solid
  //                       shape4_mat,          //its material
  //                       "Shape4_floor");           //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   pos4_floor,                    //at position
  //                   logicShape4_floor,             //its logical volume
  //                   "Shape4floor",                //its name
  //                   logicShape4_tub,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

  G4MultiUnion* solidShape4
  = new G4MultiUnion("Shape4Tunnel");

  G4RotationMatrix rot4 = G4RotationMatrix(0,0,0);
  G4ThreeVector pos4_tub = G4ThreeVector(0 * base_world_unit, 0* base_world_unit, 0 * base_world_unit);
  G4ThreeVector pos4_floor = G4ThreeVector(0 * base_world_unit, -1.25 * base_world_unit, 0 * base_world_unit);

  G4Transform3D tubetr = G4Transform3D(rot4,pos4_tub);
  G4Transform3D floortr = G4Transform3D(rot4,pos4_floor);

  solidShape4->AddNode(*solidShape4_tub,tubetr);
  solidShape4->AddNode(*solidShape4_floor,floortr);

  solidShape4->Voxelize();

  // G4UnionSolid* solidShape4
  //  = new G4UnionSolid("Shape4",
  //  solidShape4_tub,
  //  solidShape4_floor);

  G4LogicalVolume* logicShape4 =                         
    new G4LogicalVolume(solidShape4,         //its solid
                        shape4_mat,          //its material
                        "Shape4Tunnel");           //its name

  new G4PVPlacement(0,                       //no rotation
                    pos4,                    //at position
                    logicShape4,             //its logical volume
                    "Shape4Tunnel",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  
  // G4Material* shape5_mat = nist->FindOrBuildMaterial("G4_Al");
  // G4ThreeVector pos5 = G4ThreeVector(0 * base_world_unit, 0* base_world_unit + y_height, 2.8121* base_world_unit);
        
  // // Conical section shape       
  // G4double shape5_rmin =  9.4*cm, shape5_rmax = 10.67*cm;
  // G4double shape5_hz = 20.53*cm;
  // G4Tubs* solidShape5     
  //   = new G4Tubs("Shape3_cylinder",
  //               shape5_rmin,
  //               shape5_rmax,
  //               86.29*cm,
  //               startAngle,
  //               spanningAngle);
                      
  // G4LogicalVolume* logicShape5 =                         
  //   new G4LogicalVolume(solidShape5,         //its solid
  //                       shape5_mat,          //its material
  //                       "Shape5");           //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   pos5,                    //at position
  //                   logicShape5,             //its logical volume
  //                   "Shape5",                //its name
  //                   logicWorld,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking
#endif

#ifdef USE_SHIELD1

    G4Material* shape6_mat = nist->FindOrBuildMaterial("G4_PARAFFIN");
    G4ThreeVector pos6 = G4ThreeVector(0 * base_world_unit, -.3015* base_world_unit, 11 * base_world_unit + 16.66875 * cm);
          
    // Conical section shape       
    G4Box* solidShape6_1 =    
      new G4Box("Shape6_box1",
      48.5775 * cm, 69.85 * cm, 16.66875 * cm);

    // G4Box* solidShape6_2 =    
    //   new G4Box("Shape6_box2",
    //   500 * mm, 200 * mm, 120 * mm);

    // G4Tubs* solidShape3_cylinder
    //   = new G4Tubs("Shape3_cylinder",
    //               0*cm,
    //               10.67*cm,
    //               90*cm,
    //               startAngle,
    //               spanningAngle);

    // G4SubtractionSolid* solidShape6 
    //   = new G4SubtractionSolid("Shape6",
    //   solidShape6_1,
    //   solidShape6_2);

    G4LogicalVolume* logicShape6 =                         
      new G4LogicalVolume(solidShape6_1,         //its solid
                          shape6_mat,          //its material
                          "Shape6");           //its name
                
    new G4PVPlacement(0,                       //no rotation
                      pos6,                    //at position
                      logicShape6,             //its logical volume
                      "Shape6",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

#endif


#ifdef USE_SHIELD2

    G4double mcrate[3] = {24.2888 * cm, 16.6687 * cm,13.97 * cm};

    double dist = 15;

    G4Material* shape6_mat = nist->FindOrBuildMaterial("G4_PARAFFIN");
    G4ThreeVector pos6a = G4ThreeVector(0 * base_world_unit, -.3015* base_world_unit, dist * base_world_unit + mcrate[2]);
    G4ThreeVector pos6b = G4ThreeVector(-2 * mcrate[0] + 14.173 * cm, -.366626* base_world_unit + 3.3374 * cm, dist* base_world_unit);
    G4ThreeVector pos6c = G4ThreeVector(2 * mcrate[0] - 14.173 * cm, -.366626* base_world_unit  + 3.3374 * cm, dist * base_world_unit);


    G4RotationMatrix* rot6b = new G4RotationMatrix;
    G4RotationMatrix* rot6c = new G4RotationMatrix;

    rot6b->rotateY(-4*M_PI/9*rad);    
    rot6c->rotateY(4*M_PI/9*rad);


    

    

    // Conical section shape       
    G4Box* solidShape6_1 =    
      new G4Box("ShieldPlate",
      mcrate[0] * 2, mcrate[1] * 4, mcrate[2]);

      // 48.5775 * cm, 69.85 * cm, 16.66875 * cm);

    G4Box* solidShape6_2 =    
      new G4Box("ShieldWing",
      mcrate[0] * 2, mcrate[1] * 4, mcrate[2]);
    // G4Tubs* solidShape3_cylinder
    //   = new G4Tubs("Shape3_cylinder",
    //               0*cm,
    //               10.67*cm,
    //               90*cm,
    //               startAngle,
    //               spanningAngle);

    // G4SubtractionSolid* solidShape6 
    //   = new G4SubtractionSolid("Shape6",
    //   solidShape6_1,
    //   solidShape6_2);

    // G4LogicalVolume* logicShape6a=                         
    //   new G4LogicalVolume(solidShape6_1,         //its solid
    //                       shape6_mat,          //its material
    //                       "ShieldPlate");           //its name
                
    // new G4PVPlacement(0,                       //no rotation
    //                   pos6a,                    //at position
    //                   logicShape6a,             //its logical volume
    //                   "ShieldPlate",                //its name
    //                   logicWorld,                //its mother  volume
    //                   false,                   //no boolean operation
    //                   0,                       //copy number
    //                   checkOverlaps);          //overlaps checking

    G4LogicalVolume* logicShape6b =                         
      new G4LogicalVolume(solidShape6_2,         //its solid
                          shape6_mat,          //its material
                          "ShieldWingR");           //its name
                
    new G4PVPlacement(rot6b,                       //no rotation
                      pos6b,                    //at position
                      logicShape6b,             //its logical volume
                      "ShieldWingR",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

    G4LogicalVolume* logicShape6c =                         
      new G4LogicalVolume(solidShape6_2,         //its solid
                          shape6_mat,          //its material
                          "ShieldWingL");           //its name
                
    new G4PVPlacement(rot6c,                       //no rotation
                      pos6c,                    //at position
                      logicShape6c,             //its logical volume
                      "ShieldWingL",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

#endif

  /* 
  *******************************************************************************
  *******************************************************************************
  * 
  * This section contains the colors used separate parts for visualization
  * as well as any needed custom materials for the structure
  * 
  *******************************************************************************
  *******************************************************************************
  */

  G4Colour steel(0.44313725,0.4745098,0.49411765);
  G4Colour concrete(0.67058824,0.62352941,0.56862745);
  G4Colour concrete_structure(0.27058824,0.25098039,0.23137255);
  G4Colour borpoly30(0.4, 0.2, 0.6);
  G4Colour lead(0.18039216,0.18039216,0.06666667);
  G4Colour polyeth(0.85882353,0.85882353,0.80784314);
  G4Colour aluminum(0.752941176,0.752941176,0.752941176);
  G4Colour nylon(0.85882353,0.85882353,0.6);
  G4Colour generic(0.20392157,0.92156863,0.92156863);


  G4VisAttributes* steelVisAttributes = new G4VisAttributes(steel);
  G4VisAttributes* concreteVisAttributes = new G4VisAttributes(concrete);
  G4VisAttributes* concrete_structureVisAttributes = new G4VisAttributes(concrete_structure);
  G4VisAttributes* borpoly30VisAttributes = new G4VisAttributes(borpoly30);
  G4VisAttributes* leadVisAttributes = new G4VisAttributes(lead);
  G4VisAttributes* polyethVisAttributes = new G4VisAttributes(polyeth);
  G4VisAttributes* aluminumVisAttributes = new G4VisAttributes(aluminum);
  G4VisAttributes* nylonVisAttributes = new G4VisAttributes(nylon);
  G4VisAttributes* genericVisAttributes = new G4VisAttributes(generic);




  /*
  *
  * New materials
  * 
  */

  // 30% Borated Polyethelyne. Density value from: https://marshield.com/borated-polyethylene-neutron-shielding
  G4Element* BoratedPoly30_B = nist->FindOrBuildElement(5);
  G4Material* BoratedPoly30_PolyEth = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4double BoratedPoly30_B_density = 0.918 * g/cm3;
  G4Material* BoratedPoly30 = new G4Material(name="BorPoly30",density,ncomponents=2);
  BoratedPoly30->AddElement(BoratedPoly30_B, 30*perCent);
  BoratedPoly30->AddMaterial(BoratedPoly30_PolyEth, 70*perCent);

  // Concrete block (less dense than regular concrete ~2.1 g/cm^3)
  G4Material* ConcreteBlock_concrete = nist->FindOrBuildMaterial("G4_CONCRETE");
  G4Material* ConcreteBlock_air = nist->FindOrBuildMaterial("G4_AIR");
  G4double ConcreteBlock_density = 2.17909050602022 * g/cm3;
  G4Material* ConcreteBlock = new G4Material(name="ConcreteBlock",density,ncomponents=2);
  ConcreteBlock->AddMaterial(ConcreteBlock_air, (100-94.74306547914)*perCent);
  ConcreteBlock->AddMaterial(ConcreteBlock_concrete, 94.74306547914*perCent);



  /* 
  *******************************************************************************
  *******************************************************************************
  * 
  * This section contains the physical structure of the tunnel and tunnel
  * collimator areas. Does not contain any of the source room structures.
  * 
  *******************************************************************************
  *******************************************************************************
  */

  const G4double ConcreteYAdj = .5 * cm;


  // Tunnel Structure (such as walls, floor)
  {  
    auto TunnelStructureMesh = CADMesh::TessellatedMesh::FromSTL("CADModels/tunnel_structure.stl");
    TunnelStructureMesh->SetScale(1000);
    G4Material* TunnelStructureMat = nist->FindOrBuildMaterial("G4_CONCRETE");

    // TunnelStructureMesh->SetMaterial(TunnelStructureMat);
    G4VSolid* TunnelStructureMeshSolid = TunnelStructureMesh->GetTessellatedSolid();

    auto TunnelStructureMesh_pos = G4ThreeVector(-106.47*cm, (-106.1*cm + ConcreteYAdj), 2.3*m);
    auto TunnelStructureMesh_rot = new G4RotationMatrix;


    G4LogicalVolume* logicTunnelStructureMesh =                         
      new G4LogicalVolume(TunnelStructureMeshSolid,         //its solids
                          TunnelStructureMat,          //its material
                          "TunnelStructureMesh");           //its name

    logicTunnelStructureMesh->SetVisAttributes(concreteVisAttributes);

    new G4PVPlacement(TunnelStructureMesh_rot,                       //no rotation
                      TunnelStructureMesh_pos,                    //at position
                      logicTunnelStructureMesh,             //its logical volume
                      "TunnelStructureMesh",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // Fixed Front Concrete Collimator Assembly
  {
    auto FixedFrontConcreteMesh = CADMesh::TessellatedMesh::FromSTL("CADModels/fixed_front_concrete.stl");
    FixedFrontConcreteMesh->SetScale(1000);
    G4Material* FixedFrontConcreteMat = nist->FindOrBuildMaterial("G4_CONCRETE");

    // TunnelStructureMesh->SetMaterial(TunnelStructureMat);
    G4VSolid* FixedFrontConcreteMeshSolid = FixedFrontConcreteMesh->GetTessellatedSolid();

    auto FixedFrontConcreteMesh_pos = G4ThreeVector(-106.47*cm, (-106.1*cm + ConcreteYAdj), 56*cm);
    auto FixedFrontConcreteMesh_rot = new G4RotationMatrix;


    G4LogicalVolume* logicFixedFrontConcreteMeshMesh =                         
      new G4LogicalVolume(FixedFrontConcreteMeshSolid,         //its solids
                          FixedFrontConcreteMat,          //its material
                          "FixedFrontConcreteCollimator");           //its name

    logicFixedFrontConcreteMeshMesh->SetVisAttributes(concreteVisAttributes);

                
    new G4PVPlacement(FixedFrontConcreteMesh_rot,                       //no rotation
                      FixedFrontConcreteMesh_pos,                    //at position
                      logicFixedFrontConcreteMeshMesh,             //its logical volume
                      "FixedFrontConcreteCollimator",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // TOF Tunnel Entrance Door
  {
    auto TOFTunnelDoorMesh = CADMesh::TessellatedMesh::FromSTL("CADModels/sliding_door.stl");
    TOFTunnelDoorMesh->SetScale(1000);
    G4Material* TOFTunnelDoorMat = nist->FindOrBuildMaterial("G4_CONCRETE");

    // TunnelStructureMesh->SetMaterial(TunnelStructureMat);
    G4VSolid* TOFTunnelDoorMeshSolid = TOFTunnelDoorMesh->GetTessellatedSolid();

    auto TOFTunnelDoorMesh_pos = G4ThreeVector(1*cm, (-106.1*cm + ConcreteYAdj), 228.5*cm);
    auto TOFTunnelDoorMesh_rot = new G4RotationMatrix;


    G4LogicalVolume* logicTOFTunnelDoorMeshMesh =                         
      new G4LogicalVolume(TOFTunnelDoorMeshSolid,         //its solids
                          TOFTunnelDoorMat,          //its material
                          "TOFTunnelDoor");           //its name

    logicTOFTunnelDoorMeshMesh->SetVisAttributes(concreteVisAttributes);
                
    new G4PVPlacement(TOFTunnelDoorMesh_rot,                       //no rotation
                      TOFTunnelDoorMesh_pos,                    //at position
                      logicTOFTunnelDoorMeshMesh,             //its logical volume
                      "TOFTunnelDoor",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // Concrete bricks surrounding front of tunnel collimator
  {
    G4Material* TOFCollimatorBlocksMat = ConcreteBlock;

        // Slab is (x, y, z) = (61 cm, 61 cm, 1.27 cm)
    G4Box* TOFCollimatorBlocks_main =    
      new G4Box("TOFCollimatorBlocks_main",
      137.825/2 * cm, 244/2 * cm, 39.6875/2 * cm);

    auto TOFCollimatorBlocks_collimator_cutout_pos = G4ThreeVector(-12.7/2*cm, (-15.9*cm-ConcreteYAdj), 0*cm);


    G4Box* TOFCollimatorBlocks_collimator_cutout =    
      new G4Box("TOFCollimatorBlocks_collimator_cutout",
      30.5 * cm, 30.5 * cm, 50./2 * cm);

    G4SubtractionSolid* TOFCollimatorBlocks 
      = new G4SubtractionSolid("TOFCollimatorBlocks",
      TOFCollimatorBlocks_main,
      TOFCollimatorBlocks_collimator_cutout,0,TOFCollimatorBlocks_collimator_cutout_pos);


    auto TOFCollimatorBlocks_pos = G4ThreeVector(12.7/2*cm, (15.9 * cm + ConcreteYAdj), (256*cm-39.6875/2 * cm));
    auto TOFCollimatorBlocks_rot = new G4RotationMatrix;


    G4LogicalVolume* logicTOFCollimatorBlocks =                         
      new G4LogicalVolume(TOFCollimatorBlocks,         //its solids
                          TOFCollimatorBlocksMat,          //its material
                          "TOFCollimatorBlocks");           //its name

    logicTOFCollimatorBlocks->SetVisAttributes(concreteVisAttributes);
                
    new G4PVPlacement(TOFCollimatorBlocks_rot,                       //no rotation
                      TOFCollimatorBlocks_pos,                    //at position
                      logicTOFCollimatorBlocks,             //its logical volume
                      "TOFCollimatorBlocks",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }
  
  // Concrete Wall at end of tunnel
  {
    G4Material* ConcreteTunnelRearWall_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
    G4ThreeVector ConcreteTunnelRearWall_pos = G4ThreeVector(0*cm, 0*cm, (3092+5) * cm + (9.98795+9.41205)*mm);
          
    // Conical section shape       
    G4Box* ConcreteTunnelRearWall =    
      new G4Box("ConcreteTunnelRearWall",
      122*cm, 122 * cm, 5*cm);

    G4LogicalVolume* logicConcreteTunnelRearWall=                         
      new G4LogicalVolume(ConcreteTunnelRearWall,         //its solid
                          ConcreteTunnelRearWall_mat,          //its material
                          "ConcreteTunnelRearWall");           //its name

    logicConcreteTunnelRearWall->SetVisAttributes(concrete_structureVisAttributes);

    new G4PVPlacement(0,                       //no rotation
                      ConcreteTunnelRearWall_pos,                    //at position
                      logicConcreteTunnelRearWall,             //its logical volume
                      "ConcreteTunnelRearWall",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

// 
  /* 
  *******************************************************************************
  *******************************************************************************
  * 
  * This section contains the tunnel collimator pieces, including iron/borpoly
  * shield, snouts, and poly/nylon collimator sections.
  * 
  *******************************************************************************
  *******************************************************************************
  */

  // Concrete collimator hole steel liner
  {
    G4Material* CollimatorHoleSteelLiner_mat = nist->FindOrBuildMaterial("G4_Fe");
    G4ThreeVector CollimatorHoleSteelLiner_pos = G4ThreeVector(0 * cm, 0 * cm, (256. + 124./2- (2.45*11)/2 ) * cm);
          
    // Conical section shape       
      G4Cons *CollimatorHoleSteelLiner = new G4Cons("CollimatorHoleSteelLiner", 
      29.8/2 * cm, 33/2 * cm, 
      29.8/2 * cm, 33/2 * cm,
      151.94/2 * cm, 0 * deg, 360 * deg);

    G4LogicalVolume* logicCollimatorHoleSteelLiner=                         
      new G4LogicalVolume(CollimatorHoleSteelLiner,         //its solid
                          CollimatorHoleSteelLiner_mat,          //its material
                          "CollimatorHoleSteelLiner");           //its name

    logicCollimatorHoleSteelLiner->SetVisAttributes(steelVisAttributes);
                
    new G4PVPlacement(0,                       //no rotation
                      CollimatorHoleSteelLiner_pos,                    //at position
                      logicCollimatorHoleSteelLiner,             //its logical volume
                      "CollimatorHoleSteelLiner",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // Alternating BorPoly30 and Steel slabs (6 steel, 5 BorPoly30, each slab 1" thick)
  {
    G4Material* PolyIronSlabCollimator_Fe_mat = nist->FindOrBuildMaterial("G4_Fe");
    G4Material* PolyIronSlabCollimator_mat;

    // Slab is (x, y, z) = (61 cm, 61 cm, 1.27 cm)
    G4Box* PolyIronSlab_rect =    
      new G4Box("PolyIronSlab_rect",
      30.5 * cm, 30.5 * cm, 1.27 * cm);

    // Hole in slab is 33 cm in diameter
    G4Tubs* PolyIronSlab_hole
      = new G4Tubs("PolyIronSlab_hole",
                  0 * cm,
                  33./2 * cm,
                  1.3 * cm,
                  0 * deg,
                  360 * deg);

    G4SubtractionSolid* PolyIronSlabCollimator 
      = new G4SubtractionSolid("PolyIronSlab",
      PolyIronSlab_rect,
      PolyIronSlab_hole);

    for(auto i = 0; i < 11; ++i)
    {

      G4ThreeVector PolyIronSlabCollimator_pos = G4ThreeVector(0 * cm, 0 * cm, (256 - (i * 2.54) - 1.27) * cm);

      std::string CollimatorName = "PolyIronSlabCollimator_" + std::to_string(i);


      if((i+1)%2 == 0)
      {
        PolyIronSlabCollimator_mat = BoratedPoly30;
      }
      else
      {
        PolyIronSlabCollimator_mat = PolyIronSlabCollimator_Fe_mat;
      }


      G4LogicalVolume* logicPolyIronSlabCollimator=                         
        new G4LogicalVolume(PolyIronSlabCollimator,         //its solid
                            PolyIronSlabCollimator_mat,          //its material
                            CollimatorName);           //its name

      if((i+1)%2 == 0)
      {
        logicPolyIronSlabCollimator->SetVisAttributes(borpoly30VisAttributes);      
      }
      else
      {
        logicPolyIronSlabCollimator->SetVisAttributes(steelVisAttributes);
      }

                  
      new G4PVPlacement(0,                       //no rotation
                        PolyIronSlabCollimator_pos,                    //at position
                        logicPolyIronSlabCollimator,             //its logical volume
                        CollimatorName,                //its name
                        logicWorld,                //its mother  volume
                        false,                   //no boolean operation
                        0,                       //copy number
                        checkOverlaps);          //overlaps checking
    }
  }

  // Alternating Poly and Lead slabs (2 lead (4 cm thick), 1 Poly (2 cm thick))
  {
    G4Material* PolyLeadSlabCollimator_Pb_mat = nist->FindOrBuildMaterial("G4_Pb");
    G4Material* PolyLeadSlabCollimator_PE_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");

    G4Material* PolyLeadSlabCollimator_mat;



    // Lead insert is 4 cm long. PE insert is 2 cm long
    G4Tubs* PolyLeadSlab_Pb
      = new G4Tubs("PolyLeadSlab_Pb",
                  20.5/2 * cm,
                  29.8/2 * cm,
                  2 * cm,
                  0 * deg,
                  360 * deg);

    G4Tubs* PolyLeadSlab_PE
      = new G4Tubs("PolyLeadSlab_PE",
                  20.5/2 * cm,
                  29.8/2 * cm,
                  1 * cm,
                  0 * deg,
                  360 * deg);


    std::vector<G4double> PbPolyInsert_rel_pos = {0 * cm, 3 * cm, 6 * cm};

    for(auto i = 0; i < 3; ++i)
    {

      G4ThreeVector PolyLeadSlabCollimator_pos = G4ThreeVector(0 * cm, 0 * cm, (254) * cm + PbPolyInsert_rel_pos.at(i));
      // G4ThreeVector PolyLeadSlabCollimator_pos = G4ThreeVector(0 * cm, 0 * cm, (-42)* cm + PbPolyInsert_rel_pos.at(i));


      std::string CollimatorName = "PolyLeadSlabCollimator_" + std::to_string(i);


      G4LogicalVolume* logicPolyLeadSlabCollimator;

      if((i+1)%2 == 0)
      {
        G4LogicalVolume* logicPolyLeadSlabCollimator_PE=                         
        new G4LogicalVolume(PolyLeadSlab_PE,         //its solid
                            PolyLeadSlabCollimator_PE_mat,          //its material
                            CollimatorName);           //its name
        logicPolyLeadSlabCollimator = logicPolyLeadSlabCollimator_PE;

      }
      else
      {
        G4LogicalVolume* logicPolyLeadSlabCollimator_Pb=                         
        new G4LogicalVolume(PolyLeadSlab_Pb,         //its solid
                            PolyLeadSlabCollimator_Pb_mat,          //its material
                            CollimatorName);           //its name
        logicPolyLeadSlabCollimator = logicPolyLeadSlabCollimator_Pb;

      }


      // G4LogicalVolume* logicPolyLeadSlabCollimator=                         
      //   new G4LogicalVolume(PolyLeadSlab,         //its solid
      //                       PolyLeadSlabCollimator_mat,          //its material
      //                       CollimatorName);           //its name

      if((i+1)%2 == 0)
      {
        logicPolyLeadSlabCollimator->SetVisAttributes(polyethVisAttributes);      
      }
      else
      {
        logicPolyLeadSlabCollimator->SetVisAttributes(leadVisAttributes);
      }

                  
      new G4PVPlacement(0,                       //no rotation
                        PolyLeadSlabCollimator_pos,                    //at position
                        logicPolyLeadSlabCollimator,             //its logical volume
                        CollimatorName,                //its name
                        logicWorld,                //its mother  volume
                        false,                   //no boolean operation
                        0,                       //copy number
                        checkOverlaps);          //overlaps checking
    }
  }

  // Aluminum snout
  {
    G4Material* AluminumSnout_mat = nist->FindOrBuildMaterial("G4_Al");
    G4ThreeVector AluminumSnout_pos = G4ThreeVector(0 * cm, 0 * cm, (228.06 - 41./2 + 5) * cm);

    G4ThreeVector AluminumSnout_insert_pos = G4ThreeVector(0 * cm, 0 * cm, (41./2) * cm);

          
    // Conical section shape       
    G4Tubs *AluminumSnout_insert = new G4Tubs("AluminumSnout_insert",
                18.8/2 * cm,
                29.8/2 * cm,
                5 * cm,
                0 * deg,
                360 * deg);



    G4Tubs *AluminumSnout_protrusion = new G4Tubs("AluminumSnout_protrusion",
            18.8/2 * cm,
            19.8/2 * cm,
            41./2 * cm,
            0 * deg,
            360 * deg);
    
    G4UnionSolid* AluminumSnout = new G4UnionSolid("AluminumSnout",
    AluminumSnout_protrusion, AluminumSnout_insert, 0, AluminumSnout_insert_pos);


    G4LogicalVolume* logicAluminumSnout=                         
      new G4LogicalVolume(AluminumSnout,         //its solid
                          AluminumSnout_mat,          //its material
                          "AluminumSnout");           //its name

    logicAluminumSnout->SetVisAttributes(aluminumVisAttributes);
                
    new G4PVPlacement(0,                       //no rotation
                      AluminumSnout_pos,                    //at position
                      logicAluminumSnout,             //its logical volume
                      "AluminumSnout",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // Nylon collimator pieces (LTC) in collimator wall
  {
    std::vector<std::string> LTC_name = {"GA_GB_GC","G4","G5"};

    //  The extra 11 cm account for the lead/poly colimator and the insert depth for the aluminum snout
    G4double LTC_start = 2.865 * m;  

    std::ifstream template_file_stream;
    template_file_stream.open("CADModels/LTC.json");
    if(template_file_stream.is_open())
    {
      json template_file;
      template_file_stream >> template_file;    

      for(size_t i = 0; i < LTC_name.size(); ++i)
      {

        col_rmina = template_file[LTC_name.at(i)]["ID_start"];
        col_rminb = template_file[LTC_name.at(i)]["ID_end"];
        col_rmaxa = template_file[LTC_name.at(i)]["OD_start"];
        col_rmaxb = template_file[LTC_name.at(i)]["OD_end"];
        col_hz = template_file[LTC_name.at(i)]["length"];
        col_rel_pos = template_file[LTC_name.at(i)]["rel_pos"];

        G4Cons *LTC = new G4Cons(template_file[LTC_name.at(i)]["part_number"], 
        col_rmina/2 * cm, col_rmaxa/2 * cm, 
        col_rminb/2 * cm, col_rmaxb/2 * cm,
        col_hz/2 * cm, 0 * deg, 360 * deg);

        G4Material* LTC_mat = nist->FindOrBuildMaterial(template_file[LTC_name.at(i)]["material"]);
        G4ThreeVector LTC_pos = G4ThreeVector(0 * mm, 0*cm, (LTC_start + col_rel_pos * cm));

        G4LogicalVolume* logicLTC=                         
          new G4LogicalVolume(LTC,         //its solid
                              LTC_mat,          //its material
                              template_file[LTC_name.at(i)]["part_number"]);           //its name


        if(template_file[LTC_name.at(i)]["material"] ==  "G4_POLYETHYLENE")
        {
          logicLTC->SetVisAttributes(polyethVisAttributes);
        }
        else if(template_file[LTC_name.at(i)]["material"] ==  "G4_NYLON-6-6")
        {
          logicLTC->SetVisAttributes(nylonVisAttributes);
        }
        else
        {
          logicLTC->SetVisAttributes(genericVisAttributes);
        }
                    
        new G4PVPlacement(0,                       //no rotation
                          LTC_pos,                    //at position
                          logicLTC,             //its logical volume
                          template_file[LTC_name.at(i)]["part_number"],                //its name
                          logicWorld,                //its mother  volume
                          false,                   //no boolean operation
                          0,                       //copy number
                          checkOverlaps);  

      }
    }
    template_file_stream.close();
  }

  // PE collimator pieces in aluminum snout
  {
    std::vector<std::string> OTC_name = {"G14","G15","G16"};

    //  The extra 11 cm account for the lead/poly colimator and the insert depth for the aluminum snout
    G4double OTC_start = 228.06 * cm - 36 * cm + 20.3/2 * cm ;  

    std::ifstream template_file_stream;
    template_file_stream.open("CADModels/OTC.json");
    if(template_file_stream.is_open())
    {
      json template_file;
      template_file_stream >> template_file;    

      for(size_t i = 0; i < OTC_name.size(); ++i)
      {

        col_rmina = template_file[OTC_name.at(i)]["ID_start"];
        col_rminb = template_file[OTC_name.at(i)]["ID_end"];
        col_rmaxa = template_file[OTC_name.at(i)]["OD_start"];
        col_rmaxb = template_file[OTC_name.at(i)]["OD_end"];
        col_hz = template_file[OTC_name.at(i)]["length"];
        col_rel_pos = template_file[OTC_name.at(i)]["rel_pos"];

        G4Cons *OTC = new G4Cons(template_file[OTC_name.at(i)]["part_number"], 
        col_rmina/2 * cm, col_rmaxa/2 * cm, 
        col_rminb/2 * cm, col_rmaxb/2 * cm,
        col_hz/2 * cm, 0 * deg, 360 * deg);

        G4Material* OTC_mat = nist->FindOrBuildMaterial(template_file[OTC_name.at(i)]["material"]);
        G4ThreeVector OTC_pos = G4ThreeVector(0 * mm, 0*cm, (OTC_start + col_rel_pos * cm));

        G4LogicalVolume* logicOTC=                         
          new G4LogicalVolume(OTC,         //its solid
                              OTC_mat,          //its material
                              template_file[OTC_name.at(i)]["part_number"]);           //its name
        
        if(template_file[OTC_name.at(i)]["material"] ==  "G4_POLYETHYLENE")
        {
          logicOTC->SetVisAttributes(polyethVisAttributes);
        }
        else if(template_file[OTC_name.at(i)]["material"] ==  "G4_NYLON-6-6")
        {
          logicOTC->SetVisAttributes(nylonVisAttributes);
        }
        else
        {
          logicOTC->SetVisAttributes(genericVisAttributes);
        }     

        new G4PVPlacement(0,                       //no rotation
                          OTC_pos,                    //at position
                          logicOTC,             //its logical volume
                          template_file[OTC_name.at(i)]["part_number"],                //its name
                          logicWorld,                //its mother  volume
                          false,                   //no boolean operation
                          0,                       //copy number
                          checkOverlaps);  

      }
    }
    template_file_stream.close();
  }

  // Medium PE collimator pieces (MTC)
  {
    std::vector<std::string> MTC_name = {"M2","M3B","M3A","M4","M5","M6","M7"};

    //  The extra 11 cm account for the lead/poly colimator and the insert depth for the aluminum snout
    G4double MTC_start = 228.06 * cm - 36 * cm + 15.40/2 * cm;  

    std::ifstream template_file_stream;
    template_file_stream.open("CADModels/MTC.json");
    if(template_file_stream.is_open())
    {
      json template_file;
      template_file_stream >> template_file;    

      for(size_t i = 0; i < MTC_name.size(); ++i)
      {

        col_rmina = template_file[MTC_name.at(i)]["ID_start"];
        col_rminb = template_file[MTC_name.at(i)]["ID_end"];
        col_rmaxa = template_file[MTC_name.at(i)]["OD_start"];
        col_rmaxb = template_file[MTC_name.at(i)]["OD_end"];
        col_hz = template_file[MTC_name.at(i)]["length"];
        col_rel_pos = template_file[MTC_name.at(i)]["rel_pos"];

        G4Cons *MTC = new G4Cons(template_file[MTC_name.at(i)]["part_number"], 
        col_rmina/2 * cm, col_rmaxa/2 * cm, 
        col_rminb/2 * cm, col_rmaxb/2 * cm,
        col_hz/2 * cm, 0 * deg, 360 * deg);

        G4Material* MTC_mat = nist->FindOrBuildMaterial(template_file[MTC_name.at(i)]["material"]);
        G4ThreeVector MTC_pos = G4ThreeVector(0 * mm, 0*cm, (MTC_start + col_rel_pos * cm));

        G4LogicalVolume* logicMTC=                         
          new G4LogicalVolume(MTC,         //its solid
                              MTC_mat,          //its material
                              template_file[MTC_name.at(i)]["part_number"]);           //its name


        if(template_file[MTC_name.at(i)]["material"] ==  "G4_POLYETHYLENE")
        {
          logicMTC->SetVisAttributes(polyethVisAttributes);
        }
        else if(template_file[MTC_name.at(i)]["material"] ==  "G4_NYLON-6-6")
        {
          logicMTC->SetVisAttributes(nylonVisAttributes);
        }
        else
        {
          logicMTC->SetVisAttributes(genericVisAttributes);
        }
                    
        new G4PVPlacement(0,                       //no rotation
                          MTC_pos,                    //at position
                          logicMTC,             //its logical volume
                          template_file[MTC_name.at(i)]["part_number"],                //its name
                          logicWorld,                //its mother  volume
                          false,                   //no boolean operation
                          0,                       //copy number
                          checkOverlaps);  

      }
    }
    template_file_stream.close();
  }

  /* 
  *******************************************************************************
  *******************************************************************************
  * 
  * This section contains the physical structure of the source area, such as
  * walls, ceilings, and floors. Does not include any of the tunnel or
  * tunnel collimator structure
  * 
  *******************************************************************************
  *******************************************************************************
  */

  
  // Concrete Floor in source area
  {
    G4Material* ConcreteFloor_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
    G4ThreeVector ConcreteFloor_pos = G4ThreeVector((500 -(106.47+42.54*3./2.)) * cm, (-111.1 * cm + ConcreteYAdj), -2 * m+42.54*cm);
          
    // Conical section shape       
    G4Box* ConcreteFloor =    
      new G4Box("ConcreteFloor",
      6*m, 5 * cm, 6*m);

    G4LogicalVolume* logicConcreteFloor=                         
      new G4LogicalVolume(ConcreteFloor,         //its solid
                          ConcreteFloor_mat,          //its material
                          "ConcreteFloor");           //its name

    logicConcreteFloor->SetVisAttributes(concrete_structureVisAttributes);
                
    new G4PVPlacement(0,                       //no rotation
                      ConcreteFloor_pos,                    //at position
                      logicConcreteFloor,             //its logical volume
                      "ConcreteFloor",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // Concrete Moving Wall (goes to loading dock) in source area
  {
    G4Material* ConcreteMovingWall_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
    G4ThreeVector ConcreteMovingWall_pos = G4ThreeVector(-(106.47+42.54/2) * cm, (15.9 * cm + ConcreteYAdj), -94 * cm);
          
    // Conical section shape       
    G4Box* ConcreteMovingWall =    
      new G4Box("ConcreteMovingWall",
      (124/2)*cm, 122 * cm, 1.5*m);

    G4LogicalVolume* logicConcreteMovingWall=                         
      new G4LogicalVolume(ConcreteMovingWall,         //its solid
                          ConcreteMovingWall_mat,          //its material
                          "ConcreteMovingWall");           //its name

    logicConcreteMovingWall->SetVisAttributes(concrete_structureVisAttributes);
                
    new G4PVPlacement(0,                       //no rotation
                      ConcreteMovingWall_pos,                    //at position
                      logicConcreteMovingWall,             //its logical volume
                      "ConcreteMovingWall",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // Concrete Wall at entry to swinger area
  {
    G4Material* ConcreteSwingerAreaEntryWall_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
    G4ThreeVector ConcreteSwingerAreaEntryWall_pos = G4ThreeVector(-(106.47+124/2 + 42.54) * cm, (15.9 * cm + ConcreteYAdj), -450 * cm);
          
    // Conical section shape       
    G4Box* ConcreteSwingerAreaEntryWall =    
      new G4Box("ConcreteSwingerAreaEntryWall",
      (42.54/2)*cm, 122 * cm, 3*m);

    G4LogicalVolume* logicConcreteSwingerAreaEntryWall=                         
      new G4LogicalVolume(ConcreteSwingerAreaEntryWall,         //its solid
                          ConcreteSwingerAreaEntryWall_mat,          //its material
                          "ConcreteSwingerAreaEntryWall");           //its name

    logicConcreteSwingerAreaEntryWall->SetVisAttributes(concrete_structureVisAttributes);
                
    new G4PVPlacement(0,                       //no rotation
                      ConcreteSwingerAreaEntryWall_pos,                    //at position
                      logicConcreteSwingerAreaEntryWall,             //its logical volume
                      "ConcreteSwingerAreaEntryWall",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // Concrete Wall behind source
  {
    G4Material* ConcreteRearWall_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
    G4ThreeVector ConcreteRearWall_pos = G4ThreeVector(2.25 * m, (15.9 * cm + ConcreteYAdj), -4 * m);
          
    // Conical section shape       
    G4Box* ConcreteRearWall =    
      new G4Box("ConcreteRearWall",
      2.5*m, 122 * cm, (42.54/2)*cm);

    G4LogicalVolume* logicConcreteRearWall=                         
      new G4LogicalVolume(ConcreteRearWall,         //its solid
                          ConcreteRearWall_mat,          //its material
                          "ConcreteRearWall");           //its name

    logicConcreteRearWall->SetVisAttributes(concrete_structureVisAttributes);

    new G4PVPlacement(0,                       //no rotation
                      ConcreteRearWall_pos,                    //at position
                      logicConcreteRearWall,             //its logical volume
                      "ConcreteRearWall",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // Concrete Wall beside tunnel door alcove
  {
    G4Material* ConcreteFrontWall_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
    G4ThreeVector ConcreteFrontWall_pos = G4ThreeVector(5 * m - (34.325)*cm + 1.816*m, (15.9 * cm + ConcreteYAdj), 2.56 * m + (42.54/2)*cm);
          
    // Conical section shape       
    G4Box* ConcreteFrontWall =    
      new G4Box("ConcreteFrontWall",
      2.5*m, 122 * cm, (42.54/2)*cm);

    G4LogicalVolume* logicConcreteFrontWall=                         
      new G4LogicalVolume(ConcreteFrontWall,         //its solid
                          ConcreteFrontWall_mat,          //its material
                          "ConcreteFrontWall");           //its name

    logicConcreteFrontWall->SetVisAttributes(concrete_structureVisAttributes);

    new G4PVPlacement(0,                       //no rotation
                      ConcreteFrontWall_pos,                    //at position
                      logicConcreteFrontWall,             //its logical volume
                      "ConcreteFrontWall",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // Tunnel Door alcove side wall
  {
    G4Material* ConcreteAlcoveSideWall_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
    G4ThreeVector ConcreteAlcoveSideWall_pos = G4ThreeVector((2*1.816)*m + (34.325+42.54/2)*cm, (15.9 * cm + ConcreteYAdj), 2.56 * m + (42.54*2)*cm);
          
    // Conical section shape       
    G4Box* ConcreteAlcoveSideWall =    
      new G4Box("ConcreteAlcoveSideWall",
      (42.54/2)*cm, 122 * cm, 40.73*cm);

    G4LogicalVolume* logicConcreteAlcoveSideWall=                         
      new G4LogicalVolume(ConcreteAlcoveSideWall,         //its solid
                          ConcreteAlcoveSideWall_mat,          //its material
                          "ConcreteAlcoveSideWall");           //its name
    
    logicConcreteAlcoveSideWall->SetVisAttributes(concrete_structureVisAttributes);

    new G4PVPlacement(0,                       //no rotation
                      ConcreteAlcoveSideWall_pos,                    //at position
                      logicConcreteAlcoveSideWall,             //its logical volume
                      "ConcreteAlcoveSideWall",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // Tunnel Door alcove rear wall
  {
    G4Material* ConcreteAlcoveRearWall_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
    G4ThreeVector ConcreteAlcoveRearWall_pos = G4ThreeVector(1 * m - (34.325)*cm + 1.816*m +(42.54*1.5)*cm, (15.9 * cm + ConcreteYAdj), 2.56 * m + (42.54*3.5)*cm);
          
    // Conical section shape       
    G4Box* ConcreteAlcoveRearWall =    
      new G4Box("ConcreteAlcoveRearWall",
      1*m+(42.54/2)*cm, 122 * cm, (42.54/2)*cm);

    G4LogicalVolume* logicConcreteAlcoveRearWall=                         
      new G4LogicalVolume(ConcreteAlcoveRearWall,         //its solid
                          ConcreteAlcoveRearWall_mat,          //its material
                          "ConcreteAlcoveRearWall");           //its name

    logicConcreteAlcoveRearWall->SetVisAttributes(concrete_structureVisAttributes);

    new G4PVPlacement(0,                       //no rotation
                      ConcreteAlcoveRearWall_pos,                    //at position
                      logicConcreteAlcoveRearWall,             //its logical volume
                      "ConcreteAlcoveRearWall",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  /* 
  *******************************************************************************
  *******************************************************************************
  * 
  * This section contains the physical structure of the swinger, such as
  * target, target chamber, and magnet stucture
  * 
  *******************************************************************************
  *******************************************************************************
  */

  // Be9 target. ~5mm thick approx.
  {
    G4Material* Be9Target_mat = nist->FindOrBuildMaterial("G4_Be");
    G4ThreeVector Be9Target_pos = G4ThreeVector(0 * mm, 0* mm, 5 * mm);
          
    // Conical section shape       
    G4Box* Be9Target =    
      new G4Box("Be9Target",
      1.5 * cm, 1.5 * cm, 2.5 * mm);

    G4LogicalVolume* logicBe9Target =                         
      new G4LogicalVolume(Be9Target,         //its solid
                          Be9Target_mat,          //its material
                          "Be9Target");           //its name
                
    logicBe9Target->SetVisAttributes(genericVisAttributes);

    new G4PVPlacement(0,                       //no rotation
                      Be9Target_pos,                    //at position
                      logicBe9Target,             //its logical volume
                      "Be9Target",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // Large Target Chamber
  {
    G4Material* LargeTargetChamber_mat = nist->FindOrBuildMaterial("G4_Al");
    G4ThreeVector LargeTargetChamber_pos = G4ThreeVector(0 * mm, 0* mm, 0 * mm);
    auto LargeTargetChamber_rot = new G4RotationMatrix();
    LargeTargetChamber_rot->rotateX(M_PI/2*rad);    
    LargeTargetChamber_rot->rotateY(M_PI/2*rad);    


    
    // Hole in slab is 33 cm in diameter
    G4Tubs* LargeTargetChamber_casing
      = new G4Tubs("LargeTargetChamber_casing",
                  0 * cm,
                  30.48/2 * cm,
                  8 * cm,
                  0 * deg,
                  360 * deg);
    
    G4Tubs* LargeTargetChamber_void
      = new G4Tubs("LargeTargetChamber_void",
                  0 * cm,
                  29.48/2 * cm,
                  7 * cm,
                  0 * deg,
                  360 * deg);

    G4SubtractionSolid* LargeTargetChamber_hollow 
      = new G4SubtractionSolid("LargeTargetChamber_hollow",
      LargeTargetChamber_casing,
      LargeTargetChamber_void);

    G4Tubs* LargeTargetChamber_tophole
      = new G4Tubs("LargeTargetChamber_tophole",
                  0 * cm,
                  2.5/2 * cm,
                  30.48/2 * cm,
                  0 * deg,
                  360 * deg);

    G4ThreeVector LargeTargetChamber_tophole_pos = G4ThreeVector(0 * mm, +30.48/2 * cm, 0 * mm);
    auto LargeTargetChamber_tophole_rot = new G4RotationMatrix();
    LargeTargetChamber_tophole_rot->rotateX(-M_PI/2*rad); 
   
    G4SubtractionSolid* LargeTargetChamber_w_hole 
      = new G4SubtractionSolid("LargeTargetChamber_w_hole",
      LargeTargetChamber_hollow,
      LargeTargetChamber_tophole,LargeTargetChamber_tophole_rot,LargeTargetChamber_tophole_pos);

    G4Tubs* LargeTargetChamber_pipe
      = new G4Tubs("LargeTargetChamber_pipe",
                  1 * cm,
                  2.5/2 * cm,
                  22 * cm,
                  0 * deg,
                  360 * deg);

    G4ThreeVector LargeTargetChamber_pipe_pos = G4ThreeVector(0 * mm, (30.48/2) * cm - 0.5*cm+22*cm, 0 * mm);
    auto LargeTargetChamber_pipe_rot = new G4RotationMatrix();
    LargeTargetChamber_pipe_rot->rotateX(-M_PI/2*rad); 
   
    // G4VSolid* LTC_H = (G4VSolid *)LargeTargetChamber_hollow;

    G4UnionSolid* LargeTargetChamber 
      = new G4UnionSolid("LargeTargetChamber",
      LargeTargetChamber_w_hole,
      LargeTargetChamber_pipe,LargeTargetChamber_pipe_rot,LargeTargetChamber_pipe_pos); 

    G4LogicalVolume* logicLargeTargetChamber =                         
      new G4LogicalVolume(LargeTargetChamber,         //its solid
                          LargeTargetChamber_mat,          //its material
                          "LargeTargetChamber");           //its name

    logicLargeTargetChamber->SetVisAttributes(aluminumVisAttributes);
                
    new G4PVPlacement(LargeTargetChamber_rot,                       //no rotation
                      LargeTargetChamber_pos,                    //at position
                      logicLargeTargetChamber,             //its logical volume
                      "LargeTargetChamber",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // Standin for steel magnet structure
  {
    G4Material* IronMagnetStruct_mat = nist->FindOrBuildMaterial("G4_Fe");
    G4ThreeVector IronMagnetStruct_pos1 = G4ThreeVector(0 * mm, 5* cm, -45 * cm);
    G4ThreeVector IronMagnetStruct_pos2 = G4ThreeVector(0 * mm, -5* cm, -45 * cm);

    // Conical section shape       
    G4Box* IronMagnetStruct =    
      new G4Box("IronMagnetStruct",
      8 * cm, 3.5 * cm, 25 * cm);

    G4LogicalVolume* logicIronMagnetStruct_top =                         
      new G4LogicalVolume(IronMagnetStruct,         //its solid
                          IronMagnetStruct_mat,          //its material
                          "IronMagnetStruct_top");           //its name
    
    logicIronMagnetStruct_top->SetVisAttributes(steelVisAttributes);

    new G4PVPlacement(0,                       //no rotation
                  IronMagnetStruct_pos1,                    //at position
                  logicIronMagnetStruct_top,             //its logical volume
                  "IronMagnetStruct_top",                //its name
                  logicWorld,                //its mother  volume
                  false,                   //no boolean operation
                  0,                       //copy number
                  checkOverlaps);          //overlaps checking

    G4LogicalVolume* logicIronMagnetStruct_bottom =                         
      new G4LogicalVolume(IronMagnetStruct,         //its solid
                          IronMagnetStruct_mat,          //its material
                          "IronMagnetStruct_bottom");           //its name
    
    logicIronMagnetStruct_bottom->SetVisAttributes(steelVisAttributes);

    new G4PVPlacement(0,                       //no rotation
                      IronMagnetStruct_pos2,                    //at position
                      logicIronMagnetStruct_bottom,             //its logical volume
                      "IronMagnetStruct_bottom",                //its name
                      logicWorld,                //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking
  }

  // G4Material* matoriginbox = nist->FindOrBuildMaterial("G4_WATER");
  // G4ThreeVector posoriginbox = G4ThreeVector(0 * base_world_unit, 0* base_world_unit, 0 * base_world_unit);
        
  // // Conical section shape       
  // G4Box* originbox =    
  //   new G4Box("originbox",
  //   10 * cm, 10 * cm, 10 * cm);

  // // G4Box* solidShape6_2 =    
  // //   new G4Box("Shape6_box2",
  // //   500 * mm, 200 * mm, 120 * mm);

  // // G4Tubs* solidShape3_cylinder
  // //   = new G4Tubs("Shape3_cylinder",
  // //               0*cm,
  // //               10.67*cm,
  // //               90*cm,
  // //               startAngle,
  // //               spanningAngle);

  // // G4SubtractionSolid* solidShape6 
  // //   = new G4SubtractionSolid("Shape6",
  // //   solidShape6_1,
  // //   solidShape6_2);

  // G4LogicalVolume* logicOriginBox =                         
  //   new G4LogicalVolume(originbox,         //its solid
  //                       matoriginbox,          //its material
  //                       "OriginBox");           //its name
              
  // new G4PVPlacement(0,                       //no rotation
  //                   posoriginbox,                    //at position
  //                   logicOriginBox,             //its logical volume
  //                   "OriginBox",                //its name
  //                   logicWorld,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

  // Test CADMesh
  // auto TunnelStructureMesh = CADMesh::TetrahedralMesh::FromSTL("CADModels/fixed_front_concrete.stl");
  // TunnelStructureMesh->SetScale(1000);
  // G4Material* TunnelStructureMat = nist->FindOrBuildMaterial("G4_CONCRETE");

  // TunnelStructureMesh->SetMaterial(TunnelStructureMat);
  // G4VSolid* TunnelStructureMeshSolid = TunnelStructureMesh->GetSolid();

  // auto testmeshpos = G4ThreeVector(0 * base_world_unit, 0* base_world_unit, 2 * base_world_unit);
  // auto testmeshrot = new G4RotationMatrix;

  // auto tunnel_structure_assembly = TunnelStructureMesh->GetAssembly();

  // tunnel_structure_assembly->MakeImprint(logicWorld, testmeshpos, testmeshrot,0,false);



  // G4LogicalVolume* logicTunnelStructureMesh =                         
  //   new G4LogicalVolume(TunnelStructureMeshSolid,         //its solids
  //                       collimatormat,          //its material
  //                       "TestMesh");           //its name
              
  // new G4PVPlacement(0,                       //no rotation
  //                   testmeshpos,                    //at position
  //                   logicTunnelStructureMesh,             //its logical volume
  //                   "TestMesh",                //its name
  //                   logicWorld,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

  //
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
