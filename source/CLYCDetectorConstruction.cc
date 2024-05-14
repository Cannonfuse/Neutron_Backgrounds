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
#include "DetectorMessenger.hh"

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
#include "G4GenericMessenger.hh"

#include "G4Scintillation.hh"
#include "G4MaterialPropertiesTable.hh"

// for reading JSON files
#include "nlohmann/json.hpp"

// CADMesh, read STL geometry from disk and insert as tesselated solids into simulation
#include "CADMesh.hh"

using json = nlohmann::json;


// const G4double det_dist = 5 * m;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CLYCDetectorConstruction::CLYCDetectorConstruction()
: G4VUserDetectorConstruction(),
  fC6LYCPV(nullptr), fC7LYCPV(nullptr), fdummydetectorPV(nullptr)
  {
    // DefineCommands();
    fMessenger = new DetectorMessenger(this);

  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CLYCDetectorConstruction::~CLYCDetectorConstruction()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* CLYCDetectorConstruction::Construct()
{  
  // G4cout << "UseMTC = " << GetUseMTC() << ", UseLTC = " << GetUseLTC() << ", detDistance = " << GetDetDistance() << G4endl;


  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  // Detector distance; original distance was face of crystal @ 0.5 m, need to account for that with second line.
  G4double det_dist = 5 * m;
  // G4double det_dist = GetDetDistance();

  bool useStructure{false}, useDummy{false}, useBe9target{false}, useFTC{false}, useMTC{false}, useLTC{false}, useOneInchCLYC{false}, useThreeInchCLYC{false};

  if(GetUseStructure()) {useStructure = GetUseStructure();};
  if(GetUseDummy()) {useDummy = GetUseDummy();};
  if(GetUseBe9target()) {useBe9target = GetUseBe9target();};
  if(GetUseFTC()) {useFTC = GetUseFTC();};
  if(GetUseMTC()) {useMTC = GetUseMTC();};
  if(GetUseLTC()) {useLTC = GetUseLTC();};


  if(GetUseC6LYC()) {useOneInchCLYC = GetUseC6LYC();};
  if(GetUseC7LYC()) {useThreeInchCLYC = GetUseC7LYC();};

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
/*

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
*/

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
  // Get Li6 and print its properties    void SetUseC6LYC(G4bool value) {UseC6LYC = value;};
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

  // CLYC_molecule->SetMaterialPropertiesTable(CLYC_MPT);
  // CLYC_molecule->GetIonisation()->SetBirksConstant(6.95e-4 * cm / MeV);

  G4Material* C7LYC_molecule = new G4Material(name="C7LYC",density,ncomponents=4);
  C7LYC_molecule->AddElement(elCs, 20*perCent);
  C7LYC_molecule->AddElement(Li7Enh_el, 10*perCent);
  C7LYC_molecule->AddElement(elY, 10*perCent);
  C7LYC_molecule->AddElement(elCl, 60*perCent);

  // C7LYC_molecule->SetMaterialPropertiesTable(CLYC_MPT);
  // C7LYC_molecule->GetIonisation()->SetBirksConstant(6.95e-4 * cm / MeV);

  // bool use_ftc{false};
  // bool use_mtc{false};
  // bool use_ltc{false};
  // bool use_one_inch_clyc{false};
  // bool use_three_inch_clyc{false};

  // #ifdef USE_FTC
  // use_ftc=true;
  // #endif

  // #ifdef USE_MTC
  // use_mtc=true;
  // #endif

  // #ifdef USE_LTC
  // use_ltc=true;
  // #endif


  // use_mtc=CLYCDetectorConstruction::GetUseMTC();
  // use_ltc=CLYCDetectorConstruction::GetUseLTC();


    // printf("UseMTC = %d, UseLTC = %d\n, detDistance = %f",use_mtc,use_ltc,GetDetDistance());


  // use_one_inch_clyc = GetUseC6LYC();
  // use_three_inch_clyc = GetUseC7LYC();

  // #ifdef ONEINCHCLYC
  //   use_one_inch_clyc=true;
  // #endif

  // #ifdef THREEINCHCLYC
  //   use_three_inch_clyc=true;
  // #endif

  if(useOneInchCLYC)
  { 
    
      printf("\n\n\n\nC6LYC Distance = %f\n\n\n\n\n",GetC6LYCDistance()/m);
      // CLYC Crystal

      G4double innerRadius = 0.*cm;
      G4double outerRadius = 1.25*cm;
      G4double hz = (2.5/2)*cm;
      G4double startAngle = 0.*deg;
      G4double spanningAngle = 360.*deg;
      G4ThreeVector C6LYCpos = G4ThreeVector(GetC6LYC_X(), GetC6LYC_Y(), GetC6LYCDistance() + hz);

      printf("\n\n\n\nC6LYC Position = (%f,%f,%f\n\n\n\n\n",C6LYCpos.getX(),C6LYCpos.getY(),C6LYCpos.getZ());


      G4Tubs* c6lycTube
        = new G4Tubs("C6LYCcrystal",
                      innerRadius,
                      outerRadius,
                      hz,
                      startAngle,
                      spanningAngle);

      G4LogicalVolume* logicC6LYC =                         
        new G4LogicalVolume(c6lycTube,         //its solid
                            CLYC_molecule,          //its material
                            "C6LYCcrystal");           //its name

      fC6LYCPV = new G4PVPlacement(0,                       //no rotation
        C6LYCpos,                    //at position
        logicC6LYC,               //its logical volume
        "C6LYCcrystal",                //its name
        logicWorld,                //its mother  volume
        false,                   //no boolean operation
        0,                       //copy number
        checkOverlaps);          //overlaps checking  
    
    {
      // CLYC Brass Housing
      G4Material* CLYCDetectorHousing_Cu = nist->FindOrBuildMaterial("G4_Cu");
      G4Material* CLYCDetectorHousing_Zn = nist->FindOrBuildMaterial("G4_Zn");
      G4Material* CLYCDetectorHousing = new G4Material(name="CommonBrass",8.52 * g/cm3,ncomponents=2);
      CLYCDetectorHousing->AddMaterial(CLYCDetectorHousing_Cu, 60*perCent);
      CLYCDetectorHousing->AddMaterial(CLYCDetectorHousing_Zn, 40*perCent);

      G4Tubs* CLYCDetectorHousingTube
        = new G4Tubs("CLYCDetectorHousingTube",
                      (5.1/2.+.102) * cm, // 5.1 cm diameter of tube plus 2% to account for tube diameter variations
                      6.31/2. * cm, // 6.31 cm housing outer diameter
                      11.950/2. * cm,
                      0 * deg,
                      360 * deg);



      G4Cons *CLYCDetectorHousingFace = new G4Cons("CLYCDetectorHousingFace", 
        (3.480/2. - .102) * cm, 3.480/2. * cm, 
        (5.1/2.+.102) * cm, 6.31/2. * cm,
        0.945/2 * cm, 0 * deg, 360 * deg);
              

      G4ThreeVector CLYCDetectorHousingFace_pos = G4ThreeVector(0 * cm, 0 * cm, -6.4475 * cm);

      G4UnionSolid* CLYCDetectorHousingWithFront = new G4UnionSolid("CLYCDetectorHousingTubeFace",
      CLYCDetectorHousingTube, CLYCDetectorHousingFace, 0, CLYCDetectorHousingFace_pos);


      G4Tubs* CLYCDetectorHousingFaceCover
        = new G4Tubs("CLYCDetectorHousingFaceCover",
                      0 * cm, 
                      3.480/2. * cm, 
                      0.555/2. * cm,
                      0 * deg,
                      360 * deg);

      G4ThreeVector CLYCDetectorHousingFaceCover_pos = G4ThreeVector(0 * cm, 0 * cm, (-7.1975 + 0.555) * cm);

      G4UnionSolid* CLYCDetectorHousingComplete = new G4UnionSolid("CLYCDetectorHousingComplete",
      CLYCDetectorHousingWithFront, CLYCDetectorHousingFaceCover, 0, CLYCDetectorHousingFaceCover_pos);

      G4ThreeVector CLYCDetectorHousingComplete_pos = G4ThreeVector(C6LYCpos.getX(), C6LYCpos.getY(), det_dist + (13.45/2 - 1.19) * cm);


      G4LogicalVolume* logicCLYCDetectorHousingComplete =                         
        new G4LogicalVolume(CLYCDetectorHousingComplete,         //its solid
                            CLYCDetectorHousing,          //its material
                            "CLYCDetectorHousingComplete");           //its name

      new G4PVPlacement(0,                       //no rotation
                        CLYCDetectorHousingComplete_pos,                    //at position
                        logicCLYCDetectorHousingComplete,             //its logical volume
                        "CLYC_1_inch_Housing_Complete",                //its name
                        logicWorld,                //its mother  volume
                        false,                   //no boolean operation
                        0,                       //copy number
                        checkOverlaps);  

    }
    {
      // CLYC Window Holder 

      G4double innerRadiusWH= 1.3181*cm;
      G4double outerRadiusWH = 1.48955*cm;
      G4double hzWH = 0.1905*cm;
      G4double innerRadiusWH2 = 1.51495*cm;
      G4double outerRadiusWH2= 1.6229*cm;
      G4double hzWH2 = 0.127*cm;
      G4double startAngleWH = 0.*deg;
      G4double spanningAngleWH = 360.*deg;
      G4Material* CLYCWinH_mat = nist->FindOrBuildMaterial("G4_Al");


      // Front piece of window holder, holds crystal
      G4Tubs* clycWinH
        = new G4Tubs("CLYCWindowHolder1",
                      innerRadiusWH,
                      outerRadiusWH,
                      hzWH,
                      startAngleWH,
                      spanningAngleWH);
      
      // Back piece of window holder, holds quartz window
      G4Tubs* clycWinH2
        = new G4Tubs("CLYCWindowHolder2",
                      innerRadiusWH2,
                      outerRadiusWH2,
                      hzWH2,
                      startAngleWH,
                      spanningAngleWH);
      
      G4ThreeVector CLYCWinHpos2 = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit, 0.4445 * cm);


      G4UnionSolid* CLYCWindowHolder = new G4UnionSolid("CLYCWindowHolder",
      clycWinH, clycWinH2, 0, CLYCWinHpos2);

      G4LogicalVolume* logicclycWinH =                         
        new G4LogicalVolume(CLYCWindowHolder,         //its solid
                            CLYCWinH_mat,          //its material
                            "CLYCWindowHolder1");           //its name

      G4ThreeVector CLYCWinHpos = G4ThreeVector(C6LYCpos.getX(), C6LYCpos.getY(), det_dist+(2.5 - 0.2159) * cm);
           

      new G4PVPlacement(0,                       //no rotation
        CLYCWinHpos,                    //at position
        logicclycWinH,               //its logical volume
        "CLYCWindowHolder",                //its name
        logicWorld,                //its mother  volume
        false,                   //no boolean operation
        0,                       //copy number
        checkOverlaps);          //overlaps checking
    }
    {
        // CLYC Crystal Housing (Aluminum scintillator cup) 

        G4Material* CLYCScintCup_mat = nist->FindOrBuildMaterial("G4_Al");


        // Front piece of window holder, holds crystal
        G4Tubs* CLYCScintCup_Outer
          = new G4Tubs("CLYCScintCup_Outer",
                        0 * cm,
                        3.175/2. * cm,
                        3.175/2. * cm,
                        0.*deg,
                        360.*deg);
        
        // Back piece of window holder, holds quartz window
        G4Tubs* CLYCScintCup_Inner
          = new G4Tubs("CLYCScintCup_Inner",
                        0 * cm,
                        3.09372/2. * cm,
                        3.09372 * cm,
                        0.*deg,
                        360.*deg);
        
        G4ThreeVector CLYCScintCup_Inner_pos = G4ThreeVector(0 * cm, 0 * cm, 3.09372/2. * cm);


        G4SubtractionSolid* CLYCScintCup = new G4SubtractionSolid("CLYCScintCup",
        CLYCScintCup_Outer, CLYCScintCup_Inner, 0, CLYCScintCup_Inner_pos);

        G4LogicalVolume* logicCLYCScintCup =                         
          new G4LogicalVolume(CLYCScintCup,         //its solid
                              CLYCScintCup_mat,          //its material
                              "CLYCScintCup");           //its name

        G4ThreeVector CLYCScintCup_pos = G4ThreeVector(C6LYCpos.getX(), C6LYCpos.getY(), det_dist+(2.5 - (0.3375+0.635*2)) * cm);
            

        new G4PVPlacement(0,                       //no rotation
          CLYCScintCup_pos,                    //at position
          logicCLYCScintCup,               //its logical volume
          "CLYCScintCup",                //its name
          logicWorld,                //its mother  volume
          false,                   //no boolean operation
          0,                       //copy number
          checkOverlaps);          //overlaps checking
      }
    {
        // The glass tube
      G4Material* CLYCQuartzWindow_Si = nist->FindOrBuildMaterial("G4_Si");
      G4Material* CLYCQuartzWindow_O = nist->FindOrBuildMaterial("G4_O");
      G4Material* CLYCQuartzWindow_mat = new G4Material(name="Quartz",2.65 * g/cm3,ncomponents=2);
      CLYCQuartzWindow_mat->AddMaterial(CLYCQuartzWindow_Si, 33.33*perCent);
      CLYCQuartzWindow_mat->AddMaterial(CLYCQuartzWindow_O, 66.67*perCent);

        // Front piece of window holder, holds crystal
        G4Tubs* CLYCQuartzWindow
          = new G4Tubs("CLYCQuartzWindow_Outer",
                        0 * cm,
                        (1.186*2.54)/2. * cm,
                        (.060*2.54)/2. * cm,
                        0.*deg,
                        360.*deg);

        G4LogicalVolume* logicCLYCQuartzWindow =                         
          new G4LogicalVolume(CLYCQuartzWindow,         //its solid
                              CLYCQuartzWindow_mat,          //its material
                              "CLYCQuartzWindow");           //its name

        G4ThreeVector CLYCQuartzWindow_pos = G4ThreeVector(C6LYCpos.getX(), C6LYCpos.getY(), det_dist+(2.5 + ((.06+0.2775)*2.54)/2.) * cm);
            

        new G4PVPlacement(0,                       //no rotation
          CLYCQuartzWindow_pos,                    //at position
          logicCLYCQuartzWindow,               //its logical volume
          "CLYCQuartzWindow",                //its name
          logicWorld,                //its mother  volume
          false,                   //no boolean operation
          0,                       //copy number
          checkOverlaps);          //overlaps checking
    }
  }

  if(useThreeInchCLYC)  // CLYC Detector physical volume
  {

    printf("\n\n\n\nC7LYC Distance = %f\n\n\n\n\n",GetC7LYCDistance()/m);

    G4double innerRadius = 0.*cm;
    G4double outerRadius = 3.75*cm;
    G4double hz = (1./2)*cm;
    G4double startAngle = 0.*deg;
    G4double spanningAngle = 360.*deg;
    G4ThreeVector C7LYCpos = G4ThreeVector(GetC7LYC_X(), GetC7LYC_Y(), GetC7LYCDistance() + hz);

      printf("\n\n\n\nC7LYC Position = (%f,%f,%f\n\n\n\n\n",C7LYCpos.getX(),C7LYCpos.getY(),C7LYCpos.getZ());


    G4Tubs* c7lycTube
      = new G4Tubs("C7LYCcrystal",
                    innerRadius,
                    outerRadius,
                    hz,
                    startAngle,
                    spanningAngle);

    G4LogicalVolume* logicC7LYC =                         
      new G4LogicalVolume(c7lycTube,         //its solid
                          C7LYC_molecule,          //its material
                          "C7LYCcrystal");           //its name

    fC7LYCPV = new G4PVPlacement(0,                       //no rotation
      C7LYCpos,                    //at position
      logicC7LYC,               //its logical volume
      "C7LYCcrystal",                //its name
      logicWorld,                //its mother  volume
      false,                   //no boolean operation
      0,                       //copy number
      checkOverlaps);          //overlaps checking
  }

  /* 
  *******************************************************************************
  *******************************************************************************
  * 
  * This section contains a dummy detector to check the pattern emerging from the 
  * tunnel collimator
  * 
  *******************************************************************************
  *******************************************************************************
  */

if(useDummy)
{
  G4double dummy_innerRadius = 0.*cm;
  G4double dummy_outerRadius = 25*cm;
  G4double dummy_hz = 0.5 * mm;
  G4double dummy_startAngle = 0.*deg;
  G4double dummy_spanningAngle = 360.*deg;
  G4ThreeVector dummy_pos = G4ThreeVector(0 * mm, 0 * base_world_unit, 4.2*m);



  G4Tubs* dummydetector
    = new G4Tubs("dummydetector",
                  dummy_innerRadius,
                  dummy_outerRadius,
                  dummy_hz,
                  dummy_startAngle,
                  dummy_spanningAngle);

  G4LogicalVolume* logicdummydetector =                         
    new G4LogicalVolume(dummydetector,         //its solid
                        world_mat,          //its material
                        "dummydetector");           //its name

  fdummydetectorPV = new G4PVPlacement(0,                       //no rotation
    dummy_pos,                    //at position
    logicdummydetector,               //its logical volume
    "dummydetector",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking
}

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

  // New material colors

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

  // New materials

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

  // Variables
  const G4double ConcreteYAdj = .5 * cm;

  // Tunnel Structure (such as walls, floor)
  if(useStructure){  
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
  if(useStructure){
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
  if(useStructure){
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
  if(useStructure){
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
  if(useStructure){
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

  // Variables
  G4double col_rmina, col_rmaxa, col_rminb, col_rmaxb, col_hz, col_rel_pos;

  // Concrete collimator hole steel liner
  if(useStructure){
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
  if(useStructure){
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
  if(useStructure){
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
  if(useStructure){
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
  if(useLTC){
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
  if(useStructure){
    std::vector<std::string> OTC_name = {"G14","G15","G16"};

    //  The extra 11 cm account for the lead/poly colimator and the insert depth for the aluminum snout
    G4double OTC_start = 181 * cm  + 20.2/2 * cm ;//  - 36 * cm

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
  if(useMTC)
  {
    std::vector<std::string> MTC_name = {"M3B","M3A","M4","M5","M6","M7"};

    //  The extra 11 cm account for the lead/poly colimator and the insert depth for the aluminum snout
    G4double MTC_start = 181 * cm + 20.2/2 *cm;  //- 36 * cm 

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

        printf("Piece = %s, Z start = %f, Z end = %f\n",MTC_name.at(i).c_str(),MTC_pos.getZ()-MTC->GetZHalfLength(),MTC_pos.getZ()+MTC->GetZHalfLength());

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

    // Fine PE collimator pieces (FTC)
  if(useFTC)
  {
    std::vector<std::string> FTC_name = {"F15","F14","F13","F12","F11","F10","F9","F8","F7","F6","F5","F4","F3","F2","F1"};

    //  The extra 11 cm account for the lead/poly colimator and the insert depth for the aluminum snout
    G4double FTC_start = 228.06 * cm - 36 * cm + 2.01 * cm + 21.022*cm;  

    std::ifstream template_file_stream;
    template_file_stream.open("CADModels/FTC.json");
    if(template_file_stream.is_open())
    {
      json template_file;
      template_file_stream >> template_file;    

      for(size_t i = 0; i < FTC_name.size(); ++i)
      {

        col_rmina = template_file[FTC_name.at(i)]["ID_start"];
        col_rminb = template_file[FTC_name.at(i)]["ID_end"];
        col_rmaxa = template_file[FTC_name.at(i)]["OD_start"];
        col_rmaxb = template_file[FTC_name.at(i)]["OD_end"];
        col_hz = template_file[FTC_name.at(i)]["length"];
        col_rel_pos = template_file[FTC_name.at(i)]["rel_pos"];

        G4Cons *FTC = new G4Cons(template_file[FTC_name.at(i)]["part_number"], 
        col_rmina/2 * cm, col_rmaxa/2 * cm, 
        col_rminb/2 * cm, col_rmaxb/2 * cm,
        col_hz/2 * cm, 0 * deg, 360 * deg);

        G4Material* FTC_mat = nist->FindOrBuildMaterial(template_file[FTC_name.at(i)]["material"]);
        G4ThreeVector FTC_pos = G4ThreeVector(0 * mm, 0*cm, (FTC_start + col_rel_pos * cm));

        G4LogicalVolume* logicFTC=                         
          new G4LogicalVolume(FTC,         //its solid
                              FTC_mat,          //its material
                              template_file[FTC_name.at(i)]["part_number"]);           //its name


        if(template_file[FTC_name.at(i)]["material"] ==  "G4_POLYETHYLENE")
        {
          logicFTC->SetVisAttributes(polyethVisAttributes);
        }
        else if(template_file[FTC_name.at(i)]["material"] ==  "G4_NYLON-6-6")
        {
          logicFTC->SetVisAttributes(nylonVisAttributes);
        }
        else
        {
          logicFTC->SetVisAttributes(genericVisAttributes);
        }
                    
        new G4PVPlacement(0,                       //no rotation
                          FTC_pos,                    //at position
                          logicFTC,             //its logical volume
                          template_file[FTC_name.at(i)]["part_number"],                //its name
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
  if(useStructure){
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
  if(useStructure){
    G4Material* ConcreteMovingWall_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
    G4ThreeVector ConcreteMovingWall_pos = G4ThreeVector(-(106.47+42.54*3./2. + 1) * cm, (15.9 * cm + ConcreteYAdj), 25 * cm);
          
    // Conical section shape       
    G4Box* ConcreteMovingWall =    
      new G4Box("ConcreteMovingWall",
      (121.6/2)*cm, 122 * cm, 1.8*m);

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
  if(useStructure){
    G4Material* ConcreteSwingerAreaEntryWall_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
    G4ThreeVector ConcreteSwingerAreaEntryWall_pos = G4ThreeVector(-(106.47+124/2 + 2*42.54) * cm, (15.9 * cm + ConcreteYAdj), -375 * cm);
          
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
  if(useStructure){
    G4Material* ConcreteRearWall_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
    G4ThreeVector ConcreteRearWall_pos = G4ThreeVector(223.6/2 * cm + (146.6-77) *cm, (15.9 * cm + ConcreteYAdj), -(365.76 + 42.54/2) * cm);
          
    // Conical section shape       
    G4Box* ConcreteRearWall =    
      new G4Box("ConcreteRearWall",
      223.6*cm, 122 * cm, (42.54/2)*cm);

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
  if(useStructure){
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
  if(useStructure){
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
  if(useStructure){
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
  if(useBe9target){
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
  if(useStructure){
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
  if(useStructure){
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

  //
  //
  //always return the physical World
  //
  return physWorld;
}
