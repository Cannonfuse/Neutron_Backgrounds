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

#include "CLYCDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

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

  //     
  // World
  //
  G4double base_world_unit = 1*m;
  G4double world_sizeX = 1 * base_world_unit;
  G4double world_sizeY = 1 * base_world_unit;
  G4double world_sizeZ  = 1 * base_world_unit;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
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

  G4double a, density;
  G4int z, n, ncomponents, natoms;
  G4String name, symbol;
  G4int nLi_i{2}, nCl_i{2};

  // Define Li6 enhanced 
  G4double Li6_a,Li7_a;//,Li6enh_a;
  Li6_a = 6.01512288742 * g/mole;
  Li7_a = 7.01600343666 * g/mole;
  // Li6enh_a = 0.05 * Li7_a + 0.95 * Li6_a;

  G4Isotope* Li6_i = new G4Isotope(name="Li6",z=3,n=6,a=Li6_a);
  G4Isotope* Li7_i = new G4Isotope(name="Li7",z=3,n=7,a=Li7_a);
  G4Element* Li6Enh_el = new G4Element(name="Li6Enhanced",symbol="Li6Enh" ,nLi_i);
  // G4Element* Li6Enh_el = new G4Element(name="Li6Enhanced",symbol="Li6Enh" , z= 3., Li6enh_a);
  // Li6Enh_el->AddIsotope(Li6_i,95*perCent);
  // Li6Enh_el->AddIsotope(Li7_i,5*perCent);

  // G4Element* Li_i = nist->FindOrBuildElement("Li",true);
  // G4Isotope* Li6_i = const_cast<G4Isotope*>(Li_i->GetIsotope(6));
  // G4Isotope* Li7_i = const_cast<G4Isotope*>(Li_i->GetIsotope(7));

  Li6Enh_el->AddIsotope(Li6_i, 95. * perCent);
  Li6Enh_el->AddIsotope(Li7_i, 5. * perCent);


  // // Define natural chlorine 
  // G4double Cl35_a, Cl37_a;//, Clnat_a;
  // Cl35_a = 34.968852694 * g/mole;
  // Cl37_a = 36.965902584 * g/mole;
  // // Clnat_a = 0.7576 * Cl35_a + 0.2424 * Cl37_a;
  
  // G4Isotope* Cl35_i = new G4Isotope(name="Cl35",z=17,n=35,a=Cl35_a);
  // G4Isotope* Cl37_i = new G4Isotope(name="Cl37",z=17,n=37,a=Cl37_a);
  // G4Element* Clnat_el = new G4Element(name="natCl",symbol="Clnat", nCl_i);
  // Clnat_el->AddIsotope(Cl35_i,75.76*perCent);
  // Clnat_el->AddIsotope(Cl37_i,24.24*perCent);


  // // Define Caesium 
  // G4Element* elCs  = new G4Element(name="Caesium"  ,symbol="Cs" , z= 55., a = 132.90545196 * g/mole);
  
  // // Define Yttrium
  // G4Element* elY  = new G4Element(name="Yttrium"  ,symbol="Y" , z= 39., a = 88.90584 * g/mole);

  // Define Chlorine
  G4Material* elCl = nist->FindOrBuildMaterial("G4_Cl");

  // Define Caesium 
  G4Material* elCs  = nist->FindOrBuildMaterial("G4_Cs");
  
  // Define Yttrium
  G4Material* elY  = nist->FindOrBuildMaterial("G4_Y");

  // Define the final CLYC detector material
  density = 3.31 * g/cm3;
  G4Material* CLYC = new G4Material(name="CLYC",density,ncomponents=4);
  // CLYC->AddElement(elCs, natoms=2);
  // CLYC->AddElement(Li6Enh_el, natoms=1);
  // CLYC->AddElement(elY, natoms=1);
  // CLYC->AddElement(elCl, natoms=6);
  CLYC->AddMaterial(elCs, 20. * perCent);
  CLYC->AddElement(Li6Enh_el, 10. * perCent);
  CLYC->AddMaterial(elY, 10. * perCent);
  CLYC->AddMaterial(elCl, 60. * perCent);

  // CLYC Detector physical volume

  G4double innerRadius = 0.*cm;
  G4double outerRadius = 1.25*cm;
  G4double hz = 1.25*cm;
  G4double startAngle = 0.*deg;
  G4double spanningAngle = 360.*deg;
  G4ThreeVector CLYCpos = G4ThreeVector(0 * base_world_unit, 0 * base_world_unit, 0.5125 * base_world_unit);


   G4Tubs* clycTube
     = new G4Tubs("CLYCcrystal",
                  innerRadius,
                  outerRadius,
                  hz,
                  startAngle,
                  spanningAngle);

  G4LogicalVolume* logicCLYC =                         
    new G4LogicalVolume(clycTube,         //its solid
                        CLYC,          //its material
                        "CLYCcrystal");           //its name

  fCLYCPV = new G4PVPlacement(0,                       //no rotation
    CLYCpos,                    //at position
    logicCLYC,               //its logical volume
    "CLYCcrystal",                //its name
    logicWorld,                //its mother  volume
    false,                   //no boolean operation
    0,                       //copy number
    checkOverlaps);          //overlaps checking


  //     
  // Shape 1
  //  
  // G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  // G4ThreeVector pos1 = G4ThreeVector(0 * base_world_unit, 0* base_world_unit, 2.16635 * base_world_unit);
        
  // // Conical section shape       
  // G4double shape1_rmina =  2.49*cm, shape1_rmaxa = 9.4*cm;
  // G4double shape1_rminb =  1.06*cm, shape1_rmaxb = 9.4*cm;
  // G4double shape1_hz = 16.635*cm;
  // G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  // G4Cons* solidShape1 =    
  //   new G4Cons("Shape1", 
  //   shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
  //   shape1_phimin, shape1_phimax);
                      
  // G4LogicalVolume* logicShape1 =                         
  //   new G4LogicalVolume(solidShape1,         //its solid
  //                       shape1_mat,          //its material
  //                       "Shape1");           //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   pos1,                    //at position
  //                   logicShape1,             //its logical volume
  //                   "Shape1",                //its name
  //                   logicWorld,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

  // G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  // G4ThreeVector pos2 = G4ThreeVector(0 * base_world_unit, 0* base_world_unit, 2.97845 * m);
        
  // // Conical section shape       
  // G4double shape2_rmina =  .32*cm, shape2_rmaxa = 9.4*cm;
  // G4double shape2_rminb =  1.25*cm, shape2_rmaxb = 9.4*cm;
  // G4double shape2_hz = 64.575*cm;
  // G4double shape2_phimin = 0.*deg, shape2_phimax = 360.*deg;
  // G4Cons* solidShape2 =    
  //   new G4Cons("Shape2", 
  //   shape2_rmina, shape2_rmaxa, shape2_rminb, shape2_rmaxb, shape2_hz,
  //   shape2_phimin, shape2_phimax);
                      
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

  // G4Material* shape3_mat = nist->FindOrBuildMaterial("G4_CONCRETE");
  // G4ThreeVector pos3 = G4ThreeVector(0 * base_world_unit, 0* base_world_unit, 2.8121 * base_world_unit);
        
  // // Conical section shape       
  // G4Box* solidShape3_box =    
  //   new G4Box("Shape3_box",
  //   world_sizeX, world_sizeY, 0.8121 * base_world_unit);

  // G4Tubs* solidShape3_cylinder
  //   = new G4Tubs("Shape3_cylinder",
  //               0*cm,
  //               9.4*cm,
  //               90*cm,
  //               startAngle,
  //               spanningAngle);

  // G4SubtractionSolid* solidShape3 
  //   = new G4SubtractionSolid("Shape3",
  //   solidShape3_box,
  //   solidShape3_cylinder);

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

  // G4Material* shape4_mat = nist->FindOrBuildMaterial("G4_AIR");
  // G4ThreeVector pos4 = G4ThreeVector(0 * base_world_unit, 1* base_world_unit, 2.8121 * base_world_unit);
        
  // // Conical section shape       
  // G4double shape4_rmina =  0*cm, shape4_rmaxa = 9.4*cm;
  // G4double shape4_rminb =  0*cm, shape4_rmaxb = 9.4*cm;
  // G4double shape4_hz = 81.21*cm;
  // G4double shape4_phimin = 0.*deg, shape4_phimax = 360.*deg;
  // G4Cons* solidShape4 =    
  //   new G4Cons("Shape4", 
  //   shape4_rmina, shape4_rmaxa, shape4_rminb, shape4_rmaxb, shape4_hz,
  //   shape4_phimin, shape4_phimax);
                      
  // G4LogicalVolume* logicShape4 =                         
  //   new G4LogicalVolume(solidShape4,         //its solid
  //                       shape4_mat,          //its material
  //                       "Shape4");           //its name
               
  // new G4PVPlacement(0,                       //no rotation
  //                   pos4,                    //at position
  //                   logicShape4,             //its logical volume
  //                   "Shape4",                //its name
  //                   logicWorld,                //its mother  volume
  //                   false,                   //no boolean operation
  //                   0,                       //copy number
  //                   checkOverlaps);          //overlaps checking

/*
  //     
  // Envelope
  //  

  // Envelope parameters
  //
  G4double env_sizeX = 20*cm, env_sizeY = 20*cm, env_sizeZ = 1*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeX, 0.5*env_sizeY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking   
  */
  /*
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //     
  // Shape 1
  //  
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);
        
  // Conical section shape       
  G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  G4double shape1_hz = 3.*cm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  G4Cons* solidShape1 =    
    new G4Cons("Shape1", 
    shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb, shape1_hz,
    shape1_phimin, shape1_phimax);
                      
  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                        shape1_mat,          //its material
                        "Shape1");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //     
  // Shape 2
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4ThreeVector pos2 = G4ThreeVector(0, -1*cm, 7*cm);

  // Trapezoid shape       
  G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  G4double shape2_dz  = 6*cm;      
  G4Trd* solidShape2 =    
    new G4Trd("Shape2",                      //its name
              0.5*shape2_dxa, 0.5*shape2_dxb, 
              0.5*shape2_dya, 0.5*shape2_dyb, 0.5*shape2_dz); //its size
                
  G4LogicalVolume* logicShape2 =                         
    new G4LogicalVolume(solidShape2,         //its solid
                        shape2_mat,          //its material
                        "Shape2");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos2,                    //at position
                    logicShape2,             //its logical volume
                    "Shape2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                
  // Set Shape2 as scoring volume
  //
  fScoringVolume = logicShape2;
  */
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
