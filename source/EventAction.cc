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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
#include "CLYCDetectorConstruction.hh"
// #include "G4AnalysisManager.hh"
#include "g4root.hh"


#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
// #include "PhysicalConstants.hh"
#include <random>

#include <chrono>
using namespace std::chrono;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// G4double energyRes(G4double edep)
// {
//   std::random_device DEVICE;
//   std::mt19937_64 GENERATOR(DEVICE());
//   // std::default_random_engine GENERATOR;
//   double EDEP = static_cast<double> (edep);
//   EDEP *= 1000.;
//   // G4cout << edep * MeV << G4endl;
//   // G4cout << EDEP << G4endl;
   
//   std::normal_distribution<double> resA(-1.15446028e-04, 1.85495156e-05);
//   std::normal_distribution<double> resB(2.14550165e-01, 1.69963576e-02);
//   double stddev = (  resA(GENERATOR) * EDEP + resB(GENERATOR)) / (2 * std::sqrt(2. * std::log(2)));
//   stddev *= EDEP;
//   // G4cout << stddev << G4endl;

//   // std::normal_distribution<double> edepRes(0, stddev);
//   // double edepAdd = edepRes(GENERATOR);
//   // G4double modedep = edep + edepAdd * keV;

//   std::normal_distribution<double> edepRes(static_cast<double> (edep) * 1000., stddev);
//   return edepRes(GENERATOR) * keV;
// }

// G4double detEnergyResponse(G4double edep)
// {
//   std::random_device DEVICE;
//   std::mt19937_64 GENERATOR(DEVICE());
//   // std::default_random_engine GENERATOR;
//   double EDEP = static_cast<double> (edep);
//   EDEP *= 1000.;
//   std::normal_distribution<double> resA(0, 1);
//   if(EDEP <= 300.)
//   {
//     if(resA(GENERATOR) >= 0)
//     {
//       return EDEP * keV;
//     }
//     return 0 * keV;
//   }
//   return EDEP * keV;
// }

// G4double calcTime(G4double start_KE, G4double end_KE, G4double DIST)
// {
//   // G4cout << DIST << G4endl;
//   // lambda E,M: constants.c * np.sqrt(1 - np.power(1+E/M,-2))
//   G4double vstart = CLHEP::c_light * sqrt(1 - pow(1 + start_KE/(939.550*MeV),-2));
//   // G4double vend = CLHEP::c_light * sqrt(1 - pow(1 + end_KE/(939.550*MeV),-2));
//   // G4double vavg = (vstart - vend)/2;

//   // if(DIST > 0)
//   // {
//   //   return ((DIST/mm)/v);
//   // }

//   // G4cout << (DIST/mm)/vstart << ", " << (DIST/mm)/vavg  << G4endl;

//   if(DIST > 0)
//   {
//     return ((DIST/mm)/vstart);
//   }
//     return 0;
// }


EventAction::EventAction(CLYCDetectorConstruction* detConstruction)
: G4UserEventAction(),
fDetConstruction(detConstruction)
{
  // Clear the histogram vectors if we aren't using them
  EnergyAngle_dist.clear();
  EnergyZ_dist.clear();
  EnergyAngleZ_bins.clear();
  NeutronsData.clear();

  C6LYC_Slices = fDetConstruction->GetC6LYC_Slices();
  C7LYC_Slices = fDetConstruction->GetC7LYC_Slices();

  // Clear the other vectors
  ClearVectors();
  // Clear all variables
  ClearVariables();
} 

EventAction::EventAction(CLYCDetectorConstruction* detConstruction, const G4bool usedists, 
                const std::vector<std::vector<G4bool>> energyangledist,
                const std::vector<std::vector<G4bool>> energyzdist,
                const std::vector<std::vector<double>> energyangzbins)
: G4UserEventAction(), fDetConstruction(detConstruction),
generator(nullptr), En(nullptr), Ang(nullptr), Zh(nullptr), 
Ebin(nullptr), AngleBin(nullptr), ZBin(nullptr)
{
  // Clear the other vectors
  ClearVectors();
  // Clear all variables
  ClearVariables();

  SetEnergyAngle_dist(energyangledist);
  SetEnergyZ_dist(energyzdist);
  SetEnergyAngleZ_bins(energyangzbins);
  SetUseDists(usedists);

  C6LYC_Slices = fDetConstruction->GetC6LYC_Slices();
  C7LYC_Slices = fDetConstruction->GetC7LYC_Slices();

  generator = new std::minstd_rand();

  En = new std::uniform_int_distribution<int>;
  Ang = new std::uniform_int_distribution<int>;
  Zh = new std::uniform_int_distribution<int>;

  Ebin = new std::uniform_real_distribution<double>;
  AngleBin = new std::uniform_real_distribution<double>;
  ZBin = new std::uniform_real_distribution<double>;
} 

EventAction::EventAction(CLYCDetectorConstruction* detConstruction,
                         const G4bool useneutrons, 
                         const std::vector<std::vector<double>> neutronsdata)
: G4UserEventAction(), fDetConstruction(detConstruction),
generator(nullptr), neutron(nullptr), Ebin(nullptr), AngleBin(nullptr), ZBin(nullptr),
Rho(nullptr), Phi(nullptr)
{
  // Clear the other vectors
  ClearVectors();
  // Clear all variables
  ClearVariables();

  SetNeutronsData(neutronsdata);
  SetUseNeutronsData(useneutrons);

  generator = new std::minstd_rand();

  neutron = new std::uniform_int_distribution<int>;

  Ebin = new std::uniform_real_distribution<double>;
  AngleBin = new std::uniform_real_distribution<double>;
  ZBin = new std::uniform_real_distribution<double>;
  Rho = new std::uniform_real_distribution<double>;
  setRhoParam(0,(double)detConstruction->GetGasCellDiameter()/2);
  Phi = new std::uniform_real_distribution<double>;
  setPhiParam(-M_PI,M_PI);

  C6LYC_Slices = fDetConstruction->GetC6LYC_Slices();
  C7LYC_Slices = fDetConstruction->GetC7LYC_Slices();
  // printf("TEST\n");

} 



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{ 
  // auto start = high_resolution_clock::now();

  // auto analysisManager = G4AnalysisManager::Instance(); 

  // printf("EventID = %i\n", event->GetEventID());
  if(GetUseDists())
  {

    // auto eangle = GetEnergyAngle_dist();
    // auto ez = GetEnergyZ_dist();
    // auto eanglez = GetEnergyAngleZ_bins();

    //     std::uniform_int_distribution<int> En(0,eangle[0].size()-1);
    // std::uniform_int_distribution<int> Ang(0,eangle.size()-1);
    // std::uniform_int_distribution<int> Zh(0,ez.size()-1);

    // Set the bin distribution parameters
    setEnintParam(0, GetEnergySize()-1);
    setAngintParam(0, GetAngleSize()-1);
    setZhintParam(0, GetZedSize()-1);

    int en{getEnintValue()}, ang{getAngintValue()}, z{getZhintValue()};

    if(EnergyAngleValue(ang, en) == true && EnergyZValue(z, en)== true)
    {


      double Elo = GetEnergyAngleZ_Value(0, en);
      double Ehi = GetEnergyAngleZ_Value(0, en+1);
      double Anglo = GetEnergyAngleZ_Value(1, ang);
      double Anghi = GetEnergyAngleZ_Value(1, ang+1);
      double Zlo = GetEnergyAngleZ_Value(2, z);
      double Zhi = GetEnergyAngleZ_Value(2, z+1);

      setEParam(Elo,Ehi);
      setAngleParam(Anglo,Anghi);
      setZedParam(Zlo,Zhi);

      double E{getEValue()},Angle{getAngleValue()},Zed{getZedValue()};

      auto particlepos = event->GetPrimaryVertex()->GetPosition();
      auto momentum = event->GetPrimaryVertex()->GetPrimary()->GetMomentumDirection();
      particlepos.setZ(Zed*mm);
      momentum.setTheta(Angle*radian);


      event->GetPrimaryVertex()->GetPrimary()->SetKineticEnergy(E * MeV);
      event->GetPrimaryVertex()->SetPosition(particlepos.getX(),particlepos.getY(),particlepos.getZ());
      event->GetPrimaryVertex()->GetPrimary()->SetMomentumDirection(momentum);
      // printf("p_f, x = %f, y = %f, z = %f, theta = %f\n",particlepos.getX(),particlepos.getY(),particlepos.getZ(),particlepos.unit().getTheta());

      // printf("Jackpot! E = %i, Ang = %i, Z = %i\n",en,ang,z);
      // printf("Max Size: E = %f, Ang = %f, Z = %f\n",E,Angle,Zed);
    } 
    // else
    // {
    //   // printf("Well, shit. E = %i, Ang = %i, Z = %i\n",en,ang,z);
    // }
  }
  if(GetUseNeutronsData())
  {


    // auto stop0 = high_resolution_clock::now();

    setNeutronParam(0,GetNeutronsDataSize()-1);

    // see which neutron to sample
    auto theneutron = getNeutronValue();

    // return the parameters for a neutron bin
    std::vector<double> neutronbins = GetNeutronData(theneutron);
    

    setEParam(neutronbins[0],neutronbins[1]);
    setAngleParam(neutronbins[2],neutronbins[3]);
    setZedParam(neutronbins[4],neutronbins[5]);

    auto gascellposition = fDetConstruction->GetGasCellPosition();


    double Zed{getZedValue()+gascellposition};

    auto particlepos = event->GetPrimaryVertex()->GetPosition();

    // printf("Original momentum: R = %f, phi = %f, theta = %f\n",Pmag,Pphi,Ptheta);
    // printf("New momentum: R = %f, phi = %f, theta = %f\n",momentum.mag(),momentum.phi(),momentum.theta());
    // printf("Original momentum (unit): R = %f, phi = %f, theta = %f\n",momentum.unit().mag(),momentum.unit().phi(),momentum.unit().theta());
    // printf("# of vertex = %i, # of particles = %i\n",event->GetNumberOfPrimaryVertex(),event->GetPrimaryVertex()->GetNumberOfParticle());
    double rho=getRhoValue();
    double phi=getPhiValue();

    G4double xpos = rho*cos(phi)*mm;
    G4double ypos = rho*sin(phi)*mm;
    G4double zpos = Zed*mm;

    // printf("x,y,z = %f,%f,%f\n",xpos,ypos,zpos);

    event->GetPrimaryVertex()->SetPosition(xpos,ypos,zpos);

    for(auto i = 0; i < event->GetPrimaryVertex()->GetNumberOfParticle(); ++i)
    {
      double E{getEValue()}, Angle{getAngleValue()};
      auto momentum = event->GetPrimaryVertex()->GetPrimary(i)->GetMomentumDirection();
      auto Pphi = momentum.phi();
      auto Ptheta = momentum.theta();
      momentum.setTheta(Angle*radian);
      momentum.setPhi(Pphi*radian);
      event->GetPrimaryVertex()->GetPrimary(i)->SetKineticEnergy(E * MeV);
      event->GetPrimaryVertex()->GetPrimary(i)->SetMomentumDirection(momentum);
    }




    // auto stop4 = high_resolution_clock::now();

    // auto duration0 = duration_cast<microseconds>(stop0 - start);
    // auto duration1 = duration_cast<microseconds>(stop1 - stop0);
    // auto duration2 = duration_cast<microseconds>(stop2 - stop1);
    // auto duration3 = duration_cast<microseconds>(stop3 - stop2);
    // auto duration4 = duration_cast<microseconds>(stop4 - stop3);

    // printf("T0 = %i us, T1 = %i us,T2 = %i us,T3 = %i us,T4 = %i us\n",duration0,duration1,duration2,duration3,duration4);

  }

  // Clear the vectors or you will get a memory leak
  ClearVectors();
  // Clear all variables
  ClearVariables();

  // Add some stuff to the vectors descibing the initial position of the particle
  // and some of it's features
  {
    int startA = event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetAtomicMass();
    int startZ = event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetAtomicNumber();
    auto startposition = event->GetPrimaryVertex()->GetPosition();


    AddTo_pos_x(startposition.getX());
    AddTo_pos_y(startposition.getY());
    AddTo_pos_z(startposition.getZ());
    AddToEdepVector(0 * MeV);
    AddToSliceVector(-1);
    AddToAVector(startA);
    AddToZVector(startZ);
    AddToLstepVector(0 * mm);
    AddToEndKEVector(retGunEnergy());
    AddToDeltaTVector(0 * ns);


  }
  
  // Set all the variables at the beginning of the event
  {
    auto P = event->GetPrimaryVertex()->GetPrimary()->GetMomentumDirection();

    setEventID(event->GetEventID());
    SetGunEnergy(event->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy());
    SetTheta(P.getTheta());
  }




  // fGlobalTime = 0 * ns;
  // fEdep = 0.;
  // fPreDetectorEnergy = 0.;
  // fLstep = 0.;
  // fSteps = 0;
  // fDeltat = 0.;
  // inDetDeltaT = 0.;
  // inDetDeltaD = 0.;
  // A = 0;
  // Z = 0;
  // fdummyXPosition = 0;
  // fdummyYPosition = 0;
  // fdummyZPosition = 0;
  // fDummy=false;
  // hasCl35=false;
  // hasLi6=false;
  // goodPreDet=false;
  // Detector=0;
  // LstepVector.clear();
  // DeltaTVector.clear();
  // endKEVector.clear();
  // EdepVector.clear();
  // particleAVector.clear();
  // particleZVector.clear();
  // detectorSliceVector.clear();
  // pos_x.clear();
  // pos_y.clear();
  // pos_z.clear();







  // printf("GunEnergy = %f\n",(double)fGunEnergy);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{ 
    // if(LstepVector.size() > 4){
    // printf("sizeof(LstepVector) = %llu\n",sizeof(LstepVector.data()));
    // printf("sizeof(DeltaTVector) = %llu\n",sizeof(DeltaTVector.data()));
    // printf("sizeof(endKEVector) = %llu\n",sizeof(endKEVector.data()));
    // printf("sizeof(EdepVector) = %llu\n",sizeof(EdepVector.data()));
    // printf("sizeof(particleZVector) = %llu\n",sizeof(particleZVector.data()));
    // printf("sizeof(particleAVector) = %llu\n",sizeof(particleAVector.data()));
    // printf("sizeof(detectorSliceVector) = %llu\n",sizeof(detectorSliceVector.data()));
    // printf("sizeof(pos_x) = %llu\n",sizeof(pos_x.data()));
    // printf("sizeof(pos_y) = %llu\n",sizeof(pos_y.data()));
    // printf("sizeof(pos_z) = %llu\n",sizeof(pos_z.data()));
    // printf("sizeof(EnergyAngle_dist) = %llu\n",sizeof(EnergyAngle_dist.data()));
    // printf("sizeof(EnergyZ_dist) = %llu\n",sizeof(EnergyZ_dist.data()));
    // printf("sizeof(EnergyAngleZ_bins) = %llu\n",sizeof(EnergyAngleZ_bins.data()));
    // printf("sizeof(NeutronsData) = %llu\n",sizeof(NeutronsData.data()));}

    // hasCl35
    // hasLi6
    // goodPreDet
    // fEventID
    // fTheta
    // fEdep
    // fLstep
    // fDeltat
    // fGlobalTime
    // inDetDeltaT
    // inDetDeltaD
    // fGunEnergy
    // fKineticEnergy
    // particleX
    // particleY
    // particleZ
    // fXPosition
    // fYPosition
    // fZPosition
    // fdummyXPosition
    // fdummyYPosition
    // fdummyZPosition
    // fDummy
    // fPreDetectorEnergy
    // fSteps
    // A 
    // Z
    // Detector


  // This is the analysis manager
  auto analysisManager = G4AnalysisManager::Instance(); 
  auto  dataNtuple = analysisManager->GetFirstNtupleId();
  auto reacs6Ntuple = analysisManager->GetFirstNtupleId() + 1;
  auto reacs7Ntuple = analysisManager->GetFirstNtupleId() + 2;
  auto dummyNtuple = analysisManager->GetFirstNtupleId() + 3;

  // auto reacsParticlesNtuple = analysisManager->GetFirstNtupleId() + 2;

  
  // G4double new_Deltat{0};
  // G4int ntupleID = analysisManager->GetFirstNtupleId();

  // if(LstepVector.size() > 0)
  // {
  //   for(auto i = 0; i < LstepVector.size(); ++i)
  //   {
  //     new_Deltat += calcTime(startKEVector.at(i)/MeV,endKEVector.at(i)/MeV,LstepVector.at(i)/mm);
  //   }
  //   new_Deltat/=ns;
  // }

  // G4cout << fDeltat << ", " << new_Deltat << G4endl;

  // analysisManager->FillNtupleDColumn(dataNtuple, 0, fGunEnergy);
  // analysisManager->AddNtupleRow(dataNtuple);

  // analysisManager->FillNtupleDColumn(dataNtuple, 0, fEdep);

  if(false)
  {
  analysisManager->FillNtupleDColumn(dataNtuple, 0, fGunEnergy);
  analysisManager->FillNtupleDColumn(dataNtuple, 1, fGlobalTime);
  analysisManager->FillNtupleDColumn(dataNtuple, 2, fTheta);

  analysisManager->AddNtupleRow(dataNtuple);
  }
  
  /*
  analysisManager->FillNtupleDColumn(dataNtuple, 0, fEdep);
  analysisManager->FillNtupleDColumn(dataNtuple, 1, fLstep);
  analysisManager->FillNtupleDColumn(dataNtuple, 2, fDeltat);
  analysisManager->FillNtupleDColumn(dataNtuple, 3, fGunEnergy);
  analysisManager->FillNtupleDColumn(dataNtuple, 4, fPreDetectorEnergy);
  analysisManager->FillNtupleIColumn(dataNtuple, 5, Detector);

  analysisManager->FillNtupleDColumn(dataNtuple, 0, fEdep);
  analysisManager->FillNtupleDColumn(dataNtuple, 1, fLstep);
  analysisManager->FillNtupleDColumn(dataNtuple, 2, fDeltat);
  analysisManager->FillNtupleDColumn(dataNtuple, 3, fGunEnergy);
  analysisManager->FillNtupleIColumn(dataNtuple, 4, Z);
  analysisManager->FillNtupleIColumn(dataNtuple, 5, A);
  analysisManager->FillNtupleDColumn(dataNtuple, 6, fKineticEnergy);
  analysisManager->FillNtupleIColumn(dataNtuple, 7, fSteps);
  analysisManager->FillNtupleDColumn(dataNtuple, 8, fXPosition);
  analysisManager->FillNtupleDColumn(dataNtuple, 9, fYPosition);
  analysisManager->FillNtupleDColumn(dataNtuple, 10, fZPosition);
  analysisManager->FillNtupleDColumn(dataNtuple, 11, fPreDetectorEnergy);
  analysisManager->FillNtupleIColumn(dataNtuple, 12, Detector);

  analysisManager->AddNtupleRow(dataNtuple);
  */
 
  if(fDummy)
  {  
    // analysisManager->FillNtupleDColumn(dummyNtuple, 0, fdummyXPosition);
    // analysisManager->FillNtupleDColumn(dummyNtuple, 1, fdummyYPosition);
    // analysisManager->FillNtupleDColumn(dummyNtuple, 2, fdummyZPosition);


    analysisManager->AddNtupleRow(dummyNtuple);
    analysisManager->FillH2(analysisManager->GetFirstH2Id()+17,fdummyXPosition,fdummyYPosition);

  }  

  // G4cout << analysisManager->GetFirstNtupleId() << " , " << analysisManager->GetFirstNtupleId()+1 << G4endl;


  // G4cout << G4BestUnit(fLstep, "Length")<< ", " << G4BestUnit(fEdep, "Energy") << G4endl;

  // G4cout << analysisManager->GetNofNtuples() << G4endl;

  analysisManager->FillH1(0,fLstep);
  analysisManager->FillH1(1,fEdep);
  analysisManager->FillH1(2,fDeltat);
  analysisManager->FillH1(3,fGunEnergy);
  analysisManager->FillH1(16,fTheta);



  // analysisManager->FillH2(analysisManager->GetFirstH2Id(),fEdep,fLstep);
  // analysisManager->FillH2(analysisManager->GetFirstH2Id()+1,fDeltat,fEdep);
  // analysisManager->FillH2(analysisManager->GetFirstH2Id()+2,fGunEnergy,fEdep);
  // analysisManager->FillH2(analysisManager->GetFirstH2Id()+3,fGunEnergy,fLstep);
  // analysisManager->FillH2(analysisManager->GetFirstH2Id()+4,Z,A);
  // analysisManager->FillH2(analysisManager->GetFirstH2Id()+7,fSteps,fKineticEnergy);

  if(isC7LYC)
  {
    analysisManager->FillH1(analysisManager->GetFirstH1Id()+15,fGunEnergy);
    analysisManager->FillH1(analysisManager->GetFirstH1Id()+24,fPreDetectorEnergy);
    if(SecondaryNeutrons.size() > 0)
    {
      for(auto i = 0; i < SecondaryNeutrons.size(); ++i)
      {analysisManager->FillH1(analysisManager->GetFirstH1Id()+24,SecondaryNeutrons.at(i));}
    }
  }
  if(isC6LYC)
  {
    analysisManager->FillH1(analysisManager->GetFirstH1Id()+14,fGunEnergy);
    analysisManager->FillH1(analysisManager->GetFirstH1Id()+23,fPreDetectorEnergy);
    if(SecondaryNeutrons.size() > 0)
    {
      for(auto i = 0; i < SecondaryNeutrons.size(); ++i)
      {analysisManager->FillH1(analysisManager->GetFirstH1Id()+23,SecondaryNeutrons.at(i));}
    }

  }

  if(Detector == 6 && fDetConstruction->GetUseC6LYC())
  {

    // Fill the 3D histogram with the start position of the particle
    analysisManager->FillH3(analysisManager->GetFirstH3Id(),Get_pos_x_vector().at(0),Get_pos_y_vector().at(0),Get_pos_z_vector().at(0));

    if(fDetConstruction->GetUseC7LYC() && isC7LYC)
    {
      analysisManager->FillH1(analysisManager->GetFirstH1Id()+19,fGunEnergy);
    }


    analysisManager->FillH1(analysisManager->GetFirstH1Id()+10,fGlobalTime);
    // analysisManager->FillH1(analysisManager->GetFirstH1Id()+14,fGunEnergy);

    // analysisManager->FillNtupleDColumn(reacs6Ntuple, 0, fEdep);
    // analysisManager->FillNtupleDColumn(reacs6Ntuple, 1, fLstep);
    // analysisManager->FillNtupleDColumn(reacs6Ntuple, 2, fDeltat);
    // analysisManager->FillNtupleDColumn(reacs6Ntuple, 3, fGunEnergy);
    // analysisManager->FillNtupleIColumn(reacs6Ntuple, 4, Z);
    // analysisManager->FillNtupleIColumn(reacs6Ntuple, 5, A);
    // analysisManager->FillNtupleDColumn(reacs6Ntuple, 6, fKineticEnergy);
    // analysisManager->FillNtupleIColumn(reacs6Ntuple, 7, fSteps);
    // analysisManager->FillNtupleDColumn(reacs6Ntuple, 8, fXPosition);
    // analysisManager->FillNtupleDColumn(reacs6Ntuple, 9, fYPosition);
    // analysisManager->FillNtupleDColumn(reacs6Ntuple, 10, fZPosition);
    // analysisManager->FillNtupleDColumn(reacs6Ntuple, 11, fPreDetectorEnergy);
    // analysisManager->FillNtupleIColumn(reacs6Ntuple, 12, Detector);
    // analysisManager->FillNtupleDColumn(reacs6Ntuple, 13, fGlobalTime);
    // analysisManager->FillNtupleIColumn(reacs6Ntuple, 14, DetectorSlice);



    // analysisManager->AddNtupleRow(reacs6Ntuple);


    // analysisManager->FillH2(analysisManager->GetFirstH2Id()+5,inDetDeltaD,DetectorSlice);
    // analysisManager->FillH2(analysisManager->GetFirstH2Id()+7,fKineticEnergy,DetectorSlice);

    // if(fEdep >= fGunEnergy)
    // {
    analysisManager->FillH2(analysisManager->GetFirstH2Id()+15,fGlobalTime,fEdep);
    // }

    // if(issetCl35Reaction())
    // {
      // analysisManager->FillH1(analysisManager->GetFirstH1Id()+17,fGunEnergy);
      // analysisManager->FillH2(analysisManager->GetFirstH2Id(),fGunEnergy,fKineticEnergy);

    // }
    // if((hasH1 && hasS35) or (hasHe4 && hasP32))
    // {
    //   analysisManager->FillH1(analysisManager->GetFirstH1Id()+17,fGunEnergy);
    // }

    if(issethasCl35() or issethasLi6())
    // if(true)
    {

      if(hasH1 && hasS35)
      {
        // printf("Cl35(n,p)\n");
        analysisManager->FillH3(analysisManager->GetFirstH3Id(),fXPosition,fYPosition,fZPosition);
        analysisManager->FillH1(analysisManager->GetFirstH1Id()+17,fGunEnergy);
        analysisManager->FillH1(analysisManager->GetFirstH1Id()+12,fZPosition);
      }
      else if(hasHe4 and hasP32)
      {
        // printf("Cl35(n,alpha)\n");
        analysisManager->FillH3(analysisManager->GetFirstH3Id(),fXPosition,fYPosition,fZPosition);
        analysisManager->FillH1(analysisManager->GetFirstH1Id()+21,fGunEnergy);
        analysisManager->FillH1(analysisManager->GetFirstH1Id()+12,fZPosition);
      }

      // analysisManager->FillH1(analysisManager->GetFirstH1Id()+12,fZPosition);

      // analysisManager->FillH1(analysisManager->GetFirstH1Id()+17,fGunEnergy);
      analysisManager->FillH2(analysisManager->GetFirstH2Id(),fGunEnergy,fKineticEnergy);
      if(false)
      {
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 0, fEdep);
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 1, fLstep);
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 2, fDeltat);
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 3, fGunEnergy);
      analysisManager->FillNtupleIColumn(reacs6Ntuple, 4, Z);
      analysisManager->FillNtupleIColumn(reacs6Ntuple, 5, A);
      // if(fKineticEnergy > 0)
      // {
      //   analysisManager->FillNtupleDColumn(reacs6Ntuple, 6, fKineticEnergy);
      // }
      // else
      // {
      //   analysisManager->FillNtupleDColumn(reacs6Ntuple, 6, fGunEnergy);
      // } 
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 6, fKineticEnergy);
      analysisManager->FillNtupleIColumn(reacs6Ntuple, 7, fSteps);
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 8, fXPosition);
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 9, fYPosition);
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 10, fZPosition);
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 11, fPreDetectorEnergy);
      analysisManager->FillNtupleIColumn(reacs6Ntuple, 12, Detector);
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 13, fGlobalTime);
      analysisManager->FillNtupleIColumn(reacs6Ntuple, 14, DetectorSlice);
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 15, inDetDeltaD);
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 16, inDetDeltaT);
      analysisManager->FillNtupleDColumn(reacs6Ntuple, 17, fNonIonEdep);


      analysisManager->AddNtupleRow(reacs6Ntuple);
      }

      analysisManager->FillH1(analysisManager->GetFirstH1Id()+4,inDetDeltaT);
      analysisManager->FillH1(analysisManager->GetFirstH1Id()+5,inDetDeltaD);
      if(issetGoodPreDetEn())
      {
        analysisManager->FillH1(analysisManager->GetFirstH1Id()+8,fPreDetectorEnergy);
      }


      analysisManager->FillH2(analysisManager->GetFirstH2Id(),fGunEnergy,fPreDetectorEnergy);
      analysisManager->FillH2(analysisManager->GetFirstH2Id()+2,fXPosition,fYPosition);

      analysisManager->FillH2(analysisManager->GetFirstH2Id()+5,DetectorSlice,inDetDeltaD);
      analysisManager->FillH2(analysisManager->GetFirstH2Id()+7,fKineticEnergy,DetectorSlice);
      analysisManager->FillH2(analysisManager->GetFirstH2Id()+9,fGunEnergy,inDetDeltaD);
      analysisManager->FillH2(analysisManager->GetFirstH2Id()+11,fGunEnergy,inDetDeltaT);
      analysisManager->FillH2(analysisManager->GetFirstH2Id()+13,fZPosition,inDetDeltaD);



      // analysisManager->FillH2(analysisManager->GetFirstH2Id()+8,fDeltat,fGunEnergy);
      // analysisManager->FillH2(analysisManager->GetFirstH2Id()+9,fGunEnergy,inDetDeltaT);
      // analysisManager->FillH2(analysisManager->GetFirstH2Id()+10,fGunEnergy,inDetDeltaD);
      // analysisManager->FillH3(analysisManager->GetFirstH3Id(),particleX,particleY,particleZ);
    }
  }

  if(Detector == 7 && fDetConstruction->GetUseC7LYC())
    {
      // Fill the 3D histogram with the start position of the particle
      if(isC6LYC)
      {
        analysisManager->FillH1(analysisManager->GetFirstH1Id()+20,fGunEnergy);
      }

      // analysisManager->FillNtupleDColumn(reacs7Ntuple, 0, fEdep);
      // analysisManager->FillNtupleDColumn(reacs7Ntuple, 1, fLstep);
      // analysisManager->FillNtupleDColumn(reacs7Ntuple, 2, fDeltat);
      // analysisManager->FillNtupleDColumn(reacs7Ntuple, 3, fGunEnergy);
      // analysisManager->FillNtupleIColumn(reacs7Ntuple, 4, Z);
      // analysisManager->FillNtupleIColumn(reacs7Ntuple, 5, A);
      // analysisManager->FillNtupleDColumn(reacs7Ntuple, 6, fKineticEnergy);
      // analysisManager->FillNtupleIColumn(reacs7Ntuple, 7, fSteps);
      // analysisManager->FillNtupleDColumn(reacs7Ntuple, 8, fXPosition);
      // analysisManager->FillNtupleDColumn(reacs7Ntuple, 9, fYPosition);
      // analysisManager->FillNtupleDColumn(reacs7Ntuple, 10, fZPosition);
      // analysisManager->FillNtupleDColumn(reacs7Ntuple, 11, fPreDetectorEnergy);
      // analysisManager->FillNtupleIColumn(reacs7Ntuple, 12, Detector);
      // analysisManager->FillNtupleDColumn(reacs7Ntuple, 13, fGlobalTime);
      // analysisManager->FillNtupleIColumn(reacs7Ntuple, 14, DetectorSlice);


      // analysisManager->AddNtupleRow(reacs7Ntuple);

      analysisManager->FillH1(analysisManager->GetFirstH1Id()+11,fGlobalTime);
      // analysisManager->FillH1(analysisManager->GetFirstH1Id()+15,fGunEnergy);

      analysisManager->FillH2(analysisManager->GetFirstH2Id()+4,fGunEnergy,fTheta);

      // analysisManager->FillH2(analysisManager->GetFirstH2Id()+6,inDetDeltaD,DetectorSlice);
      // analysisManager->FillH2(analysisManager->GetFirstH2Id()+8,fKineticEnergy,DetectorSlice);
      
      // if(fEdep >= fGunEnergy)
      // {
      analysisManager->FillH2(analysisManager->GetFirstH2Id()+16,fGlobalTime,fEdep);
      // }

      // if(issetCl35Reaction())
      // {
        // analysisManager->FillH1(analysisManager->GetFirstH1Id()+18,fGunEnergy);
        // analysisManager->FillH2(analysisManager->GetFirstH2Id()+1,fGunEnergy,fKineticEnergy);

      // }

      // if((hasH1 && hasS35))
      // {
      //   analysisManager->FillH3(analysisManager->GetFirstH3Id(),fXPosition,fYPosition,fZPosition);
      //   analysisManager->FillH1(analysisManager->GetFirstH1Id()+18,fGunEnergy);
      //   analysisManager->FillH1(analysisManager->GetFirstH1Id()+13,fZPosition);
      // }
      // else if((hasHe4 && hasP32))
      // {
      //   analysisManager->FillH3(analysisManager->GetFirstH3Id(),fXPosition,fYPosition,fZPosition);
      //   analysisManager->FillH1(analysisManager->GetFirstH1Id()+22,fGunEnergy);
      //   analysisManager->FillH1(analysisManager->GetFirstH1Id()+13,fZPosition);
      // }


      // if(fEdep >= (0.1 * fGunEnergy))
      if(issethasCl35() or issethasLi6())
      // if(true)
      {

        if((hasH1 && hasS35))
        {
          analysisManager->FillH3(analysisManager->GetFirstH3Id(),fXPosition,fYPosition,fZPosition);
          analysisManager->FillH1(analysisManager->GetFirstH1Id()+18,fGunEnergy);
          analysisManager->FillH1(analysisManager->GetFirstH1Id()+13,fZPosition);
        }
        else if((hasHe4 && hasP32))
        {
          analysisManager->FillH3(analysisManager->GetFirstH3Id(),fXPosition,fYPosition,fZPosition);
          analysisManager->FillH1(analysisManager->GetFirstH1Id()+22,fGunEnergy);
          analysisManager->FillH1(analysisManager->GetFirstH1Id()+13,fZPosition);
        }

        // printf("total in det = %f\n",(double)inDetDeltaD);
        // analysisManager->FillH1(analysisManager->GetFirstH1Id()+13,fZPosition);


        // analysisManager->FillH1(analysisManager->GetFirstH1Id()+18,fGunEnergy);
        analysisManager->FillH2(analysisManager->GetFirstH2Id()+1,fGunEnergy,fKineticEnergy);



        if(false)
        {
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 0, fEdep);
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 1, fLstep);
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 2, fDeltat);
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 3, fGunEnergy);
        analysisManager->FillNtupleIColumn(reacs7Ntuple, 4, Z);
        analysisManager->FillNtupleIColumn(reacs7Ntuple, 5, A);
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 6, fKineticEnergy);
        analysisManager->FillNtupleIColumn(reacs7Ntuple, 7, fSteps);
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 8, fXPosition);
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 9, fYPosition);
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 10, fZPosition);
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 11, fPreDetectorEnergy);
        analysisManager->FillNtupleIColumn(reacs7Ntuple, 12, Detector);
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 13, fGlobalTime);
        analysisManager->FillNtupleIColumn(reacs7Ntuple, 14, DetectorSlice);
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 15, inDetDeltaD);
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 16, inDetDeltaT);
        analysisManager->FillNtupleDColumn(reacs7Ntuple, 17, fNonIonEdep);


        analysisManager->AddNtupleRow(reacs7Ntuple);
        }

        analysisManager->FillH1(analysisManager->GetFirstH1Id()+6,inDetDeltaT);
        analysisManager->FillH1(analysisManager->GetFirstH1Id()+7,inDetDeltaD);
        if(issetGoodPreDetEn())
        {
          analysisManager->FillH1(analysisManager->GetFirstH1Id()+9,fPreDetectorEnergy);
        }


        // analysisManager->FillH2(analysisManager->GetFirstH2Id()+1,fGunEnergy,fPreDetectorEnergy);
        analysisManager->FillH2(analysisManager->GetFirstH2Id()+3,fXPosition,fYPosition);

        analysisManager->FillH2(analysisManager->GetFirstH2Id()+6,DetectorSlice,inDetDeltaD);
        analysisManager->FillH2(analysisManager->GetFirstH2Id()+8,fKineticEnergy,DetectorSlice);

        if(!(hasS35 && hasH1) and !(hasP32 && hasHe4))
        {
          analysisManager->FillH2(analysisManager->GetFirstH2Id()+10,fGunEnergy,inDetDeltaD);
          analysisManager->FillH2(analysisManager->GetFirstH2Id()+12,fGunEnergy,inDetDeltaT);
        }
        analysisManager->FillH2(analysisManager->GetFirstH2Id()+14,fZPosition,inDetDeltaD);


        // auto momentum = event->GetPrimaryVertex()->GetPrimary()->GetMomentumDirection();
        // if(fGenEnergy < 12)
        // {        
        // }

        // analysisManager->FillH2(analysisManager->GetFirstH2Id()+8,fDeltat,fGunEnergy);
        // analysisManager->FillH2(analysisManager->GetFirstH2Id()+9,fGunEnergy,inDetDeltaT);
        // analysisManager->FillH2(analysisManager->GetFirstH2Id()+10,fGunEnergy,inDetDeltaD);
        // analysisManager->FillH3(analysisManager->GetFirstH3Id(),particleX,particleY,particleZ);
      }
    }

  // else
  // {
  //     analysisManager->SetNtupleActivation(analysisManager->GetFirstNtupleId()+1,false);
  //     analysisManager->SetNtupleActivation(analysisManager->GetFirstNtupleId(),true);

  //     analysisManager->FillNtupleDColumn(0, fEdep);
  //     analysisManager->FillNtupleDColumn(1, fLstep);
  //     analysisManager->FillNtupleDColumn(2, new_Deltat);
  //     analysisManager->FillNtupleDColumn(3, fGunEnergy);
  //     analysisManager->FillNtupleDColumn(4, Z);
  //     analysisManager->FillNtupleDColumn(5, A);
  //     analysisManager->FillNtupleDColumn(6, fKineticEnergy);
  //     analysisManager->FillNtupleIColumn(7, fSteps);
  // }
  // accumulate statistics in run action

}

bool EventAction::ClearVectors()
{
  LstepVector.clear();
  DeltaTVector.clear();
  endKEVector.clear();
  EdepVector.clear();
  particleZVector.clear();
  particleAVector.clear();
  detectorSliceVector.clear();
  pos_x.clear();
  pos_y.clear();
  pos_z.clear();
  SecondaryNeutrons.clear();
  
  return true;
}

bool EventAction::ClearVariables()
{
    isC7LYC = false;
    isC6LYC = false;
    hasP32 = false;
    hasH1 = false;
    hasH3 = false;
    hasHe4 = false;
    hasS35 = false;
    hasCl35 = false;
    hasLi6 = false;
    Cl35reac = false;
    goodPreDet = false;
    fEventID = 0;
    fTheta = 0;
    fEdep = 0;
    fNonIonEdep = 0;
    fLstep = 0;
    fDeltat = 0;
    fGlobalTime = 0;
    inDetDeltaT = 0;
    inDetDeltaD = 0;
    fGunEnergy = 0;
    fKineticEnergy = 0;
    particleX = 0;
    particleY = 0;
    particleZ = 0;
    fXPosition = 0;
    fYPosition = 0;
    fZPosition = 0;
    fdummyXPosition = 0;
    fdummyYPosition = 0;
    fdummyZPosition = 0;
    fDummy = false;
    fPreDetectorEnergy = 0;
    fSteps = 0;
    A = 0; 
    Z = 0;
    Detector = 0;
    DetectorSlice = -1;

    return true;
}

std::vector<double> EventAction::GetNeutronData(int neutron)
{
  return NeutronsData.at(neutron);
}

void EventAction::setNeutronParam(int low, int high)
{
  neutron->param(std::uniform_int_distribution<int>::param_type(low,high));
  return;
}

int EventAction::getNeutronValue()
{
  int neut = neutron->operator()(*generator);
  return neut;
}

void EventAction::setEParam(double low, double high)
{
  Ebin->param(std::uniform_real_distribution<double>::param_type(low,high));
  return;
}

double EventAction::getEValue()
{
  double E = Ebin->operator()(*generator);
  return E;
}

void EventAction::setAngleParam(double low, double high)
{
  AngleBin->param(std::uniform_real_distribution<double>::param_type(low,high));
  return;
}

double EventAction::getAngleValue()
{
  double Ang = AngleBin->operator()(*generator);
  return Ang;
}

void EventAction::setZedParam(double low, double high)
{
  ZBin->param(std::uniform_real_distribution<double>::param_type(low,high));
  return;
}

double EventAction::getZedValue()
{
  double Z = ZBin->operator()(*generator);
  return Z;
}

void EventAction::setRhoParam(double low, double high)
{
  Rho->param(std::uniform_real_distribution<double>::param_type(low,high));
  return;
}

double EventAction::getRhoValue()
{
  double rho = Rho->operator()(*generator);
  return rho;
}

void EventAction::setPhiParam(double low, double high)
{
  Phi->param(std::uniform_real_distribution<double>::param_type(low,high));
  return;
}

double EventAction::getPhiValue()
{
  double phi = Phi->operator()(*generator);
  return phi;
}

void EventAction::setEnintParam(int low, int high)
{
  En->param(std::uniform_int_distribution<int>::param_type(low,high));
  return;
}

int EventAction::getEnintValue()
{
  int en = En->operator()(*generator);
  return en;
}

void EventAction::setAngintParam(int low, int high)
{
  Ang->param(std::uniform_int_distribution<int>::param_type(low,high));
  return;
}

int EventAction::getAngintValue()
{
  int ang = Ang->operator()(*generator);
  return ang;
}

void EventAction::setZhintParam(int low, int high)
{
  Zh->param(std::uniform_int_distribution<int>::param_type(low,high));
  return;
}

int EventAction::getZhintValue()
{
  int zh = Zh->operator()(*generator);
  return zh;
}

double EventAction::GetEnergyAngleZ_Value(int row, int bin)
{
  return EnergyAngleZ_bins.at(row).at(bin);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
