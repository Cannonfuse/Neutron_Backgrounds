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
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <random>
#include <vector>
#include <cmath>

class CLYCDetectorConstruction;

/// Event action class
///

class EventAction : public G4UserEventAction
{
  public:
    EventAction(CLYCDetectorConstruction* detConstruction);
    EventAction(CLYCDetectorConstruction* detConstruction, const G4bool usedists, 
                const std::vector<std::vector<G4bool>> energyangledist,
                const std::vector<std::vector<G4bool>> energyzdist,
                const std::vector<std::vector<double>> energyangzbins);
    EventAction(CLYCDetectorConstruction* detConstruction, const G4bool useneutrons, 
                const std::vector<std::vector<double>> neutronsdata);
    virtual ~EventAction();

    bool ClearVectors();
    bool ClearVariables();

    void SetEnergyAngle_dist(std::vector<std::vector<G4bool>> dist) {EnergyAngle_dist = dist;};
    void SetEnergyZ_dist(std::vector<std::vector<G4bool>> dist) {EnergyZ_dist = dist;};
    void SetEnergyAngleZ_bins(std::vector<std::vector<double>> dist) {EnergyAngleZ_bins = dist;};
    void SetNeutronsData(std::vector<std::vector<double>> dist) {NeutronsData = dist;};
    void SetUseDists(G4bool usedists) {UseDists = usedists;};
    void SetUseNeutronsData(G4bool useneutronsdata) {UseNeutronsData = useneutronsdata;};


    std::vector<std::vector<G4bool>>& GetEnergyAngle_dist() {return EnergyAngle_dist;};
    G4bool EnergyAngleValue(int enint, int angint) {return EnergyAngle_dist[angint][enint];};
    int GetEnergySize() {return EnergyAngle_dist.at(0).size();};
    int GetAngleSize() {return EnergyAngle_dist.size();};
    std::vector<std::vector<G4bool>>& GetEnergyZ_dist() {return EnergyZ_dist;};
    G4bool EnergyZValue(int enint, int zhint) {return EnergyZ_dist[zhint][enint];};
    int GetZedSize() {return EnergyZ_dist.size();};
    std::vector<std::vector<double>>& GetEnergyAngleZ_bins() {return EnergyAngleZ_bins;};
    double GetEnergyAngleZ_Value(int row, int bin);

    std::vector<std::vector<double>>& GetNeutronsData() {return NeutronsData;};
    int GetNeutronsDataSize() {return NeutronsData.size();};
    std::vector<double> GetNeutronData(int neutron);
    G4bool GetUseDists() {return UseDists;};
    G4bool GetUseNeutronsData() {return UseNeutronsData;};

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddEdep(G4double edep) {fEdep += edep;};
    void AddNonIonEdep(G4double edep) {fNonIonEdep += edep;};

    void setSKE(G4double ske) {fStartKineticEnergy = ske;};
    void setFKE(G4double fke) {fKineticEnergy = fke;};
    G4double getFKE() {return fKineticEnergy;};
    void AddLstep(G4double lstep);
    G4double retLstep();
    // void AdddeltaE(G4double deltae);
    void AdddeltaT(G4double deltat);
    void SetGlobalTime(G4double atime) {fGlobalTime = atime;}
    void AddInDetDeltaT(G4double deltat) {inDetDeltaT += deltat;};
    void AddInDetDeltaD(G4double deltad) {inDetDeltaD += deltad;};
    void SetParticleX(G4double xray) {particleX = xray;};
    void SetParticleY(G4double yankee) {particleY = yankee;};
    void SetParticleZ(G4double zulu) {particleZ = zulu;};
    void SetGunEnergy(G4double gamma) {fGunEnergy = gamma;};
    void SetTheta(G4double tango) {fTheta = tango;};
    void SetSlice(G4int sigma) {DetectorSlice = sigma;};
    void AddStep();
    void SetA(G4int e_a);
    void SetZ(G4int e_z);
    void SetDetector(G4int detectornumber);
    void SetPosX(G4double posX);
    void SetPosY(G4double posY);
    void SetPosZ(G4double posZ);
    void SetDummyPosX(G4double posX);
    void SetDummyPosY(G4double posY);
    void SetDummyPosZ(G4double posZ);
    
    G4int GetDetector() {return Detector;};

    void SetDummy();
    bool GetDummy();

    void SetPreDetectorEnergy(G4double PreDetE);
    G4double retGunEnergy() {return fGunEnergy;};
    G4double retPreDetEnergy() {return fPreDetectorEnergy;};

    std::vector<G4double>& GetLstepVector();
    void AddToLstepVector(G4double lstep);
    std::vector<G4double>& GetDeltaTVector();
    void AddToDeltaTVector(G4double DT);
    std::vector<G4double>& GetEndKEVector();
    void AddToEndKEVector(G4double KE);
    std::vector<G4double>& GetEdepVector();
    void AddToEdepVector(G4double Edep);

    std::vector<int>& GetAVector();
    void AddToAVector(int alpha);
    std::vector<int>& GetZVector();
    void AddToZVector(int zulu);
    std::vector<int>& GetSliceVector();
    void AddToSliceVector(int stingray);


    std::vector<G4double>& Get_pos_x_vector();
    void AddTo_pos_x(G4double x);
    std::vector<G4double>& Get_pos_y_vector();
    void AddTo_pos_y(G4double y);
    std::vector<G4double>& Get_pos_z_vector();
    void AddTo_pos_z(G4double z);


    void sethasCl35() {hasCl35 = true;};
    void clearsethasCl35() {hasCl35 = false;};
    G4bool issethasCl35() {return hasCl35;};

    void setCl35Reaction() {Cl35reac = true;};
    G4bool issetCl35Reaction() {return Cl35reac;};


    void sethasLi6() {hasLi6 = true;}
    void clearsethasLi6() {hasLi6 = false;}
    G4bool issethasLi6() {return hasLi6;}

    void sethasH1() {hasH1 = true;}
    void sethasHe4() {hasHe4 = true;}
    void sethasS35() {hasS35 = true;}
    void sethasP32() {hasP32 = true;}

    void setGoodPreDetEn() {goodPreDet = true;};
    G4bool issetGoodPreDetEn() {return goodPreDet;};

    void setNeutronParam(int low, int high);
    int getNeutronValue();

    void setEParam(double low, double high);
    double getEValue();
    void setAngleParam(double low, double high);
    double getAngleValue();
    void setZedParam(double low, double high);
    double getZedValue();
    void setRhoParam(double low, double high);
    double getRhoValue();
    void setPhiParam(double low, double high);
    double getPhiValue();

    void setEnintParam(int low, int high);
    int getEnintValue();
    void setAngintParam(int low, int high);
    int getAngintValue();
    void setZhintParam(int low, int high);
    int getZhintValue();

    G4int getEventID() {return fEventID;};
    void setEventID(G4int value) {fEventID = value;};

    G4int getC6LYC_Slices() {return C6LYC_Slices;};
    G4int getC7LYC_Slices() {return C7LYC_Slices;};
    
    void setC6LYC() {isC6LYC = true;};
    void setC7LYC() {isC7LYC = true;};



    // CLYCDetectorConstruction* retDetConstruction() const {return fDetConstruction;};

    // friend G4double energyRes(G4double edep);
    // friend G4double detEnergyResponse(G4double edep);
    // friend G4double calcTime(G4double start_KE, G4double end_KE, G4double DIST);


  private:
  
    CLYCDetectorConstruction* fDetConstruction;

    bool isC7LYC;
    bool isC6LYC;
    bool hasCl35{false};
    bool hasH1{false};
    bool hasHe4{false};
    bool hasS35{false};
    bool hasP32{false};

    bool hasLi6{false};
    bool Cl35reac{false};
    bool Li6reac{false};
    bool goodPreDet{false};
    G4int     fEventID;
    G4double     fTheta{0};
    G4double     fEdep{0};
    G4double     fNonIonEdep{0};
    G4double     fLstep{0};
    // G4double     fDeltae;
    G4double     fDeltat{0};
    G4double      fGlobalTime{0};
    G4double     inDetDeltaT{0};
    G4double     inDetDeltaD{0};
    G4double     fGunEnergy{0};
    G4double     fKineticEnergy{0};
    G4double     fStartKineticEnergy{0};
    G4double      particleX{0};
    G4double      particleY{0};
    G4double      particleZ{0};
    G4double      fXPosition{0};
    G4double      fYPosition{0};
    G4double      fZPosition{0};
    G4double      fdummyXPosition{0};
    G4double      fdummyYPosition{0};
    G4double      fdummyZPosition{0};
    G4bool        fDummy{false};
    G4double      fPreDetectorEnergy{0};
    G4int fSteps{0};
    G4int A{0}; 
    G4int Z{0};
    G4int Detector{0};
    G4int DetectorSlice{0};
    G4int C6LYC_Slices{1};
    G4int C7LYC_Slices{1};
    std::vector<G4double> LstepVector{NULL};
    std::vector<G4double> DeltaTVector{NULL};
    std::vector<G4double> endKEVector{NULL};
    std::vector<G4double> EdepVector{NULL};
    std::vector<int> particleZVector{NULL};
    std::vector<int> particleAVector{NULL};
    std::vector<int> detectorSliceVector{NULL};
    std::vector<G4double> pos_x{NULL};
    std::vector<G4double> pos_y{NULL};
    std::vector<G4double> pos_z{NULL};

    
    std::vector<std::vector<G4bool>> EnergyAngle_dist;
    std::vector<std::vector<G4bool>> EnergyZ_dist;
    std::vector<std::vector<double>> EnergyAngleZ_bins;
    std::vector<std::vector<double>> NeutronsData;
    G4bool UseDists{false};
    G4bool UseNeutronsData{false};

    std::minstd_rand *generator;

    std::uniform_real_distribution<double> *Ebin;
    std::uniform_real_distribution<double> *AngleBin;
    std::uniform_real_distribution<double> *ZBin;
    std::uniform_real_distribution<double> *Rho;
    std::uniform_real_distribution<double> *Phi;


    std::uniform_int_distribution<int> *neutron;

    std::uniform_int_distribution<int> *En;
    std::uniform_int_distribution<int> *Ang;
    std::uniform_int_distribution<int> *Zh;









};

// inline functions

// inline void EventAction::setNeutronParam(int low, int high)
// {
//   neutron.param(std::uniform_int_distribution<int>::param_type(low,high));
// }

// inline int EventAction::getNeutronValue()
// {
//   return neutron(*generator);
// }

// inline void EventAction::AddEdep(G4double edep)
// {
//   fEdep += edep; 
// }

// inline void EventAction::setFKE(G4double fke)
// {
//   fKineticEnergy = fke; 
// }

inline void EventAction::AddLstep(G4double lstep)
{
  fLstep += lstep;
}

// inline void EventAction::AdddeltaE(G4double deltae)
// {
//   fDeltae += deltae; 
// }

inline void EventAction::AdddeltaT(G4double deltat)
{
  fDeltat += deltat;
}

inline void EventAction::SetPosX(G4double posX)
{
  fXPosition = posX;
}

inline void EventAction::SetPosY(G4double posY)
{
  fYPosition = posY;
}

inline void EventAction::SetPosZ(G4double posZ)
{
  fZPosition = posZ;
}

inline void EventAction::SetDummyPosX(G4double posX)
{
  fdummyXPosition = posX;
}

inline void EventAction::SetDummyPosY(G4double posY)
{
  fdummyYPosition = posY;
}

inline void EventAction::SetDummyPosZ(G4double posZ)
{
  fdummyZPosition = posZ;
}

inline void EventAction::SetDummy()
{
  fDummy = true;
}

inline bool EventAction::GetDummy()
{
  return fDummy;
}

inline void EventAction::SetPreDetectorEnergy(G4double PreDetE)
{
  fPreDetectorEnergy = PreDetE;
}

inline void EventAction::SetA(G4int e_a)
{
  A = e_a;
}

inline void EventAction::SetZ(G4int e_z)
{
  Z = e_z;
}

inline void EventAction::SetDetector(G4int detectornumber)
{
  Detector = detectornumber;  
}


inline void EventAction::AddStep()
{
  fSteps += 1;
}

// inline G4double EventAction::retGE()
// {
//   return fGunEnergy; 
// }

inline G4double EventAction::retLstep()
{
  return fLstep; 
}

inline std::vector<G4double>& EventAction::GetLstepVector()
{
  return LstepVector;
}

inline void EventAction::AddToLstepVector(G4double lstep)
{
  LstepVector.push_back(lstep);
}

inline std::vector<G4double>& EventAction::GetDeltaTVector()
{
  return DeltaTVector;
}

inline void EventAction::AddToDeltaTVector(G4double DT)
{
  DeltaTVector.push_back(DT);
}

inline std::vector<G4double>& EventAction::GetEndKEVector()
{
  return endKEVector;
}

inline void EventAction::AddToEndKEVector(G4double KE)
{
  endKEVector.push_back(KE);
}

inline std::vector<G4double>& EventAction::GetEdepVector()
{
  return EdepVector;
}

inline void EventAction::AddToEdepVector(G4double Edep)
{
  EdepVector.push_back(Edep);
}

inline std::vector<int>& EventAction::GetAVector()
{
  return particleAVector;
}

inline void EventAction::AddToAVector(int alpha)
{
  particleAVector.push_back(alpha);
}

inline std::vector<int>& EventAction::GetZVector()
{
  return particleZVector;
}

inline void EventAction::AddToZVector(int zulu)
{
  particleZVector.push_back(zulu);
}

inline std::vector<int>& EventAction::GetSliceVector()
{
  return detectorSliceVector;
}

inline void EventAction::AddToSliceVector(int stingray)
{
  detectorSliceVector.push_back(stingray);
}

inline std::vector<G4double>& EventAction::Get_pos_x_vector()
{
  return pos_x;
}

inline void EventAction::AddTo_pos_x(G4double x)
{
  pos_x.push_back(x);
}

inline std::vector<G4double>& EventAction::Get_pos_y_vector()
{
  return pos_y;
}

inline void EventAction::AddTo_pos_y(G4double y)
{
  pos_y.push_back(y);
}

inline std::vector<G4double>& EventAction::Get_pos_z_vector()
{
  return pos_z;
}

inline void EventAction::AddTo_pos_z(G4double z)
{
  pos_z.push_back(z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

    
