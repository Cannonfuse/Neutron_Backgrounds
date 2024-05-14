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
#include <cmath>

// class RunAction;

/// Event action class
///

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddEdep(G4double edep);
    void setFKE(G4double fke);
    void AddLstep(G4double lstep);
    G4double retLstep();
    // void AdddeltaE(G4double deltae);
    void AdddeltaT(G4double deltat);
    void AddInDetDeltaT(G4double deltat) {inDetDeltaT += deltat;};
    void AddInDetDeltaD(G4double deltad) {inDetDeltaD += deltad;};
    void SetParticleX(G4double X) {particleX = X;};
    void SetParticleY(G4double Y) {particleY = Y;};
    void SetParticleZ(G4double Z) {particleZ = Z;};
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
    void AddToAVector(int A);
    std::vector<int>& GetZVector();
    void AddToZVector(int Z);

    std::vector<G4double>& Get_pos_x_vector();
    void AddTo_pos_x(G4double x);
    std::vector<G4double>& Get_pos_y_vector();
    void AddTo_pos_y(G4double y);
    std::vector<G4double>& Get_pos_z_vector();
    void AddTo_pos_z(G4double z);


    void sethasCl35() {hasCl35 = true;}
    G4bool issethasCl35() {return hasCl35;}

    void sethasLi6() {hasLi6 = true;}
    G4bool issethasLi6() {return hasLi6;}

    void setGoodPreDetEn() {goodPreDet = true;};
    G4bool issetGoodPreDetEn() {return goodPreDet;};


    friend G4double energyRes(G4double edep);
    friend G4double detEnergyResponse(G4double edep);
    friend G4double calcTime(G4double start_KE, G4double end_KE, G4double DIST);


  private:
    bool hasCl35{false};
    bool hasLi6{false};
    bool goodPreDet{false};
    G4double     fEdep{0};
    G4double     fLstep{0};
    // G4double     fDeltae;
    G4double     fDeltat{0};
    G4double     inDetDeltaT{0};
    G4double     inDetDeltaD{0};
    G4double     fGunEnergy{0};
    G4double     fKineticEnergy{0};
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
    std::vector<G4double> LstepVector{NULL};
    std::vector<G4double> DeltaTVector{NULL};
    std::vector<G4double> endKEVector{NULL};
    std::vector<G4double> EdepVector{NULL};
    std::vector<int> particleZVector{NULL};
    std::vector<int> particleAVector{NULL};
    std::vector<G4double> pos_x{NULL};
    std::vector<G4double> pos_y{NULL};
    std::vector<G4double> pos_z{NULL};







};

// inline functions

inline void EventAction::AddEdep(G4double edep)
{
  fEdep += edep; 
}

inline void EventAction::setFKE(G4double fke)
{
  fKineticEnergy = fke; 
}

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

inline void EventAction::AddToAVector(int A)
{
  particleAVector.push_back(A);
}

inline std::vector<int>& EventAction::GetZVector()
{
  return particleZVector;
}

inline void EventAction::AddToZVector(int Z)
{
  particleZVector.push_back(Z);
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

    
