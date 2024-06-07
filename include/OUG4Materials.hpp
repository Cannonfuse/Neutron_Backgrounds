#ifndef OUG4MATERIALS
#define OUG4MATERIALS 1

#include "globals.hh"

class G4Material;
class G4NistManager;

class OUG4Materials
{
    public:

    // Constructors
    OUG4Materials();
    OUG4Materials(G4NistManager* thenist);

    // Destructor
    ~OUG4Materials();

    // Material Get functions
    // Defined in place - they only return the material named
    G4Material* GetC6LYCMaterial_95() {return C6LYC_95;};
    G4Material* GetC7LYCMaterial_99() {return C7LYC_99;};
    G4Material* GetHAVARMaterial() {return HAVAR;};
    G4Material* Get_MuMetal_UNS14080Material() {return MuMetal_UNS14080;};
    G4Material* GetQuartzMaterial() {return Quartz;};

    // Set the G4NistManager
    void SetNistManager(G4NistManager* thenist);

    // Set the materials if G4NistManager has changed/updated
    void SetMaterials();

    private:

    // Nist Manager
    G4NistManager* NistManager{nullptr};

    // Material Set functions
    void SetC6LYCMaterial_95();
    void SetC7LYCMaterial_99();
    void SetHAVARMaterial();
    void Set_MuMetal_UNS14080Material();
    void SetQuartzMaterial();

    // Materials
    G4Material* MuMetal_UNS14080;
    G4Material* HAVAR;
    G4Material* C7LYC_99;
    G4Material* C6LYC_95;
    G4Material* Quartz;
    



    
};


#endif