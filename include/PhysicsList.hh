//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//Original code adapted from Hadrontherapy example 
//Which was developed developed by:
//R. Calcagno(a), G.A.P. Cirrone(a)*, G.Cuttone(a), L. Pandola(a)
// F.Romano(a)*, A.Varisano(a)
//(a) Laboratori Nazionali del Sud of the INFN, Catania, Italy
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4EmConfigurator.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class G4StepLimiter;
class PhysicsListMessenger;

class PhysicsList: public G4VModularPhysicsList
{
public:

  PhysicsList();
  virtual ~PhysicsList();

  void ConstructParticle();

  void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);

  void AddPhysicsList(const G4String& name);
  void ConstructProcess();

  void AddStepMax();
  G4StepLimiter* GetStepMaxProcess() {return fStepMaxProcess;};
  void AddPackage(const G4String& name);

private:
  G4bool hadronic; 
  G4EmConfigurator em_config;

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;


  G4String                             emName;
  G4VPhysicsConstructor*               emPhysicsList;
  G4VPhysicsConstructor*               decPhysicsList;
  G4VPhysicsConstructor*               raddecayList;
  std::vector<G4VPhysicsConstructor*>  hadronPhys;

  G4StepLimiter* fStepMaxProcess;

  PhysicsListMessenger* pMessenger;
};

#endif
