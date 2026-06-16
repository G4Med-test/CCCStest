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
#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4PhysListFactory.hh"
#include "G4VPhysicsConstructor.hh"

#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4ProcessManager.hh"
#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4StepLimiter.hh"

#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsINCLXX.hh"

#include "G4IonElasticPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonQMDPhysics.hh"
#include "G4IonINCLXXPhysics.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

/////////////////////////////////////////////////////////////////////////////
PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 10.3*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  hadronic = false;
  fStepMaxProcess=nullptr;

  pMessenger = new PhysicsListMessenger(this);

  SetVerboseLevel(1);
  
  // EM physics
  emPhysicsList = new G4EmStandardPhysics_option3(1);
  emName = G4String("emstandard_opt3");

  // Radioactive decay
  //raddecayList = new G4RadioactiveDecayPhysics();

  // Decay physics and all particles
  decPhysicsList = new G4DecayPhysics();
  
}

/////////////////////////////////////////////////////////////////////////////
PhysicsList::~PhysicsList()
{
  delete pMessenger;
  delete emPhysicsList;
  delete decPhysicsList; 
  delete raddecayList;
  for(size_t i=0; i<hadronPhys.size(); i++) {delete hadronPhys[i];}
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::AddPackage(const G4String& name)
{

  G4PhysListFactory factory;
  G4VModularPhysicsList* phys =factory.GetReferencePhysList(name);
  G4int i=0;
  const G4VPhysicsConstructor* elem= phys->GetPhysics(i);
  G4VPhysicsConstructor* tmp = const_cast<G4VPhysicsConstructor*> (elem);
  while (elem !=0)
	{
	  RegisterPhysics(tmp);
	  elem= phys->GetPhysics(++i) ;
	  tmp = const_cast<G4VPhysicsConstructor*> (elem);
	}
 
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::ConstructParticle()
{
  decPhysicsList->ConstructParticle();
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::ConstructProcess()
{

  // transportation
  //
  AddTransportation();

  // electromagnetic physics list
  //
  //emPhysicsList->ConstructProcess();
  //em_config.AddModels();

  // decay physics list
  //
  //decPhysicsList->ConstructProcess();
  //raddecayList->ConstructProcess();

  // hadronic physics lists
  for(size_t i=0; i<hadronPhys.size(); ++i) {
    hadronPhys[i]->ConstructProcess();
  }

  // step limitation (as a full process)
  // It is switched off
  //AddStepMax();
  
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::AddPhysicsList(const G4String& name)
{ 
  /////////////////////////////////////////////////////////////////////////////
  //   ELECTROMAGNETIC MODELS
  /////////////////////////////////////////////////////////////////////////////
  if(verboseLevel > 0) {
    G4cout << "### AddPhysicsList: <" << name << ">" << G4endl;
  }
  if (name == "standard_opt3" && name != emName) {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option3();

  } else if (name == "LowE_Livermore" && name != emName) {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmLivermorePhysics();

  } else if (name == "LowE_Penelope" && name != emName) {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmPenelopePhysics();

  } else if (name == "QGSP_BIC_EMY") {
    AddPhysicsList("standard_opt3");
    hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC());
    hadronPhys.push_back( new G4EmExtraPhysics());
    hadronPhys.push_back( new G4HadronElasticPhysics());
    hadronPhys.push_back( new G4StoppingPhysics());
    hadronPhys.push_back( new G4IonPhysics());
    hadronPhys.push_back( new G4IonElasticPhysics());

  } else if (name == "QGSP_BIC_HP_EMY") {
    AddPhysicsList("standard_opt3");
    hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC_HP());
    hadronPhys.push_back( new G4EmExtraPhysics());
    hadronPhys.push_back( new G4HadronElasticPhysicsHP());
    hadronPhys.push_back( new G4StoppingPhysics());
    hadronPhys.push_back( new G4IonPhysics());
    hadronPhys.push_back( new G4IonElasticPhysics());

  } else if (name == "QMD") {
    AddPhysicsList("standard_opt3");

    hadronPhys.push_back( new G4HadronElasticPhysicsHP());
    hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC_HP());
    hadronPhys.push_back( new G4EmExtraPhysics());
    hadronPhys.push_back( new G4StoppingPhysics());
    hadronPhys.push_back( new G4IonQMDPhysics());
    hadronPhys.push_back( new G4IonElasticPhysics());// in a second stage


  } else if (name == "INCL") {
    emName = name;
    hadronPhys.push_back( new G4HadronPhysicsINCLXX());
    hadronPhys.push_back( new G4EmExtraPhysics());
    hadronPhys.push_back( new G4HadronElasticPhysics());
    hadronPhys.push_back( new G4StoppingPhysics());
    hadronPhys.push_back( new G4IonINCLXXPhysics());
    hadronPhys.push_back( new G4IonElasticPhysics());// in a second stage

  } else if (name == "QBBC") {
    emName = name;
    hadronPhys.push_back( new G4HadronInelasticQBBC());
    hadronPhys.push_back( new G4EmExtraPhysics());
    hadronPhys.push_back( new G4HadronElasticPhysicsXS());
    hadronPhys.push_back( new G4StoppingPhysics());
    hadronPhys.push_back( new G4IonBinaryCascadePhysics());
    hadronPhys.push_back( new G4IonElasticPhysics());// in a second stage

  } else { 
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
	   << " is not defined"
	   << G4endl;
  }
}

/////////////////////////////////////////////////////////////////////////////

void PhysicsList::AddStepMax()
{
  // Step limitation seen as a process

  fStepMaxProcess = new G4StepLimiter();

  auto partIt = GetParticleIterator();
  partIt->reset();
  while ((*partIt)()){
    G4ParticleDefinition* particle = partIt->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (fStepMaxProcess->IsApplicable(*particle) && pmanager)
      {
	pmanager ->AddDiscreteProcess(fStepMaxProcess);
      }
  }
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::SetCuts()
{

  if (verboseLevel >0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  G4double lowLimit = 1. * keV;
  G4double highLimit = 100. * GeV;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  
  G4String regionName = "DetectorLog";
  G4Region* region = G4RegionStore::GetInstance()->GetRegion(regionName);
  G4ProductionCuts* cuts = new G4ProductionCuts ;
  G4double regionCut = 10.*mm;
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("gamma"));
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("e-"));
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("e+"));
  region -> SetProductionCuts(cuts);

  if (verboseLevel>0) DumpCutValuesTable();
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}


