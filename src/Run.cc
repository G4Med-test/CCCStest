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
/// \file Run.cc
/// \brief Implementation of the Run class
//
// $Id: Run.cc 71376 2013-06-14 07:44:50Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4ProcessTable.hh"
#include "G4HadronicProcessStore.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
: G4Run(),
  fDetector(det), fParticle(0), fEkin(0.),
  fTotalCount(0), fGammaCount(0),
  fSumTrack(0.), fSumTrack2(0.),
  cccsTotalCnt(0)
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Run::CountProcesses(G4VProcess* process) 
{
  G4String procName = process->GetProcessName();
  std::map<G4String,G4int>::iterator it = fProcCounter.find(procName);
  if ( it == fProcCounter.end()) {
    fProcCounter[procName] = 1;
  }
  else {
    fProcCounter[procName]++; 
  }
}                 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumTrack(G4double trackl)
{
  fTotalCount++;
  fSumTrack += trackl; fSumTrack2 += trackl*trackl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::ParticleCount(G4String name, G4double Ekin)
{
  std::map<G4String, ParticleData>::iterator it = fParticleDataMap.find(name);
  if ( it == fParticleDataMap.end()) {
    fParticleDataMap[name] = ParticleData(1, Ekin, Ekin, Ekin);
  }
  else {
    ParticleData& data = it->second;
    data.fCount++;
    data.fEmean += Ekin;
    //update min max
    G4double emin = data.fEmin;
    if (Ekin < emin) data.fEmin = Ekin;
    G4double emax = data.fEmax;
    if (Ekin > emax) data.fEmax = Ekin; 
  }   
}

void Run::ParticleCountCCCS(G4String name, G4double Ekin)
{
  std::map<G4String, ParticleData>::iterator it = fParticleDataMapCCCS.find(name);
  if ( it == fParticleDataMapCCCS.end()) {
    fParticleDataMapCCCS[name] = ParticleData(1, Ekin, Ekin, Ekin);
  }
  else {
    ParticleData& data = it->second;
    data.fCount++;
    data.fEmean += Ekin;
    //update min max
    G4double emin = data.fEmin;
    if (Ekin < emin) data.fEmin = Ekin;
    G4double emax = data.fEmax;
    if (Ekin > emax) data.fEmax = Ekin; 
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun(G4bool print) 
{
  G4int prec = 5, wid = prec + 2;  
  G4int dfprec = G4cout.precision(prec);
  
  //run condition
  //
  G4Material* material = fDetector->GetMaterial();
  G4double density = material->GetDensity();
   
  G4String Particle = fParticle->GetParticleName();    
  G4cout << "\n The run is " << numberOfEvent << " "<< Particle << " of "
         << G4BestUnit(fEkin,"Energy") << " through " 
         << G4BestUnit(fDetector->GetSize(),"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;

  if (numberOfEvent == 0) { G4cout.precision(dfprec);   return;}
             
  //frequency of processes
  //
 // G4cout << "\n Process calls frequency:" << G4endl;  
  G4int survive = 0;
  std::map<G4String,G4int>::iterator it;
  G4int total=0;
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
     G4String procName = it->first;
     G4int    count    = it->second;
     total+=count;
      // Uncomment for debugging purposes
     //G4cout << std::setw(15) << procName << "= " << std::setw(8) << count << "evts "  
     //       << "CrossSection=" << std::setw(10) << 1./(fSumTrack/count*density)/mm2*g << " mm2/g "
     //       << "(MeanFreePath=" << std::setw(10) << (fSumTrack/count*density)/g*cm2 << " g/cm2)" << G4endl;
     if (procName == "Transportation") survive = count;
  }
  
 /*
  G4cout << std::setw(15) << "total" << "= " << std::setw(8) << total << "evts "  
         << "CrossSection=" << std::setw(10) << 1./(fSumTrack/total*density)/mm2*g << " mm2/g "
         << "(MeanFreePath=" << std::setw(10) << (fSumTrack/total*density)/g*cm2 << " g/cm2)" << G4endl;
  G4cout << G4endl;
   
   */   
  if (survive > 0) {
    G4cout << "\n Nb of incident particles surviving after "
           << G4BestUnit(fDetector->GetSize(),"Length") << " of "
           << material->GetName() << " : " << survive << G4endl;
  }
  
  if (fTotalCount == 0) fTotalCount = 1;   //force printing anyway
  
  //
  //check cross section from G4HadronicProcessStore
  //
  //G4cout << " Verification: "
  //       << "crossSections from G4HadronicProcessStore:";
  if (print)
  {
  G4ProcessTable* processTable  = G4ProcessTable::GetProcessTable();
  G4HadronicProcessStore* store = G4HadronicProcessStore::Instance();
  G4double sumc1 = 0.0, sumc2 = 0.0; 
  if (material->GetNumberOfElements() == 1) {
    const G4Element* element = material->GetElement(0);
    for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
      G4String procName = it->first;
      G4VProcess* process = processTable->FindProcess(procName, fParticle);
      G4double xs1 =
      store->GetCrossSectionPerVolume(fParticle,fEkin,process,material);
      G4double massSigma = xs1/density;
      sumc1 += massSigma;      
      G4double xs2 =
      store->GetCrossSectionPerAtom(fParticle,fEkin,process,element,material);
      sumc2 += xs2;
      G4cout << "\n" << std::setw(20) << procName << "= "
             << G4BestUnit(massSigma, "Surface/Mass") << "\t"
             << G4BestUnit(xs2, "Surface");
      
    }             
    G4cout << "\n" << std::setw(20) << "total" << "= "
           << G4BestUnit(sumc1, "Surface/Mass") << "\t" 
           << G4BestUnit(sumc2, "Surface") << G4endl;  
  } else {
    for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
      G4String procName = it->first;
      G4VProcess* process = processTable->FindProcess(procName, fParticle);
      G4double xs =
      store->GetCrossSectionPerVolume(fParticle,fEkin,process,material);
      G4double massSigma = xs/density;
      sumc1 += massSigma;
      G4cout << "\n" << std::setw(15)  << procName << "  "
             << std::setw(26) << "CrossSection=" 
             << std::setw(10) << G4BestUnit(massSigma, "Surface/Mass");
    }             
    G4cout << "\n" << std::setw(15) << "total" << "  " 
           << std::setw(26) << "CrossSection="
           << std::setw(10) << G4BestUnit(sumc1, "Surface/Mass") << G4endl;  
  }
  G4cout << G4endl;
 }
 /*G4cout << "test for retrieve cross section" << G4endl;
 G4ParticleTable *ptable = G4ParticleTable::GetParticleTable();
 G4String pName[5]={"B11","Be9","C10","C12","C14"};
 for(int npart=0; npart<5; npart++){
   G4ParticleDefinition* lparticle = ptable->FindParticle(pName[npart]);
   for(int a=100; a<800; a+=50){
     double xs=store->GetCrossSectionPerVolume(lparticle,
                       a,
                       processTable->FindProcess("ionInelastic", lparticle),
                       material);
     G4cout << lparticle->GetParticleName() << " ionInelastic " << a << " " << xs << G4endl;
   }
 }*/

 // charge changing cross section
 std::map<G4String, G4int> CountsByIon;
 std::map<G4String,ParticleData>::iterator itn;
 //int totalCounts=0;
 for (itn = fParticleDataMapCCCS.begin(); itn != fParticleDataMapCCCS.end(); itn++) { 
    G4String name = itn->first;
    ParticleData data = itn->second;
    G4int count = data.fCount;
    G4String IonName;
    if(name.size()<2) IonName = name;
    else if(name.size()>4) 
          {
           if(name == "alpha" || name == "proton" || name == "neutron" || name =="deuteron" || name =="triton")  IonName=name;
           	else if (std::isdigit(name[1])) IonName=name[0]; 
           	  else if (std::isdigit(name[2])) IonName=name.substr(0,2);  
          	         else IonName = name;
         }
    else if(std::isdigit(name[1])) IonName=name[0];  
      else if (std::isdigit(name[2])) IonName=name.substr(0,2);  
          else IonName = name;
    
    CountsByIon[IonName]+=count;

    if (print)         
    G4cout << std::setw(8) << IonName << std::setw(13) << name << ": " << std::setw(7) << count
           << " MeanFreePath=" << std::setw(wid) << G4BestUnit(fSumTrack/count,"Length")
           << " CrossSection=" << std::setw(wid) << 1./(fSumTrack/count*(1*6.02e23/(18.01528*1e28))) << " [millibarn]"
           << G4endl;
 }

 G4cout << G4endl;
 G4cout << std::setw(8) << "Total"
        << std::setw(8) << cccsTotalCnt
        << std::setw(8) << 1./(fSumTrack/cccsTotalCnt*(1*6.02e23/(18.01528*1e28))) << " [millibarn]" << G4endl;
 G4cout << "-------------------------------------------------------" << G4endl;
 for(auto itrCntIon=CountsByIon.begin(); itrCntIon!=CountsByIon.end(); itrCntIon++){
   G4cout << std::setw(8) << itrCntIon->first  
          << std::setw(8) << itrCntIon->second 
          << std::setw(8) << 1./(fSumTrack/itrCntIon->second*(1*6.02e23/(18.01528*1e28))) << " [millibarn]" << G4endl;
 }
 
           
  //remove all contents in fProcCounter, fCount 
  fProcCounter.clear();
  fParticleDataMap.clear();
  fParticleDataMapCCCS.clear();
                          
  //restore default format         
  G4cout.precision(dfprec);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
