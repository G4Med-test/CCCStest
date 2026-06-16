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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 105744 2017-08-16 13:13:16Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4HadronicProcess.hh"
                           
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
: G4UserSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
 Run* run 
   = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
         
  // count processes
  // 
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  G4VProcess* process   = 
                   const_cast<G4VProcess*>(endPoint->GetProcessDefinedStep());
  run->CountProcesses(process);
  
  // check that an real interaction occured (eg. not a transportation)
  G4StepStatus stepStatus = endPoint->GetStepStatus();
  G4bool transmit = (stepStatus==fGeomBoundary || stepStatus==fWorldBoundary);
  if (transmit) return;
                      
  //real processes : sum track length
  //
  G4double stepLength = aStep->GetStepLength();
  run->SumTrack(stepLength);
  
  //
  //secondaries
  //
  
  const G4StepPoint* ePoint = aStep->GetPostStepPoint();
  const std::vector<const G4Track*>* secondary 
                                    = aStep->GetSecondaryInCurrentStep();
  bool isTotalCCCS=true;
  for (size_t lp=0; lp<(*secondary).size(); lp++) {
    G4ParticleDefinition* particle = (*secondary)[lp]->GetDefinition(); 
    G4String name   = particle->GetParticleName();
    G4String type   = particle->GetParticleType();      
    G4double energy = (*secondary)[lp]->GetKineticEnergy();
    run->ParticleCount(name,energy);
    if(particle->GetAtomicMass()==0) continue;
    if(((*secondary)[lp]->GetMomentum().x()/particle->GetAtomicMass())>600){
      run->ParticleCountCCCS(name,energy);
      if(name[0]=='C') isTotalCCCS=false;
    }
  }
  if(isTotalCCCS){
    if((ePoint->GetMomentum().x()/aStep->GetTrack()->GetDefinition()->GetAtomicMass())<600) run->AddCCCSTotalCount();
  }  

  // kill event after first interaction
  //
  if((ePoint->GetMomentum().x()/aStep->GetTrack()->GetDefinition()->GetAtomicMass())>600) return;
  G4RunManager::GetRunManager()->AbortEvent();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

