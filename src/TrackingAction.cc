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
/// \file electromagnetic/TestEm1/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class

#include "TrackingAction.hh"
using namespace std;

TrackingAction::TrackingAction( B4PrimaryGeneratorAction* prim )
  :G4UserTrackingAction(),
   fPrimary(prim)
{ }

void TrackingAction::PreUserTrackingAction(const G4Track*)
{
/*
  //debug  
  const G4DynamicParticle* dynamic = aTrack->GetDynamicParticle();
  G4double dynamCharge = dynamic->GetCharge();
  G4int occup          = dynamic->GetTotalOccupancy();
  G4double   pdgMass   = dynamic->GetParticleDefinition()->GetPDGMass();    
  G4double invarMass   = dynamic->Get4Momentum().m();  
  G4double dynamMass   = dynamic->GetMass();
  G4double dif1 = invarMass - pdgMass;
  G4double dif2 = dynamMass - pdgMass;
  
  G4cout
    << "\n  Begin of track :" 
    << "\n    charge= " <<  dynamCharge << "  occupancy= " << occup
    << "\n   pdgMass= " << G4BestUnit (pdgMass  , "Energy")    
///    << "\n invarMass= " << G4BestUnit (invarMass, "Energy")
///    << "   invar-pdg= " << G4BestUnit (dif1     , "Energy")
    << "\n dynamMass= " << G4BestUnit (dynamMass, "Energy")
    << "   dynam-pdg= " << G4BestUnit (dif2     , "Energy")
    << G4endl;          
*/             
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  //increase nb of processed tracks 
  //count nb of steps of this track
  G4int   nbSteps = aTrack->GetCurrentStepNumber();
  G4double Trleng = aTrack->GetTrackLength();
  
  // Run* run = static_cast<Run*>(
  //            G4RunManager::GetRunManager()->GetNonConstCurrentRun()); 
    
  // if (aTrack->GetDefinition()->GetPDGCharge() == 0.) {
  //   run->CountTraks0(1); 
  //   run->CountSteps0(nbSteps);
  
  // } else {
  //   run->CountTraks1(1); 
  //   run->CountSteps1(nbSteps);
  // }

  const G4Step* step = aTrack->GetStep();
  const G4DynamicParticle* particle = aTrack->GetDynamicParticle();
  const G4ParticleDefinition* particle_def = aTrack->GetParticleDefinition();
  
  //true and projected ranges for primary particle
  if (aTrack->GetTrackID() == 1) {
    //run->AddTrueRange(Trleng);
    G4ThreeVector vertex = fPrimary->GetParticleGun()->GetParticlePosition();
    G4ThreeVector position = aTrack->GetPosition() - vertex;      
    //run->AddProjRange(position.x());
    //run->AddTransvDev(position.y());
    //run->AddTransvDev(position.z());
      
    //    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4double beam_energy = aTrack->GetVertexKineticEnergy();
    //    cout << "Beam energy: " << beam_energy << endl;
    B4cEventAction::beam_energy = beam_energy;
    //    analysisManager->FillNtupleDColumn( 0, beam_energy ); // 0: Beam
    //analysisManager->FillH1(1,Trleng);
    //analysisManager->FillH1(2,(float)nbSteps);
    //    std::cout << beam_energy << endl;
  }
  //  Print( aTrack );

}

void TrackingAction::Print( const G4Track* aTrack )
{
  G4int   nbSteps = aTrack->GetCurrentStepNumber();
  G4double Trleng = aTrack->GetTrackLength();

  const G4Step* step = aTrack->GetStep();
  const G4DynamicParticle* particle = aTrack->GetDynamicParticle();
  const G4ParticleDefinition* particle_def = aTrack->GetParticleDefinition();

    //  G4ThreeVector vertex = fPrimary->GetParticleGun()->GetParticlePosition();
  // G4ThreeVector position = aTrack->GetPosition() - vertex;      
  cout << setw(5)  << aTrack->GetTrackID() << "->" // good
       << setw(10) << left << particle_def->GetParticleName() << "@" // good
       << "( " 
       << setw(5)  << aTrack->GetVertexPosition().X << ", " // useless
       << setw(5)  << aTrack->GetVertexPosition().Y << ", " // useless
       << setw(5)  << aTrack->GetVertexPosition().Z << ") " // useless
    
       // << "( " 
       // << setw(5)  << aTrack->GetPosition().X << ", " // useless
       // << setw(5)  << aTrack->GetPosition().Y << ", " // useless
       // << setw(5)  << aTrack->GetPosition().Z << ")  " // useless
    
       << setw(5)  << aTrack->GetParentID() << "  " // good
    //       << setw(5)  << nbSteps << "  "
    //       << setw(10) << setprecision(2) << Trleng << "  "
    //       << setw(5)  << particle->GetMass() << " "
    //       << setw(3)  << particle->GetCharge() << "  "
    //       << setw(5)  << particle->GetTotalEnergy() << ",  " // useless
       << setw(5)  << aTrack->GetVertexKineticEnergy() << "  " // good
       // << "( "
       // << setw(5)  << aTrack->GetVertexMomentumDirection().X << ", "  // useless
       // << setw(5)  << aTrack->GetVertexMomentumDirection().Y << ", "  // useless
       // << setw(5)  << aTrack->GetVertexMomentumDirection().Z << ") "  // useless
    //       << setw(5)  << particle->GetTotalMomentum() << "  " // useless
    //       << setw(5)  << particle->GetKineticEnergy() << "  |  " // useless
    //       << setw(5)  << step->GetTotalEnergyDeposit()
       << endl;

  /*
  //debug  
  const G4DynamicParticle* dynamic = aTrack->GetDynamicParticle();
  G4double dynamCharge = dynamic->GetCharge();
  G4int occup          = dynamic->GetTotalOccupancy();
  G4double   pdgMass   = dynamic->GetParticleDefinition()->GetPDGMass();    
  G4double invarMass   = dynamic->Get4Momentum().m();  
  G4double dynamMass   = dynamic->GetMass();
  G4double dif1 = invarMass - pdgMass;
  G4double dif2 = dynamMass - pdgMass;
  
  G4cout
    << "\n  End of track :"    
    << "\n    charge= " <<  dynamCharge << "  occupancy= " << occup
    << "\n   pdgMass= " << G4BestUnit (pdgMass  , "Energy")    
///    << "\n invarMass= " << G4BestUnit (invarMass, "Energy")
///    << "   invar-pdg= " << G4BestUnit (dif1     , "Energy")
    << "\n dynamMass= " << G4BestUnit (dynamMass, "Energy")
    << "   dynam-pdg= " << G4BestUnit (dif2     , "Energy")
    << G4endl;          
*/                       

}
