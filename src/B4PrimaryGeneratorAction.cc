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
/// \file B4PrimaryGeneratorAction.cc
/// \brief Implementation of the B4PrimaryGeneratorAction class

#include "B4PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

B4PrimaryGeneratorAction::B4PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(nullptr),
   fMessenger()
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  DefineCommands();

  //  auto UImanager = G4UImanager::GetUIpointer();
  
  // default particle kinematic
  auto particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle( "e-" );
  fParticleGun->SetParticleDefinition( particleDefinition );
  fParticleGun->SetParticleMomentumDirection( G4ThreeVector( 0., 0., 1. ) );
  //  fParticleGun->SetParticleEnergy(50.*MeV);
  fParticleGun->SetParticleEnergy( 1 * GeV );
}

B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void B4PrimaryGeneratorAction::DefineCommands()
{
  fMessenger
    = new G4GenericMessenger(this, "/cal/beam/", "Commands for beam in this application");

  //////////////////////////////////////////////////////////////////////////////
  // Switch for the momentum spread
  G4GenericMessenger::Command& setMomentumSpread
    = fMessenger->DeclareProperty( "isMomentumSpread", isMomentumSpread );
  setMomentumSpread.SetGuidance( "Switch the surroundings around the calorimeters ON(true) or OFF(false)." );
  setMomentumSpread.SetParameterName( "IsMomentumSpread", false ); // (name, is_omittable)
  setMomentumSpread.SetDefaultValue( "false" );

  //////////////////////////////////////////////////////////////////////////////
  // Beam energy
  // G4GenericMessenger::Command& setBeamEnergy
  //   = fMessenger->DeclarePropertyWithUnit( "beamEnergy", "GeV", beam_energy );
  // setBeamEnergy.SetGuidance( "Beam energy in the units of GeV" );
  // setBeamEnergy.SetParameterName( "BeamEnergy", false ); // (name, is_omittable)
  // setBeamEnergy.SetDefaultValue( 1.0 );

}

void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume 
  // from G4LogicalVolumeStore
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast < G4Box* > ( worldLV->GetSolid() );
  }

  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();  
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  } 
  
  // Set gun position
  fParticleGun
    ->SetParticlePosition(G4ThreeVector( 1 * mm, 1 * mm , -worldZHalfLength));
    //    ->SetParticlePosition(G4ThreeVector( pos, pos , -worldZHalfLength));
    //    ->SetParticlePosition(G4ThreeVector( 0, 0 , -worldZHalfLength));

  //std::cout << CLHEP::RandGauss::shoot(1, 0.01 ) << "\t" << CLHEP::RandGauss::shoot(1, 0.01 ) << std::endl;
  //auto UImanager = G4UImanager::GetUIpointer();

  if( isMomentumSpread )
    {
      //  G4double beam_energy = UImanager->GetCurrentDoubleValue( "/gun/energy" );
      G4double factor = CLHEP::RandGauss::shoot( 1.0, 0.02 );
      //  std::cout << B4RunAction::beam_energy << "\t"
      //<< factor << "\t"
      //<< B4RunAction::beam_energy * factor << std::endl;
      //  std::cout << factor << std::endl;
      
      //    	    << UImanager->GetCurrentDoubleValue( "/gun/energy" ) << "\t";
      
      // std::cout << B4RunAction::beam_energy << " | "
      // 	    << B4RunAction::beam_energy * factor * MeV << " (*MeV)\t"
      // 	    << B4RunAction::beam_energy * factor * GeV << " (*GeV)\t"
      // 	    << B4RunAction::beam_energy * factor / GeV << " (/GeV)"
      // 	    << std::endl;
      
      //  fParticleGun->SetParticleEnergy( B4RunAction::beam_energy * factor * GeV);
      //std::cout << "Momentum spread ON, " <<  B4RunAction::beam_energy * factor << std::endl;;
      fParticleGun->SetParticleEnergy( B4RunAction::beam_energy * factor);
      
    }

  // not used
  //  else 
  //    fParticleGun->SetParticleEnergy( beam_energy );

  fParticleGun->GeneratePrimaryVertex( anEvent );
}

