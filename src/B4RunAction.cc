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
/// \file B4RunAction.cc
/// \brief Implementation of the B4RunAction class

#include "B4RunAction.hh"

G4double B4RunAction::beam_energy;
std::vector < G4double > B4RunAction::edep_absorber_tower_ecal( num_of_towers_ecal * num_of_towers_ecal );
std::vector < G4double > B4RunAction::edep_gap_tower_ecal( num_of_towers_ecal * num_of_towers_ecal );
std::vector < G4double > B4RunAction::edep_tower_ecal( num_of_towers_ecal * num_of_towers_ecal );
std::vector < G4double > B4RunAction::edep_absorber_tower_hcal( num_of_towers_hcal * num_of_towers_hcal );
std::vector < G4double > B4RunAction::edep_gap_tower_hcal( num_of_towers_hcal * num_of_towers_hcal );
std::vector < G4double > B4RunAction::edep_tower_hcal( num_of_towers_hcal * num_of_towers_hcal );

B4RunAction::B4RunAction( B4PrimaryGeneratorAction* pga_arg ) : G4UserRunAction()
//B4RunAction::B4RunAction() : G4UserRunAction()
{
  // set printing event number per each event
  //  G4RunManager::GetRunManager()->SetPrintProgress(1);
  DefineCommands();

  pga = pga_arg;
  // Create analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  //analysisManager->SetVerboseLevel( 9999 ); // usesless
  analysisManager->SetNtupleMerging(true);


  // Creating ntuple
  analysisManager->CreateNtuple("tr", "Edep");

  // sum of energy deposit to the detector

  // Beam energy in GeV
  analysisManager->CreateNtupleDColumn( "Beam" );

  // sum of energy deposit to all sensitive detector
  analysisManager->CreateNtupleDColumn( "Total" );

  // energy deposit to the Ecal
  analysisManager->CreateNtupleDColumn( "Ecal" );
  analysisManager->CreateNtupleDColumn( "EcalAbs" );
  analysisManager->CreateNtupleDColumn( "EcalGap" );

  // energy deposit to the Hcal
  analysisManager->CreateNtupleDColumn( "Hcal" );
  analysisManager->CreateNtupleDColumn( "HcalAbs" );
  analysisManager->CreateNtupleDColumn( "HcalGap" );

  // energy deposit to the wall
  analysisManager->CreateNtupleDColumn( "Wall" );

  // energy deposit to surrounding blocks around the Ecal
  analysisManager->CreateNtupleDColumn( "LeakEcal" );
  
  // energy deposit to surrounding blocks around the Hcal
  analysisManager->CreateNtupleDColumn( "LeakHcal" );

  // energy deposit to the lid on the upstream side
  analysisManager->CreateNtupleDColumn( "LeakUpstream" );

  // energy deposit to the lid on the downstream side
  analysisManager->CreateNtupleDColumn( "LeakDownstream" );

  // sum of energy deposit to absorbers, gaps, or absorbers+gaps of the tower of the Ecal
  analysisManager->CreateNtupleDColumn( "EcalTower"	, edep_tower_ecal );
  analysisManager->CreateNtupleDColumn( "EcalTowerAbs"	, edep_absorber_tower_ecal );
  analysisManager->CreateNtupleDColumn( "EcalTowerGap"	, edep_gap_tower_ecal );

  // those to the Hcal
  analysisManager->CreateNtupleDColumn( "HcalTower"	, edep_tower_hcal );
  analysisManager->CreateNtupleDColumn( "HcalTowerAbs"	, edep_absorber_tower_hcal );
  analysisManager->CreateNtupleDColumn( "HcalTowerGap"	, edep_gap_tower_hcal );

  analysisManager->FinishNtuple();

}


B4RunAction::~B4RunAction()
{
  delete G4AnalysisManager::Instance();
}

////////////////////////////////////////////////////////////////////////////
// Private function                                                       //
////////////////////////////////////////////////////////////////////////////
void B4RunAction::DefineCommands()
{
  fMessenger
    = new G4GenericMessenger(this, "/cal/run/", "Commands for run in this application");

  //////////////////////////////////////////////////////////////////////////////
  // Tag for the output file to indentify or something
  G4GenericMessenger::Command& setTag
    = fMessenger->DeclareProperty( "tag", tag );
  setTag.SetGuidance( "Tag for the output file." );
  setTag.SetParameterName( "tag", true ); // (name, is_omittable)
  setTag.SetDefaultValue( "" );

  //////////////////////////////////////////////////////////////////////////////
  // Path to a directory for the output file
  G4GenericMessenger::Command& setDir
    = fMessenger->DeclareProperty( "dir", output_dir );
  setDir.SetGuidance( "A path to a direcory to be output" );
  setDir.SetParameterName( "dir", true ); // (name, is_omittable)
  setDir.SetDefaultValue( "" );

  //////////////////////////////////////////////////////////////////////////////
  // Beam energy
  // G4GenericMessenger::Command& setBeamEnergy
  //   = fMessenger->DeclarePropertyWithUnit( "beamEnergy", "GeV", beam_energy );
  // setBeamEnergy.SetGuidance( "Beam energy in the units of GeV" );
  // setBeamEnergy.SetParameterName( "BeamEnergy", false ); // (name, is_omittable)
  // setBeamEnergy.SetDefaultValue( 1.0 );

}

////////////////////////////////////////////////////////////////////////////
// Public functions                                                       //
////////////////////////////////////////////////////////////////////////////
void B4RunAction::BeginOfRunAction(const G4Run* run/*run*/)
{
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();  
  auto UImanager = G4UImanager::GetUIpointer();

  if( isFirst ){
    beam_energy =  UImanager->GetCurrentDoubleValue( "/gun/energy" ) * GeV; // GeV is used whatever I give? changed to MeV

    // UImanager->ApplyCommand( "/gun/energy 7 MeV" );
    // std::cout << "UI manager, /gun/energy: "
    // 	      << UImanager->GetCurrentStringValue( "/gun/energy" , 1 ) << " "
    // 	      << UImanager->GetCurrentStringValue( "/gun/energy" , 2 ) << " "
    // 	      << std::endl;
    // std::cout << beam_energy << "\t" << beam_energy / GeV << std::endl;
    isFirst = false;
  }

  // define output name
  G4String particle = UImanager->GetCurrentStringValue( "/gun/particle" );
  std::stringstream ss;
  ss << beam_energy / GeV;
  G4String energy   = ss.str() + "GeV";
  //  G4String event_num = UImanager->GetCurrentStringValue( "/gun/beamOn") + "events";

  // if energy is less than 1, it should be in MeV
  if( beam_energy < 1 ){
    std::stringstream ss2;
    ss2 << beam_energy << "MeV";
    energy =  ss.str();
  }

  // Open an output file
  G4String fileName = particle + "_" + energy;// + "_" + event_num;
  if( pga->IsMomentumSpread() )
    fileName += "_EbeamFracON";
  else
    fileName += "_EbeamFracOFF";

  if( tag != "" )
    fileName += "_" + tag;
  
  // In the case of precision of the beam energy is better than 1 GeV, for example 1.1 GeV,
  // a suffix ".root" should be added
  if( fileName.contains( ".root" ) == false )
    fileName += ".root";

  if( output_dir != "" )
    fileName = output_dir + "/" + fileName;
  
  analysisManager->OpenFile(fileName);
}


void B4RunAction::EndOfRunAction(const G4Run* run)
{
  auto analysisManager = G4AnalysisManager::Instance();

  analysisManager->Write();
  analysisManager->CloseFile();

}
