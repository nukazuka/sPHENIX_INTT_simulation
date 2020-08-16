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
/// \file B4cEventAction.cc
/// \brief Implementation of the B4cEventAction class

#include "B4cEventAction.hh"
G4double B4cEventAction::beam_energy;

/*B4cEventAction::B4cEventAction()
 : G4UserEventAction(),
 fEcalAbsHCID(-1),
 fEcalGapHCID(-1),
 fHcalAbsHCID(-1),
 fHcalGapHCID(-1)
 {}
*/
B4cEventAction::B4cEventAction()
 : G4UserEventAction(),
   hc_id_si_0(-1),
   hc_id_si_1(-1)
{}

B4cEventAction::~B4cEventAction()
{}

// definition: using B4cCalorHitsCollection = G4THitsCollection<B4cCalorHit>;
B4cCalorHitsCollection*
B4cEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection
    = static_cast<B4cCalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));

  if ( ! hitsCollection )
    {
      G4ExceptionDescription msg;
      msg << "Cannot access hitsCollection ID " << hcID;
      G4Exception("B4cEventAction::GetHitsCollection()",
		  "MyCode0003", FatalException, msg);
    }
  
  return hitsCollection;
}

void B4cEventAction::PrintEventStatistics(
                              G4double absoEdep, G4double absoTrackLength,
                              G4double gapEdep, G4double gapTrackLength) const
{
  /*
  // print event statistics
  G4cout
     << "   Absorber: total energy: "
     << std::setw(7) << G4BestUnit(absoEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
     << G4endl
     << "        Gap: total energy: "
     << std::setw(7) << G4BestUnit(gapEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
     << G4endl;
     */
}

void B4cEventAction::BeginOfEventAction(const G4Event* /*event*/)
{
  std::cout << std::string( 100 , '-' ) << std::endl;
  std::cout << " Start event"  << std::endl;
  
}

void B4cEventAction::GetHitsCollectionIDs()
{

  /*
  // skip everything if IDs were already got
  if ( fEcalAbsHCID != -1 )
    return;

  //////////////////////////////////////////////////////////////////////////////
  // Ecal                                                                     //
  //////////////////////////////////////////////////////////////////////////////
  for( int i=0; i<num_of_towers_ecal; i++ )
    {      
      for( int j=0; j<num_of_towers_ecal; j++ )
	{
	      
	  std::stringstream num;
	  num  << std::setw(2) <<  std::setfill( '0' ) << i << std::setw(2) << std::setfill('0') << j;

	  // for absorber
	  std::stringstream absorber_EM_HC;
	  absorber_EM_HC << "AbsorberEMHitsCollection" << num.str();
	  fEcalAbsHCID = G4SDManager::GetSDMpointer()->GetCollectionID( absorber_EM_HC.str().c_str() );
	  ID_abs_tower_ecal.push_back( fEcalAbsHCID );

	  // for gap
	  std::stringstream gap_EM_HC;
	  gap_EM_HC << "GapEMHitsCollection" << num.str();
	  fEcalAbsHCID = G4SDManager::GetSDMpointer()->GetCollectionID( gap_EM_HC.str().c_str() );
	  ID_gap_tower_ecal.push_back( fEcalAbsHCID );

	}
    }

  //////////////////////////////////////////////////////////////////////////////
  // Hcal                                                                     //
  //////////////////////////////////////////////////////////////////////////////
  for( int i=0; i<num_of_towers_hcal; i++ )
    {      
      for( int j=0; j<num_of_towers_hcal; j++ )
	{
	      
	  std::stringstream num;
	  num  << std::setw(2) <<  std::setfill( '0' ) << i << std::setw(2) << std::setfill('0') << j;

	  // for absorber
	  std::stringstream absorber_EM_HC;
	  absorber_EM_HC << "AbsorberHitsCollection" << num.str();
	  fHcalAbsHCID = G4SDManager::GetSDMpointer()->GetCollectionID( absorber_EM_HC.str().c_str() );
	  ID_abs_tower_hcal.push_back( fHcalAbsHCID );

	  // for gap
	  std::stringstream gap_EM_HC;
	  gap_EM_HC << "GapHitsCollection" << num.str();
	  fHcalAbsHCID = G4SDManager::GetSDMpointer()->GetCollectionID( gap_EM_HC.str().c_str() );
	  ID_gap_tower_hcal.push_back( fHcalAbsHCID );

	}
    }

  //////////////////////////////////////////////////////////////////////////////
  // Tha wall                                                                 //
  //////////////////////////////////////////////////////////////////////////////
  // Give the collect name of the hit collection of the wall
  // when you implement it
  ID_wall = G4SDManager::GetSDMpointer()->GetCollectionID( "Wall_HitsCollection" );
  
  //////////////////////////////////////////////////////////////////////////////
  // Surroundings                                                             //
  //////////////////////////////////////////////////////////////////////////////

  G4String directions[4] = { "top", "bottom", "left", "right" };

  // loop over directions
  for( int i=0; i<4; i++ )
    {

      // loop over calorimeters
      G4String calos[2] = { "ecal" , "hcal" };
      for( int j=0; j<2; j++ )
	{

	  // make name of the hits collection
	  std::stringstream cal_HC;
	  cal_HC << "Surrounding_" << calos[j] << "_" << directions[i] << "_HitsCollection";

	  // set the pointer of the id of the block
	  G4int* id;
	  if( j == 0 )
	    id = ID_surroundings_ecal;
	  else 
	    id = ID_surroundings_hcal;

	  id[i] = G4SDManager::GetSDMpointer()->GetCollectionID( cal_HC.str().c_str() );

	  //std::cout << id[i] << std::endl;
	}
    }

  ID_lid_upstream = G4SDManager::GetSDMpointer()->GetCollectionID( "Lid_upstream_HitsCollection" );
  ID_lid_downstream = G4SDManager::GetSDMpointer()->GetCollectionID( "Lid_downstream_HitsCollection" );
  */
}

void B4cEventAction::EndOfEventAction(const G4Event* event)
{
  std::cout << " Processing for the end of this event" << std::endl;

  // Get hits collections IDs (only once)
  if ( hc_id_si_0 == -1 ) {
    hc_id_si_0
      = G4SDManager::GetSDMpointer()->GetCollectionID( "SiStrip0HitsCollection" );
    hc_id_si_1
      = G4SDManager::GetSDMpointer()->GetCollectionID( "SiStrip1HitsCollection" );
  }

  //  GetHitsCollectionIDs();
  // Get hits collections
  B4cCalorHitsCollection* hc_si_0 = GetHitsCollection( hc_id_si_0, event );
  B4cCalorHitsCollection* hc_si_1 = GetHitsCollection( hc_id_si_0, event );

  for( int i=0; i<hc_si_0->entries()-1; i++ )
    {
      auto hc = (*hc_si_0)[i];
      std::cout << i << "\t" << hc->GetEdep() << std::endl;
    }
  
  // Print per event (modulo n)
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  /*
  //////////////////////////////////////////////////////////////////////////////
  // Ecal                                                                     //
  //////////////////////////////////////////////////////////////////////////////
  // loop over all Ecal towers to get energy deposit
  G4double sum = 0.0;
  for( int i=0; i<ID_abs_tower_ecal.size(); i++ )
    {

      // for absorbers
      auto test = GetHitsCollection( ID_abs_tower_ecal[i], event );
      auto testHit = (*test)[test->entries()-1];
      B4RunAction::edep_absorber_tower_ecal[i] = testHit->GetEdep() / GeV;
      //      B4RunAction::GetVector( "absorber_tower_ecal" ).push_back( testHit->GetEdep() / GeV );
      //run_action->GetVector( "absorber_tower_ecal" ).push_back( testHit->GetEdep() / GeV );

      // for gaps
      test = GetHitsCollection( ID_gap_tower_ecal[i], event );
      testHit = (*test)[test->entries()-1];
      B4RunAction::edep_gap_tower_ecal[i] = testHit->GetEdep() / GeV;

      B4RunAction::edep_tower_ecal[i] = B4RunAction::edep_absorber_tower_ecal[i] + B4RunAction::edep_gap_tower_ecal[i];
      sum += B4RunAction::edep_absorber_tower_ecal[i];
      sum += B4RunAction::edep_gap_tower_ecal[i];
    }

  // get sum of energy deposit to absorbers of Ecal
  // std::accumulate( iterator of the first element, iterator of the last element, inital value )
  // This function returns: (inital value + 1st + 2nd + ... + last)
  G4double edep_ecal_abs = std::accumulate( B4RunAction::edep_absorber_tower_ecal.begin(),
					    B4RunAction::edep_absorber_tower_ecal.end(), 0.0 );
  
  G4double edep_ecal_gap = std::accumulate( B4RunAction::edep_gap_tower_ecal.begin(),
					    B4RunAction::edep_gap_tower_ecal.end(), 0.0 );
  
  G4double edep_ecal = edep_ecal_abs + edep_ecal_gap;

  //////////////////////////////////////////////////////////////////////////////
  // Hcal                                                                     //
  //////////////////////////////////////////////////////////////////////////////
  // loop over all Hcal towers to get energy deposit
  for( int i=0; i<ID_abs_tower_hcal.size(); i++ )
    {

      // for absorbers
      auto test = GetHitsCollection( ID_abs_tower_hcal[i], event );

      auto testHit = (*test)[test->entries()-1];
      B4RunAction::edep_absorber_tower_hcal[i] = testHit->GetEdep() / GeV;

      // for gaps
      test = GetHitsCollection( ID_gap_tower_hcal[i], event );
      testHit = (*test)[test->entries()-1];
      B4RunAction::edep_gap_tower_hcal[i] = testHit->GetEdep() / GeV;
      B4RunAction::edep_tower_hcal[i] = B4RunAction::edep_absorber_tower_hcal[i] + B4RunAction::edep_gap_tower_hcal[i];      
    }

  // get sum of energy deposit to absorbers of Hcal
  G4double edep_hcal_abs = std::accumulate( B4RunAction::edep_absorber_tower_hcal.begin(),
					    B4RunAction::edep_absorber_tower_hcal.end(), 0.0 );
  G4double edep_hcal_gap = std::accumulate( B4RunAction::edep_gap_tower_hcal.begin(),
					    B4RunAction::edep_gap_tower_hcal.end(), 0.0 );
  G4double edep_hcal = edep_hcal_abs + edep_hcal_gap;


  //////////////////////////////////////////////////////////////////////////////
  // Wall                                                                     //
  //////////////////////////////////////////////////////////////////////////////
  // loop over all Wall towers (should be one-tower) to get energy deposit
  G4double edep_wall = 0.0;

  // get energy deposit to the wall if the wall is defined
  if( ID_wall != -1 )
    {
      auto test_wall = GetHitsCollection( ID_wall, event );  
      auto testHit_wall = (*test_wall)[test_wall->entries()-1];
      edep_wall = testHit_wall->GetEdep() / GeV;

    }


  //////////////////////////////////////////////////////////////////////////////
  // Surroundings                                                             //
  //////////////////////////////////////////////////////////////////////////////
  G4double edep_surroundings_ecal = 0.0;
  G4double edep_surroundings_hcal = 0.0;
  G4double edep_lid_upstream = 0.0;
  G4double edep_lid_downstream = 0.0;

  // get infomation only if the blocks are on
  if( ID_surroundings_ecal[0] != -1 )
    {
      // loop over blocks around ecal
      for( int i=0; i<4; i++ )
	{
	  auto test = GetHitsCollection( ID_surroundings_ecal[i], event );
	  auto testHit = (*test)[test->entries()-1];
	  edep_surroundings_ecal += testHit->GetEdep() / GeV;
	}
      
      // loop over blocks around hcal
      for( int i=0; i<4; i++ )
	{
	  auto test = GetHitsCollection( ID_surroundings_hcal[i], event );
	  auto testHit = (*test)[test->entries()-1];
	  edep_surroundings_hcal += testHit->GetEdep() / GeV;
	}
      
      auto test = GetHitsCollection( ID_lid_upstream, event );
      auto testHit = (*test)[test->entries()-1];
      edep_lid_upstream = testHit->GetEdep() / GeV;

      test = GetHitsCollection( ID_lid_downstream, event );
      testHit = (*test)[test->entries()-1];
      edep_lid_downstream = testHit->GetEdep() / GeV;
     
    }
  */
  
  //  auto analysisManager = G4AnalysisManager::Instance();
  /*
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;

    // PrintEventStatistics(
    //   absoHit->GetEdep(), absoHit->GetTrackLength(),
    //   gapHit->GetEdep(), gapHit->GetTrackLength());
  }
  */

  /*
  // fill ntuple, if the wall and/or the surrounding blocks are off, it's OK to use the variables
  // because their default values are 0.0
  G4double edep_total = edep_ecal + edep_hcal
    + edep_wall
    + edep_surroundings_ecal
    + edep_surroundings_hcal
    + edep_lid_upstream
    + edep_lid_downstream;

  //  std::cout << "EventAction, beam energy = " << beam_energy << " is filled" << std::endl;
  analysisManager->FillNtupleDColumn( 0, beam_energy / GeV );      // 0: Beam
  analysisManager->FillNtupleDColumn( 1, edep_total);              // 1 : Total
  analysisManager->FillNtupleDColumn( 2, edep_ecal );              // 2 : Ecal
  analysisManager->FillNtupleDColumn( 3, edep_ecal_abs );          // 3 : EcalAbs
  analysisManager->FillNtupleDColumn( 4, edep_ecal_gap );           // 4 : EcalGap
  analysisManager->FillNtupleDColumn( 5, edep_hcal );              // 5 : Hcal
  analysisManager->FillNtupleDColumn( 6, edep_hcal_abs );          // 6 : HcalAbs
  analysisManager->FillNtupleDColumn( 7, edep_hcal_gap );          // 7 : HcalGap
  analysisManager->FillNtupleDColumn( 8, edep_wall );              // 8 : Wall
  analysisManager->FillNtupleDColumn( 9, edep_surroundings_ecal ); // 9 : LeakEcal
  analysisManager->FillNtupleDColumn(10, edep_surroundings_hcal ); // 10: LeakHcal
  analysisManager->FillNtupleDColumn(11, edep_lid_upstream);       // 11: LeakUpstream
  analysisManager->FillNtupleDColumn(12, edep_lid_downstream);     // 12: LeakDownsstream
  */
  
  // Write vector values
  //analysisManager->AddNtupleRow();
}
