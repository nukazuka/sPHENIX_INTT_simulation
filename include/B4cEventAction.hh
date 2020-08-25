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
/// \file B4cEventAction.hh
/// \brief Definition of the B4cEventAction class

#ifndef B4cEventAction_h
#define B4cEventAction_h 1

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4UserEventAction.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

#include "B4cCalorHit.hh"
#include "B4RunAction.hh"
#include "B4cCalorimeterSD.hh"
#include "B4Analysis.hh"

#include "globals.hh"

#include <sstream>
#include <iomanip>
#include <vector>
#include <numeric> // for accumulate
#include <iomanip>

#include "Randomize.hh"

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy
/// deposit and track lengths of charged particles in Absober and Gap layers
/// stored in the hits collections.

class B4cEventAction : public G4UserEventAction
{
public:
  B4cEventAction();
  virtual ~B4cEventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);

  void AddWhatever( G4double de, G4double dl );
  
  static G4double beam_energy;
  
private:
  // methods
  B4cCalorHitsCollection* GetHitsCollection(G4int hcID, const G4Event* event) const;
  void PrintEventStatistics(G4double absoEdep, G4double absoTrackLength,
                            G4double gapEdep, G4double gapTrackLength) const;

  void GetHitsCollectionIDs();
  
  // data members
  // Ecal
  G4int  fEcalGapHCID;
  G4int  fEcalAbsHCID;
  std::vector < G4int > ID_abs_ecal;

  std::vector < G4int > ID_abs_tower_ecal;
  std::vector < G4int > ID_gap_tower_ecal;

  std::vector < G4int > ID_abs_tower_hcal;
  std::vector < G4int > ID_gap_tower_hcal;

  std::vector < G4int > ID_abs_tower_wall;
  std::vector < G4int > ID_gap_tower_wall;

  G4int ID_wall = -1;
  G4int ID_surroundings_ecal[4] = { -1 };
  G4int ID_surroundings_hcal[4] = { -1 };
  G4int ID_lid_upstream = -1;
  G4int ID_lid_downstream = -1;

  // Hcal
  G4int  fHcalGapHCID;
  G4int  fHcalAbsHCID;

  G4int  hc_id_si_0;
  G4int  hc_id_si_1;  
};

#endif
