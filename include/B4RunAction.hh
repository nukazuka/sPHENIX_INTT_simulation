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
/// \file B4RunAction.hh
/// \brief Definition of the B4RunAction class

#ifndef B4RunAction_h
#define B4RunAction_h

#include <sstream>

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4UserRunAction.hh"
#include "globals.hh"
#include "Constants.hh" // should be included after headers starting from "G4" otherwise many difinitions are missed

#include "B4Analysis.hh"
#include "B4PrimaryGeneratorAction.hh"

class G4Run;
class B4PrimaryGeneratorAction;
/// Run action class
///
/// It accumulates statistic and computes dispersion of the energy deposit 
/// and track lengths of charged particles with use of analysis tools:
/// H1D histograms are created in BeginOfRunAction() for the following 
/// physics quantities:
/// - Edep in absorber
/// - Edep in gap
/// - Track length in absorber
/// - Track length in gap
/// The same values are also saved in the ntuple.
/// The histograms and ntuple are saved in the output file in a format
/// accoring to a selected technology in B4Analysis.hh.
///
/// In EndOfRunAction(), the accumulated statistic and computed 
/// dispersion is printed.
///

class B4RunAction : public G4UserRunAction
{
public:               
  B4RunAction( B4PrimaryGeneratorAction* pga_arg );
  //B4RunAction();
  virtual ~B4RunAction();
  
  virtual void BeginOfRunAction(const G4Run*);
  virtual void   EndOfRunAction(const G4Run*);
  
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //  
  // column 1  2   3   4   5   6   layer num    1     2     3   ...   N       //
  // row +-----------------------+           +-----+-----+-----+---+-----+    //
  //  0  | 1 | 2 | 3 | 4 | 5 | 6 |           |A1 G1 A2 G2 A3 G4 ... AX GX|    //
  //  1  | 7 | 8 | 9 |10 |11 |12 |           |A1 G1 A2 G2 A3 G4 ... AX GX|    //
  //  2  |13 |14 |15 |16 |17 |18 |           |A1 G1 A2 G2 A3 G4 ... AX GX|    //
  //  3  |19 |20 |21 |22 |23 |24 |           |A1 G1 A2 G2 A3 G4 ... AX GX|    //
  //  4  |25 |26 |27 |28 |29 |30 |           |A1 G1 A2 G2 A3 G4 ... AX GX|    //
  //  5  |31 |32 |33 |34 |35 |36 |           |A1 G1 A2 G2 A3 G4 ... AX GX|    //
  //     +-----------------------+           +-----+-----+-----+---+-----+    // 
  //               ----> x-axis                     ----> z-axis              //
  //    View from upstream                View from side                      //
  //                                                                          //
  //  where A1, ... AX mean n-th layer of the absorber                        //
  //  and G1, ... AGX mean n-th layer of the gap.                             //
  //  Numbers in the view from upstream are ID of towers.                     //
  //  Energy deposit is recoreded for each absorber and gap layers.           //
  //  Variables of vector of vector contain the values.                       //
  //  For example, val_absorber[layer + tower* tower row]                     //
  //  Due to technical reasons, 1 dimentional vector is used.                 //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

  //  std::vector < G4double >& GetVector( std::string type );

  //! a flag to 
  G4bool isFirst = true;
  //! beam energy at the beginning of the run, 
  static G4double beam_energy;
  //////////////////////////////////////////////////////////////////////////////
  // variables for Ecal
  // each tower
  static std::vector < G4double > edep_absorber_tower_ecal;
  static std::vector < G4double > edep_gap_tower_ecal;
  // each tower for both absorebers and gaps
  static std::vector < G4double > edep_tower_ecal;

  //////////////////////////////////////////////////////////////////////////////
  // variables for Hcal
  // each tower
  static std::vector < G4double > edep_absorber_tower_hcal;
  static std::vector < G4double > edep_gap_tower_hcal;
  // each tower for both abosorbers and gaps
  static std::vector < G4double > edep_tower_hcal;

private:
  void DefineCommands();

  G4String tag;
  G4String output_dir;
  B4PrimaryGeneratorAction* pga;
  G4GenericMessenger* fMessenger;
};


#endif

