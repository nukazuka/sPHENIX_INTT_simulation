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
/// \file B4cActionInitialization.cc
/// \brief Implementation of the B4cActionInitialization class

#include "B4cActionInitialization.hh"

B4cActionInitialization::B4cActionInitialization(B4cDetectorConstruction* det_const_arg)
  : G4VUserActionInitialization(),
    detector_constraction_( det_const_arg )
{}

B4cActionInitialization::~B4cActionInitialization()
{}

void B4cActionInitialization::BuildForMaster() const
{
  SetUserAction( new B4RunAction( new B4PrimaryGeneratorAction() ) );
}

void B4cActionInitialization::Build() const
{
  auto pga = new B4PrimaryGeneratorAction();
  SetUserAction( pga );
  //  SetUserAction(new B4PrimaryGeneratorAction);
  
  SetUserAction( new B4RunAction( pga ) );
  B4cEventAction* event_action = new B4cEventAction;
  SetUserAction( event_action );
  SetUserAction( new TrackingAction( pga ) );
  SetUserAction( new SteppingAction(detector_constraction_, event_action)  );
  //SetUserAction(new SteppingAction(fDetConstruction,eventAction));
}  
