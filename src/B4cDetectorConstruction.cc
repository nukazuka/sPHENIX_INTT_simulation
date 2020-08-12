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
// $Id$
//
/// \file B4cDetectorConstruction.cc
/// \brief Implementation of the B4cDetectorConstruction class

#include "B4cDetectorConstruction.hh"

//G4GenericMessenger* B4cDetectorConstruction::fMessenger;

B4cDetectorConstruction::B4cDetectorConstruction()
  : G4VUserDetectorConstruction(),
    fMessenger(0),
    //   fMagField(0),
    fCheckOverlaps(true),
    isSurroundings(false),
    wallThickness(0.0)
{

  DefineCommands();
  
  layerThickness_ecal	= absoThickness_ecal + gapThickness_ecal;
  layerThickness_hcal	= absoThickness_hcal + gapThickness_hcal;
  calorThickness_ecal	= nofLayers * layerThickness_ecal;
  calorThickness_hcal	= nofLayers2 * layerThickness_hcal;

  worldSizeXY		= 10.0 * ( 5.0 * calorSizeXY );
  worldSizeZ		= 5.0 * ( 4 * calorThickness_hcal + calorThickness_ecal );

  ////////////////////////////////////////////////////
  // I want to set following 2 parameters in macro. //
  // But initialization of geometry starts before   //
  // executing macros. How can I do this?           //
  //  isSurroundings = true;                        //
  wallThickness = 0.0 * cm;                      //
  wallThickness = 1.5 * cm;                       //
  //wallThickness = 10.0 * cm;                      //
  ////////////////////////////////////////////////////  

  isSurroundings = true;
  //isSurroundings = false;
  surroundingsThickness = 50 * cm;
  //surroundingsThickness = 2 * m;

  // as you want
  // lidThickness = surroundingsThickness;
  lidThickness = 30 * cm;


  DefineVisAttributes();
}

B4cDetectorConstruction::~B4cDetectorConstruction()
{
  //  delete fMagField;
  delete fMessenger;
}

G4VPhysicalVolume* B4cDetectorConstruction::Construct()
{
  std::cout << "\n\nConstruction start!\n\n" << std::endl;
  std::cout << wallThickness << std::endl;
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();

}

void B4cDetectorConstruction::DefineCommands()
{
  // Define /B4/det commands using generic messenger class
  fMessenger
    = new G4GenericMessenger(this, "/cal/detector/", "Commands for detectors in this application");

  //////////////////////////////////////////////////////////////////////////////
  // Switch for the surroundings
  G4GenericMessenger::Command& setSurroundings
    = fMessenger->DeclareProperty( "isSurroundings", isSurroundings );
  setSurroundings.SetGuidance( "Switch the surroundings around the calorimeters ON(true) or OFF(false)." );
  setSurroundings.SetParameterName( "IsSurroundings", false ); // (name, is_omittable)
  //  setSurroundings.SetDefaultValue( "false" );
  setSurroundings.SetDefaultValue( "true" );

  //////////////////////////////////////////////////////////////////////////////
  // Thickness of the wall
  G4GenericMessenger::Command& setWallThickness
    //    = fMessenger->DeclarePropertyWithUnit( "wallThickness", "cm", wallThickness );
    = fMessenger->DeclareMethodWithUnit( "wallThickness", "cm",
					 &B4cDetectorConstruction::SetWallThickness,
					 "Set thickness of the wall in cm" );
  setWallThickness.SetGuidance( "Give the thickness of the wall. Negative value and 0 remove the wall" );
  setWallThickness.SetParameterName( "WallThickness", false ); // (name, is_omittable)
  setWallThickness.SetDefaultValue( 0.0 );
  
  // Define /B4/det/setMagField command
  //  G4GenericMessenger::Command& setMagFieldCmd
  //    = fMessenger->DeclareMethod("setMagField",
  //                                &B4cDetectorConstruction::SetMagField,
  //                                "Define magnetic field value (in X direction");
  //  setMagFieldCmd.SetUnitCategory("Magnetic flux density");


  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//////////////////////////////////////////////////////////////////////////////
//   Geometry in y-z plane (surrounding blocks not drawn)                   //
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//      calorThickness_ecal                                                 //
//      |<---------->|                                                      //
//      |            | wallThickness                                        //
//      |            |<-->|                                                 //
//      |            |    |  calorThickness_hcal                            //
//      |            |    |<--------------------->|                         //
//      |            |    |                       |                         //
//      |             ---  -----------------------                          //
//      |            |   ||                       |                         //
//       ----------- |   ||                       |                         //
// Beam |           || W ||                       |                         //
// ---> |    Ecal   || a ||         Hcal          |                         //
//      |           || l ||                       |                         //
//       -----+----- | l ||                       |                         //
//            |      |   ||                       |                         //
//            |       -+-  -----------+-----------+                         //
//            |        |              |           |                         //
//            |        |              |           |                         //
// -----------|--------|--------------|-----------|-----> z                 //
// center     |        |              |           0                         //
// positions  |        |              |                                     //
//            |        |              -calorThickness_hcal/2                //
//            |        |                                                    //
//            |        - calorThickness_hcal - wallThickness/2              //
//            |                                                             //
// 	      -calorThickness_ecal/2 - calorThickness_hcal - wallThickness  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

void B4cDetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  G4bool fromIsotopes = false;
  nistManager->FindOrBuildMaterial( "G4_Pb", fromIsotopes );
  nistManager->FindOrBuildMaterial( "G4_Fe", fromIsotopes );
  nistManager->FindOrBuildMaterial( "G4_POLYSTYRENE", fromIsotopes );
  nistManager->FindOrBuildMaterial( "G4_He", fromIsotopes );
  nistManager->FindOrBuildMaterial( "G4_AIR", fromIsotopes );
  nistManager->FindOrBuildMaterial( "G4_Ag", fromIsotopes );
  nistManager->FindOrBuildMaterial( "G4_WATER", fromIsotopes );
  nistManager->FindOrBuildMaterial( "G4_H", fromIsotopes );
  nistManager->FindOrBuildMaterial( "G4_lH2", fromIsotopes );
  nistManager->FindOrBuildMaterial( "G4_O", fromIsotopes );
  nistManager->FindOrBuildMaterial( "G4_lO2", fromIsotopes );

  // Vacuum
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;

  new G4Material("Galactic", z = 1., a = 1.01 * g / mole, density = universe_mean_density,
		 kStateGas, 2.73 * kelvin, 3.e-18 * pascal);

  new G4Material("Stopper", z = 1., a = 1.01 * g / mole, density = 1.0e10 * g / m3,
		 kStateSolid );


  /////////////////////////////////////////////////////////////////
  // some material to know their radiation length

  G4int ncomponets, natoms;
  G4Material* epoxy = new G4Material( "epoxy", 1.2 * g / cm3, ncomponets = 2 );
  epoxy->AddElement( G4Element::GetElement( "H" ), natoms = 2 );
  epoxy->AddElement( G4Element::GetElement( "C" ), natoms = 2 );

  G4double fractionmass;
  G4Material* ag_epoxy_glue = new G4Material( "ag_epoxy_glue", density = 3.2 * g / cm3, ncomponets = 2 );
  ag_epoxy_glue->AddMaterial( epoxy, fractionmass = 0.79 );
  ag_epoxy_glue->AddMaterial( G4Material::GetMaterial( "G4_Ag"), fractionmass = 0.21 );

  std::cout << std::string(100, '=' ) << std::endl;
  std::cout << "Water: "    << G4Material::GetMaterial( "G4_WATER" )->GetRadlen() << " mm" << std::endl;
  std::cout << "H: "        << G4Material::GetMaterial( "G4_H" )->GetRadlen()     << " mm" << std::endl;
  std::cout << "lH2: "      << G4Material::GetMaterial( "G4_lH2" )->GetRadlen()   << " mm" << std::endl;
  std::cout << "O: "        << G4Material::GetMaterial( "G4_O" )->GetRadlen()     << " mm" << std::endl;
  std::cout << "lO2: "      << G4Material::GetMaterial( "G4_lO2" )->GetRadlen()   << " mm" << std::endl;
  std::cout << "Ag: "       << G4Material::GetMaterial( "G4_Ag" )->GetRadlen()    << " mm" << std::endl;
  std::cout << "Epoxy: "    << epoxy->GetRadlen()                                 << " mm" << std::endl;
  std::cout << "The glue: " << ag_epoxy_glue->GetRadlen()                         << " mm" << std::endl;
  std::cout << std::string(100, '=' ) << std::endl;
  
  // Print materials
  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  defaultMaterial	= G4Material::GetMaterial( "Galactic" );
  absorberMaterial	= G4Material::GetMaterial( "G4_Pb" );
  absorberMaterial2	= G4Material::GetMaterial( "G4_Fe" );
  gapMaterial   	= G4Material::GetMaterial( "G4_POLYSTYRENE" );
  airMaterial   	= G4Material::GetMaterial( "G4_AIR" );
  stopperMaterial	= G4Material::GetMaterial( "Stopper" );

  // if definition of materials failed...
  if ( ! defaultMaterial
       || ! absorberMaterial
       || ! absorberMaterial2
       || ! gapMaterial
       || ! airMaterial
       || ! stopperMaterial )
    {
      G4cerr << "Cannot retrieve materials already defined. " << G4endl;
      G4cerr << "Exiting application " << G4endl;
      exit(1);
    }

}

G4VPhysicalVolume* B4cDetectorConstruction::DefineVolumes()
{

  // World
  G4VSolid* worldS = new G4Box("World", worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // name and size
  G4LogicalVolume* worldLV = new G4LogicalVolume( worldS, defaultMaterial, "World"); // solid, material, and name

  // if you want to fill the world with air,
  // solid, material, and name
  //G4LogicalVolume* worldLV = new G4LogicalVolume( worldS, airMaterial, "World");

  // Constructor of G4VPhysicalVolume
  // 1: no rotation, 2: position at (0,0,0), 3: logical volume,
  // 4: name,        5: mother  volume,      6: no boolean operation,
  // 7: copy number, 8: checking overlaps
  G4VPhysicalVolume* worldPV = new G4PVPlacement( 0, G4ThreeVector(), worldLV,
						  "World", 0, false,
						  0, fCheckOverlaps);

  // Visualization attributes
  worldLV->SetVisAttributes( G4VisAttributes::Invisible );

  // for Ecal
  DefineVolumes_Ecal( worldLV );

  // for Hcal
  DefineVolumes_Hcal( worldLV );

  // for the wall if needed
  if( wallThickness > 0.1 )
    DefineVolumes_Wall( worldLV );

  // for surrounding blocks to measure energy leakage
  if( isSurroundings == true )
    DefineVolumes_Surroundings( worldLV );

  // not used, to be deleted
  //DefineVolumes_NotUsed( worldLV );
  
  // Always return the physical World
  return worldPV;
}

void B4cDetectorConstruction::DefineVolumes_Ecal( G4LogicalVolume* worldLV )
{
  
  G4LogicalVolume* calorEMLV[num_of_towers_ecal][num_of_towers_ecal];
  G4LogicalVolume* layerEMLV[num_of_towers_ecal][num_of_towers_ecal];
  G4LogicalVolume* absorberEMLV[num_of_towers_ecal][num_of_towers_ecal];
  G4LogicalVolume* gapEMLV[num_of_towers_ecal][num_of_towers_ecal];

  for( int i=0; i<num_of_towers; i++ )
    {
      for( int j=0; j<num_of_towers; j++ )
	{
	  std::stringstream num;
	  num  << std::setw(2) <<  std::setfill( '0' ) << i
	       << std::setw(2) << std::setfill('0') << j;

	  ////////////////////////////////////////////
	  // Calorimeter                            //
	  ////////////////////////////////////////////
	  std::stringstream calorEMLV_name;
	  calorEMLV_name << "Ecal" << num.str();

	  // name and size
	  G4VSolid* calorimeterEM
	    = new G4Box( "Ecal", emcell / 2, emcell / 2, calorThickness_ecal / 2);
	  calorEMLV[i][j] = new G4LogicalVolume( calorimeterEM, defaultMaterial, calorEMLV_name.str().c_str() );
	  calorEMLV[i][j]->SetVisAttributes( simpleBoxVisAtt );

	  int copy_num = i * num_of_towers + j;
	  G4ThreeVector pos_cal
	    = G4ThreeVector( num_of_towers * emcell / 2 - emcell / 2 - j * emcell,
			     num_of_towers * emcell / 2 - emcell / 2 - i * emcell,
			     -calorThickness_ecal/2 - calorThickness_hcal - wallThickness  );
	  
	  new G4PVPlacement( 0, pos_cal,
			     calorEMLV[i][j], "Ecal", worldLV,
			     false, copy_num, fCheckOverlaps );

	  ////////////////////////////////////////////
	  // Layer                                  //
	  ////////////////////////////////////////////
	  std::stringstream layerEMLV_name;
	  layerEMLV_name << "LayerEcal" << num.str();

	  G4VSolid* layerEM
	    = new G4Box("Layer", emcell / 2, emcell / 2, layerThickness_ecal / 2); // name and size
	  layerEMLV[i][j] = new G4LogicalVolume( layerEM, defaultMaterial, layerEMLV_name.str().c_str() );
	  layerEMLV[i][j]->SetVisAttributes( simpleBoxVisAtt );

	  new G4PVReplica( "Layer", layerEMLV[i][j], calorEMLV[i][j],
			   kZAxis, nofLayers, layerThickness_ecal );

	  ////////////////////////////////////////////
	  // Absorber                               //
	  ////////////////////////////////////////////
	  std::stringstream absorberEMLV_name;
	  absorberEMLV_name << "AbsorberEcal" << num.str();

	  G4VSolid* absorberEM = new G4Box("Absorber", emcell / 2, emcell / 2, absoThickness_ecal / 2); 
	  
	  absorberEMLV[i][j] = new G4LogicalVolume( absorberEM, absorberMaterial, absorberEMLV_name.str().c_str() );

	  absorberEMLV[i][j]->SetVisAttributes( ecalAbsorberCol );
	  
	  // Constructor of G4PVPlacement:
	  // 1: rotation = no rotation, 2: position,       3: logical volume,
	  // 4: name                    5: mother  volume, 6: boolean operation = no,
	  // 7: copy number,            8: checking overlaps
	  new G4PVPlacement( 0, G4ThreeVector(0.,0.,-gapThickness_ecal/2), absorberEMLV[i][j],
			     "Absorber", layerEMLV[i][j], false,
			     copy_num, fCheckOverlaps );

	  ////////////////////////////////////////////
	  // Gap                                    //
	  ////////////////////////////////////////////
	  std::stringstream gapEMLV_name;
	  gapEMLV_name << "GapEcal" << num.str();
	  G4VSolid* gapEM = new G4Box("Gap", emcell / 2, emcell / 2, gapThickness_ecal / 2); // name, size
	  gapEMLV[i][j] = new G4LogicalVolume( gapEM, gapMaterial, gapEMLV_name.str().c_str() );
	  gapEMLV[i][j]->SetVisAttributes( gapCol );

	  // Constructor of G4PVPlacement:
	  // 1: rotation = no rotation, 2: position,       3: logical volume,
	  // 4: name                    5: mother  volume, 6: boolean operation = no,
	  // 7: copy number,            8: checking overlaps
	  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, absoThickness_ecal / 2),
			     gapEMLV[i][j], "Gap", layerEMLV[i][j],
			     false, copy_num, fCheckOverlaps );
	}
    }
}

void B4cDetectorConstruction::DefineVolumes_Hcal( G4LogicalVolume* worldLV )
{
  G4LogicalVolume* calorLV[num_of_towers_hcal][num_of_towers_hcal];
  G4LogicalVolume* layerLV[num_of_towers_hcal][num_of_towers_hcal];
  G4LogicalVolume* absorberLV[num_of_towers_hcal][num_of_towers_hcal];
  G4LogicalVolume* gapLV[num_of_towers_hcal][num_of_towers_hcal];

  for( int i=0; i<num_of_towers_hcal; i++ )
    {
      for( int j=0; j<num_of_towers_hcal; j++ )
	{
	  std::stringstream num;
	  num  << std::setw(2) <<  std::setfill( '0' ) << i << std::setw(2) << std::setfill('0') << j;
	  
	  ////////////////////////////////////////////
	  // Calorimeter                            //
	  ////////////////////////////////////////////
	  std::stringstream calorLV_name;
	  calorLV_name << "Hcal" << num.str();

	  G4VSolid* calorimeterS
	    = new G4Box("Hcal",hadcell/2, hadcell/2, calorThickness_hcal/2); // name and size
	  calorLV[i][j]  = new G4LogicalVolume( calorimeterS, defaultMaterial, calorLV_name.str().c_str() ); // solid, material, and name
	  calorLV[i][j]->SetVisAttributes( simpleBoxVisAtt );

	  int copy_num = i * num_of_towers + j;
	  G4ThreeVector pos_cal
	    = G4ThreeVector( num_of_towers_hcal * hadcell / 2 - hadcell / 2 - j * hadcell,
			     num_of_towers_hcal * hadcell / 2 - hadcell / 2 - i * hadcell,
			     -calorThickness_hcal/2 );  // downstream side is at z=0 cm

	  new G4PVPlacement( 0, pos_cal,                // no rotation, postion
			     calorLV[i][j], "Hcal", worldLV, // logical volume, nameand mother volume
			     false, copy_num, fCheckOverlaps ); // no boolean operation, copy number, and checking overlaps


	  ////////////////////////////////////////////
	  // Layer                                  //
	  ////////////////////////////////////////////
	  std::stringstream layerLV_name;
	  layerLV_name << "LayerHcal" << num.str();

	  G4VSolid* layerS
	    = new G4Box("Layer", hadcell/2, hadcell/2, layerThickness_hcal/2); // name and size

	  layerLV[i][j]  = new G4LogicalVolume( layerS, defaultMaterial, layerLV_name.str().c_str()  ); // solid, material, and name
	  layerLV[i][j]->SetVisAttributes( simpleBoxVisAtt );
	  
	  new G4PVReplica( "Layer", layerLV[i][j], calorLV[i][j], // name, logical volume,and mother
			  kZAxis, nofLayers2, layerThickness_hcal );  // axis of replication, number of replica, and  witdth of replica

	  ////////////////////////////////////////////
	  // Absorber                               //
	  ////////////////////////////////////////////
	  std::stringstream absorberLV_name;
	  absorberLV_name << "AbsorberHcal" << num.str();

	  // name and size
	  G4VSolid* absorberS = new G4Box("Absorber", hadcell / 2, hadcell / 2, absoThickness_hcal / 2); 

	  absorberLV[i][j]  = new G4LogicalVolume( absorberS, absorberMaterial2, absorberLV_name.str().c_str() );  // solid, material, and name 
	  absorberLV[i][j]->SetVisAttributes( hcalAbsorberCol );

	  // Constructor of G4PVPlacement:
	  // 1: rotation = no rotation, 2: position,       3: logical volume,
	  // 4: name                    5: mother  volume, 6: boolean operation = no,
	  // 7: copy number,            8: checking overlaps
	  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, -gapThickness_hcal / 2), absorberLV[i][j],
			     "Absorber", layerLV[i][j], false,
			     copy_num, fCheckOverlaps );


	  ////////////////////////////////////////////
	  // Gap                                    //
	  ////////////////////////////////////////////
	  std::stringstream gapLV_name;
	  gapLV_name << "GapHcal" << num.str();

	  // name and size
	  G4VSolid* gapS = new G4Box("Gap", hadcell / 2, hadcell / 2, gapThickness_hcal / 2);

	  // solid, material, and name
	  gapLV[i][j]  = new G4LogicalVolume( gapS, gapMaterial, gapLV_name.str().c_str() );
	  gapLV[i][j]->SetVisAttributes( gapCol );


	  // Constructor of G4PVPlacement:
	  // 1: rotation = no rotation, 2: position,       3: logical volume,
	  // 4: name                    5: mother  volume, 6: boolean operation = no,
	  // 7: copy number,            8: checking overlaps
	  new G4PVPlacement( 0, G4ThreeVector(0.0,  0.0,  absoThickness_hcal / 2), gapLV[i][j],
			     "Gap", layerLV[i][j], false,
			     copy_num, fCheckOverlaps );
	  
	}
    }

}

void B4cDetectorConstruction::DefineVolumes_Wall( G4LogicalVolume* worldLV )
{
  // size of cell, the naming is bit strange because there is always only 1 cell for the wall
  G4double wallcell = num_of_towers_hcal * hadcell;

  // name and size
  G4VSolid* wallS
    = new G4Box("Wall", wallcell / 2, wallcell/2, wallThickness / 2);

  // solid, material, and name
  G4LogicalVolume* wallLV = new G4LogicalVolume( wallS, absorberMaterial2, "Wall" );
  wallLV->SetVisAttributes( wallCol );

  G4ThreeVector pos_wall
    = G4ThreeVector( 0, 0, -calorThickness_hcal - wallThickness / 2 );
  
  // Constructor of G4VPhysicalVolume
  // 1: no rotation, 2: position at (0,0,0), 3: logical volume,
  // 4: name,        5: mother  volume,      6: no boolean operation,
  // 7: copy number, 8: checking overlaps
  new G4PVPlacement( 0, pos_wall, wallLV,
		     "Wall", worldLV, false,
		     0, fCheckOverlaps );
  
}

void B4cDetectorConstruction::DefineVolumes_Surroundings( G4LogicalVolume* worldLV )
{
  // top-bottom: y-axis
  // left-right: x-axis

  //////////////////////////////////////////////////////////////////////////////////
  // Surrounding Hcal                                                             //
  //////////////////////////////////////////////////////////////////////////////////
  G4double transverse_hcal = num_of_towers_hcal * hadcell;

  G4VSolid* surrounding_hcalS
    = new G4Box("Surrounding_hcal",
		transverse_hcal / 2 + surroundingsThickness / 2, // transverse size and some to cover fully
		surroundingsThickness / 2, // thickness
		calorThickness_hcal / 2); // length

  // make logical volumes
  G4LogicalVolume* surroundings_hcal_topLV    = new G4LogicalVolume( surrounding_hcalS, stopperMaterial, "Surrounding_hcal_top"    );
  G4LogicalVolume* surroundings_hcal_bottomLV = new G4LogicalVolume( surrounding_hcalS, stopperMaterial, "Surrounding_hcal_bottom" );
  G4LogicalVolume* surroundings_hcal_leftLV   = new G4LogicalVolume( surrounding_hcalS, stopperMaterial, "Surrounding_hcal_left"   );
  G4LogicalVolume* surroundings_hcal_rightLV  = new G4LogicalVolume( surrounding_hcalS, stopperMaterial, "Surrounding_hcal_right"  );

  surroundings_hcal_topLV   ->SetVisAttributes( surroundingsCol );
  surroundings_hcal_bottomLV->SetVisAttributes( surroundingsCol );
  surroundings_hcal_leftLV  ->SetVisAttributes( surroundingsCol );
  surroundings_hcal_rightLV ->SetVisAttributes( surroundingsCol );

  // make vector for position
  G4ThreeVector pos_surrounding_hcal_top
    = G4ThreeVector( -surroundingsThickness / 2,
		     transverse_hcal / 2 + surroundingsThickness / 2,
		     -calorThickness_hcal / 2 );

  G4ThreeVector pos_surrounding_hcal_bottom
    = G4ThreeVector( surroundingsThickness / 2,
		     -transverse_hcal / 2 - surroundingsThickness / 2,
		     -calorThickness_hcal / 2 );

  G4ThreeVector pos_surrounding_hcal_left
    = G4ThreeVector( transverse_hcal / 2 + surroundingsThickness / 2,
		     surroundingsThickness / 2 , 
		     -calorThickness_hcal / 2 );

  G4ThreeVector pos_surrounding_hcal_right
    = G4ThreeVector( -transverse_hcal / 2 - surroundingsThickness / 2,
		     -surroundingsThickness / 2 , 
		     -calorThickness_hcal / 2 );

  // put the block!
  // top
  new G4PVPlacement( 0, pos_surrounding_hcal_top, surroundings_hcal_topLV,
		     "Surrounding_hcal_top", worldLV, false,
		     0, fCheckOverlaps );


  // bottom
  new G4PVPlacement( 0, pos_surrounding_hcal_bottom, surroundings_hcal_bottomLV,
		     "Surrounding_hcal_bottom", worldLV, false,
		     0, fCheckOverlaps );

  // Rotate the box for the top and bottom by 90 degree
  // so that the box can be used for the right and left sides!
  G4RotationMatrix* rot = new G4RotationMatrix();
  *rot = G4RotationMatrix();
  rot->rotateZ( 90 * deg );

  // left
  new G4PVPlacement( rot, pos_surrounding_hcal_left, surroundings_hcal_topLV,
		     "Surrounding_hcal_left", worldLV, false,
		     0, fCheckOverlaps );

  // right
  new G4PVPlacement( rot, pos_surrounding_hcal_right, surroundings_hcal_topLV,
		     "Surrounding_hcal_right", worldLV, false,
		     0, fCheckOverlaps );
  

  //////////////////////////////////////////////////////////////////////////////////
  // Surrounding Ecal                                                             //
  //////////////////////////////////////////////////////////////////////////////////
  G4double transverse_ecal = num_of_towers_ecal * hadcell;
  G4VSolid* surrounding_ecalS
    = new G4Box("Surrounding_ecal",
		transverse_ecal / 2 + surroundingsThickness / 2, // transverse size and some to cover fully
		surroundingsThickness / 2, // thickness
		calorThickness_ecal / 2  + wallThickness / 2); // length

  // make logical volumes
  G4LogicalVolume* surroundings_ecal_topLV
    = new G4LogicalVolume( surrounding_ecalS, stopperMaterial, "Surrounding_ecal_top"    );

  G4LogicalVolume* surroundings_ecal_bottomLV
    = new G4LogicalVolume( surrounding_ecalS, stopperMaterial, "Surrounding_ecal_bottom" );

  G4LogicalVolume* surroundings_ecal_leftLV
    = new G4LogicalVolume( surrounding_ecalS, stopperMaterial, "Surrounding_ecal_left"   );

  G4LogicalVolume* surroundings_ecal_rightLV
    = new G4LogicalVolume( surrounding_ecalS, stopperMaterial, "Surrounding_ecal_right"  );

  surroundings_ecal_topLV   ->SetVisAttributes( surroundingsCol );
  surroundings_ecal_bottomLV->SetVisAttributes( surroundingsCol );
  surroundings_ecal_leftLV  ->SetVisAttributes( surroundingsCol );
  surroundings_ecal_rightLV ->SetVisAttributes( surroundingsCol );

  // position in z-axis, depend on thickness of the wall
  G4double z_pos =  -calorThickness_hcal - calorThickness_ecal / 2  - wallThickness / 2;

  // make vector for position
  G4ThreeVector pos_surrounding_ecal_top
    = G4ThreeVector( -surroundingsThickness / 2,
		     transverse_ecal / 2 + surroundingsThickness / 2,
		     z_pos);

  G4ThreeVector pos_surrounding_ecal_bottom
    = G4ThreeVector( surroundingsThickness / 2,
		     -transverse_ecal / 2 - surroundingsThickness / 2,
		     z_pos );

  G4ThreeVector pos_surrounding_ecal_left
    = G4ThreeVector( transverse_ecal / 2 + surroundingsThickness / 2,
		     surroundingsThickness / 2, 
		     z_pos );
  
  G4ThreeVector pos_surrounding_ecal_right
    = G4ThreeVector( -transverse_ecal / 2 - surroundingsThickness / 2,
		     -surroundingsThickness / 2 , 
		     z_pos );  // downstream side is at z=0 cm

  // put the blocks!
  // top
  new G4PVPlacement( 0, pos_surrounding_ecal_top, surroundings_ecal_topLV,
		     "Surrounding_ecal_top", worldLV, false,
		     0, fCheckOverlaps );


  // bottom
  new G4PVPlacement( 0, pos_surrounding_ecal_bottom, surroundings_ecal_bottomLV,
		     "Surrounding_ecal_bottom", worldLV, false,
		     0, fCheckOverlaps );

  // left, the same rotation as ones for Hcal used
  new G4PVPlacement( rot, pos_surrounding_ecal_left, surroundings_ecal_topLV,
		     "Surrounding_ecal_left", worldLV, false,
		     0, fCheckOverlaps );

  // right, the same rotation as ones for Hcal used
  new G4PVPlacement( rot, pos_surrounding_ecal_right, surroundings_ecal_topLV,
		     "Surrounding_ecal_right", worldLV, false,
		     0, fCheckOverlaps );

  //////////////////////////////////////////////////////////////////////////////////
  // Lids                                                                         //
  ///////////////////////////////
  //   Geometry of the lids    //
  ///////////////////////////////
  //                           //
  //  --------      --------   //
  // |   --   |    |        |  //   
  // |  |  |  |    |        |  //
  // |   --   |    |        |  //
  //  --------      --------   //
  //  Upstream     Downstream  //
  //                           //
  //  Seen along z-axis        //
  ///////////////////////////////
  
  // normal box for the downstream side
  G4VSolid* lid_boxS
    = new G4Box("Lid_box",
		transverse_hcal / 2 + surroundingsThickness, // transverse size and some to cover fully
		transverse_hcal / 2 + surroundingsThickness, // transverse size and some to cover fully
		lidThickness / 2); // length
  
  G4double hole_size = 1 * cm;
  G4VSolid* removalS
    = new G4Box("Removal", hole_size, hole_size, lidThickness / 2); // length

  // for the upstream side
  G4VSolid* lidS
    = new G4SubtractionSolid( "Lid", lid_boxS, removalS, 0, G4ThreeVector( 0.0, 0.0, 0.0 ) );
  
  G4LogicalVolume* lid_upstreamLV   = new G4LogicalVolume( lidS, stopperMaterial, "Lid_upstream"    );
  G4LogicalVolume* lid_downstreamLV = new G4LogicalVolume( lid_boxS, stopperMaterial, "Lid_downstream"  );

  lid_upstreamLV   ->SetVisAttributes( surroundingsCol );
  lid_downstreamLV ->SetVisAttributes( surroundingsCol );

  G4double lid_upstream_z = -calorThickness_hcal - calorThickness_ecal - wallThickness - lidThickness / 2;
  G4double lid_downstream_z = lidThickness / 2;

  G4ThreeVector pos_lid_upstream    = G4ThreeVector( 0.0, 0.0, lid_upstream_z   );
  G4ThreeVector pos_lid_downstream  = G4ThreeVector( 0.0, 0.0, lid_downstream_z );

  new G4PVPlacement( 0, pos_lid_upstream, lid_upstreamLV,
		     "Lid_upstream", worldLV, false,
		     0, fCheckOverlaps );
  
  new G4PVPlacement( 0, pos_lid_downstream, lid_downstreamLV,
		     "Lid_downstream", worldLV, false,
		     0, fCheckOverlaps );
  
}

/*
void B4cDetectorConstruction::DefineVolumes_NotUsed( G4LogicalVolume* worldLV )
{
  // Calorimeter
  G4LogicalVolume* calor2LV[num_of_towers_hcal][num_of_towers_hcal];
  char bu[200];
  char bu2[200];
  char bu3[200];

  for( int i=0; i<num_of_towers_hcal; i++ )
    {
      for( int j=0; j<num_of_towers_hcal; j++ )
	{
	  sprintf( bu, "Calorimeter1%02d%02d", i, j );
	  //sprintf( bu2, "calor%02d%02d", i, j );
	  sprintf( bu3, "Calorimeter2%02d%02d", i, j);

	  G4VSolid* calorimeterS
	    = new G4Box("Calorimeter",hadcell/2, hadcell/2, calorThickness_hcal/2); // name and size
	  calor2LV[i][j] = new G4LogicalVolume( calorimeterS, defaultMaterial, bu3);

	  G4ThreeVector pos_cal2 = G4ThreeVector( -num_of_towers_hcal * hadcell / 2 + hadcell / 2 + i * hadcell,
						  num_of_towers_hcal * hadcell / 2 - hadcell / 2 - j * hadcell,
						  //						  3 * calorThickness_hcal / 2 + calorThickness_ecal / 2 );
						  calorThickness_hcal / 2 );

	  new G4PVPlacement(0, pos_cal2, 
	   		    calor2LV[i][j], bu3, worldLV,
	   		    false, i+j, fCheckOverlaps );
	  
	} // end of j-loop
    } // end of i-loop


  // Layer
  G4LogicalVolume* layer2LV[num_of_towers_hcal][num_of_towers_hcal];

  for( int i=0; i<num_of_towers_hcal; i++ )
    {

      for( int j=0; j<num_of_towers_hcal; j++ )
	{
	  sprintf( bu , "Layer1%02d%02d" , i, j );
	  sprintf( bu2, "Layer2%02d%02d", i, j );
	  sprintf( bu3, "Layer3%02d%02d", i, j );
	  G4VSolid* layerS
	    = new G4Box("Layer", hadcell/2, hadcell/2, layerThickness_hcal/2); // name and size

	  layer2LV[i][j] = new G4LogicalVolume( layerS, defaultMaterial, bu2 );
	  
	  new G4PVReplica("Layer2", layer2LV[i][j], calor2LV[i][j],
			  kZAxis, nofLayers2, layerThickness_hcal );
	}
    }


  // Absorber
  G4LogicalVolume* absorberLV2[num_of_towers_hcal][num_of_towers_hcal];

  for( int i=0; i<num_of_towers_hcal; i++ )
    {
      for( int j=0; j<num_of_towers_hcal; j++ )
	{
	  sprintf(bu,"Abso1%02d%02d",i,j);
	  sprintf(bu2,"Abso2%02d%02d",i,j);
	  sprintf(bu3,"Abso2%02d%02d",i,j);

	  G4VSolid* absorberS = new G4Box("Abso", hadcell / 2, hadcell / 2, absoThickness_hcal / 2); // name and size

	  absorberLV2[i][j] = new G4LogicalVolume( absorberS, absorberMaterial2, bu2 );

	  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, -gapThickness_hcal / 2),
			     absorberLV2[i][j], "Abso2", layer2LV[i][j],
			     false, i * num_of_towers_hcal + j, fCheckOverlaps );

	}
    }

  // Gap
  char buffer [200];
  G4LogicalVolume* gapLV2[num_of_towers_hcal][num_of_towers_hcal];

  for( int i=0; i<num_of_towers_hcal; i++ )
    {
      for( int j=0; j<num_of_towers_hcal; j++ )
	{
	  sprintf( buffer, "gap1%02d%02d" , i, j );
	  sprintf( bu2   , "gap2%02d%02d", i, j );
	  sprintf( bu3   , "gap3%02d%02d", i, j );

	  G4VSolid* gapS = new G4Box("Gap", hadcell / 2, hadcell / 2, gapThickness_hcal / 2); // name and size

	  gapLV2[i][j] = new G4LogicalVolume( gapS, gapMaterial, bu2 );

	  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, absoThickness_hcal / 2),
			     gapLV2[i][j], bu2, layer2LV[i][j],
			     false, i * num_of_towers_hcal + j, fCheckOverlaps );
	}
    }


  
  for( int i=0; i<num_of_towers_hcal; i++ )
    {
      for( int j=0; j<num_of_towers_hcal; j++ )
	{
	  calor2LV[i][j]->SetVisAttributes( simpleBoxVisAtt );

	  layer2LV[i][j]->SetVisAttributes( simpleBoxVisAtt );

	  absorberLV2[i][j]->SetVisAttributes( hcalAbsorberCol2 );

	  gapLV2[i][j]->SetVisAttributes( gapCol );
	}
    }

}
*/

void B4cDetectorConstruction::DefineVisAttributes()
{

  simpleBoxVisAtt= new G4VisAttributes( G4Colour(1,0,0, 0) );
  simpleBoxVisAtt->SetVisibility( false );

  ecalAbsorberCol = new G4VisAttributes( true , G4Colour( 0.725, 0.403, 0.780, 0.7 ) ); // Pb
  ecalAbsorberCol->SetForceSolid( true );

  hcalAbsorberCol = new G4VisAttributes( true , G4Colour( 0.238, 0.732, 0.827, 0.9 ) ); // Fe
  hcalAbsorberCol->SetForceSolid( true );

  gapCol = new G4VisAttributes( true , G4Colour( 1, 1, 1, 1 ) ); // polystyrene
  gapCol->SetForceSolid( true );
  //  gapCol->SetForceWireframe( true );

  wallCol = new G4VisAttributes( true , G4Colour( 1, 1, 1, 1 ) ); // polystyrene
  wallCol->SetForceSolid( true );

  //  surroundingsCol = new G4VisAttributes( true , G4Colour( 1, 1, 1, 1 ) ); // white
  surroundingsCol = new G4VisAttributes( true , G4Colour( 0.0, 0.0, 0.0, 0.9 ) ); 
  //  surroundingsCol = new G4VisAttributes( true , G4Colour( 0.4, 0.4, 0.4, 1 ) ); 
  surroundingsCol->SetForceWireframe( true );
  //  surroundingsCol->SetForceSolid( true );

}

//void B4cDetectorConstruction::SetMagField(G4double fieldValue)
void B4cDetectorConstruction::ConstructSDandField()
{

  // Sensitive detectors
  B4cCalorimeterSD* absoSD[num_of_towers_hcal][num_of_towers_hcal];
  B4cCalorimeterSD* gapSD[num_of_towers_hcal][num_of_towers_hcal];
  B4cCalorimeterSD* abso2SD[num_of_towers_hcal][num_of_towers_hcal];
  B4cCalorimeterSD* gap2SD[num_of_towers_hcal][num_of_towers_hcal];
  //  B4cCalorimeterSD* abso3SD[num_of_towers_hcal][num_of_towers_hcal];
  //  B4cCalorimeterSD* gap3SD[num_of_towers_hcal][num_of_towers_hcal];

  ////////////////////////////////////////////
  // Hcal                                   //
  ////////////////////////////////////////////
  for( int i=0; i<num_of_towers_hcal; i++ )
    {
      for( int j=0; j<num_of_towers_hcal; j++ )
	{
	  std::stringstream num;
	  num  << std::setw(2) <<  std::setfill( '0' ) << i << std::setw(2) << std::setfill('0') << j;

	  std::stringstream absorber_LV, absorber_HC, absorber_SD, gap_LV, gap_HC, gap_SD;	      
	  absorber_LV << "AbsorberHcal" << num.str();
	  absorber_HC << "Absorber" << "HitsCollection" << num.str();
	  absorber_SD << "Absorber" << "SD" << num.str();
	      	      
	  absoSD[i][j] = new B4cCalorimeterSD( absorber_SD.str().c_str(), absorber_HC.str().c_str(), nofLayers2 );
	  G4SDManager::GetSDMpointer()->AddNewDetector( absoSD[i][j] );
	  SetSensitiveDetector( absorber_LV.str().c_str(), absoSD[i][j] );

	  gap_LV << "GapHcal" << num.str();
	  gap_HC << "Gap" << "HitsCollection" << num.str();
	  gap_SD << "Gap" << "SD" << num.str();
	  
	  gapSD[i][j] = new B4cCalorimeterSD( gap_SD.str().c_str(), gap_HC.str().c_str(), nofLayers2 );
	  G4SDManager::GetSDMpointer()->AddNewDetector( gapSD[i][j] );
	  SetSensitiveDetector( gap_LV.str().c_str(), gapSD[i][j] );
	      
	} // end of the loop j
    } // end of the loop i
  
  ////////////////////////////////////////////
  // Ecal                                   //
  ////////////////////////////////////////////
  B4cCalorimeterSD* absoEMSD[num_of_towers_ecal][num_of_towers_ecal];
  B4cCalorimeterSD* gapEMSD[num_of_towers_ecal][num_of_towers_ecal];

  for( int i=0; i<num_of_towers_ecal; i++ )
    {
      for( int j=0; j<num_of_towers_ecal; j++ )
	{
	  std::stringstream num;
	  num  << std::setw(2) <<  std::setfill( '0' ) << i << std::setw(2) << std::setfill('0') << j;

	  std::stringstream absorber_LV, absorber_HC, absorber_SD, gap_LV, gap_HC, gap_SD;

	  absorber_LV << "AbsorberEcal" << num.str();
	  absorber_HC << "AbsorberEMHitsCollection" << num.str();
	  absorber_SD << "AbsorberEMSD" << num.str();

	  gap_LV << "GapEcal" << num.str();
	  gap_HC << "GapEMHitsCollection" << num.str();
	  gap_SD << "GapEMSD" << num.str();

	  // for absorbers
	  absoEMSD[i][j] = new B4cCalorimeterSD( absorber_SD.str().c_str(), absorber_HC.str().c_str(), nofLayers ); // name, hits coolection name, no of cells
	  G4SDManager::GetSDMpointer()->AddNewDetector( absoEMSD[i][j] );
	  SetSensitiveDetector( absorber_LV.str().c_str(), absoEMSD[i][j] );

	  // for gaps
	  gapEMSD[i][j] = new B4cCalorimeterSD( gap_SD.str().c_str(), gap_HC.str().c_str(), nofLayers );
	  G4SDManager::GetSDMpointer()->AddNewDetector( gapEMSD[i][j] );
	  SetSensitiveDetector( gap_LV.str().c_str() , gapEMSD[i][j]);
	}

    }

  ////////////////////////////////////////////
  // Wall                                   //
  ////////////////////////////////////////////
  if( wallThickness > 0 )
    {
      B4cCalorimeterSD* wall_SD = new B4cCalorimeterSD( "Wall_SD", "Wall_HitsCollection", 1 );
      G4SDManager::GetSDMpointer()->AddNewDetector( wall_SD );
      SetSensitiveDetector( "Wall", wall_SD );
    }
  
  ////////////////////////////////////////////
  // Surroundings                           //
  ////////////////////////////////////////////
  // if the blocks are off, here is the end
  if( isSurroundings == false )
    return;

  const G4int direction_num = 4;
  G4String directions[ direction_num ] = { "top", "bottom", "left", "right" };  
  B4cCalorimeterSD* surroundings_ecalSD[ direction_num ];
  B4cCalorimeterSD* surroundings_hcalSD[ direction_num ];

  // loop over 4 directions
  for( int i=0; i<direction_num; i++ )
    {

      // loop over 2 calorimeters
      G4String calos[2] = { "ecal", "hcal" };
      for( int j=0; j<2; j++ )
	{

	  // make name of logical volume, sensitive detector, and hits collection
	  std::stringstream cal_LV, cal_SD, cal_HC;
	  cal_LV << "Surrounding_" << calos[j] << "_" << directions[i];
	  cal_SD << cal_LV.str() << "_SD";
	  cal_HC << cal_LV.str() << "_HitsCollection";

	  // set the pointer of the block depending on j
	  B4cCalorimeterSD* calSD;
	  if( j == 0 )
	    calSD = surroundings_ecalSD[i];
	  else
	    calSD = surroundings_hcalSD[i];

	  // make sensitive detector
	  calSD = new B4cCalorimeterSD( cal_SD.str().c_str(), cal_HC.str().c_str(), 1 );	  
	  G4SDManager::GetSDMpointer()->AddNewDetector( calSD );
	  SetSensitiveDetector( cal_LV.str().c_str(), calSD );
	} // end of loop over detectors
    } // end of loop over directions

  B4cCalorimeterSD* lid_upstream_SD = new B4cCalorimeterSD( "Lid_upstream_SD", "Lid_upstream_HitsCollection", 1 );
  G4SDManager::GetSDMpointer()->AddNewDetector( lid_upstream_SD );
  SetSensitiveDetector( "Lid_upstream", lid_upstream_SD );

  B4cCalorimeterSD* lid_downstream_SD = new B4cCalorimeterSD( "Lid_downstream_SD", "Lid_downstream_HitsCollection", 1 );
  G4SDManager::GetSDMpointer()->AddNewDetector( lid_downstream_SD );
  SetSensitiveDetector( "Lid_downstream", lid_downstream_SD );


  /*
  // Apply a global uniform magnetic field along X axis
  G4FieldManager* fieldManager
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  // Delete the existing magnetic field
  if ( fMagField )  delete fMagField;

  if ( fieldValue != 0. )
    {
      // create a new one if not null
      fMagField
	= new G4UniformMagField(G4ThreeVector(fieldValue, 0., 0.));
      
      fieldManager->SetDetectorField(fMagField);
      fieldManager->CreateChordFinder(fMagField);
    }
  else
    {
      fMagField = 0;
      fieldManager->SetDetectorField(fMagField);
    }
  */
}


void B4cDetectorConstruction::SetWallThickness( G4double thickness )
{
  std::cout << std::string( 100, '-' ) << std::endl;
  std::cout << wallThickness / cm << std::endl;
  wallThickness = thickness;
  G4RunManager::GetRunManager()->GeometryHasBeenModified( true );

  std::cout << "modified?" << std::endl;
  std::cout << wallThickness / cm << std::endl;
  std::cout << std::string( 100, '-' ) << std::endl;
  //  Construct();
}
