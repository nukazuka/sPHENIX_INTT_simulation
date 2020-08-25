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
    wallThickness(0.0),
    whatever( nullptr )
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
  // G4GenericMessenger::Command& setWallThickness
  //   //    = fMessenger->DeclarePropertyWithUnit( "wallThickness", "cm", wallThickness );
  //   = fMessenger->DeclareMethodWithUnit( "wallThickness", "cm",
  // 					 &B4cDetectorConstruction::SetWallThickness,
  // 					 "Set thickness of the wall in cm" );
  // setWallThickness.SetGuidance( "Give the thickness of the wall. Negative value and 0 remove the wall" );
  // setWallThickness.SetParameterName( "WallThickness", false ); // (name, is_omittable)
  // setWallThickness.SetDefaultValue( 0.0 );
  
  // Define /B4/det/setMagField command
  //  G4GenericMessenger::Command& setMagFieldCmd
  //    = fMessenger->DeclareMethod("setMagField",
  //                                &B4cDetectorConstruction::SetMagField,
  //                                "Define magnetic field value (in X direction");
  //  setMagFieldCmd.SetUnitCategory("Magnetic flux density");


  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

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
  nistManager->FindOrBuildMaterial( "G4_Cu", fromIsotopes );
  nistManager->FindOrBuildMaterial( "G4_Si", fromIsotopes );
  nistManager->FindOrBuildMaterial( "G4_KAPTON", fromIsotopes );

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
  //  surroundingsCol->SetForceWireframe( true );
  //  surroundingsCol->SetForceSolid( true );


  wireframe_vis_att = new G4VisAttributes( true, G4Colour(0, 0, 0, 1) );
  wireframe_vis_att->SetForceWireframe( true );

  cu_vis_att = new G4VisAttributes( true, G4Colour( 0.7, 0.4, 0, 1) );
  cu_vis_att->SetForceSolid( true );

  kapton_vis_att = new G4VisAttributes( true, G4Colour(0.0, 0.590, 1.0, 0.5 ) ); // blue
  kapton_vis_att->SetForceSolid( true );

  si_vis_att = new G4VisAttributes( true, G4Colour(0, 0, 1, 0.5) ); // transparent blue
  si_vis_att->SetForceSolid( true );

  fphx_vis_att = new G4VisAttributes( true, G4Colour(1.0, 0.843, 0.0, 0.5) ); // HTML gold
  fphx_vis_att->SetForceSolid( true );  
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

  // Coordinate :  coresoftware  | this
  // x : thickness, vertical     | beam
  // y : transverse              | transverse
  // z : beam                    | vertical
  
  // int param
  const int kNstrips_phi_cell =   256;
  const int kNstrips_phi_sensor = 256;
  //const int kNstrips_z_sensor_0 =   8;
  const int kNstrips_x_sensor_0 =   8;
  //const int kNstrips_z_sensor_1 =   5;
  const int kNstrips_x_sensor_1 =   5;

  // double param
  //const double kFphx_x		= 0.032  * cm; // coresoftware
  const double kFphx_z		= 0.032  * cm;
  const double kFphx_y		= 0.27   * cm;
  //const double kFphx_z		= 0.91   * cm;
  const double kFphx_x		= 0.91   * cm;
  //const double kFphx_offset_z	= 0.005  * cm;
  const double kFphx_offset_x	= 0.005  * cm;
  const double kGap_sensor_fphx = 0.1    * cm;  

  //  const double kSi_glue_x	= 0.0014 * cm;
  const double kSi_glue_z	= 0.0014 * cm;
  //const double kFphx_glue_x	= 0.005  * cm;
  const double kFphx_glue_z	= 0.005  * cm;

  //const double kHdi_copper_x	= 0.0037 * cm;
  const double kHdi_copper_z	= 0.0037 * cm;
  //const double kHdi_edge_z	= 0.0    * cm;
  const double kHdi_edge_x	= 0.0    * cm;
  //  const double kHdi_kapton_x	= 0.038  * cm;
  const double kHdi_kapton_z	= 0.038  * cm;
  const double kHdi_y		= 3.8    * cm;

  const double kSensor_edge_phi = 0.13   * cm;
  //const double kSensor_edge_z	= 0.1    * cm;
  const double kSensor_edge_x	= 0.1    * cm;
  const double kSensor_offset_y = 0.     * cm;  

  //const double kStrip_x		= 0.032  * cm;
  const double kStrip_z		= 0.032  * cm;
  const double kStrip_y		= 0.007  * cm;
  //const double kStrip_z_0	= 1.6    * cm;
  const double kStrip_x_0	= 1.6    * cm;
  //const double kStrip_z_1	= 2.0    * cm;
  const double kStrip_x_1	= 2.0    * cm;

  // Ratation matrix to change the coordinate in the coresoftware repo to this codes
  //G4RotationMatrix* rot_core2this = new G4RotationMatrix();
  //rot_core2this->rotateY( 90 * deg );

    /*
  // SEGMENTATION_PHI //////////////////////////////////////
  
  const double kHalfladder_z =        40.00 ;
  const double kHalfladder_inside_z = 23.9622;

  const double kStave_straight_cooler_x = 0.03  ;
  const double kStave_straight_cooler_y = 1.47684;
  const double kStave_slant_cooler_y =   0.6322614829;
  const double kStave_straight_outer_y =    0.3322;
  const double kStave_straight_rohacell_y = 0.5884;
  */  


  //     [][][][]                          [][][][]
  // [][][]                                      [][][]
  // []                                              []
  // []              BUILD EACH ELEMENT              []
  // []                                              []
  // [] [][] []     [][]     []     [][]     [] [][] []
  // [] [][] [][][]   [][] [][][] [][]   [][][] [][] []
    
  //////////////////////////////////////////////////////////
  // Silicon sensors                                      //
  //////////////////////////////////////////////////////////
  // Si strip
  // +-------------+
  // +-------------+
  const double kSi_strip_0_x = kStrip_x_0;
  const double kSi_strip_1_x = kStrip_x_1;
  const double kSi_strip_0_y = kStrip_y;// * kNstrips_phi_sensor / 2;
  const double kSi_strip_0_z = kStrip_z;  // thickness

  G4VSolid* si_strip_0  = new G4Box( "Si_strip_0", kSi_strip_0_x / 2, kSi_strip_0_y / 2, kSi_strip_0_z / 2 );
  G4LogicalVolume* si_strip_0_lv = new G4LogicalVolume( si_strip_0, G4Material::GetMaterial( "G4_Si"), "Si_strip_0" );
  si_strip_0_lv->SetVisAttributes( si_vis_att );

  G4VSolid* si_strip_1  = new G4Box( "Si_strip_1", kSi_strip_1_x / 2, kSi_strip_0_y / 2, kSi_strip_0_z / 2 );
  G4LogicalVolume* si_strip_1_lv = new G4LogicalVolume( si_strip_1, G4Material::GetMaterial( "G4_Si" ), "Si_strip_1" );
  si_strip_1_lv->SetVisAttributes( si_vis_att );
  
  // Si cell, 2 types 
  //
  // +-------------+
  // |- - - - - - -| 
  // |- - - - - - -| 
  // |- - - - - - -| 
  // |- - - - - - -| 
  // +-------------+
  //
  const double kSi_cell_0_x = kStrip_x_0;
  const double kSi_cell_1_x = kStrip_x_1;
  const double kSi_cell_0_y = kStrip_y * kNstrips_phi_sensor / 2;
  const double kSi_cell_0_z = kStrip_z;  // thickness

  G4VSolid* si_cell_0  = new G4Box( "Si_cell_0", kSi_cell_0_x / 2, kSi_cell_0_y / 2, kSi_cell_0_z / 2 );
  G4LogicalVolume* si_cell_0_lv = new G4LogicalVolume( si_cell_0, defaultMaterial, "Si_cell_0" );
  si_cell_0_lv->SetVisAttributes( wireframe_vis_att );

  G4VSolid* si_cell_1  = new G4Box( "Si_cell_1", kSi_cell_1_x / 2, kSi_cell_0_y / 2, kSi_cell_0_z / 2 );
  G4LogicalVolume* si_cell_1_lv = new G4LogicalVolume( si_cell_1, defaultMaterial, "Si_cell_1" );
  si_cell_1_lv->SetVisAttributes( wireframe_vis_att );

  // Si area
  // +-------------+-------------+-------------+-------------+-------------+
  // |- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|
  // |- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|
  // |- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|
  // |- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|
  // +-------------+-------------+-------------+-------------+-------------+
  // |- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|
  // |- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|
  // |- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|
  // |- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|- - - - - - -|
  // +-------------+-------------+-------------+-------------+-------------+
  const double kSi_area_0_x = kNstrips_x_sensor_0 * kSi_cell_0_x;
  const double kSi_area_1_x = kNstrips_x_sensor_1 * kSi_cell_1_x;
  const double kSi_area_0_y = kSi_cell_0_y;
  const double kSi_area_1_y = kSi_area_0_y;;
  const double kSi_area_0_z = kSi_cell_0_z;
  const double kSi_area_1_z = kSi_area_0_z;
  
  G4VSolid* si_area_0  = new G4Box( "Si_area_0", kSi_area_0_x / 2, kSi_area_0_y / 2, kSi_area_0_z / 2 );				    
  G4LogicalVolume* si_area_0_lv = new G4LogicalVolume( si_area_0, defaultMaterial, "Si_area_0" );
  si_area_0_lv->SetVisAttributes( wireframe_vis_att );

  G4VSolid* si_area_1  = new G4Box( "Si_area_1", kSi_area_1_x / 2, kSi_area_1_y / 2, kSi_area_1_z / 2 );				    
  G4LogicalVolume* si_area_1_lv = new G4LogicalVolume( si_area_1, defaultMaterial, "Si_area_1" );
  si_area_1_lv->SetVisAttributes( wireframe_vis_att );
  
  // Fill si strips into the Si cell
  new G4PVReplica( "Si_strip_replica_0", si_strip_0_lv, si_cell_0_lv, kYAxis, kNstrips_phi_sensor / 2, kStrip_y );
  new G4PVReplica( "Si_strip_replica_1", si_strip_1_lv, si_cell_1_lv, kYAxis, kNstrips_phi_sensor / 2, kStrip_y );

  // Fill si cells into the Si area
  new G4PVReplica( "Si_cell_replica_0", si_cell_0_lv, si_area_0_lv, kXAxis, kNstrips_x_sensor_0, kStrip_x_0 );
  new G4PVReplica( "Si_cell_replica_1", si_cell_1_lv, si_area_1_lv, kXAxis, kNstrips_x_sensor_1, kStrip_x_1 );
  
  //////////////////////////////////////////////////////////
  // HDI                                                  //
  //////////////////////////////////////////////////////////
  
  // HDI: one for half ladder
  //const double kHdi_y = kStrip_y * kNstrips_phi_sensor; // 1.9968 cm
  const double kHdi_x = kStrip_x_0 * kNstrips_x_sensor_0 + 2 * kSensor_edge_x
    + kStrip_x_1 * kNstrips_x_sensor_1 + 2 * kSensor_edge_x; // 23.2 cm

  const double kHdi_0_x = kStrip_x_0 * kNstrips_x_sensor_0 + 2 * kSensor_edge_x;
  const double kHdi_1_x = kStrip_x_1 * kNstrips_x_sensor_1 + 2 * kSensor_edge_x;
      
  // std::cout << std::string( 100, '=' ) << std::endl;
  // std::cout << "kHdi_x = " << kHdi_x << std::endl;
  // std::cout << "kHdi_y = " << kHdi_y << std::endl;
  // std::cout << std::string( 100, '=' ) << std::endl;

  // type-1
  // Kapton layer of HDI
  G4VSolid* hdi_0_kapton = new G4Box( "HDI_0_kapton", kHdi_0_x / 2, kHdi_y / 2, kHdi_kapton_z / 2 );
  G4LogicalVolume* hdi_0_kapton_lv = new G4LogicalVolume( hdi_0_kapton, G4Material::GetMaterial( "G4_KAPTON" ), "HDI_0_Kapton" );
  hdi_0_kapton_lv->SetVisAttributes( kapton_vis_att );

  // Copper layer of HDI
  G4VSolid* hdi_0_copper = new G4Box( "HDI_0_copper", kHdi_0_x / 2, kHdi_y / 2, kHdi_copper_z / 2 );
  G4LogicalVolume* hdi_0_copper_lv = new G4LogicalVolume( hdi_0_copper, G4Material::GetMaterial( "G4_Cu" ), "HDI_0_Cu" );
  hdi_0_copper_lv->SetVisAttributes( cu_vis_att );

  // type-1
  // Kapton layer of HDI
  G4VSolid* hdi_1_kapton = new G4Box( "HDI_1_kapton", kHdi_1_x / 2, kHdi_y / 2, kHdi_kapton_z / 2 );
  G4LogicalVolume* hdi_1_kapton_lv = new G4LogicalVolume( hdi_1_kapton, G4Material::GetMaterial( "G4_KAPTON" ), "HDI_1_Kapton" );
  hdi_1_kapton_lv->SetVisAttributes( kapton_vis_att );

  // Copper layer of HDI
  G4VSolid* hdi_1_copper = new G4Box( "HDI_1_copper", kHdi_1_x / 2, kHdi_y / 2, kHdi_copper_z / 2 );
  G4LogicalVolume* hdi_1_copper_lv = new G4LogicalVolume( hdi_1_copper, G4Material::GetMaterial( "G4_Cu" ), "HDI_1_Cu" );
  hdi_1_copper_lv->SetVisAttributes( cu_vis_att );

  //////////////////////////////////////////////////////////
  // FPHX chips                                           //
  //////////////////////////////////////////////////////////
  G4VSolid* fphx = new G4Box( "FPHX", kFphx_x / 2, kFphx_y / 2, kFphx_z / 2 );
  G4LogicalVolume* fphx_lv = new G4LogicalVolume( fphx, G4Material::GetMaterial( "G4_Si" ), "FPHX" );
  fphx_lv->SetVisAttributes( fphx_vis_att );

  //     [][][][]                          [][][][]
  // [][][]                                      [][][]
  // []                                              []
  // []          MAKE  PART OF HALF-LADDER           []
  // []                                              []
  // [] [][] []     [][]     []     [][]     [] [][] []
  // [] [][] [][][]   [][] [][][] [][]   [][][] [][] []

  /////////////////////////////////////////////////////
  // Make volume for half ladder
  // type-0
  double half_ladder_0_thickness = kSi_area_0_z + kHdi_kapton_z + kHdi_copper_z;
  double half_ladder_1_thickness = half_ladder_0_thickness;

  // type-0
  G4VSolid* half_ladder_0 = new G4Box( "half_ladder_0", kHdi_0_x / 2, kHdi_y / 2,
				       half_ladder_0_thickness / 2 );
  G4LogicalVolume* half_ladder_0_lv = new G4LogicalVolume( half_ladder_0, defaultMaterial, "Half_Ladder_0" );
  half_ladder_0_lv->SetVisAttributes( wireframe_vis_att );

  // type-1
  G4VSolid* half_ladder_1 = new G4Box( "half_ladder_1", kHdi_1_x / 2, kHdi_y / 2,
				       (kSi_area_0_z + kHdi_kapton_z + kHdi_copper_z) / 2 );
  G4LogicalVolume* half_ladder_1_lv = new G4LogicalVolume( half_ladder_1, defaultMaterial, "Half_Ladder_1" );
  half_ladder_1_lv->SetVisAttributes( wireframe_vis_att );

  /////////////////////////////////////////////////////
  // Assemble parts 
  // ------------> z 
  // |<------->| half_ladder_0(1)_thickness
  // +----+----+    
  // |Si|Cu|Kap|    - Most downstream layer: Kapton of HDI,  +half_ladder_0(1)_thickness/2 - kapton thickness / 2
  // |Si|Cu|Kap|    - Central layer        : Copper of HDI,  - kapton thickness /2 - copper thickness / 2 from kapton center
  // |Si|Cu|Kap|    - Most upstream layer  : Silicon sensor, - copper thickness / 2 - silicon thickness / 2 from copper center  
  // |Si|Cu|Kap|    -                      : FPHX chips,     same as silicon center
  // |Si|Cu|Kap|    
  // |Si|Cu|Kap|    
  // |Si|Cu|Kap|    
  // +----+----+    
  // |    |    |    
  // |    |    +half_ladder_0(1)_thickness/2
  // |    center
  // -half_ladder_0(1)_thickness/2
  
  // type-0
  G4ThreeVector pos_hdi_kapton_half_ladder = G4ThreeVector( 0, 0,  half_ladder_0_thickness / 2 - kHdi_kapton_z / 2 );
  new G4PVPlacement( 0, pos_hdi_kapton_half_ladder, hdi_0_kapton_lv, "HDI_Kapton", half_ladder_0_lv, false, 0, fCheckOverlaps );
  
  G4ThreeVector pos_hdi_copper_half_ladder = pos_hdi_kapton_half_ladder + G4ThreeVector( 0, 0,  -kHdi_kapton_z / 2 - kHdi_copper_z / 2 );
  new G4PVPlacement( 0, pos_hdi_copper_half_ladder, hdi_0_copper_lv, "HDI_Copper", half_ladder_0_lv, false, 0, fCheckOverlaps );
  
  G4ThreeVector pos_si_area_0_up = pos_hdi_copper_half_ladder + G4ThreeVector( 0, kSi_area_0_y / 2,  - kHdi_copper_z / 2 - kSi_area_0_z / 2);
  new G4PVPlacement( 0, pos_si_area_0_up, si_area_0_lv, "Si_up", half_ladder_0_lv, false, 0, fCheckOverlaps );

  G4ThreeVector pos_si_area_0_down = pos_si_area_0_up + G4ThreeVector( 0, - kSi_area_0_y, 0 );
  new G4PVPlacement( 0, pos_si_area_0_down, si_area_0_lv, "Si_down", half_ladder_0_lv, false, 0, fCheckOverlaps );
  
  // type-1  
  G4ThreeVector pos_hdi_kapton_half_ladder_1 = G4ThreeVector( 0, 0,  half_ladder_1_thickness / 2 - kHdi_kapton_z / 2 );
  new G4PVPlacement( 0, pos_hdi_kapton_half_ladder_1, hdi_1_kapton_lv, "HDI_Kapton", half_ladder_1_lv, false, 0, fCheckOverlaps );
  
  G4ThreeVector pos_hdi_copper_half_ladder_1 = pos_hdi_kapton_half_ladder + G4ThreeVector( 0, 0,  -kHdi_kapton_z / 2 - kHdi_copper_z / 2 );
  new G4PVPlacement( 0, pos_hdi_copper_half_ladder_1, hdi_1_copper_lv, "HDI_Copper", half_ladder_1_lv, false, 0, fCheckOverlaps );
  
  G4ThreeVector pos_si_area_1_up = pos_hdi_copper_half_ladder + G4ThreeVector( 0, kSi_area_1_y / 2,  - kHdi_copper_z / 2 - kSi_area_1_z / 2);
  new G4PVPlacement( 0, pos_si_area_1_up, si_area_1_lv, "Si_up", half_ladder_1_lv, false, 0, fCheckOverlaps );

  G4ThreeVector pos_si_area_1_down = pos_si_area_1_up + G4ThreeVector( 0, - kSi_area_1_y, 0 );
  new G4PVPlacement( 0, pos_si_area_1_down, si_area_1_lv, "Si_down", half_ladder_1_lv, false, 0, fCheckOverlaps );


  //     [][][][]                          [][][][]
  // [][][]                                      [][][]
  // []                                              []
  // []               MAKE HALF-LADDER               []
  // []                                              []
  // [] [][] []     [][]     []     [][]     [] [][] []
  // [] [][] [][][]   [][] [][][] [][]   [][][] [][] []

  /////////////////////////////////////////////////////
  // Make volume for half ladder
  double half_ladder_width = kHdi_0_x + kHdi_1_x;
  double half_ladder_thickness = half_ladder_0_thickness;
  G4VSolid* half_ladder = new G4Box( "half_ladder", half_ladder_width / 2, kHdi_y / 2,
				     half_ladder_thickness / 2 );
  G4LogicalVolume* half_ladder_lv = new G4LogicalVolume( half_ladder_0, defaultMaterial, "HalfLadder" );
  half_ladder_lv->SetVisAttributes( wireframe_vis_att );

  G4ThreeVector pos_half_ladder_0_in_half_ladder = G4ThreeVector( half_ladder_width / 2 - kHdi_0_x / 2, 0, 0 );
  new G4PVPlacement( 0, pos_half_ladder_0_in_half_ladder, half_ladder_0_lv, "HalfLadder_0", half_ladder_lv, false, 0, fCheckOverlaps );

  G4ThreeVector pos_half_ladder_1_in_half_ladder = G4ThreeVector( -half_ladder_width / 2 + kHdi_1_x / 2, 0, 0 );
  new G4PVPlacement( 0, pos_half_ladder_1_in_half_ladder, half_ladder_1_lv, "HalfLadder_1", half_ladder_lv, false, 0, fCheckOverlaps );

  //     [][][][]                          [][][][]
  // [][][]                                      [][][]
  // []                                              []
  // []        PUT ELEMENTS IN MOTHER VOLUME         []
  // []                                              []
  // [] [][] []     [][]     []     [][]     [] [][] []
  // [] [][] [][][]   [][] [][][] [][]   [][][] [][] []

  //  G4VSolid* si_strip_0 = new G4Box( "Si_strip_0", kStrip_x / 2 , kStrip_y / 2, kStrip_z_0 / 2 );
  whatever = new G4PVPlacement( 0, G4ThreeVector(0, 0, 0), half_ladder_lv, "HalfLadder", worldLV, false, 0, fCheckOverlaps );
  
  /*G4LogicalVolume*
      // FPHX
      G4VSolid *fphx_box = new G4Box((boost::format("fphx_box_%d_%d") % inttlayer % itype).str(), fphx_x / 2., fphx_y / 2., fphx_z / 2.);
      G4LogicalVolume *fphx_volume = new G4LogicalVolume(fphx_box, G4Material::GetMaterial("G4_Si"),
                                                         (boost::format("fphx_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      if ((m_IsAbsorberActiveMap.find(inttlayer))->second > 0)
      {
        m_PassiveVolumeTuple.insert(make_pair(fphx_volume, make_tuple(inttlayer, PHG4InttDefs::FPHX)));
      }
      m_DisplayAction->AddVolume(fphx_volume, "FPHX");

      const double gap_sensor_fphx = params->get_double_param("gap_sensor_fphx") * cm;

      //  FPHX Container
      // make a container for the FPHX chips needed for this sensor, and  then place them in the container
      G4VSolid *fphxcontainer_box = new G4Box((boost::format("fphxcontainer_box_%d_%d") % inttlayer % itype).str(),
                                              fphx_x / 2., fphx_y / 2., hdi_z / 2.);
      G4LogicalVolume *fphxcontainer_volume = new G4LogicalVolume(fphxcontainer_box, G4Material::GetMaterial("G4_AIR"),
                                                                  (boost::format("fphxcontainer_volume_%d_%d") % inttlayer % itype).str(), 0, 0, 0);
      m_DisplayAction->AddVolume(fphxcontainer_volume, "FPHXContainer");

      // Install multiple FPHX volumes in the FPHX container volume
      // one FPHX chip per cell - each cell is 128 channels
      const double fphx_offsetx = 0.;
      const double fphx_offsety = 0.;      
      int ncopy;
      double offsetz, cell_length_z;

      if (laddertype == PHG4InttDefs::SEGMENTATION_Z)  // vertical strips
      {
        // For laddertype 0, we have 5 cells per sensor, but the strips are vertical, so we have to treat it specially
        ncopy = nstrips_z_sensor / 128.0;
      }
      else if (laddertype == PHG4InttDefs::SEGMENTATION_PHI)
      {
        ncopy = nstrips_z_sensor;
      }
      else
      {
        cout << PHWHERE << "invalid laddertype " << laddertype << endl;
        gSystem->Exit(1);
        // this is just to make the optimizer happy which otherwise complains about possibly
        // uninitialized variables. It doesn't know gSystem->Exit(1) quits,
        // this exit here terminates the program for it
        exit(1);
      }
      cell_length_z = strip_z * nstrips_z_sensor / ncopy;
      offsetz = (ncopy % 2 == 0) ? -2. * cell_length_z / 2. * double(ncopy / 2) + cell_length_z / 2. + fphx_offset_z : -2. * cell_length_z / 2. * double(ncopy / 2) + fphx_offset_z;

      G4VPVParameterisation *fphxparam = new PHG4InttFPHXParameterisation(fphx_offsetx, +fphx_offsety, offsetz, 2. * cell_length_z / 2., ncopy);
      new G4PVParameterised((boost::format("fphxcontainer_%d_%d") % inttlayer % itype).str(),
                            fphx_volume, fphxcontainer_volume, kZAxis, ncopy, fphxparam, OverlapCheck());
  */  
  
  // Always return the physical World
  return worldPV;
}

//void B4cDetectorConstruction::SetMagField(G4double fieldValue)
void B4cDetectorConstruction::ConstructSDandField()
{
  // Silicon type-0
  B4cCalorimeterSD* si_strip_0_sd = new B4cCalorimeterSD("SiStrip0SD", "SiStrip0HitsCollection", 256 );
  G4SDManager::GetSDMpointer()->AddNewDetector( si_strip_0_sd );
  SetSensitiveDetector("Si_strip_0", si_strip_0_sd );

  // Silicon type-1
  B4cCalorimeterSD* si_strip_1_sd = new B4cCalorimeterSD("SiStrip1SD", "SiStrip1HitsCollection", 128 );
  G4SDManager::GetSDMpointer()->AddNewDetector( si_strip_1_sd );
  SetSensitiveDetector("Si_strip_1", si_strip_1_sd );

  /* minimum 
     B4cCalorimeterSD* absoSD[num_of_towers_hcal][num_of_towers_hcal];
     absoSD[i][j] = new B4cCalorimeterSD( absorber_SD.str().c_str(), absorber_HC.str().c_str(), nofLayers2 );
     G4SDManager::GetSDMpointer()->AddNewDetector( absoSD[i][j] );
     SetSensitiveDetector( absorber_LV.str().c_str(), absoSD[i][j] );     
  */

  /*
 280   │   //
 281   │   auto absoSD
 282   │     = new B4cCalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
 283   │   G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
 284   │   SetSensitiveDetector("AbsoLV",absoSD);
*/
  
  
  /*
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
  */

}



const G4VPhysicalVolume* B4cDetectorConstruction::GetPV() const
{ 
  return whatever;
}
