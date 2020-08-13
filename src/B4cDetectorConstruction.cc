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
  nistManager->FindOrBuildMaterial( "G4_Cu", fromIsotopes );
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

  // HDI: one for half ladder
  //const double kHdi_y = kStrip_y * kNstrips_phi_sensor; // 1.9968 cm
  const double kHdi_x = kStrip_x_0 * kNstrips_x_sensor_0 + 2 * kSensor_edge_x
    + kStrip_x_1 * kNstrips_x_sensor_1 + 2 * kSensor_edge_x; // 23.2 cm
  
  std::cout << std::string( 100, '=' ) << std::endl;
  std::cout << "kHdi_x = " << kHdi_x << std::endl;
  std::cout << "kHdi_y = " << kHdi_y << std::endl;
  std::cout << std::string( 100, '=' ) << std::endl;

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

  double pos_z = 0.0;
  //////////////////////////////////////////////////////////
  // HDI                                                  //
  //////////////////////////////////////////////////////////

  // Kapton layer of HDI
  G4VSolid* hdi_kapton = new G4Box( "HDI_kapton", kHdi_x / 2, kHdi_y / 2, kHdi_kapton_z / 2 );
  G4LogicalVolume* hdi_kapton_lv = new G4LogicalVolume( hdi_kapton, G4Material::GetMaterial( "G4_KAPTON" ), "HDI_Kapton" );
  hdi_kapton_lv->SetVisAttributes( kapton_vis_att );

  new G4PVPlacement( 0, G4ThreeVector(0, 0, pos_z ), hdi_kapton_lv, "hdi_kapton", worldLV, false, 0, fCheckOverlaps );

  // Copper layer of HDI
  G4VSolid* hdi_copper = new G4Box( "HDI_copper", kHdi_x / 2, kHdi_y / 2, kHdi_copper_z / 2 );
  G4LogicalVolume* hdi_copper_lv = new G4LogicalVolume( hdi_copper, G4Material::GetMaterial( "G4_Cu" ), "HDI_Cu" );
  hdi_copper_lv->SetVisAttributes( cu_vis_att );

  pos_z -= kHdi_kapton_z / 2 + kHdi_copper_z / 2;
  G4ThreeVector pos_hdi_copper = G4ThreeVector( 0, 0, pos_z );
  new G4PVPlacement( 0, pos_hdi_copper, hdi_copper_lv, "hdi_copper", worldLV, false, 0, fCheckOverlaps );

  //////////////////////////////////////////////////////////
  // FPHX chips                                           //
  //////////////////////////////////////////////////////////
  G4VSolid* fphx = new G4Box( "FPHX", kFphx_x / 2, kFphx_y / 2, kFphx_z / 2 );
  //  G4LogicalVolume* fphx_lv[];


  //////////////////////////////////////////////////////////
  // Silicon sensors                                      //
  //////////////////////////////////////////////////////////
  // Si strip
  // +-------------+
  // +-------------+
  
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
  const double kSi_area_0_z = kSi_cell_0_z;
  
  G4VSolid* si_area_0  = new G4Box( "Si_area_0", kSi_area_0_x / 2, kSi_area_0_y / 2, kSi_area_0_z / 2 );				    
  G4LogicalVolume* si_area_0_lv = new G4LogicalVolume( si_area_0, defaultMaterial, "Si_area_0" );
  si_area_0_lv->SetVisAttributes( wireframe_vis_att );

  G4VSolid* si_area_1  = new G4Box( "Si_area_1", kSi_area_1_x / 2, kSi_area_0_y / 2, kSi_area_0_z / 2 );				    
  G4LogicalVolume* si_area_1_lv = new G4LogicalVolume( si_area_1, defaultMaterial, "Si_area_1" );
  si_area_1_lv->SetVisAttributes( wireframe_vis_att );


  // Placement
  new G4PVReplica( "name", si_cell_0_lv, si_area_0_lv, kXAxis, kNstrips_x_sensor_0, kStrip_x_0 );
  new G4PVReplica( "name", si_cell_1_lv, si_area_1_lv, kXAxis, kNstrips_x_sensor_1, kStrip_x_1 );

  //  double pos_x = kSi_area_0_x / 2;
  double hdi_margin_x = (kHdi_x - kSi_area_0_x - kSi_area_1_x ) / 2;
  double pos_x = kSi_area_0_x / 2 -( kHdi_x / 2 - kSi_area_1_x - hdi_margin_x );
  pos_z -= kHdi_copper_z / 2 + kStrip_z / 2;
  G4ThreeVector pos_si_0_up = G4ThreeVector( pos_x, kSi_area_0_y / 2,  pos_z );
  new G4PVPlacement( 0, pos_si_0_up, si_area_0_lv, "si_area_0", worldLV, false, 0, fCheckOverlaps );

  G4ThreeVector pos_si_0_down = G4ThreeVector(pos_x, -kSi_area_0_y / 2,  pos_z );
  new G4PVPlacement( 0, pos_si_0_down, si_area_0_lv, "si_area_0", worldLV, false, 0, fCheckOverlaps );

  // for shorter type
  //pos_x = -kSi_area_1_x / 2;
  //pos_x = -1 * (pos_x - kSi_area_0_x / 2 + kSi_area_1_x / 2 );
  pos_x = -1 * (kSi_area_1_x / 2 - ( kHdi_x / 2 - kSi_area_1_x - hdi_margin_x ) );
  G4ThreeVector pos_si_1_up = G4ThreeVector( pos_x, kSi_area_0_y / 2,  pos_z );
  new G4PVPlacement( 0, pos_si_1_up, si_area_1_lv, "si_area_1", worldLV, false, 0, fCheckOverlaps );

  G4ThreeVector pos_si_1_down = G4ThreeVector( pos_x, -kSi_area_0_y / 2,  pos_z );
  new G4PVPlacement( 0, pos_si_1_down, si_area_1_lv, "si_area_1", worldLV, false, 0, fCheckOverlaps );

  //  G4VSolid* si_strip_0 = new G4Box( "Si_strip_0", kStrip_x / 2 , kStrip_y / 2, kStrip_z_0 / 2 );

  
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



  
  /*
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
  */
  
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

//void B4cDetectorConstruction::SetMagField(G4double fieldValue)
void B4cDetectorConstruction::ConstructSDandField()
{

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
