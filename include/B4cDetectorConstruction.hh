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
/// \file B4cDetectorConstruction.hh
/// \brief Definition of the B4cDetectorConstruction class

#ifndef B4cDetectorConstruction_h
#define B4cDetectorConstruction_h 1

#include <cstdio>
#include <sstream>
#include <iomanip> // for setprecision

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4UniformMagField.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4GenericMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4SubtractionSolid.hh"
#include "G4RunManager.hh"

#include "B4cCalorimeterSD.hh"
//#include "B4RunManager.hh"

#include "globals.hh"
#include "Constants.hh"

//class G4GlobalMagFieldMessenger;
class G4VPhysicalVolume;

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In ConstructSDandField() sensitive detectors of B4cCalorimeterSD type
/// are created and associated with the Absorber and Gap volumes.
/// In addition a transverse uniform magnetic field is defined
/// via G4GlobalMagFieldMessenger class.


class B4cDetectorConstruction : public G4VUserDetectorConstruction
{
private:
  // methods
  void DefineCommands();
  void DefineMaterials();
  void DefineVisAttributes();
  G4VPhysicalVolume* DefineVolumes();
  void DefineVolumes_Ecal( G4LogicalVolume* worldLV );
  void DefineVolumes_Hcal( G4LogicalVolume* worldLV );
  void DefineVolumes_Wall( G4LogicalVolume* worldLV );
  void DefineVolumes_Surroundings( G4LogicalVolume* worldLV );
  void DefineVolumes_NotUsed( G4LogicalVolume* worldLV );

  //void SetWallThickness( G4double thickness );
  // data members
  //  static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
  // magnetic field messenger
  //
  G4Material* defaultMaterial;
  G4Material* absorberMaterial;
  G4Material* absorberMaterial2;
  G4Material* gapMaterial;
  G4Material* airMaterial;
  G4Material* stopperMaterial;

  G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
  //  static G4GenericMessenger* fMessenger;
  G4GenericMessenger* fMessenger;
  /*
  G4int num_of_towers_hcal	= 6;
  G4int num_of_towers	= 12;

  G4double calorSizeXY		= 600.0 * mm;

  G4int   nofLayers_ecal = 66;     // number of layers
  G4int   nofLayers_hcal = 36;     // number of layers

  G4int    nofLayers		= 66;//66
  G4int    nofLayers2		= 36;

  G4double hadcell		= 10.0 * cm;
  G4double emcell		= 5.0 * cm;

  G4double absoThickness_ecal	= 1.5 * mm;
  G4double gapThickness_ecal		= 4.0 * mm;

  G4double absoThickness_hcal	= 20. * mm;
  G4double gapThickness_hcal	= 6. * mm;//used to be 3.0 mm
  */
  
  G4double layerThickness_ecal;
  G4double layerThickness_hcal;
  G4double calorThickness_ecal;
  G4double calorThickness_hcal;

  G4double worldSizeXY;
  G4double worldSizeZ;

  G4double wallThickness;
  G4double surroundingsThickness;
  G4double lidThickness;
  
  G4VisAttributes* simpleBoxVisAtt;
  G4VisAttributes* ecalAbsorberCol;
  G4VisAttributes* hcalAbsorberCol;
  G4VisAttributes* gapCol;
  G4VisAttributes* wallCol;
  G4VisAttributes* surroundingsCol;


  G4VisAttributes* wireframe_vis_att;
  G4VisAttributes* cu_vis_att;
  G4VisAttributes* kapton_vis_att;
  G4VisAttributes* si_vis_att;
  G4VisAttributes* fphx_vis_att;
  
  G4bool isSurroundings = true;

  G4VPhysicalVolume* whatever;

  
public:
  B4cDetectorConstruction();
  virtual ~B4cDetectorConstruction();
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
  //  void SetMagField( G4double fieldValue );

  // number of towers
  G4int GetTowerNumEcal(){ return num_of_towers;};
  G4int GetTowerNumHcal(){ return num_of_towers_hcal;};

  // number of layers of a tower
  G4int GetLayerNumEcal(){ return nofLayers;};
  G4int GetLayerNumHcal(){ return nofLayers;}; // which is Hcal?

  // size of cells in x-y dimensions
  G4double GetCellSizeEcal(){ return emcell;};
  G4double GetCellSizeHcal(){ return hadcell;};

  // thickness of absorbers
  G4double GetThicknessAbsorberEcal(){ return absoThickness_ecal;};
  G4double GetThicknessAbsorberHcal(){ return absoThickness_hcal;};

  // thickness of gaps
  G4double GetThicknessGapEcal(){ return gapThickness_ecal;};
  G4double GetThicknessGapHcal(){ return gapThickness_hcal;};


  const G4VPhysicalVolume* GetPV() const;
  
};


#endif
