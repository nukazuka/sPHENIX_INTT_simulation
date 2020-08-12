#pragma once

#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

// the actual number of towers should be num_of_towers * num_of_towers...
const G4int num_of_towers_hcal	= 4;
const G4int num_of_towers	= 4;

const G4int num_of_towers_ecal	= num_of_towers;

const G4double calorSizeXY	= 600.0 * mm;

const G4int   nofLayers_ecal	= 66;     // number of layers
const G4int   nofLayers_hcal	= 36;     // number of layers

const G4int    nofLayers	= 66;//66
const G4int    nofLayers2	= 36;

const G4double hadcell		= 10.0 * cm;
const G4double emcell		= 5.0 * cm;

const G4double absoThickness_ecal	= 1.5 * mm;
const G4double gapThickness_ecal	= 4.0 * mm;

const G4double absoThickness_hcal	= 20. * mm;
const G4double gapThickness_hcal	= 6. * mm;//used to be 3.0 mm
