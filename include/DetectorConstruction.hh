#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include <iostream>
#include <string>

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "SensitiveDetector.hh"

class DetectorConstruction : public G4VUserDetectorConstruction {

   public:
      DetectorConstruction();
      ~DetectorConstruction();
	
   virtual G4VPhysicalVolume *Construct();

   private: 
	G4LogicalVolume *logicStripR1;
	virtual void ConstructSDandField();
};

#endif
