#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4MuonMinus.hh" 
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
//#include "G4VisAttributes.hh"
//#include "G4Tubs.hh"
#include "Randomize.hh"

#include <cmath>
#include <math.h>

class DetectorConstruction; //¿?

/// Action initialization class

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

  public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);
    
  private:
    G4ParticleGun *fParticleGun;
    
};

#endif
