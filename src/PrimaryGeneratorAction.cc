#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(){
   fParticleGun = new G4ParticleGun(1);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction(){
   delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {

   G4ParticleDefinition* particle = G4MuonMinus::Definition();
    	
   G4ThreeVector pos(6.0025*m, 0., 0.);
   G4double r = G4UniformRand() * 100;
   G4double theta = G4UniformRand() * CLHEP::pi;
   G4double phi = G4UniformRand() * 2.0 * CLHEP::pi;
   G4ThreeVector mom(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
	
   fParticleGun->SetParticlePosition(pos);
   fParticleGun->SetParticleMomentumDirection(mom);
   fParticleGun->SetParticleMomentum(r*GeV);
   fParticleGun->SetParticleDefinition(particle);
   fParticleGun->GeneratePrimaryVertex(anEvent);
   
}
