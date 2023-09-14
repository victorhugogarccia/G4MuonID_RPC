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
/// \file MuonID_sim.cc
/// \brief Main program of the geometry of the MuonID detector for ALICE3 with Resistive Rlate Chambers (RPC)
/// \author: Héctor David Regules, rhectord@cern.ch
///          Víctor Hugo García

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
//#include "G4VisAttributes.hh"
//#include "G4Colour.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "Physics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int main(int argc,char** argv) {

   G4RunManager *runManager = new G4RunManager();
   
   runManager->SetUserInitialization(new DetectorConstruction());
   runManager->SetUserInitialization(new PhysicsList());
   runManager->SetUserInitialization(new ActionInitialization());
   
   runManager->Initialize();
	
   G4UIExecutive *ui = new G4UIExecutive(argc, argv);
	
   G4VisManager *visManager = new G4VisExecutive();
   visManager->Initialize();
   
   G4UImanager *UImanager = G4UImanager::GetUIpointer();
   
   /*UImanager->ApplyCommand("/vis/open OGL");
   UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
   UImanager->ApplyCommand("/vis/viewer/set/autoRefresh true");*/
   
   UImanager->ApplyCommand("/vis/open OGL");
   UImanager->ApplyCommand("/vis/viewer/set/viewpointVector 1 1 1");
   UImanager->ApplyCommand("/vis/drawVolume");
   UImanager->ApplyCommand("/vis/viewer/update");
   UImanager->ApplyCommand("/vis/viewer/set/autoRefresh true");
   UImanager->ApplyCommand("/vis/scene/add/trajectories smooth");
   UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
   
   ui->SessionStart();
	
   return 0;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

