#include "SensitiveDetector.hh"

SensitiveDetector::SensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{}

SensitiveDetector::~SensitiveDetector()
{}

G4bool SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{
	G4Track *track = aStep->GetTrack();

	//track->SetTrackStatus(fStopAndKill);
	
	G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
	G4StepPoint *postStepPoint = aStep->GetPostStepPoint();//duda sobre particulas cargadas pues menciona que varia de acuerdo a la particula si tiene o no carga 
	
	G4ThreeVector posPhoton = preStepPoint->GetPosition();
	
	//G4cout << "Muon position: " << posPhoton << G4endl; //Este comando me dice la posicion en la que el foton interactuo con el detector
	
	const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
	
	G4int copyNo = touchable->GetCopyNumber();
	
	G4cout << "Copy number: " << copyNo << G4endl; //Este sirve para saber el numero del volumen logico del detector en el que se encuentra el  foton
	
	G4VPhysicalVolume *physVol = touchable->GetVolume();
	G4ThreeVector posDetector = physVol->GetTranslation(); //Este ultimo sirve para saber la posicion del detector que fue atravesado por un foton
	
	G4cout << "Detector position: " << posDetector << G4endl;
	
	//G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	
	//G4AnalysisManager *man = G4AnalysisManager::Instance();
	//man->FillNtupleIColumn(0,evt);
	//man->FillNtupleDColumn(1, posDetector[0]);
	//man->FillNtupleDColumn(2, posDetector[1]);
	//man->FillNtupleDColumn(3, posDetector[2]);
	//man->AddNtupleRow(0);
	
}
