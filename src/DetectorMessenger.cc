#include "DetectorMessenger.hh"
#include "DetConstr.hh"
#include "OptPropManager.hh"

#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>
#include <G4ParticleTable.hh>
#include <G4UIdirectory.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWith3Vector.hh>
#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWithABool.hh>
#include <G4Tokenizer.hh>
#include <G4ios.hh>

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>


DetectorMessenger::DetectorMessenger(DetConstrOptPh *pDetector)
:fDetector(pDetector)
{
	std::stringstream txt;
	
	fDetectorDir = new G4UIdirectory("/argoncube/detector/");
	fDetectorDir->SetGuidance("ArgonCube detector geometry and material properties control.");
	
	fDetectorOptDir = new G4UIdirectory("/argoncube/detector/optical/");
	fDetectorOptDir->SetGuidance("ArgonCube optical properties setup controls.");
	
	fTpbThicknCmd = new G4UIcmdWithADoubleAndUnit("/argoncube/detector/setTpbThick",this);
	fTpbThicknCmd->SetGuidance("Set the thickness of the arcLight TPB layer.");
	fTpbThicknCmd->SetParameterName("TpbThick",false);
	fTpbThicknCmd->AvailableForStates(G4State_PreInit);
	fTpbThicknCmd->SetDefaultUnit("mm");
	fTpbThicknCmd->SetUnitCandidates("mm cm");
	
	fDetConstrVerb = new G4UIcmdWithAnInteger("/argoncube/detector/verbosity", this);
	fDetConstrVerb->SetGuidance("Set the verbosity for the detector constructor.");
	fDetConstrVerb->SetParameterName("DetVerb",false);
	fDetConstrVerb->AvailableForStates(G4State_PreInit, G4State_Idle);
	txt.str(""); txt << "DetVerb>=0 && DetVerb<=" << static_cast<G4int>(DetVerbosity::kDebug);
	fDetConstrVerb->SetRange( txt.str().c_str() );
	
	fStepTrackLimCmd = new G4UIcommand("/argoncube/detector/setMaxTrackStep", this);
	fStepTrackLimCmd->SetGuidance("Set the maximum tracking step length for a given logical volume.");
	fStepTrackLimCmd->SetGuidance("usage: /argoncube/detector/setMaxTrackStep LogVol step unit");
	fStepTrackLimCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	G4UIparameter *param;
	param = new G4UIparameter("LogVol", 's', false);
	fStepTrackLimCmd->SetParameter(param);
	param = new G4UIparameter("Step", 'd', false);
	fStepTrackLimCmd->SetParameter(param);
	param = new G4UIparameter("unit", 's', false);
	fStepTrackLimCmd->SetParameter(param);
	
	fPhysVolCoordCmd = new G4UIcmdWithAString("/argoncube/detector/PhysVolCoord", this);
	fPhysVolCoordCmd->SetParameterName("physvol", false);
	fPhysVolCoordCmd->AvailableForStates(G4State_Idle);
	
	fPhysVolList = new G4UIcmdWithoutParameter("/argoncube/detector/PhysVolList", this);
	fPhysVolList->AvailableForStates(G4State_Idle);
	
	fPhysVolInfoCmd = new G4UIcmdWithAString("/argoncube/detector/PhysVolInfo", this);
	//fPhysVolInfoCmd->SetParameterName("physvol", false);
	fPhysVolInfoCmd->AvailableForStates(G4State_Idle);
	
	fLoadOpticalSettingsFile = new G4UIcmdWithAString("/argoncube/detector/optical/loadOptSett", this);
	fLoadOpticalSettingsFile->SetGuidance("Load the json file with all the optical settings used");
	fLoadOpticalSettingsFile->SetParameterName("SettFile",false);
	fLoadOpticalSettingsFile->AvailableForStates(G4State_Idle);
	
	fOpticalSettingsVerb = new G4UIcmdWithAnInteger("/argoncube/detector/optical/verbosity", this);
	fOpticalSettingsVerb->SetGuidance("Set the verbosity for the manager of the optical properties settings.");
	fOpticalSettingsVerb->SetParameterName("OptVerb",false);
	fOpticalSettingsVerb->AvailableForStates(G4State_PreInit, G4State_Idle);
	txt.str(""); txt << "OptVerb>=0 && OptVerb<=" << OptPropManager::kDebug;
	fOpticalSettingsVerb->SetRange( txt.str().c_str() );
}

DetectorMessenger::~DetectorMessenger()
{
	//delete fLArAbsorbtionLengthCmd;
	//delete fLArRayScatterLengthCmd;
	
	delete fTpbThicknCmd;
	delete fDetConstrVerb;
	delete fOpticalSettingsVerb;
	delete fPhysVolCoordCmd;
	delete fPhysVolInfoCmd;
	delete fPhysVolList;
	delete fLoadOpticalSettingsFile;
	delete fDetectorOptDir;
	delete fDetectorDir;
	delete fStepTrackLimCmd;
	
}

void DetectorMessenger::SetNewValue(G4UIcommand *pUIcommand, G4String hNewValue)
{
	if(pUIcommand == fDetConstrVerb){
		G4cout << "Info --> DetectorMessenger::SetNewValue(...): called command fDetConstrVerb" << G4endl;
		fDetector->SetVerbosity( static_cast<DetVerbosity>(std::stoi(hNewValue)) );
	}
	
	if(pUIcommand == fTpbThicknCmd){
		G4cout << "Info --> DetectorMessenger::SetNewValue(...): called command fTpbThicknCmd" << G4endl;
		fDetector->SetTpbThickness( fTpbThicknCmd->GetNewDoubleValue(hNewValue) );
	}
	
	if(pUIcommand == fStepTrackLimCmd){
		G4cout << "\nInfo --> MediPixDetMessenger::SetNewValue: called command \"fStepTrackLimCmd\"" << G4endl;
		
		ProcessTrackLimCmd(fStepTrackLimCmd, hNewValue);
		return;
	}
	
	if(pUIcommand == fPhysVolCoordCmd){
		G4cout << "Info --> DetectorMessenger::SetNewValue(...): called command fPhysVolCoordCmd" << G4endl;
		fDetector->PrintVolumeCoordinates( hNewValue );
	}
	
	if(pUIcommand == fPhysVolInfoCmd){
		G4cout << "Info --> DetectorMessenger::SetNewValue(...): called command fPhysVolInfoCmd" << G4endl;
		fDetector->PrintVolumeInfo( hNewValue );
	}
	
	if(pUIcommand == fPhysVolList){
		G4cout << "Info --> DetectorMessenger::SetNewValue(...): called command fPhysVolList" << G4endl;
		fDetector->PrintListOfPhysVols();
	}
	
	if(pUIcommand == fLoadOpticalSettingsFile){
		G4cout << "Info --> DetectorMessenger::SetNewValue(...): called command fLoadOpticalSettingsFile" << G4endl;
		OptPropManager::GetInstance()->ProcessJsonFile( hNewValue );
	}
	
	if(pUIcommand == fOpticalSettingsVerb){
		G4cout << "Info --> DetectorMessenger::SetNewValue(...): called command fOpticalSettingsVerb" << G4endl;
		OptPropManager::GetInstance()->SetVerbosity( (OptPropManager::verbosity)std::stoi(hNewValue) );
	}
}


void DetectorMessenger::ProcessTrackLimCmd(const G4UIcommand *cmd, const G4String& newValues)
{
	G4Tokenizer next(newValues);

	// check argument
	G4String LogVol = next();
	G4double stepLim = StoD(next());
	G4String unit = next();
	
	G4cout << "Info --> DetectorMessenger::ProcessTrackLimCmd: Setting track length limit for <" << LogVol << "> logical volume to " << stepLim << " " << unit << "." << G4endl;
	
	fDetector->SetStepTrackLimit(LogVol, stepLim*cmd->ValueOf(unit.c_str()));
}

