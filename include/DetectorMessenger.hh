#ifndef __DET_MESSANGER_HH__
#define __DET_MESSANGER_HH__

#include <G4UImessenger.hh>
#include <globals.hh>


class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;

class DetConstrOptPh;

class DetectorMessenger: public G4UImessenger
{
public:
	DetectorMessenger(DetConstrOptPh *pDetector);
	~DetectorMessenger();

	void SetNewValue(G4UIcommand *pUIcommand, G4String hString);

private:
	DetConstrOptPh* fDetector;
	
	G4UIdirectory *fDetectorDir, *fDetectorOptDir;
	
	//G4UIcmdWithADoubleAndUnit *fLArAbsorbtionLengthCmd;
	//G4UIcmdWithADoubleAndUnit *fLArRayScatterLengthCmd;
	
	G4UIcmdWithAnInteger *fDetConstrVerb;
	
	G4UIcmdWithADoubleAndUnit *fTpbThicknCmd;
	
	G4UIcmdWithAString *fPhysVolCoordCmd;
	
	G4UIcmdWithAString *fPhysVolInfoCmd;
	
	G4UIcmdWithoutParameter *fPhysVolList;
	
	G4UIcmdWithAString *fLoadOpticalSettingsFile;
	
	G4UIcmdWithAnInteger *fOpticalSettingsVerb;
	
};

#endif

