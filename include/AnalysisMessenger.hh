#ifndef __ANALISYS_MESSENGER_OPT_PH_HH__
#define __ANALISYS_MESSENGER_OPT_PH_HH__


#include "G4UIcommand.hh"
#include "G4ThreeVector.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4Tokenizer.hh"

#include "G4ios.hh"

#include "G4UImessenger.hh"
#include "globals.hh"




class AnalysisManagerOptPh;

class AnalysisOptPhMessenger: public G4UImessenger
{
public:
	AnalysisOptPhMessenger(AnalysisManagerOptPh *pAnManager);
	virtual ~AnalysisOptPhMessenger();

	void SetNewValue(G4UIcommand *pUIcommand, G4String hNewValue);

private:
	AnalysisManagerOptPh* fAnManager;

	G4UIdirectory *fAnalysisDir;
	
	G4UIcmdWithAnInteger *fRandSeedCmd;
	G4UIcmdWithAnInteger *fVerboseCmd;
	G4UIcmdWithAnInteger *fPrintModuloCmd;
	G4UIcmdWithABool *fStepDebugCmd;
	G4UIcmdWithAnInteger *fSaveDataCmd;
	G4UIcmdWithAString *fFileNameCmd;
	G4UIcmdWithAString *fDefOptSDCmd;
	
	G4UIcmdWithAnInteger *fAutoFlushCmd;
	G4UIcmdWithAnInteger *fAutoSaveCmd;
	
	
	
};

#endif

