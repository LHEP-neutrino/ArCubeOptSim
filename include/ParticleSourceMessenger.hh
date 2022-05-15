#ifndef __PARTICLE_SOURCE_OPT_PH_MESSENGER_HH__
#define __PARTICLE_SOURCE_OPT_PH_MESSENGER_HH__

#include "G4UImessenger.hh"
#include "globals.hh"


class ParticleSourceOptPh;
class G4ParticleTable;
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

class ParticleSourceOptPhMessenger: public G4UImessenger
{
public:
	ParticleSourceOptPhMessenger(ParticleSourceOptPh *pParticleSource);
	virtual ~ParticleSourceOptPhMessenger();
  
	void SetNewValue(G4UIcommand *pCommand, G4String hNewValues);
 
private:
	ParticleSourceOptPh *fParticleSource;
	G4ParticleTable *fParticleTable;
    
     G4UIdirectory              *fDirectory;
	 
	 G4UIcmdWithAnInteger       *fPrimNbCmd;
     G4UIcmdWithAString         *fTypeCmd;
     G4UIcmdWithAString         *fShapeCmd;
     G4UIcmdWith3VectorAndUnit  *fCenterCmd;
     G4UIcmdWithADoubleAndUnit  *fHalfxCmd;
     G4UIcmdWithADoubleAndUnit  *fHalfyCmd;
     G4UIcmdWithADoubleAndUnit  *fHalfzCmd;
     G4UIcmdWithADoubleAndUnit  *fRadiusCmd;
     G4UIcmdWithAString         *fConfineCmd;
     G4UIcmdWithAString         *fAngTypeCmd;
     G4UIcmdWithAnInteger       *fVerbosityCmd;
     G4UIcmdWithAString         *fParticleCmd;
     G4UIcmdWith3VectorAndUnit  *fPositionCmd;
     G4UIcmdWith3Vector         *fDirectionCmd;
     G4UIcmdWith3Vector         *fPolarCmd;
     G4UIcmdWithADoubleAndUnit  *fEnergyCmd;
     G4UIcmdWithoutParameter    *fGetDirectCmd;
	 G4UIcmdWithoutParameter    *fGetPolarCmd;
	 G4UIcmdWithoutParameter    *fGetPartCmd;
	 G4UIcmdWithAString         *fEneryDistCmd;
	 G4UIcmdWithAString         *fSpectrumFileCmd;
};

#endif
