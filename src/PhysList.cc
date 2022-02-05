/**
 * AUTHOR: 
 * CONTACT: 
 * FIRST SUBMISSION: 12-28-2010
 *
 * REVISION:
 *
 * mm-dd-yyyy, What is changed, Whoami
 *
 * 12-28-2010, the switching from Geant4.9.2 to Geant4.9.3
 *             is finished starting from e+, Xiang
 */

// ---------------------------------------------------------------------------

#include "PhysListMessenger.hh"
#include "PhysList.hh"

#include <iomanip>

#include "globals.hh"

#include "G4PhysicsListHelper.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4ParticleWithCuts.hh"
#include "G4UserLimits.hh"
#include "G4ios.hh"

#include "G4StepLimiter.hh"
#include "G4OpAbsorption.hh"
#include "G4OpWLS.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4SystemOfUnits.hh"

//Instance and initialize all the static members
//G4ThreadLocal G4int PhysListOptPh::fVerboseLevel = 0;
//G4ThreadLocal G4int PhysListOptPh::fOpVerbLevel = 0;
G4ThreadLocal G4OpAbsorption* PhysListOptPh::fAbsorptionProcess = nullptr;
G4ThreadLocal G4OpWLS* PhysListOptPh::fWLSProcess = nullptr;
G4ThreadLocal G4OpRayleigh* PhysListOptPh::fRayleighScatteringProcess = nullptr;
//G4ThreadLocal G4OpMieHG* PhysListOptPh::fMieHGScatteringProcess = nullptr;
G4ThreadLocal G4OpBoundaryProcess* PhysListOptPh::fBoundaryProcess = nullptr;
G4ThreadLocal G4StepLimiter *PhysListOptPh::fStepLimiter = nullptr;


PhysListOptPh::PhysListOptPh():
G4VUserPhysicsList(),
fVerboseLevel(PhysVerbosity::kSilent),
fDefaultCutValue(1. * mm)
{
	
	fMessenger = new PhysListOptPhMessenger(this);
	
	/*
	cutForGamma = defaultCutValue;
	cutForElectron = defaultCutValue;
	cutForPositron = defaultCutValue;
	VerboseLevel = 0;
	OpVerbLevel = 0;
	*/
	SetVerboseLevel( static_cast<G4int>(PhysListOptPh::fVerboseLevel) );

}


PhysListOptPh::~PhysListOptPh()
{
	if(fMessenger) delete fMessenger;
}


void PhysListOptPh::ConstructParticle()
{
	G4Geantino::GeantinoDefinition();
	//G4ChargedGeantino::ChargedGeantinoDefinition();
	
	G4OpticalPhoton::OpticalPhotonDefinition();
}


void PhysListOptPh::ConstructProcess()
{
	AddTransportation();
	ConstructOptical();
}





//#include "XeMaxTimeCuts.hh"
//#include "XeMinEkineCuts.hh"





// Optical Processes ////////////////////////////////////////////////////////

void PhysListOptPh::ConstructOptical()
{
	//I don't need all the scintillation and chernkov processses
	
	// optical processes
	fStepLimiter = new G4StepLimiter();
	fAbsorptionProcess = new G4OpAbsorption();
	fWLSProcess = new G4OpWLS();
	fRayleighScatteringProcess = new G4OpRayleigh();
	//fMieHGScatteringProcess = new G4OpMieHG();
	fBoundaryProcess = new G4OpBoundaryProcess();
	
	//  theAbsorptionProcess->DumpPhysicsTable();
	//  theWLSProcess->DumpPhysicsTable();
	//  theRayleighScatteringProcess->DumpPhysicsTable();
	
	
	
	fStepLimiter->SetVerboseLevel(static_cast<G4int>(PhysListOptPh::fVerboseLevel));
	fAbsorptionProcess->SetVerboseLevel(static_cast<G4int>(PhysListOptPh::fVerboseLevel));
	fWLSProcess->SetVerboseLevel(static_cast<G4int>(PhysListOptPh::fVerboseLevel));
	fRayleighScatteringProcess->SetVerboseLevel(static_cast<G4int>(PhysListOptPh::fVerboseLevel));
	//fMieHGScatteringProcess->SetVerboseLevel(static_cast<G4int>(PhysListOptPh::fVerboseLevel));
	fBoundaryProcess->SetVerboseLevel(static_cast<G4int>(PhysListOptPh::fVerboseLevel));
	
	
	//G4OpticalSurfaceModel themodel = unified;
	
	//The unified model for the G4OpBoundaryProcess is the standard one
	//fBoundaryProcess->SetModel( (G4OpticalSurfaceModel)unified ); //G4OpticalSurfaceModel is an enumerator!
	
	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
	
	G4ParticleTable::G4PTblDicIterator* theParticleIterator = GetParticleIterator();
	theParticleIterator->reset();
	
	while( (*theParticleIterator)() ){
		G4ParticleDefinition *particle = theParticleIterator->value();
		G4ProcessManager *pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		
		
		
		if(particleName == "opticalphoton")
		{
			if(fVerboseLevel>=PhysVerbosity::kInfo){
				G4cout << "Info --> PhysListOptPh::ConstructOptical: Adding DiscreteProcesses for \"opticalphoton\" particle" << G4endl;
			}
			ph->RegisterProcess(fStepLimiter, particle);
			
			//pmanager->AddDiscreteProcess(fAbsorptionProcess);
			//pmanager->AddDiscreteProcess(fWLSProcess);
			//pmanager->AddDiscreteProcess(fRayleighScatteringProcess);
			//pmanager->AddDiscreteProcess(fMieHGScatteringProcess);
			//pmanager->AddDiscreteProcess(fBoundaryProcess);
			
			ph->RegisterProcess(fAbsorptionProcess, particle);
			ph->RegisterProcess(fWLSProcess, particle);
			ph->RegisterProcess(fRayleighScatteringProcess, particle);
			ph->RegisterProcess(fBoundaryProcess, particle);
			
		}
	}
}


// Cuts /////////////////////////////////////////////////////////////////////
void PhysListOptPh::SetCuts()
{
	/*
	if(verboseLevel > 1) G4cout << "PhysList::SetCuts:";

	
	SetCutsWithDefault();

	if(verboseLevel>0) DumpCutValuesTable();
	*/
}

