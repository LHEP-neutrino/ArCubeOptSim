#ifndef __ANALYSIS_MANAGER_OPT_PH_HH__
#define __ANALYSIS_MANAGER_OPT_PH_HH__

#include "globals.hh"
#include "G4ProcessTable.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TDirectory.h"


#include <set>
#include <vector>
#include <map>


class G4Run;
class G4Event;
class G4Step;
class G4SteppingManager;

class EventDataOptPh;
class PrimaryGeneratorActionOptPh;
class AnalysisOptPhMessenger;
class G4Navigator;

class AnalysisManagerOptPh
{
public:
	enum datasave
	{
		kOff, //No saving
		kHits, //Only stuff useful for LUT (hits, and primary info)
		kHitsExt, //Hits extended informations
		kSdSteps, //All the steps in the physical volumes defined by the user
		kAll //All the steps everywhere in any physical volume: same data as kSdSteps but everywhere
	};
	
	enum verbosity{
		kSilent,
		kInfo,
		kDetails,
		kDebug
	};
	
	AnalysisManagerOptPh(PrimaryGeneratorActionOptPh *pPrimaryGeneratorAction);
	virtual ~AnalysisManagerOptPh();
	void BeginOfRun(const G4Run *pRun); 
	void EndOfRun(const G4Run *pRun); 
	void BeginOfEvent(const G4Event *pEvent); 
	void EndOfEvent(const G4Event *pEvent); 
	void Step(const G4Step *pStep, const G4SteppingManager* pStepMan);	
	
	void SetRanSeed(G4int hRanSeed) { fRanSeed = hRanSeed; };
	
	void DefineOptPhSensDet(G4String volList);
	void SetSaveData(AnalysisManagerOptPh::datasave save){ fSave=save; };
	void SetStepsDebug(G4bool dbg=true){ fStepsDebug=dbg; };
	void SetDataFilename(const G4String &hFilename) { fDataFilename = hFilename; };
	void SetNbEventsToSimulate(G4int iNbEventsToSimulate) { fNbEventsToSimulate = iNbEventsToSimulate; };
	void SetNbOfPrimariesPerEvent(G4int nPrim){ fNbPrim = nPrim; };
	virtual void SetVerbosity(G4int verb){fVerbose=verb;};
	virtual void SetPrintModulo(G4int modulo){fPrintModulo=modulo;};
	
	void SetAutoFlush(Long64_t autof){fAutoFlushEvs=autof;};
	void SetAutoSave(Long64_t autos){fAutoSaveEvs=autos;};
	
	G4String GetDataFilename() const { return fDataFilename; };
	G4int GetNbEventsToSimulate() const { return fNbEventsToSimulate; };
	G4int GetNbOfPrimariesPerEvent() const { return fNbPrim; };
	AnalysisManagerOptPh::datasave GetSaveStatus() const {return fSave;};
	
	
private:
	
	//This function takes care of making the look-up tables for physical volumes at the start of the run
	void MakeVolMaps();
	
	//This function tries to make physical volumes with consecutive indexing if at the same tree level
	//It behaves very differently from the same function in the detector construction
	void ScanVols(const G4LogicalVolume* LogVol, int& volindex);
	
	int FindVolumeIndex( const G4VPhysicalVolume* aVolume );
	
	//This method should determine the unique touchable copy corresponding to a physical volume by means of the touchable history (the geometry tree)
	int FindVolId(G4TouchableHandle& touch);
	
	//This actually returns the process position in the table of processes (it is a vector)
	int FindProcessIndex( const G4VProcess* aProcess );
	
	//Makes the list of all processes and produces a json file like dictionary usable in analysis phase. Method executed at the start of the run.
	std::string BuildProcsDict();
	
	//Makes the list of all physical volumes and produces a json file like dictionary usable in analysis phase. Method executed at the start of the run. Only the volumes defined by the user are taken into account if the saving level is high enough. In the case that kAll is the defined saving mode the entire list of physical volumes is considered
	std::string BuildPysVolDict();
	
	//The same as before but only for the volumes that are sensitive volumes
	std::string BuildSDvolDict();
	
	
private:
	
	//Structure to hold data of primary particles at their first step (momentum energy polarization)
	typedef struct priminfo{
		G4double en;
		G4ParticleMomentum momVec;
		G4ThreeVector polVec;
	}priminfo;
	
	
	G4int fRanSeed;
	
	PrimaryGeneratorActionOptPh *fPrimaryGeneratorAction;
	AnalysisOptPhMessenger *fMessenger;
	EventDataOptPh *fEventData;
	
	G4Navigator* fNav;
	//G4TouchableHandle fTouchableHandle;
	
	TFile *fTreeFile;
	TTree *fTree;
	
	AnalysisManagerOptPh::datasave fSave;
	G4bool fStepsDebug;
	G4int fVerbose;
	G4int fPrintModulo;
	
	G4bool fWasAtBoundary;
	
	G4int fCurrentEvent;
	
	//G4int fHitsCollectionID;
	
	G4String fDataFilename;
	G4int fNbEventsToSimulate, fNbPrim;
	
	G4ProcessTable *fProcTable; //This should NEVER be instanced or removed
	G4ProcessVector *fProcVec; //This should NEVER be instanced or removed
	G4int fNprocs;
	
	
	Long64_t fAutoSaveEvs, fAutoFlushEvs;
	
	Long64_t fNbTotHits;
	
	
	std::map<G4VPhysicalVolume*, int> fPhysVolMap; //Map of all the existing physics volumes
	std::map<int, G4String> fPhysVolNamesMap; //Map of all the existing physics volumes: link their number to their name
	std::map<G4VPhysicalVolume*, int> fPhysVolUniqueMap; //This is to save in the TTree a unique index corresponding to a phys vol name
	std::map<int, G4String> fPhysVolUniqueNamesMap; //This is to make the correspondence of the volume index saved in the TTree to the name of phys volume in the json dictionary
	
	std::map<G4VPhysicalVolume*, int> fPhysVolCpnmMap; //Map of all the existing physics volumes by copy number
	std::map<G4VPhysicalVolume*, int> fOptPhSenDetVolPtrsMap; //Index of each physical volume
	
	std::set<G4VPhysicalVolume*> fOptPhSenDetVolPtrs;
	
	
	
	std::map<int, G4String> fOptPhPhysVolsMap;
	
	std::map<int, G4String> fOptPhProcessesMap;
	
	
	
	std::set<int> fTrackIDs; //Storage of the tracks IDs seen during an event when they are saved
	std::map<int, int> fTrackParentIDsMap; //Here for each track ID (first entry, key) there is the ID of the parent track
	std::map<int, int> fTrackGenerationMap; //Here for each track ID (first entry, key) there is the track generation level
	std::map<int, int> fFirstParentIDMap; //Here for each track ID (first entry, key) there is the track ID of the first parent (0 in case it is a primary track)
	std::map<int, int> fTrackCreatProc; //Here for each track ID (first entry, key) there is the process ID of the creation process (0 in casse of primary track)
	
	
	//Stuff used to reduce the number of calls at stepping level
	int fLastTrackId; //Track ID of the last step
	G4VPhysicalVolume* fLastPhysVol; //Pointer of the physical volume seen at the last step.
	int fLastVolIdx;
	int fLastCopyNum;
	int fLastVolId;
};

#endif

