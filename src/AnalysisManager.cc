#include "DetConstr.hh"
#include "PrimGenAction.hh"
#include "AnalysisManager.hh"
#include "AnalysisMessenger.hh"
//#include "OptPhHit.hh"
#include "EventData.hh"

#include "TNamed.h"

#include "G4String.hh"
#include "Randomize.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4TrackStatus.hh"
#include "G4StepStatus.hh"
#include "G4SteppingManager.hh"
#include "G4VProcess.hh"
#include "G4PhysicalVolumeStore.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"


#include "nlohmann/json.hpp"


#include <sys/time.h>
#include <numeric>
#include <fstream>
#include <sstream>


using std::vector;
using std::stringstream;
using std::set;
using std::ofstream;

using namespace CLHEP;

using json = nlohmann::json;


AnalysisManagerOptPh::AnalysisManagerOptPh(PrimaryGeneratorActionOptPh *pPrimaryGeneratorAction):
fPrimaryGeneratorAction(pPrimaryGeneratorAction),
fPartSrc(fPrimaryGeneratorAction->GetParticleSource()),
fRunSeed(0),
fNav(nullptr),
fVerbose(AnalysisVerbosity::kSilent),
fPrintModulo(0),
fCurrentEvent(-1),
fCurrentTrackId(-1),
fCurrentTrack(nullptr),
fSave(DatasaveLevel::kOff),
fStepsDebug(false),
fWasAtBoundary(false),
fProcTable(nullptr),
fProcVec(nullptr),
fNprocs(0),
fParticlesTable(nullptr),
fDataFilename("events.root"),
fTreeFile(nullptr),
fTree(nullptr),
fAutoSaveEvs(100),
fAutoFlushEvs(100),
fLastTrackId(-1),
fLastPhysVol(nullptr),
fLastPrePhysVol(nullptr)
{
	fMessenger = new AnalysisOptPhMessenger(this);
	
	fEventData = new EventDataOptPh();
}

AnalysisManagerOptPh::~AnalysisManagerOptPh()
{
	if(fEventData) delete fEventData;
	if(fMessenger) delete fMessenger;
}



void AnalysisManagerOptPh::BeginOfRun(const G4Run *pRun)
{
	//G4cout << "\nEntering in AnalysisManagerOptPh::BeginOfRun(...)" << G4endl;
	G4int randseed;
	
	if(fRunSeed > 0){
		CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
		CLHEP::HepRandom::setTheSeed(fRunSeed);
		randseed = fRunSeed;
	}else{
		// initialize with time.....
		struct timeval hTimeValue;
		gettimeofday(&hTimeValue, nullptr);
		CLHEP::HepRandom::setTheSeed(hTimeValue.tv_usec);
		randseed = hTimeValue.tv_usec;
	}
	
	if(fVerbose>=AnalysisVerbosity::kInfo) std::cout << "\nInfo --> AnalysisManager::BeginOfRun: Random numbers generator initialized with seed " << randseed << std::endl;
	
	
	fProcTable = G4ProcessTable::GetProcessTable();
	
	fCurrentEvent =- 1;
	
	fParticlesTable = G4ParticleTable::GetParticleTable();
	
	MakeVolMaps();
	//MakeParticlesMaps();
	
	std::string vol_dict = BuildPysVolDict();
	std::string sd_vol_dict = BuildSDvolDict();
	std::string proc_dict = BuildProcsDict();
	//std::string partdef_dict = BuildParticlesDict();
	
#ifndef NDEBUG	
	//if(fVerbose>=AnalysisVerbosity::kDebug){
		//G4cout << "\nDebug --> AnalysisManager::BeginOfRun: Content of the \"partdef_dict\":" <<  G4endl;
		//G4cout << "          " << partdef_dict << G4endl;
	//}
#endif
	
	
	if(fSave<=DatasaveLevel::kOff){
		G4cout << "\nWARNING --> AnalysisManagerOptPh::BeginOfRun(...): Data will not be saved." << G4endl;
		return;
	}
	
	fTreeFile = TFile::Open(fDataFilename.c_str(), "RECREATE", "File containing event data of optical photon simulations of ArgonCube");
	
	
	if(!fTreeFile){//Here there is a problem
		fSave = DatasaveLevel::kOff; //Do this to save from application crashing
		return;
	}
	
	
	
	TNamed *tn_vol_dict = new TNamed("vol_dict", vol_dict.c_str());
	
	TNamed *tn_sdvol_dict = nullptr;
	if(fOptPhSenDetVolPtrs.size()>0) tn_sdvol_dict = new TNamed("sdvol_dict", sd_vol_dict.c_str());
	
	TNamed *tn_proc_dict = new TNamed("proc_dict", proc_dict.c_str());
	
	TParameter<int>* ptTParNbEventsToSimulate = new TParameter<int>("EventsNb", fNbEventsToSimulate);
	
	fNbPrim = fPrimaryGeneratorAction->GetPrimNb();
	TParameter<int>* ptTParNbPrimariesPerEvent = new TParameter<int>("PrimNb", fNbPrim);
	
	TParameter<G4int>* ptTParRandSeed = new TParameter<G4int>("RandSeed", randseed);
	
	if(fTreeFile){
		fTreeFile->WriteTObject(tn_vol_dict, 0, "overwrite");
		if(tn_sdvol_dict) fTreeFile->WriteTObject(tn_sdvol_dict, 0, "overwrite");
		fTreeFile->WriteTObject(tn_proc_dict, 0, "overwrite");
		//fTreeFile->WriteTObject(tn_partdef_dict, 0, "overwrite");
		fTreeFile->WriteTObject(ptTParNbEventsToSimulate, 0, "overwrite");
		fTreeFile->WriteTObject(ptTParNbPrimariesPerEvent, 0, "overwrite");
		fTreeFile->WriteTObject(ptTParRandSeed, 0, "overwrite");
	}
	
	
	delete tn_vol_dict;
	if(tn_sdvol_dict) delete tn_sdvol_dict;
	delete tn_proc_dict;
	//delete tn_partdef_dict;
	delete ptTParNbEventsToSimulate;
	delete ptTParNbPrimariesPerEvent;
	delete ptTParRandSeed;
	
	
	
	fTree = new TTree("t1", "Tree containing event data for ArgonCube optical photon simulations.");

	gROOT->ProcessLine("#include <vector>");
	gROOT->ProcessLine("#include <string>");

	fTree->Branch("EvId", &fEventData->fEventId, "eventid/I");
	
//Primary particle savings
	//fTree->Branch("prim_type", "vector<Int_t>", &fEventData->fPrimPartType);
	
	fTree->Branch("prim_vol_index", &fEventData->fPrimaryVolumeIndex, "prim_vol_idx/I");//Fill at start of tracking stage only once
	fTree->Branch("prim_vol_cpnm", &fEventData->fPrimaryVolumeCopyNum, "prim_vol_cpnm/I");//Fill at start of tracking stage only once
	fTree->Branch("prim_vol_globcp", "string", &fEventData->fPrimaryVolumeGlobCp);//Fill at start of tracking stage only once
	fTree->Branch("prim_Xpos", &fEventData->fPrimary_Xpos, "prim_Xpos/D");//Fill at start of tracking stage only once
	fTree->Branch("prim_Ypos", &fEventData->fPrimary_Ypos, "prim_Ypos/D");//Fill at start of tracking stage only once
	fTree->Branch("prim_Zpos", &fEventData->fPrimary_Zpos, "prim_Zpos/D");//Fill at start of tracking stage only once
	fTree->Branch("prim_id", "vector<Int_t>", &fEventData->fPrimary_Id);//Fill at start of tracking stage only once
	fTree->Branch("prim_wavelen", "vector<Double_t>", &fEventData->fPrimWaveLength);//This is the wavelength of the optical photon (in nm units), taken from the Ekin. //Fill start of tracking stage only once
	if(fSave > DatasaveLevel::kLUT) fTree->Branch("prim_Xmom", "vector<Double_t>", &fEventData->fPrimary_Xmom);//Fill at start of tracking stage only once
	if(fSave > DatasaveLevel::kLUT) fTree->Branch("prim_Ymom", "vector<Double_t>", &fEventData->fPrimary_Ymom);//Fill at start of tracking stage only once
	if(fSave > DatasaveLevel::kLUT) fTree->Branch("prim_Zmom", "vector<Double_t>", &fEventData->fPrimary_Zmom);//Fill at start of tracking stage only once
	if(fSave > DatasaveLevel::kLUT) fTree->Branch("prim_Xpol", "vector<Double_t>", &fEventData->fPrimary_Xpol);//Fill at start of tracking stage only once
	if(fSave > DatasaveLevel::kLUT) fTree->Branch("prim_Ypol", "vector<Double_t>", &fEventData->fPrimary_Ypol);//Fill at start of tracking stage only once
	if(fSave > DatasaveLevel::kLUT) fTree->Branch("prim_Zpol", "vector<Double_t>", &fEventData->fPrimary_Zpol);//Fill at start of tracking stage only once
	
	
	
	//Hits related data
	if(fSave < DatasaveLevel::kSdSteps) fTree->Branch("totalhits", &fEventData->fNbTotHits, "totalhits/L");
	//if(fSave < DatasaveLevel::kSdSteps) fTree->Branch("hit_part_type", "vector<Int_t>", &fEventData->fPartType);
	if(fSave < DatasaveLevel::kSdSteps) fTree->Branch("hit_vol_index", "vector<Int_t>", &fEventData->fVolIndex);//ID of the physical volume (it is a whish!). //Fill at step stage
	if(fSave < DatasaveLevel::kSdSteps) fTree->Branch("hit_vol_copy", "vector<Int_t>", &fEventData->fHitVolCopyNum);//This is the copy number of a specific physics volume. //Fill at step stage
	if(fSave < DatasaveLevel::kSdSteps) fTree->Branch("hit_vol_globcp", "vector<string>", &fEventData->fHitVolGlobCp);//This MUST become the unique ID of the touchable volume. //Fill at step stage
	if(fSave < DatasaveLevel::kSdSteps) fTree->Branch("hit_time", "vector<Double_t>", &fEventData->fTime);//Fill at step stage
	//if(fSave < DatasaveLevel::kSdSteps) fTree->Branch("hit_trackid", "vector<Int_t>", &fEventData->fTrackId);//Fill at step stage
	if(fSave < DatasaveLevel::kSdSteps) fTree->Branch("hit_firstparentid", "vector<Int_t>", &fEventData->fFirstParentId); //Fill at step stage
	if(fSave < DatasaveLevel::kSdSteps) fTree->Branch("hit_phot_wavelen", "vector<Double_t>", &fEventData->fPhWaveLength); //This is the wavelength of the optical photon (in nm units), taken from the Ekin. //Fill at step stage 
	
	
	//Extended information of hits related data
	if(fSave == DatasaveLevel::kHitsExt) fTree->Branch("hit_trackid", "vector<Int_t>", &fEventData->fTrackId);//Fill at step stage
	if(fSave == DatasaveLevel::kHitsExt) fTree->Branch("hit_partgen", "vector<Int_t>", &fEventData->fPartGener);//Fill at step stage
	//if(fSave == DatasaveLevel::kHitsExt) fTree->Branch("hit_phot_wvlen", "vector<Double_t>", &fEventData->fPhWaveLength); //This is the wavelength of the optical photon (in nm units), taken from the Ekin. /Fill at step stage.
	if(fSave == DatasaveLevel::kHitsExt) fTree->Branch("hit_xpos", "vector<Double_t>", &fEventData->fXpos);//Fill at step stage
	if(fSave == DatasaveLevel::kHitsExt) fTree->Branch("hit_ypos", "vector<Double_t>", &fEventData->fYpos);//Fill at step stage
	if(fSave == DatasaveLevel::kHitsExt) fTree->Branch("hit_zpos", "vector<Double_t>", &fEventData->fZpos);//Fill at step stage
	
	if(fSave == DatasaveLevel::kHitsExt) fTree->Branch("hit_xmom", "vector<Double_t>", &fEventData->fXmom);//Fill at step stage
	if(fSave == DatasaveLevel::kHitsExt) fTree->Branch("hit_ymom", "vector<Double_t>", &fEventData->fYmom);//Fill at step stage
	if(fSave == DatasaveLevel::kHitsExt) fTree->Branch("hit_zmom", "vector<Double_t>", &fEventData->fZmom);//Fill at step stage
	
	if(fSave == DatasaveLevel::kHitsExt) fTree->Branch("hit_xpol", "vector<Double_t>", &fEventData->fXpol);//Fill at step stage
	if(fSave == DatasaveLevel::kHitsExt) fTree->Branch("hit_ypol", "vector<Double_t>", &fEventData->fYpol);//Fill at step stage
	if(fSave == DatasaveLevel::kHitsExt) fTree->Branch("hit_zpol", "vector<Double_t>", &fEventData->fZpol);//Fill at step stage
	
	
	//Full step mode
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("totsteps", &fEventData->fNbTotHits, "totsteps/I");//Fill at step stage
	//if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("part_type", "vector<Int_t>", &fEventData->fPartType);
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("vol_index", "vector<Int_t>", &fEventData->fVolIndex);//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("vol_copy", "vector<Int_t>", &fEventData->fHitVolCopyNum);//ID of the touchable volume//Fill at step stage
	//if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("vol_id", "vector<Long64_t>", &fEventData->fHitVolId);//ID of the touchable volume//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("vol_globcp", "vector<string>", &fEventData->fHitVolGlobCp);//ID of the touchable volume//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("time", "vector<Double_t>", &fEventData->fTime);//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("trackid", "vector<Int_t>", &fEventData->fTrackId);//Fill at end of Event
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("partgener", "vector<Int_t>", &fEventData->fPartGener);//Fill at end of Event
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("parentid", "vector<Int_t>", &fEventData->fParentId);//Fill at end of Event
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("firstparentid", "vector<Int_t>", &fEventData->fFirstParentId);//Fill at end of Event
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("creatproc", "vector<Int_t>", &fEventData->fCreatProc);//Fill at end of Event
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("deposproc", "vector<Int_t>", &fEventData->fDepProc);//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("wavelen", "vector<Double_t>", &fEventData->fPhWaveLength);//This is the wavelength of the optical photon (in nm units), taken from the Ekin. //Fill at step stage
	
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("xpos", "vector<Double_t>", &fEventData->fXpos);//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("ypos", "vector<Double_t>", &fEventData->fYpos);//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("zpos", "vector<Double_t>", &fEventData->fZpos);//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("xmom", "vector<Double_t>", &fEventData->fXmom);//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("ymom", "vector<Double_t>", &fEventData->fYmom);//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("zmom", "vector<Double_t>", &fEventData->fZmom);//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("xpol", "vector<Double_t>", &fEventData->fXpol);//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("ypol", "vector<Double_t>", &fEventData->fYpol);//Fill at step stage
	if(fSave >= DatasaveLevel::kSdSteps) fTree->Branch("zpol", "vector<Double_t>", &fEventData->fZpol);//Fill at step stage
	
	
	fEventData->fPrimaryVolumeIndex = -1;
	
	//These assignment are for preallocation of memory.
	//The vectors will be resized to 0 at start of event
	fEventData->fPrimWaveLength->assign(fNbPrim,0);
	
	fEventData->fPrimary_Xmom->assign(fNbPrim,0);
	fEventData->fPrimary_Ymom->assign(fNbPrim,0);
	fEventData->fPrimary_Zmom->assign(fNbPrim,0);
	
	fEventData->fPrimary_Xpol->assign(fNbPrim,0);
	fEventData->fPrimary_Ypol->assign(fNbPrim,0);
	fEventData->fPrimary_Zpol->assign(fNbPrim,0);
	
	//fTree->SetMaxTreeSize((int)1e6);
	fTree->SetAutoFlush(fAutoFlushEvs);
	fTree->SetAutoSave(fAutoSaveEvs);
	
	fTrackIDs.clear();
	fTrackParentIDsMap.clear();
	fTrackGenerationMap.clear();
	fFirstParentIDMap.clear();
	
	
#ifndef NDEBUG
	if(fVerbose>=AnalysisVerbosity::kDebug) G4cout << "\nDebug --> Exiting from AnalysisManager::BeginOfRun" << G4endl;
#endif
}


void AnalysisManagerOptPh::EndOfRun(const G4Run *pRun)
{
	fProcTable = nullptr;
	fProcVec = nullptr;
	
	if(fSave==DatasaveLevel::kOff) return;
	
	if(fTreeFile){
		if(fTree){
			fTreeFile->WriteTObject(fTree, 0, "overwrite");
			delete fTreeFile; //This deletes also the TTree owned by the TFile
			fTree=nullptr;
		}
		fTreeFile=nullptr;
	}
}


void AnalysisManagerOptPh::BeginOfEvent(const G4Event *pEvent)
{
	// grab event ID
	fCurrentEvent = pEvent->GetEventID();
	
	// print this information event by event (modulo n)  	
	if( (fVerbose>=AnalysisVerbosity::kInfo) && (fPrintModulo>0) && (fCurrentEvent%fPrintModulo == 0) ){
		G4cout << "\nInfo --> AnalysisManagerOptPh::BeginOfEvent(...): Begin of event: " << fCurrentEvent  << G4endl;
	}
	
	
	fWasAtBoundary = false;
	
	//These initialisations are needed at the start of event
	fLastTrackId = -1;
	fLastPhysVol = nullptr;
	fLastVolIdx = -1;
	fLastCopyNum = -1;
	fLastVolGlobalCopy = std::string("");

	fEventData->Reset();
	
	
	fTrackIDs.clear();
	fTrackParentIDsMap.clear();
	fTrackGenerationMap.clear(); 
	fFirstParentIDMap.clear();
	
	
	//Primary particles information
	//Volume id.......
	G4ThreeVector posVec = fPrimaryGeneratorAction->GetPrimPos();
	
	//Here I can still use the tracking navigator as the tracking has not started yet
	if(!fNav) fNav = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
	G4VPhysicalVolume *fPrimVol = fNav->LocateGlobalPointAndSetup(posVec, nullptr, true);
	G4TouchableHandle touch = fNav->CreateTouchableHistory();
	
	
	if(fSave>DatasaveLevel::kOff){
		fEventData->fEventId = fCurrentEvent;
		
		fEventData->fPrimaryVolumeIndex = FindVolumeIndex(fPrimVol);
		fEventData->fPrimaryVolumeCopyNum = touch->GetCopyNumber();
		(*fEventData->fPrimaryVolumeGlobCp) = FindVolGlobalCopy(touch);
		
		fEventData->fPrimary_Xpos = posVec.x();
		fEventData->fPrimary_Ypos = posVec.y();
		fEventData->fPrimary_Zpos = posVec.z();
	}
	
	
#ifndef NDEBUG
	if(fVerbose>=AnalysisVerbosity::kDebug){
		G4cout << "Debug --> AnalysisManagerOptPh::BeginOfEvent(...): EventID: " << fCurrentEvent << "; Primary volume: " << fPrimVol->GetName() << "; Copy number: " << touch->GetCopyNumber() << G4endl;
	}
#endif
}


void AnalysisManagerOptPh::EndOfEvent(const G4Event *pEvent)
{
	if(fSave>DatasaveLevel::kOff){
		if(fTree) fTree->Fill();
	}
}


///////////////////////////////////////////////////////

// G4Mutex mutex_an_man_tracking = G4MUTEX_INITIALIZER;
void AnalysisManagerOptPh::PreUserTrackingAction(const G4Track* pTrack){
	
	//G4AutoLock l(&mutex_an_man_tracking);
	//l.lock();
	
	fCurrentTrack = (G4Track*)pTrack;
	G4int fCurrentTrackId = fCurrentTrack->GetTrackID();
	G4int parentid = fCurrentTrack->GetParentID();
	G4bool isPrimary = (!fCurrentTrack->GetCreatorProcess()) || (parentid<=0);
	
	//Check if the track is known and if it is a primary track
	//This book-keeping code must be executed whatever the saving options are
	if( fTrackIDs.find(fCurrentTrackId)==fTrackIDs.end() ){
		fTrackParentIDsMap[fCurrentTrackId] = parentid;
		if(isPrimary){ //It is a primary track at its very first step (in the primary volume)
			fTrackGenerationMap[fCurrentTrackId] = 1;
			fFirstParentIDMap[fCurrentTrackId] = fCurrentTrackId;
			
			if(fSave>DatasaveLevel::kOff){
				//If data will be saved the status of the primary track must be saved
				/*
				if(fCurrentTrack->GetParticleDefinition()->IsGeneralIon()){
					//fEventData->fPrimPartType->push_back(fGenIonId);
				}else{
					fEventData->fPrimPartType->push_back(fParticlesDefsMap[fCurrentTrack->GetParticleDefinition()]);
				}
				*/
				
				fEventData->fPrimary_Id->push_back(fCurrentTrackId);
				
				fEventData->fPrimWaveLength->push_back( h_Planck*c_light/(fCurrentTrack->GetVertexKineticEnergy())/nm );
				
				const G4ThreeVector& momdir = fCurrentTrack->GetVertexMomentumDirection();
				
				fEventData->fPrimary_Xmom->push_back( momdir.x() );
				fEventData->fPrimary_Ymom->push_back( momdir.y() );
				fEventData->fPrimary_Zmom->push_back( momdir.z() );
				
				//There is no a G4Track::GetVertexPolarizationDirection() interface and must use this way
				if(fCurrentTrack->GetCurrentStepNumber()==0){
					//This is a redundant control but useful anyway
					const G4ThreeVector& poldir = fCurrentTrack->GetPolarization();
					
					fEventData->fPrimary_Xpol->push_back( poldir.x() );
					fEventData->fPrimary_Ypol->push_back( poldir.y() );
					fEventData->fPrimary_Zpol->push_back( poldir.z() );
				}
				
			}
		}else{
			fTrackGenerationMap[fCurrentTrackId] = fTrackGenerationMap[parentid]+1;
			fFirstParentIDMap[fCurrentTrackId] = fFirstParentIDMap[parentid];
		}
		
		
		if( fProcTable && (fSave>=DatasaveLevel::kSdSteps) ){
			if(isPrimary){//It's a primary track
				fTrackCreatProc[fCurrentTrackId] = 0;
			}else{
				int retcode = FindProcessIndex(fCurrentTrack->GetCreatorProcess()); //This is the index if the process is found (>=0), otherwise a negative return code is returned
				if(retcode>=0){
					fTrackCreatProc[fCurrentTrackId] = retcode+1; //Add 1 to the index as 0 is reserved for primary tracks (no creator process)
				}else{
					fTrackCreatProc[fCurrentTrackId] = retcode; //This is a negative number and indicates for problems
				}
			}
		}
		
		fTrackIDs.insert(fCurrentTrackId);
	}
	//l.unlock();
}


///////////////////////////////////////////////////////
void AnalysisManagerOptPh::Step(const G4Step *pStep, const G4SteppingManager* pStepMan)
{
#ifndef NDEBUG
	if(fVerbose>=AnalysisVerbosity::kDebug) G4cout << "\nDebug ---> AnalysisManagerOptPh::Step(...): Entering in AnalysisManagerOptPh::Step\n" << G4endl;
#endif
	
	//if( (fVerbose<3) && (fOptPhSenDetVolNames.empty()) )return;
	if(!pStep) return; //This just avoids problems, although would be a big problem to be fixed
	
	//G4AutoLock l(&mutex_an_man);
	//l.lock();
	
	if(fCurrentTrack!=pStep->GetTrack()){
			G4cerr << "\nERROR --> AnalysisManagerOptPh::Step: The track pointer from the G4Step is different from the \"fCurrentTrack\" pointer for Event ID: " << fCurrentEvent << "." << G4endl;
	}
	
	
	G4StepPoint *preStepPoint = pStep->GetPreStepPoint();
	
	G4TouchableHandle pre_touch = preStepPoint->GetTouchableHandle();
	G4VPhysicalVolume *pre_Vol = preStepPoint->GetTouchableHandle()->GetVolume();
	
	
	G4StepPoint *postStepPoint = pStep->GetPostStepPoint();
	
	G4TouchableHandle post_touch = postStepPoint->GetTouchableHandle();
	G4VPhysicalVolume *post_Vol = postStepPoint->GetTouchableHandle()->GetVolume();
	
	
	//G4String partName = fCurrentTrack->GetParticleDefinition()->GetParticleName();
	
	//This is the step point from where the stuff is saved. It changes to preStepPoint only in saving mode and in the the case the photon is absorbed in the physical volume soon after it went through the boundary surface.
	G4StepPoint *saveStepPoint = postStepPoint;
	G4VPhysicalVolume *saveVol = post_Vol; //This is from where I get the current quantities
	G4TouchableHandle *saveTouch = &post_touch; //This defines where the step occurred
	
	
	if(!post_Vol) return; //Most likely it is at the world boundary (the track will be killed)
	
	if(fOptPhSenDetVolPtrs.find(pre_Vol)!=fOptPhSenDetVolPtrs.end()){
		
		//Kill the optical photon track in hit saving mode (I want to save the track only the point of entry)
		if( pStep->IsFirstStepInVolume() && ((fSave>DatasaveLevel::kOff) && (fSave<DatasaveLevel::kSdSteps))){
			fCurrentTrack->SetTrackStatus(fStopAndKill);
			saveStepPoint = preStepPoint;
			saveVol = pre_Vol; 
			saveTouch = &pre_touch;
		}else if(fSave>=DatasaveLevel::kSdSteps){
			//Here I am in stepping mode and I save the informations of the post step point and the track is not killed
			saveStepPoint = postStepPoint;
			saveVol = post_Vol; 
			saveTouch = &post_touch;
		}
	}
	
	/*
	//Old section to be removed
	if(fWasAtBoundary && (fOptPhSenDetVolPtrs.find(pre_Vol)!=fOptPhSenDetVolPtrs.end()) && pStep->IsFirstStepInVolume() ){
		//This is the first step inside the new volume where the optical photon gets absorbed
		if( (fSave>DatasaveLevel::kOff) && (fSave<DatasaveLevel::kSdSteps) ){
			fCurrentTrack->SetTrackStatus(fStopAndKill);
			
			saveStepPoint = preStepPoint;
			saveVol = pre_Vol; 
			saveTouch = &pre_touch;
		}
	}
	*/
	
	//If the volume is defined as an absorption volume I kill the optical track anyway
	if(fOptPhAbsVolPtrs.find(pre_Vol)!=fOptPhAbsVolPtrs.end()){
		fCurrentTrack->SetTrackStatus(fStopAndKill);
	}
	
	G4TrackStatus trstatus = fCurrentTrack->GetTrackStatus();
	
	
#ifndef NDEBUG	
	//Volume printouts
	if(fVerbose>=AnalysisVerbosity::kDebug || fStepsDebug){
		
		G4String TrackStat = "";
		
		
		if(trstatus==fAlive) TrackStat = "Alive";
		if(trstatus==fStopButAlive) TrackStat = "StopButAlive";
		if(trstatus==fStopAndKill) TrackStat = "StopAndKill";
		if(trstatus==fKillTrackAndSecondaries) TrackStat = "KillTrackAndSecondaries";
		if(trstatus==fSuspend) TrackStat = "Suspend";
		if(trstatus==fPostponeToNextEvent) TrackStat = "PostponeToNextEvent";
		
		
		if(fStepsDebug){ //Boundary related printouts
			if((postStepPoint->GetStepStatus()==fGeomBoundary) || fWasAtBoundary){
				if(!fWasAtBoundary){ //Printout at post step point boundary
					G4cout << "\nStepDebug --> Event " << fCurrentEvent << ", trackID: " << fCurrentTrack->GetTrackID() << ". Optical photon at volumes boundary" << G4endl;
					G4cout << "               Volume 1: <" << pre_Vol->GetName() << ">, copy num: " << pre_touch->GetCopyNumber() << G4endl; 
					G4cout << "               Volume 2: <" << post_Vol->GetName() << ">, copy num:" << post_touch->GetCopyNumber() << G4endl;
					G4cout << "           Track status: " << TrackStat << G4endl;
					G4cout << "            At boundary: " << (postStepPoint->GetStepStatus()==fGeomBoundary) << G4endl;
					G4cout << "               Sel proc: " << postStepPoint->GetProcessDefinedStep()->GetProcessName() << G4endl;
					G4cout << "   First step in volume: " << pStep->IsFirstStepInVolume() << G4endl;
					G4cout << "    Last step in volume: " << pStep->IsLastStepInVolume() << G4endl;
					G4cout << "       Step length (mm): " << pStep->GetStepLength()/mm << G4endl;
				}else{ //Printout after boundary post step point
					G4cout << "\nStepDebug --> Event " << fCurrentEvent << ", trackID: " << fCurrentTrack->GetTrackID() << ". Optical photon after volumes boundary" << G4endl;
					G4cout << "               Volume 1: <" << pre_Vol->GetName() << ">, copy num: " << pre_touch->GetCopyNumber() << G4endl; 
					G4cout << "               Volume 2: <" << post_Vol->GetName() << ">, copy num:" << post_touch->GetCopyNumber() << G4endl;
					G4cout << "           Track status: " << TrackStat << G4endl;
					G4cout << "            At boundary: " << (postStepPoint->GetStepStatus()==fGeomBoundary) << G4endl;
					G4cout << "               Sel proc: " << postStepPoint->GetProcessDefinedStep()->GetProcessName() << G4endl;
					G4cout << "   First step in volume: " << pStep->IsFirstStepInVolume() << G4endl;
					G4cout << "    Last step in volume: " << pStep->IsLastStepInVolume() << G4endl;
					G4cout << "       Step length (mm): " << pStep->GetStepLength()/mm << G4endl;
				}
			}
		}else{
			//Print out for any condition
			//if( (fSave==DatasaveLevel::kAll) || (fOptPhSenDetVolPtrs.find(saveVol)!=fOptPhSenDetVolPtrs.end()) ){
			
			if(true){
				G4String creatProc("Primary");
				if(fCurrentTrack->GetCreatorProcess()) creatProc = fCurrentTrack->GetCreatorProcess()->GetProcessName();
				
				G4double preGlobTime = 0.;
				if(preStepPoint) preGlobTime = preStepPoint->GetGlobalTime()/CLHEP::ns;
				
				G4cout << "Debug ---> AnalysisManagerOptPh::Step:\n" << G4endl;
				G4cout << "Event " << fCurrentEvent << " | trackID: " << fCurrentTrack->GetTrackID() << " | Parent ID: " << fCurrentTrack->GetParentID() << " | Creat. Proc: " << creatProc << G4endl;
				G4cout << "               Step num: " << fCurrentTrack->GetCurrentStepNumber() << G4endl;
				G4cout << "               Volume 1: <" << pre_Vol->GetName() << ">, copy num: " << pre_touch->GetCopyNumber() << G4endl; 
				G4cout << "               Volume 2: <" << post_Vol->GetName() << ">, copy num:" << post_touch->GetCopyNumber() << G4endl;
				G4cout << "           Track status: " << TrackStat << G4endl;
				G4cout << "            At boundary: " << (postStepPoint->GetStepStatus()==fGeomBoundary) << G4endl;
				G4cout << "               Sel proc: " << postStepPoint->GetProcessDefinedStep()->GetProcessName() << G4endl;
				G4cout << "   First step in volume: " << pStep->IsFirstStepInVolume() << G4endl;
				G4cout << "    Last step in volume: " << pStep->IsLastStepInVolume() << G4endl;
				G4cout << "       Step length (mm): " << pStep->GetStepLength()/mm << G4endl;
				G4cout << "        Wavelength (nm): " << h_Planck*c_light/postStepPoint->GetKineticEnergy()/nm << G4endl;
				if(preStepPoint){//Avoids out of world steps
					G4cout << "     Track time P1 (ns): " << preGlobTime << G4endl;
					G4cout << "     Track time P2 (ns): " << postStepPoint->GetGlobalTime()/CLHEP::ns << G4endl;
				}
				
			}
		}
	} //Closes if(fVerbose>=AnalysisVerbosity::kDebug || fStepsDebug)...
#endif	
	
	
	if((postStepPoint->GetStepStatus()==fGeomBoundary) && (!fWasAtBoundary)){
		fWasAtBoundary = true;
		fLastPrePhysVol = pre_Vol;
	}else{
		//For optical photons it might happen that in successive steps they go from a boundary to another
		if(postStepPoint->GetStepStatus()!=fGeomBoundary) fWasAtBoundary = false;
	}
	
	
	if( fSave==DatasaveLevel::kOff ){
#ifndef NDEBUG
		if(fVerbose>=AnalysisVerbosity::kDebug) G4cout << "\nDebug --> AnalysisManagerOptPh::Step: Exiting from method without saving.\n" << G4endl;
#endif
		return;
	}
	
	
	//For saves modes lower than kAll check whether the particle is in one of the sensitive volumes defined by the user
	if( (fSave<DatasaveLevel::kAll) ){
		//Here the mode is either "SD stepping mode" or one of the "hits modes"
		
		if( fOptPhSenDetVolPtrs.find(saveVol)==fOptPhSenDetVolPtrs.end() ){
			//The step occurred in a sens vol if the pre_step point physical volume is a sensitive volume
			//In all other saving modes I am interested only in hits or steps in specific physical volumes (sensitive volumes)
#ifndef NDEBUG
			if(fVerbose>=AnalysisVerbosity::kDebug){
				G4cout << "\nDebug --> AnalysisManagerOptPh::Step: Hit in <" << pre_Vol->GetName() << "> volume. Exiting the function." << G4endl;
			}
#endif
			return;
		}
		
		if( (fSave!=DatasaveLevel::kSdSteps) && (trstatus!=fStopAndKill) ){
			
			//When in hit mode the hit is saved only if the optical photon is going to be absorbed (killed) in one of the sensitive volumes
			return;
		}
	}
	
	
	//Here start to get stuff to be saved
	fEventData->fNbTotHits += 1;//This is the number of the total recorded steps in "stepping mode" or the number absorption (track stopped and kiled) for an optical photon is in a SD volume when in "hits mode"
	
	
	//Recalculate the volume id (recursive process) only if the volume pointer is different from before
	if( (pre_Vol!=fLastPhysVol) || (fCurrentTrackId!=fLastTrackId) ){
		fLastTrackId = fCurrentTrackId;
		/*
		if(fCurrentTrack->GetParticleDefinition()->IsGeneralIon()){
			fLastPartType = fGenIonId;
		}else{
			fLastPartType = fParticlesDefsMap[fCurrentTrack->GetParticleDefinition()];
		}
		*/
		fLastPhysVol = pre_Vol;
		fLastVolIdx = fPhysVolUniqueMap[pre_Vol];
		fLastCopyNum = pre_touch->GetCopyNumber();
		fLastVolGlobalCopy = FindVolGlobalCopy(pre_touch);
	}
	
	
	
	fEventData->fVolIndex->push_back( fLastVolIdx ); //ID of the physical volume (from a std::map)
	fEventData->fHitVolCopyNum->push_back( fLastCopyNum ); //Copy number of the physical volume
	fEventData->fHitVolGlobCp->push_back( fLastVolGlobalCopy );
	fEventData->fFirstParentId->push_back( fFirstParentIDMap[fLastTrackId] );
	fEventData->fTime->push_back( saveStepPoint->GetGlobalTime() );
	fEventData->fPhWaveLength->push_back( h_Planck*c_light/(saveStepPoint->GetKineticEnergy())/nm );
	//fEventData->fPartType->push_back(fLastPartType);
	
	if(fSave>=DatasaveLevel::kHitsExt){
		
		fEventData->fTrackId->push_back( fLastTrackId );
		
		//fEventData->fPhWaveLength->push_back( saveStepPoint->GetKineticEnergy() );
		fEventData->fXpos->push_back( (saveStepPoint->GetPosition()).x() );
		fEventData->fYpos->push_back( (saveStepPoint->GetPosition()).y() );
		fEventData->fZpos->push_back( (saveStepPoint->GetPosition()).z() );
		fEventData->fXmom->push_back( (saveStepPoint->GetMomentumDirection()).x() );
		fEventData->fYmom->push_back( (saveStepPoint->GetMomentumDirection()).y() );
		fEventData->fZmom->push_back( (saveStepPoint->GetMomentumDirection()).z() );
		fEventData->fXpol->push_back( (saveStepPoint->GetPolarization()).x() );
		fEventData->fYpol->push_back( (saveStepPoint->GetPolarization()).y() );
		fEventData->fZpol->push_back( (saveStepPoint->GetPolarization()).z() );
		
		if(fSave>=DatasaveLevel::kSdSteps){
			if(fProcTable){
				
				if(!saveStepPoint->GetProcessDefinedStep()){
					fEventData->fDepProc->push_back( 0 );//This is a primary track!
				}else{
					int retcode = FindProcessIndex(saveStepPoint->GetProcessDefinedStep()); //This is the index if the process is found (>=0), otherwise a negative return code is returned
					if(retcode>=0){
						fEventData->fDepProc->push_back( retcode+1 );//Add 1 to the index as 0 is reserved for primary tracks
					}else{
						fEventData->fDepProc->push_back( retcode );//This is a negative number and indicates there are problems
					}
				}
			}
		} // if(fSave>=kSdSteps)...
	} // if(fSave>=kHitsExt)...
	
#ifndef NDEBUG
	if(fVerbose>=AnalysisVerbosity::kDebug || fStepsDebug) G4cout << "\nDebug --> AnalysisManagerOptPh::Step: Exiting from method.\n" << G4endl;
#endif
}



//Service methods
void AnalysisManagerOptPh::DefineOptPhSensDet(G4String volList)
{
	fOptPhSenDetVolPtrs.clear();
	
	if(volList == G4String("NULL")){
		return;
	}
	
	G4PhysicalVolumeStore *pPhysVolStore = G4PhysicalVolumeStore::GetInstance();
	
	G4int nVols = pPhysVolStore->size();
	
	if(nVols <= 0) return;
	
	stringstream hStream;
	hStream.str(volList);
	G4String hVolumeName;
	
	
	// store all the volume names
	std::set<G4String> candidatevolnames;
	while(!hStream.eof()){
		hStream >> hVolumeName;
		candidatevolnames.insert(hVolumeName);
		if(!hStream) continue;
	}
	
	
	for(set<G4String>::iterator pIt = candidatevolnames.begin(); pIt != candidatevolnames.end(); pIt++){
		G4String hRequiredVolumeName = *pIt;
		G4bool bMatch = (hRequiredVolumeName.last('*') != std::string::npos);
		
		if(bMatch) hRequiredVolumeName = hRequiredVolumeName.strip(G4String::trailing, '*');
		
		for(G4int iVol=0; iVol<nVols; iVol++){
			G4String hName = pPhysVolStore->at(iVol)->GetName();
			
			if( (hName == hRequiredVolumeName) || (bMatch && (hName.substr(0, hRequiredVolumeName.size())) == hRequiredVolumeName) ){
				fOptPhSenDetVolPtrs.insert(pPhysVolStore->at(iVol));
			}
		}
	}
}


void AnalysisManagerOptPh::DefineOptPhAbsVols(G4String volList)
{
	fOptPhAbsVolPtrs.clear();
	
	
	if(volList == G4String("NULL")){
		return;
	}
	
	G4PhysicalVolumeStore *pPhysVolStore = G4PhysicalVolumeStore::GetInstance();
	
	G4int nVols = pPhysVolStore->size();
	
	if(nVols <= 0) return;
	
	
	stringstream hStream;
	hStream.str(volList);
	G4String hVolumeName;
	
	
	// store all the volume names
	std::set<G4String> candidatevolnames;
	while(!hStream.eof()){
		hStream >> hVolumeName;
		candidatevolnames.insert(hVolumeName);
		if(!hStream) continue;
	}
	
	
	for(set<G4String>::iterator pIt = candidatevolnames.begin(); pIt != candidatevolnames.end(); pIt++){
		G4String hRequiredVolumeName = *pIt;
		G4bool bMatch = (hRequiredVolumeName.last('*') != std::string::npos);
		
		if(bMatch) hRequiredVolumeName = hRequiredVolumeName.strip(G4String::trailing, '*');
		
		for(G4int iVol=0; iVol<nVols; iVol++){
			G4String hName = pPhysVolStore->at(iVol)->GetName();
			
			if( (hName == hRequiredVolumeName) || (bMatch && (hName.substr(0, hRequiredVolumeName.size())) == hRequiredVolumeName) ){
				fOptPhAbsVolPtrs.insert(pPhysVolStore->at(iVol));
			}
		}
	}
}


///////////////////////////////////////////////////////

void AnalysisManagerOptPh::MakeVolMaps()
{
	fPhysVolMap.clear();
	fPhysVolNamesMap.clear();
	fPhysVolUniqueNamesMap.clear();
	fPhysVolCpnmMap.clear();
	
	const G4VPhysicalVolume* worldPhysVol = (dynamic_cast<const DetConstrOptPh*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction()) )->GetWorldVolume();
	
	if(worldPhysVol->GetMotherLogical()){
		//This is not the world volume!!!
		return;
	}
	
	
	int volindex = 0;
	fPhysVolMap[(G4VPhysicalVolume*)worldPhysVol] = volindex;
	fPhysVolNamesMap[volindex] = worldPhysVol->GetName();
	fPhysVolCpnmMap[(G4VPhysicalVolume*)worldPhysVol] = worldPhysVol->GetCopyNo();
	
	ScanVols(worldPhysVol->GetLogicalVolume(), volindex);
	
	
	int idx=0;
	std::map<G4VPhysicalVolume*, int>::iterator it;
	std::map<G4String, int> names_map;
	for(it=fPhysVolMap.begin(); it!=fPhysVolMap.end(); it++){
		if( names_map.find(it->first->GetName())==names_map.end() ){
			fPhysVolUniqueNamesMap[idx] = it->first->GetName();
			names_map[it->first->GetName()] = idx;
			idx++;
		}
		fPhysVolUniqueMap[it->first] = names_map[it->first->GetName()];
	}
}


///////////////////////////////////////////////////////

void AnalysisManagerOptPh::ScanVols(const G4LogicalVolume* LogVol, int& volindex)
{
	G4int nDaught= LogVol->GetNoDaughters();
	
	G4VPhysicalVolume* PhysVol;
	
	for(G4int iDtr=0; iDtr<nDaught; iDtr++){
		PhysVol = LogVol->GetDaughter(iDtr);
		if( fPhysVolMap.find( PhysVol )==fPhysVolMap.end() ){
			volindex++;
			fPhysVolMap[PhysVol] = volindex;
			fPhysVolNamesMap[volindex] = PhysVol->GetName();
			fPhysVolCpnmMap[PhysVol] = PhysVol->GetCopyNo();
		}
	}
	
	
	for(G4int iDtr=0; iDtr<nDaught; iDtr++){
		PhysVol = LogVol->GetDaughter(iDtr);
		ScanVols(PhysVol->GetLogicalVolume(), volindex);
	}
}


///////////////////////////////////////////////////////

int AnalysisManagerOptPh::FindProcessIndex( const G4VProcess* aProcess )
{
	if(!aProcess) return -1;//This is a return code
	
	if(!fProcVec){
		if(!fProcTable){
			fProcTable = G4ProcessTable::GetProcessTable();
			if(!fProcTable) return -2; //This is a problem!
		}
		fProcVec = fProcTable->FindProcesses();
		if(!fProcVec) return -3;
		fNprocs = fProcVec->size();
	}
	
	for (int iProc = 0; iProc<fNprocs; iProc++) {
		if((*fProcVec)(iProc)==aProcess){
			return iProc;
		}
	}
	
	return -4;//This should not happen at this stage as the process is not found in the list of all processes
	
}


///////////////////////////////////////////////////////

#include "nlohmann/json.hpp"

int AnalysisManagerOptPh::FindVolumeIndex( const G4VPhysicalVolume* aVolume )
{
	if(!aVolume) return -1;//This is an error return code
	
	if(fPhysVolUniqueMap.size()==0){
		return -2;
	}
	
	return fPhysVolUniqueMap[(G4VPhysicalVolume*)aVolume];
}


///////////////////////////////////////////////////////

std::string AnalysisManagerOptPh::BuildProcsDict()
{
	fOptPhProcessesMap.clear();
	
	std::string dictstr("");
	
	if(!fProcVec){
		if(!fProcTable){
			fProcTable = G4ProcessTable::GetProcessTable();
			if(!fProcTable) return dictstr; //This is a problem!
		}
		fProcVec = fProcTable->FindProcesses();
		if(!fProcVec) return dictstr;
	}
	
	fNprocs = fProcVec->size();
	
	for(int iProc=0; iProc<fNprocs; iProc++) {
		fOptPhProcessesMap[iProc] = (*fProcVec)(iProc)->GetProcessName();
	}
	
	if(fOptPhProcessesMap.size()!=0){
		//Make the json dict
		json obj;
		
		stringstream ss_tmp; ss_tmp.str("");
		std::map<int, G4String>::iterator it;
		for(it=fOptPhProcessesMap.begin(); it!=fOptPhProcessesMap.end(); it++){
			ss_tmp << it->first;
			//obj[ std::stoi(it->first).c_str() ] = it->second;
			obj[ ss_tmp.str().c_str() ] = it->second;
			ss_tmp.str("");
		}
		
		dictstr = obj.dump();
	}
}


///////////////////////////////////////////////////////

std::string AnalysisManagerOptPh::BuildPysVolDict()
{
	std::string dictstr("");
	
	if(fPhysVolUniqueNamesMap.size()>0){
		//Make the json dict
		json obj;
		
		stringstream ss_tmp; ss_tmp.str("");
		std::map<int, G4String>::iterator it;
		for(it=fPhysVolUniqueNamesMap.begin(); it!=fPhysVolUniqueNamesMap.end(); it++){
			ss_tmp << it->first;
			//obj[ std::stoi(it->first).c_str() ] = it->second;
			obj[ ss_tmp.str().c_str() ] = it->second;
			ss_tmp.str("");
		}
		
		dictstr = obj.dump();
	}
	
	return dictstr;
}


///////////////////////////////////////////////////////

std::string AnalysisManagerOptPh::BuildSDvolDict()
{
	std::string dictstr("");
	
	if(fOptPhSenDetVolPtrs.size()>0){
		//Make the json dict
		json obj;
		
		stringstream ss_tmp; ss_tmp.str("");
		std::set<G4VPhysicalVolume*>::iterator it;
		for(it=fOptPhSenDetVolPtrs.begin(); it!=fOptPhSenDetVolPtrs.end(); it++){
			ss_tmp << fPhysVolUniqueMap[*it];
			//obj[ std::stoi(it->first).c_str() ] = it->second;
			obj[ ss_tmp.str().c_str() ] = fPhysVolUniqueNamesMap[fPhysVolUniqueMap[*it]];
			ss_tmp.str("");
		}
		
		dictstr = obj.dump();
	}
	
	return dictstr;
}


///////////////////////////////////////////////////////

std::string AnalysisManagerOptPh::FindVolGlobalCopy(const G4TouchableHandle& touch)
{
	G4int nLevs = touch->GetHistoryDepth();
	
	vector<Long64_t> copyNumVec(nLevs+1);
	
	for(int iLev=0; iLev<=nLevs; iLev++){
		//pCpNumVec->at(nLevs-iLev) = touch->GetCopyNumber(iLev);
		copyNumVec.at(nLevs-iLev) = touch->GetCopyNumber(iLev);
		
		if((touch->GetCopyNumber(iLev))<0){
			std::cout << " copynumber: " << touch->GetCopyNumber(iLev) << ", volname: " << touch->GetVolume(iLev)->GetName() << std::endl;
		}
		
	}
	
	stringstream ss_tmp; ss_tmp.str("");
	for(unsigned iCp=0; iCp<copyNumVec.size(); iCp++){
		if(iCp==0){
			ss_tmp << copyNumVec.at(iCp);
		}else{
			ss_tmp << "." << copyNumVec.at(iCp);
		}
	}
	
	return ss_tmp.str();
}
