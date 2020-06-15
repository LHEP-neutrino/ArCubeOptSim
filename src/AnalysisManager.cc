#include "DetConstr.hh"
#include "PrimGenAction.hh"
#include "AnalysisManager.hh"
#include "AnalysisMessenger.hh"
//#include "OptPhHit.hh"
#include "EventData.hh"

#include "TNamed.h"

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

#include "nlohmann/json.hpp"


#include <sys/time.h>
#include <numeric>
#include <fstream>
#include <sstream>


using std::vector;
using std::stringstream;
using std::set;
using std::ofstream;

using json = nlohmann::json;

AnalysisManagerOptPh::AnalysisManagerOptPh(PrimaryGeneratorActionOptPh *pPrimaryGeneratorAction)
{
	fMessenger = new AnalysisOptPhMessenger(this);
	
	fRanSeed = 0;
	
	fNav = NULL;
	
	fVerbose=AnalysisManagerOptPh::kSilent;
	fPrintModulo=0;
	
	fCurrentEvent=-1;
	
	fPrimaryGeneratorAction = pPrimaryGeneratorAction;
	
	fSave = AnalysisManagerOptPh::kOff;
	
	fStepsDebug = false;
	fWasAtBoundary = false;
	
	fProcTable = NULL;
	fProcVec = NULL;
	fNprocs = 0;
	
	//fPmtHitsCollectionID = -1;
	
	fDataFilename = "events.root";
	
	fEventData = new EventDataOptPh();
	
	fTreeFile=NULL;
	fTree=NULL;
	
	fAutoSaveEvs = 100;
	fAutoFlushEvs = 100;
	
	fLastTrackId = -1;
	fLastPhysVol = NULL;
}

AnalysisManagerOptPh::~AnalysisManagerOptPh()
{
	if(fEventData) delete fEventData;
	delete fMessenger;
}



void AnalysisManagerOptPh::BeginOfRun(const G4Run *pRun)
{
	//G4cout << "\nEntering in AnalysisManagerOptPh::BeginOfRun(...)" << G4endl;
	G4int randseed;
	
	if(fRanSeed > 0){
		CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
		CLHEP::HepRandom::setTheSeed(fRanSeed);
		randseed = fRanSeed;
	}else{
		// initialize with time.....
		struct timeval hTimeValue;
		gettimeofday(&hTimeValue, NULL);
		CLHEP::HepRandom::setTheSeed(hTimeValue.tv_usec);
		randseed = hTimeValue.tv_usec;
	}
	
	if(fVerbose>=kInfo) std::cout << "Info --> AnalysisManagerOptPh::BeginOfRun(...): Initialize random numbers with seed = " << randseed << std::endl;
	
	
	fProcTable = G4ProcessTable::GetProcessTable();
	
	fCurrentEvent =- 1;
	
	
	MakeVolMaps();
	
	std::string vol_dict = BuildPysVolDict();
	std::string sd_vol_dict = BuildSDvolDict();
	std::string proc_dict = BuildProcsDict();
	
	
	
	
	if(fSave<=AnalysisManagerOptPh::kOff){
		G4cout << "\nWARNING --> AnalysisManagerOptPh::BeginOfRun(...): Data will not be saved." << G4endl;
		return;
	}
	
	fTreeFile = TFile::Open(fDataFilename.c_str(), "RECREATE", "File containing event data of optical photon simulations of ArgonCube");
	
	
	if(!fTreeFile){//Here there is a problem
		fSave = AnalysisManagerOptPh::kOff; //Do this to save from application crashing
		return;
	}
	
	
	
	TNamed *tn_vol_dict = new TNamed("vol_dict", vol_dict.c_str());
	
	TNamed *tn_sdvol_dict = NULL;
	if(fOptPhSenDetVolPtrs.size()>0) tn_sdvol_dict = new TNamed("sdvol_dict", sd_vol_dict.c_str());
	
	TNamed *tn_proc_dict = new TNamed("proc_dict", proc_dict.c_str());
	
	TParameter<int>* ptTParNbEventsToSimulate = new TParameter<int>("EventsNb", fNbEventsToSimulate);
	
	fNbPrim = fPrimaryGeneratorAction->GetPrimNb();
	TParameter<int>* ptTParNbPrimaryPhotPerEvent = new TParameter<int>("PrimNb", fNbPrim);
	
	TParameter<G4int>* ptTParRandSeed = new TParameter<G4int>("RandSeed", randseed);
	
	if(fTreeFile){
		fTreeFile->WriteTObject(tn_vol_dict, 0, "overwrite");
		if(tn_sdvol_dict) fTreeFile->WriteTObject(tn_sdvol_dict, 0, "overwrite");
		fTreeFile->WriteTObject(tn_proc_dict, 0, "overwrite");
		fTreeFile->WriteTObject(ptTParNbEventsToSimulate, 0, "overwrite");
		fTreeFile->WriteTObject(ptTParNbPrimaryPhotPerEvent, 0, "overwrite");
		fTreeFile->WriteTObject(ptTParRandSeed, 0, "overwrite");
	}
	
	
	delete tn_vol_dict;
	if(tn_sdvol_dict) delete tn_sdvol_dict;
	delete tn_proc_dict;
	delete ptTParNbEventsToSimulate;
	delete ptTParNbPrimaryPhotPerEvent;
	delete ptTParRandSeed;
	
	
	
	fTree = new TTree("t1", "Tree containing event data for ArgonCube optical photon simulations.");

	gROOT->ProcessLine("#include <vector>");
	gROOT->ProcessLine("#include <string>");

	fTree->Branch("EvId", &fEventData->fEventId, "eventid/I");
	
	//Primary particle savings
	//if(!fLite) fTree->Branch("type_pri", "vector<string>", &fEventData->m_pPrimaryParticleType);
	
	fTree->Branch("prim_vol_index", &fEventData->fPrimaryVolumeIndex, "prim_vol_idx/I");//Fill at start of Event
	fTree->Branch("prim_vol_cpnm", &fEventData->fPrimaryVolumeCopyNum, "prim_vol_cpnm/I");//Fill at start of Event
	fTree->Branch("prim_vol_id", &fEventData->fPrimaryVolumeId, "prim_vol_id/I");//Fill at start of Event
	fTree->Branch("prim_Xpos", &fEventData->fPrimary_Xpos, "prim_Xpos/D");//Fill at start of Event
	fTree->Branch("prim_Ypos", &fEventData->fPrimary_Ypos, "prim_Ypos/D");//Fill at start of Event
	fTree->Branch("prim_Zpos", &fEventData->fPrimary_Zpos, "prim_Zpos/D");//Fill at start of Event
	fTree->Branch("prim_id", "vector<Int_t>", &fEventData->fPrimary_Id);//Fill at step stage only once
	fTree->Branch("prim_en", "vector<Double_t>", &fEventData->fPrimEn);//Fill at step stage only once
	fTree->Branch("prim_Xmom", "vector<Double_t>", &fEventData->fPrimary_Xmom);//Fill at step stage only once
	fTree->Branch("prim_Ymom", "vector<Double_t>", &fEventData->fPrimary_Ymom);//Fill at step stage only once
	fTree->Branch("prim_Zmom", "vector<Double_t>", &fEventData->fPrimary_Zmom);//Fill at step stage only once
	fTree->Branch("prim_Xpol", "vector<Double_t>", &fEventData->fPrimary_Xpol);//Fill at step stage only once
	fTree->Branch("prim_Ypol", "vector<Double_t>", &fEventData->fPrimary_Ypol);//Fill at step stage only once
	fTree->Branch("prim_Zpol", "vector<Double_t>", &fEventData->fPrimary_Zpol);//Fill at step stage only once
	
	
	
	//Hits related data
	if(fSave < AnalysisManagerOptPh::kSdSteps) fTree->Branch("totalhits", &fEventData->fNbTotHits, "totalhits/L");
	if(fSave < AnalysisManagerOptPh::kSdSteps) fTree->Branch("hit_vol_index", "vector<Int_t>", &fEventData->fVolIndex);//ID of the touchable volume (it is a whish!). //Fill at step stage
	if(fSave < AnalysisManagerOptPh::kSdSteps) fTree->Branch("hit_vol_copy", "vector<Int_t>", &fEventData->fHitVolCopyNum);//This is the copy number of a specific physics volume. //Fill at step stage
	if(fSave < AnalysisManagerOptPh::kSdSteps) fTree->Branch("hit_vol_id", "vector<Int_t>", &fEventData->fHitVolId);//This MUST become the unique ID of the touchable volume. //Fill at step stage
	if(fSave < AnalysisManagerOptPh::kSdSteps) fTree->Branch("hit_time", "vector<Double_t>", &fEventData->fTime);//Fill at step stage
	if(fSave < AnalysisManagerOptPh::kSdSteps) fTree->Branch("hit_firstparentid", "vector<Int_t>", &fEventData->fFirstParentId); //Fill at step stage
	
	//Extended information of hits related data
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_trackid", "vector<Int_t>", &fEventData->fTrackId);//Fill at step stage
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_partgen", "vector<Int_t>", &fEventData->fPartGener);//Fill at step stage
	
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_xpos", "vector<Double_t>", &fEventData->fXpos);//Fill at step stage
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_ypos", "vector<Double_t>", &fEventData->fYpos);//Fill at step stage
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_zpos", "vector<Double_t>", &fEventData->fZpos);//Fill at step stage
	
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_xmom", "vector<Double_t>", &fEventData->fXmom);//Fill at step stage
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_ymom", "vector<Double_t>", &fEventData->fYmom);//Fill at step stage
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_zmom", "vector<Double_t>", &fEventData->fZmom);//Fill at step stage
	
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_xpol", "vector<Double_t>", &fEventData->fXpol);//Fill at step stage
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_ypol", "vector<Double_t>", &fEventData->fYpol);//Fill at step stage
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_zpol", "vector<Double_t>", &fEventData->fZpol);//Fill at step stage
	
	//Full step mode
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("totsteps", &fEventData->fNbTotHits, "totsteps/I");//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("vol_index", "vector<Int_t>", &fEventData->fVolIndex);//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("vol_copy", "vector<Int_t>", &fEventData->fHitVolCopyNum);//ID of the touchable volume//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("vol_id", "vector<Int_t>", &fEventData->fHitVolId);//ID of the touchable volume//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("time", "vector<Double_t>", &fEventData->fTime);//Fill at step stage
	
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("trackid", "vector<Int_t>", &fEventData->fTrackId);//Fill at end of Event
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("partgener", "vector<Int_t>", &fEventData->fPartGener);//Fill at end of Event
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("parentid", "vector<Int_t>", &fEventData->fParentId);//Fill at end of Event
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("firstparentid", "vector<Int_t>", &fEventData->fFirstParentId);//Fill at end of Event
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("creatproc", "vector<Int_t>", &fEventData->fCreatProc);//Fill at end of Event
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("deposproc", "vector<Int_t>", &fEventData->fDepProc);//Fill at step stage
	
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("xpos", "vector<Double_t>", &fEventData->fXpos);//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("ypos", "vector<Double_t>", &fEventData->fYpos);//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("zpos", "vector<Double_t>", &fEventData->fZpos);//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("xmom", "vector<Double_t>", &fEventData->fXmom);//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("ymom", "vector<Double_t>", &fEventData->fYmom);//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("zmom", "vector<Double_t>", &fEventData->fZmom);//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("xpol", "vector<Double_t>", &fEventData->fXpol);//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("ypol", "vector<Double_t>", &fEventData->fYpol);//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("zpol", "vector<Double_t>", &fEventData->fZpol);//Fill at step stage
	
	
	fEventData->fPrimaryVolumeIndex = -1;
	
	//These assignment are for preallocating memory.
	//The vectors will be resized to 0 at start of event
	fEventData->fPrimEn->assign(fNbPrim,0);
	
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
	
	
	
	
	//G4cout << "\nExiting from AnalysisManagerOptPh::BeginOfRun(...)" << G4endl;
}


void AnalysisManagerOptPh::EndOfRun(const G4Run *pRun)
{
	fProcTable = NULL;
	fProcVec = NULL;
	
	if(fSave==AnalysisManagerOptPh::kOff) return;
	
	if(fTreeFile){
		if(fTree){
			fTreeFile->WriteTObject(fTree, 0, "overwrite");
			delete fTreeFile; //This deletes also the TTree owned by the TFile
			fTree=NULL;
		}
		fTreeFile=NULL;
	}
}


void AnalysisManagerOptPh::BeginOfEvent(const G4Event *pEvent)
{
	// grab event ID
	fCurrentEvent = pEvent->GetEventID();
	
	// print this information event by event (modulo n)  	
	if( (fVerbose>=kInfo) && (fPrintModulo>0) && (fCurrentEvent%fPrintModulo == 0) ){
		G4cout << "\nInfo --> AnalysisManagerOptPh::BeginOfEvent(...): Begin of event: " << fCurrentEvent  << G4endl;
	}
	
	
	fWasAtBoundary = false;
	
	//These initialisation are needed at the start of event
	fLastTrackId = -1;
	fLastPhysVol = NULL;
	fLastVolIdx = -1;
	fLastCopyNum = -1;
	fLastVolId = -1;
	
	fEventData->Reset();
	
	fTrackIDs.clear();
	fTrackParentIDsMap.clear();
	fTrackGenerationMap.clear(); 
	fFirstParentIDMap.clear();
	
	
	//Primary particles information
	//Volume id.......
	G4ThreeVector posVec = fPrimaryGeneratorAction->GetPrimPos();
	
	//Here I can still use the tracking navigator as the tracking have not start yet
	if(!fNav) fNav = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
	G4VPhysicalVolume *fPrimVol = fNav->LocateGlobalPointAndSetup(posVec, NULL, true);
	G4TouchableHandle touch = fNav->CreateTouchableHistory();
	
	
	if(fSave>AnalysisManagerOptPh::kOff){
		
		fEventData->fEventId = fCurrentEvent;
		
		fEventData->fPrimaryVolumeIndex = FindVolumeIndex(fPrimVol);
		fEventData->fPrimaryVolumeCopyNum = touch->GetCopyNumber();
		fEventData->fPrimaryVolumeId = FindVolId(touch);
		
		fEventData->fPrimary_Xpos = posVec.x();
		fEventData->fPrimary_Ypos = posVec.y();
		fEventData->fPrimary_Zpos = posVec.z();
		
	}
	
	
	if(fVerbose>=kDebug){
		G4cout << "Debug --> AnalysisManagerOptPh::BeginOfEvent(...): EventID: " << fCurrentEvent << "; Primary volume: " << fPrimVol->GetName() << "; Copy number: " << touch->GetCopyNumber() << G4endl;
	}
	
}


void AnalysisManagerOptPh::EndOfEvent(const G4Event *pEvent)
{
	if(fSave>AnalysisManagerOptPh::kOff){
		if(fTree) fTree->Fill();
	}
}


void AnalysisManagerOptPh::Step(const G4Step *pStep, const G4SteppingManager* pStepMan)
{
	if(fVerbose>=kDebug) G4cout << "\nDebug ---> AnalysisManagerOptPh::Step(...): Entering in AnalysisManagerOptPh::Step\n" << G4endl;
	
	//if( (fVerbose<3) && (fOptPhSenDetVolNames.empty()) )return;
	if(!pStep) return; //This just avoids problems, although would be a big problem to be fixed
	
	G4Track *track = pStep->GetTrack();
	
	G4int trackid = track->GetTrackID();
	
	G4TrackStatus trstatus = track->GetTrackStatus();
	
	G4StepPoint *preStepPoint = pStep->GetPreStepPoint();
	
	G4StepPoint *postStepPoint = pStep->GetPostStepPoint();
	
	G4TouchableHandle touch = postStepPoint->GetTouchableHandle();
	
	G4VPhysicalVolume *Vol = touch->GetVolume();
	
	if(!Vol) return;
	
	
	//Volume printouts
	if(fVerbose>=kDebug || fStepsDebug){
		
		G4VPhysicalVolume *Vol_pre = preStepPoint->GetTouchableHandle()->GetVolume();
		
		G4TouchableHandle touch_pre = preStepPoint->GetTouchableHandle();
		
		
		G4String TrackStat = "";
		
		
		if(trstatus==fAlive) TrackStat = "Alive";
		if(trstatus==fStopButAlive) TrackStat = "StopButAlive";
		if(trstatus==fStopAndKill) TrackStat = "StopAndKill";
		if(trstatus==fKillTrackAndSecondaries) TrackStat = "KillTrackAndSecondaries";
		if(trstatus==fSuspend) TrackStat = "Suspend";
		if(trstatus==fPostponeToNextEvent) TrackStat = "PostponeToNextEvent";
		
		if((postStepPoint->GetStepStatus()==fGeomBoundary) || fWasAtBoundary || (fOptPhSenDetVolPtrs.find(Vol)!=fOptPhSenDetVolPtrs.end())) {
			if(fStepsDebug){
				if(fWasAtBoundary && (postStepPoint->GetStepStatus()!=fGeomBoundary)){
					G4cout << "\nStepDebug --> " << "Event " << fCurrentEvent << ", trackID: " << track->GetTrackID() << ". Optical photon after volumes boundary" << G4endl;
					G4cout << "              Volume 1: <" << Vol_pre->GetName() << ">, copy num: " << touch_pre->GetCopyNumber() << G4endl; 
					G4cout << "              Volume 2: <" << Vol->GetName() << ">, copy num:" << touch->GetCopyNumber() << G4endl;
					G4cout << "              " << "Track status: " << TrackStat << G4endl;
					G4cout << "              " << "Sel proc: " << postStepPoint->GetProcessDefinedStep()->GetProcessName() << "\n" << G4endl;
					fWasAtBoundary = false;
				}else if(!fWasAtBoundary && fOptPhSenDetVolPtrs.find(Vol)==fOptPhSenDetVolPtrs.end()){
					G4cout << "\nStepDebug --> " << "Event " << fCurrentEvent << ", trackID: " << track->GetTrackID() << ". Optical photon at volumes boundary" << G4endl;
					G4cout << "              Volume 1: <" << Vol_pre->GetName() << ">, copy num: " << touch_pre->GetCopyNumber() << G4endl; 
					G4cout << "              Volume 2: <" << Vol->GetName() << ">, copy num:" << touch->GetCopyNumber() << G4endl;
					G4cout << "              " << "Track status: " << TrackStat << G4endl;
					G4cout << "              " << "Sel proc: " << postStepPoint->GetProcessDefinedStep()->GetProcessName() << "\n" << G4endl;
					fWasAtBoundary = true;
				}
        if (fOptPhSenDetVolPtrs.find(Vol)!=fOptPhSenDetVolPtrs.end()){
					G4cout << "\nStepDebug --> " << "Event " << fCurrentEvent << ", trackID: " << track->GetTrackID() << ":" << G4endl;
					G4cout << "              Volume: <" << Vol->GetName() << ">, copy num: " << touch->GetCopyNumber() << G4endl;
					G4cout << "              " << "Track status: " << TrackStat << G4endl;
					G4cout << "              " << "Sel proc: " << postStepPoint->GetProcessDefinedStep()->GetProcessName() << "\n" << G4endl;
					G4cout << "              " << "Energy: " << postStepPoint->GetKineticEnergy() << "\n" << G4endl;
				}
				//std::string foo;
				//G4cout << "Press a enter to continue..."; std::cin >> foo;
			}else{
				
				if( (fSave==kAll) || (fOptPhSenDetVolPtrs.find(Vol)!=fOptPhSenDetVolPtrs.end()) ){
					G4cout << "Debug ---> AnalysisManagerOptPh::Step(...):\n" << "     Event " << fCurrentEvent << ", trackID: " << track->GetTrackID() << ". Optical photon at volumes boundary"<< G4endl;
					G4cout << "     Volume 1: <" << Vol_pre->GetName() << ">, copy num: " << touch_pre->GetCopyNumber() << G4endl; 
					G4cout << "     Volume 2: <" << Vol->GetName() << ">, copy num:" << touch->GetCopyNumber() << G4endl;
					G4cout << "     Track status: " << TrackStat << G4endl;
					G4cout << "     Sel proc: " << postStepPoint->GetProcessDefinedStep()->GetProcessName() << G4endl;
				}
				
			}
		}
	}
	
	
	if( fSave==kOff ){
		return;
	}
	
	//Check if the track is known and if it is a primary track
	//This book-keeping code must be executed whatever the saving options are
	if( (fTrackIDs.insert(trackid)).second ){
		//All the code here inside is executed only once in the same event
		
		int parentid = track->GetParentID();
		fTrackParentIDsMap[trackid] = parentid;
		if(parentid==0){ //It is a primary track at its very first step (in the primary volume)
			fTrackGenerationMap[trackid] = 1;
			fFirstParentIDMap[trackid] = trackid;
			
			if(fSave>kOff){
				//If data will be saved the status of the primary track must be saved
				fEventData->fPrimary_Id->push_back(trackid);
				
				fEventData->fPrimEn->push_back( preStepPoint->GetKineticEnergy() );
				
				fEventData->fPrimary_Xmom->push_back( (preStepPoint->GetMomentumDirection()).x() );
				fEventData->fPrimary_Ymom->push_back( (preStepPoint->GetMomentumDirection()).y() );
				fEventData->fPrimary_Zmom->push_back( (preStepPoint->GetMomentumDirection()).z() );
				
				fEventData->fPrimary_Xpol->push_back( (preStepPoint->GetPolarization()).x() );
				fEventData->fPrimary_Ypol->push_back( (preStepPoint->GetPolarization()).y() );
				fEventData->fPrimary_Zpol->push_back( (preStepPoint->GetPolarization()).z() );
			}
		}else{
			fTrackGenerationMap[trackid] = fTrackGenerationMap[parentid]+1;
			fFirstParentIDMap[trackid] = fFirstParentIDMap[parentid];
		}
		
		
		if( fProcTable && (fSave>=kSdSteps) ){
			if(!track->GetCreatorProcess()){//It's a primary track
				fTrackCreatProc[trackid] = 0;
			}else{
				int retcode = FindProcessIndex(track->GetCreatorProcess()); //This is the index if the process is found (>=0), otherwise a negative return code is returned
				if(retcode>=0){
					fTrackCreatProc[trackid] = retcode+1; //Add 1 to the index as 0 is reserved for primary tracks (no creator process)
				}else{
					fTrackCreatProc[trackid] = retcode; //This is a negative number and indicates for problems
				}
			}
		}
	}
	
	
	//For saves modes lower than kAll check whether the particle is in one of the sensitive volumes defined by the user
	if( (fSave<kAll) ){
		if( fOptPhSenDetVolPtrs.find(Vol)==fOptPhSenDetVolPtrs.end() ){
			//In all other saving modes I am interested only in hits or steps in specific physical volumes (sensitive volumes)
			return;
		}
		
		if( (fSave!=kSdSteps) && (trstatus!=fStopAndKill) ){
			//Here the mode is either "SD stepping mode" or one of the "hits modes"
			//When in hit mode the hit is saved only if the optical photon is going to be absorbed (killed) in one of the sensitive volumes
			return;
		}
	}
	
	//Here start to get stuff to be saved
	fEventData->fNbTotHits += 1;//This is the number of the total recorded steps in "stepping mode" or the number absorption (track stopped and kiled) for an optical photon is in a SD volume when in "hits mode"
	
	
	//This should be replaced with a unique ID of the physical volume
	
	if( (trackid!=fLastTrackId) || (Vol!=fLastPhysVol) ){
		//Recalculate the volume id (recursive process) only if the volume pointer is different from before
		fLastTrackId = trackid;
		fLastPhysVol = Vol;
		fLastVolIdx = fPhysVolUniqueMap[Vol];
		fLastCopyNum = touch->GetCopyNumber();
		fLastVolId = FindVolId(touch);
	}
	
	
	fEventData->fVolIndex->push_back( fLastVolIdx );//Index of the physical volume (from a std::map)
	fEventData->fHitVolCopyNum->push_back( fLastCopyNum );//Copy number of the physical volume
	fEventData->fHitVolId->push_back( fLastVolId );
	fEventData->fFirstParentId->push_back( fFirstParentIDMap[fLastTrackId] );
	fEventData->fTime->push_back( postStepPoint->GetGlobalTime() );
	
	if(fSave>=kHitsExt){
		
		fEventData->fTrackId->push_back( fLastTrackId );
		
		fEventData->fXpos->push_back( (postStepPoint->GetPosition()).x() );
		fEventData->fYpos->push_back( (postStepPoint->GetPosition()).y() );
		fEventData->fZpos->push_back( (postStepPoint->GetPosition()).z() );
		fEventData->fXmom->push_back( (postStepPoint->GetMomentumDirection()).x() );
		fEventData->fYmom->push_back( (postStepPoint->GetMomentumDirection()).y() );
		fEventData->fZmom->push_back( (postStepPoint->GetMomentumDirection()).z() );
		fEventData->fXpol->push_back( (postStepPoint->GetPolarization()).x() );
		fEventData->fYpol->push_back( (postStepPoint->GetPolarization()).y() );
		fEventData->fZpol->push_back( (postStepPoint->GetPolarization()).z() );
		
		if(fSave>=kSdSteps){
			if(fProcTable){
				
				if(!postStepPoint->GetProcessDefinedStep()){
					fEventData->fDepProc->push_back( 0 );//This is a primary track!
				}else{
					int retcode = FindProcessIndex(postStepPoint->GetProcessDefinedStep()); //This is the index if the process is found (>=0), otherwise a negative return code is returned
					if(retcode>=0){
						fEventData->fDepProc->push_back( retcode+1 );//Add 1 to the index as 0 is reserved for primary tracks
					}else{
						fEventData->fDepProc->push_back( retcode );//This is a negative number and indicates for problems
					}
				}
			}
		} // if(fSave>=kSdSteps)...
	} // if(fSave>=kHitsExt)...
	
	if(fVerbose>=kDebug || fStepsDebug) G4cout << "Debug ---> AnalysisManagerOptPh::Step(...): Exiting from AnalysisManagerOptPh::Step(...)\n" << G4endl;
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
		G4bool bMatch = false;
		
		if(bMatch = (hRequiredVolumeName.last('*') != std::string::npos)) hRequiredVolumeName = hRequiredVolumeName.strip(G4String::trailing, '*');
		
		for(G4int iVol=0; iVol<nVols; iVol++){
			G4String hName = pPhysVolStore->at(iVol)->GetName();
			
			if( (hName == hRequiredVolumeName) || (bMatch && (hName.substr(0, hRequiredVolumeName.size())) == hRequiredVolumeName) ){
				fOptPhSenDetVolPtrs.insert(pPhysVolStore->at(iVol));
			}
		}
	}
}


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


#include "nlohmann/json.hpp"


int AnalysisManagerOptPh::FindVolumeIndex( const G4VPhysicalVolume* aVolume )
{
	if(!aVolume) return -1;//This is an error return code
	
	if(fPhysVolUniqueMap.size()==0){
		return -2;
	}
	
	return fPhysVolUniqueMap[(G4VPhysicalVolume*)aVolume];
}


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
	
	for (int iProc = 0; iProc<fNprocs; iProc++) {
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


int AnalysisManagerOptPh::FindVolId(G4TouchableHandle& touch)
{
	G4int nLevs = touch->GetHistoryDepth();
	
	int volId = 1;
	for(int iLev=0; iLev<=nLevs; iLev++){
		volId = volId*(touch->GetVolume(iLev)->GetCopyNo()+1);
	}
	
	return volId;
}

