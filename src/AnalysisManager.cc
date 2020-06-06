#include "DetConstr.hh"
#include "PrimGenAction.hh"
#include "AnalysisManager.hh"
#include "AnalysisMessenger.hh"
//#include "OptPhHit.hh"
#include "EventData.hh"

#include "TNamed.h"

#include "G4SDManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4TrackStatus.hh"
#include "G4StepStatus.hh"
#include "G4VProcess.hh"
#include "G4PhysicalVolumeStore.hh"

#include "nlohmann/json.hpp"

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
}

AnalysisManagerOptPh::~AnalysisManagerOptPh()
{
	if(fEventData) delete fEventData;
	delete fMessenger;
}



void AnalysisManagerOptPh::BeginOfRun(const G4Run *pRun)
{
	//G4cout << "\nEntering in AnalysisManagerOptPh::BeginOfRun(...)" << G4endl;
	
	fProcTable = G4ProcessTable::GetProcessTable();
	
	fCurrentEvent =- 1;
	
	
	std::string vol_dict = BuildPysVolDict();
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
	
	TNamed *tn_proc_dict = new TNamed("proc_dict", proc_dict.c_str());
	
	TParameter<int>* ptTParNbEventsToSimulate = new TParameter<int>("EventsNb", fNbEventsToSimulate);
	
	fNbPrim = fPrimaryGeneratorAction->GetPrimNb();
	TParameter<int>* ptParNbPrimaryPhotPerEvent = new TParameter<int>("PrimNb", fNbPrim);
	
	if(fTreeFile) fTreeFile->WriteTObject(tn_vol_dict, 0, "overwrite");
	if(fTreeFile) fTreeFile->WriteTObject(tn_proc_dict, 0, "overwrite");
	if(fTreeFile) fTreeFile->WriteTObject(ptTParNbEventsToSimulate, 0, "overwrite");
	if(fTreeFile) fTreeFile->WriteTObject(ptParNbPrimaryPhotPerEvent, 0, "overwrite");
	
	
	delete tn_vol_dict;
	delete tn_proc_dict;
	delete ptTParNbEventsToSimulate;
	delete ptParNbPrimaryPhotPerEvent;
	
	
	
	fTree = new TTree("t1", "Tree containing event data for ArgonCube optical photon simulations.");

	gROOT->ProcessLine("#include <vector>");
	gROOT->ProcessLine("#include <string>");

	fTree->Branch("EvId", &fEventData->fEventId, "eventid/I");
	
	//Primary particle savings
	//if(!fLite) fTree->Branch("type_pri", "vector<string>", &fEventData->m_pPrimaryParticleType);
	//fTree->Branch("prim_vol_index", "vector<long>", &fEventData->fPrimaryVolumeIndex);
	fTree->Branch("prim_vol_index", &fEventData->fPrimaryVolumeIndex, "prim_vol_index/L");//Fill at start of Event
	fTree->Branch("prim_vol_cpnm", &fEventData->fPrimaryVolumeCopyNum, "prim_vol_cpnm/L");//Fill at start of Event
	fTree->Branch("prim_en", "vector<Double_t>", &fEventData->fPrimEn);//Fill at start of Event
	fTree->Branch("prim_Xpos", &fEventData->fPrimary_Xpos, "prim_Xpos/D");//Fill at start of Event
	fTree->Branch("prim_Ypos", &fEventData->fPrimary_Ypos, "prim_Ypos/D");//Fill at start of Event
	fTree->Branch("prim_Zpos", &fEventData->fPrimary_Zpos, "prim_Zpos/D");//Fill at start of Event
	fTree->Branch("prim_Xmom", "vector<Double_t>", &fEventData->fPrimary_Xmom);//Fill at start of Event
	fTree->Branch("prim_Ymom", "vector<Double_t>", &fEventData->fPrimary_Ymom);//Fill at start of Event
	fTree->Branch("prim_Zmom", "vector<Double_t>", &fEventData->fPrimary_Zmom);//Fill at start of Event
	fTree->Branch("prim_Xpol", "vector<Double_t>", &fEventData->fPrimary_Xpol);//Fill at start of Event
	fTree->Branch("prim_Ypol", "vector<Double_t>", &fEventData->fPrimary_Ypol);//Fill at start of Event
	fTree->Branch("prim_Zpol", "vector<Double_t>", &fEventData->fPrimary_Zpol);//Fill at start of Event
	
	//Hits related data
	if(fSave < AnalysisManagerOptPh::kSdSteps) fTree->Branch("totalhits", &fEventData->fNbTotHits, "totalhits/L");
	if(fSave < AnalysisManagerOptPh::kSdSteps) fTree->Branch("hit_vol_index", "vector<Long64_t>", &fEventData->fVolIndex);//ID of the touchable volume (it is a whish!). //Fill at step stage
	if(fSave < AnalysisManagerOptPh::kSdSteps) fTree->Branch("hit_vol_copy", "vector<Long64_t>", &fEventData->fHitVolId);//This MUST become the unique ID of the touchable volume. //Fill at step stage
	if(fSave < AnalysisManagerOptPh::kSdSteps) fTree->Branch("hit_time", "vector<Double_t>", &fEventData->fTime);//Fill at step stage
	
	//Extended information of hits related data
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_trackid", "vector<Long64_t>", &fEventData->fTrackId);//Fill at step stage
	if(fSave == AnalysisManagerOptPh::kHitsExt) fTree->Branch("hit_partgen", "vector<Long64_t>", &fEventData->fPartGener);//Fill at step stage
	
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
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("vol_index", "vector<Long64_t>", &fEventData->fVolIndex);//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("vol_copy", "vector<Long64_t>", &fEventData->fHitVolId);//ID of the touchable volume//Fill at step stage
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("time", "vector<Double_t>", &fEventData->fTime);//Fill at step stage
	
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("trackid", "vector<Long64_t>", &fEventData->fTrackId);//Fill at end of Event
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("partgener", "vector<Long64_t>", &fEventData->fPartGener);//Fill at end of Event
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("parentid", "vector<Long64_t>", &fEventData->fParentId);//Fill at end of Event
	if(fSave >= AnalysisManagerOptPh::kSdSteps) fTree->Branch("firstparentid", "vector<Long64_t>", &fEventData->fFirstParentId);//Fill at end of Event
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
	
	/*
	if(fHitsCollectionID == -1){
		G4SDManager *pSDManager = G4SDManager::GetSDMpointer();
		fHitsCollectionID = pSDManager->GetCollectionID("OptPhHitsCollection");
	}
	*/
	
	fEventData->Reset(fNbPrim);
	
	//Primary particles information
	//Volume id.......
	G4ThreeVector posVec = fPrimaryGeneratorAction->GetPrimPos();
	
	//Here I can still use the tracking navigator as the tracking have not start yet
	if(!fNav) fNav = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
	G4TouchableHandle touch = fNav->CreateTouchableHistory();
	G4VPhysicalVolume *fPrimVol = fNav->LocateGlobalPointAndSetup(posVec, NULL, true);
	//G4VPhysicalVolume *fPrimVol = fPrimaryGeneratorAction->GetPrimVol();
	
	vector<G4ParticleMomentum> momVec = fPrimaryGeneratorAction->GetPrimMom();
	vector<G4ThreeVector> polVec = fPrimaryGeneratorAction->GetPrimPol();
	
	if(fSave>AnalysisManagerOptPh::kOff){
		
		fEventData->fEventId = fCurrentEvent;
		
		fEventData->fPrimaryVolumeIndex = FindVolumeIndex(fPrimVol);
	
		fEventData->fPrimary_Xpos = posVec.x();
		fEventData->fPrimary_Ypos = posVec.y();
		fEventData->fPrimary_Zpos = posVec.z();
	
		for(unsigned iPrim=0; iPrim<momVec.size(); iPrim++){
			fEventData->fPrimary_Xmom->at(iPrim) = momVec.at(iPrim).x();
			fEventData->fPrimary_Ymom->at(iPrim) = momVec.at(iPrim).y();
			fEventData->fPrimary_Zmom->at(iPrim) = momVec.at(iPrim).z();
			
			fEventData->fPrimary_Xpol->at(iPrim) = polVec.at(iPrim).x();
			fEventData->fPrimary_Ypol->at(iPrim) = polVec.at(iPrim).y();
			fEventData->fPrimary_Zpol->at(iPrim) = polVec.at(iPrim).z();
		}
	}
	
	
	if(fVerbose>=kDebug){
		G4cout << "Debug --> AnalysisManagerOptPh::BeginOfEvent(...): EventID: " << fCurrentEvent << "; Primary volume: " << fPrimVol->GetName() << "; Copy number: " << touch->GetCopyNumber() << G4endl;
	}
	
}


void AnalysisManagerOptPh::EndOfEvent(const G4Event *pEvent)
{
	if(fSave>AnalysisManagerOptPh::kOff){
		std::map<int, int>::iterator iT;
		int nTracks = fTrackParentIDsMap.size();
		fEventData->fPartGener->resize(nTracks);
		fEventData->fFirstParentId->resize(nTracks);
		fEventData->fCreatProc->resize(nTracks);
		
		if(nTracks>0){
			int iTrack=0;
			
			for(iT=fTrackParentIDsMap.begin(); iT!=fTrackParentIDsMap.end(); ++iT){
				int trackid = iT->first;
				int parentid = iT->second;
				if(parentid!=0){
					fFirstParentIDMap[trackid] = fFirstParentIDMap[parentid];
					fTrackGenerationMap[trackid] = fTrackGenerationMap[parentid]+1;
				}
			
				fEventData->fPartGener->at( iTrack ) = fTrackGenerationMap[trackid];
				fEventData->fFirstParentId->at( iTrack ) = fFirstParentIDMap[trackid];
				fEventData->fCreatProc->at( iTrack ) = fTrackCreatProc[trackid];
			
				iTrack++;
			}
		}
		
		//fProcTable = NULL;
	
		if(fTree) fTree->Fill();
	}
}


void AnalysisManagerOptPh::Step(const G4Step *pStep)
{
	if(fVerbose>=kDebug) G4cout << "\nDebug ---> AnalysisManagerOptPh::Step(...): Entering in AnalysisManagerOptPh::Step\n" << G4endl;
	
	//if( (fVerbose<3) && (fOptPhSenDetVolNames.empty()) )return;
	if(!pStep) return; //This just avoids problems, although would be a big problem to be fixed
	
	G4Track *track = pStep->GetTrack();
	
	G4StepPoint *postStepPoint = pStep->GetPostStepPoint();
	
	G4TouchableHandle touch = postStepPoint->GetTouchableHandle();
	
	G4VPhysicalVolume *Vol = touch->GetVolume();
	
	if(!Vol) return;//This would be a big problem.... implement an error message
	
	
	//Volume printouts
	if(fVerbose>=kDebug || fStepsDebug){
		
		
		G4StepPoint *preStepPoint = pStep->GetPreStepPoint();
	
		G4VPhysicalVolume *Vol_pre = preStepPoint->GetTouchableHandle()->GetVolume();
	
		G4TouchableHandle touch_pre = preStepPoint->GetTouchableHandle();
		
		
		G4String TrackStat = "";
	
		G4TrackStatus trstatus = track->GetTrackStatus();
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
	
	
	//Check whether the particle is in one of the sensitive volumes defined by the user
	if( (fSave<kAll) && (fOptPhSenDetVolPtrs.find(Vol)==fOptPhSenDetVolPtrs.end()) ){
		return; //There isn't anything to save here
	}
	
	
	//Here start to get stuff to be saved
	fEventData->fNbTotHits += 1; //This is the number of the total recorded steps in "stepping analisys mode"
	
	//This should be replaced with a unique ID of the physical volume (or the touchable)... how?
	//fEventData->fVolName->push_back( vol->GetName() );
	
	fEventData->fVolIndex->push_back( fOptPhSenDetVolPtrsMap[Vol] );//ID of the physical volume (from a std::map)
	fEventData->fHitVolId->push_back( touch->GetCopyNumber() );//Copy number of the physical volume
	//fEventData->fHitVolId->push_back( touch->GetId() ); //This doesn't exists of course
	
	fEventData->fTime->push_back( postStepPoint->GetGlobalTime() );
	
	if(fSave>=kHitsExt){
		int trackid = track->GetTrackID();
		if( fProcTable && (fTrackParentIDsMap.find(trackid) == fTrackParentIDsMap.end()) ){
			
			//This code is executed only the first time that the specific track is found
			int parentid = track->GetParentID();
			fTrackParentIDsMap[trackid] = parentid;
			
			if(parentid==0){//This is a primary track
				fTrackGenerationMap[trackid] = 1;
				fFirstParentIDMap[trackid] = trackid;
			}
			
			if(fSave>=kSdSteps){
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
		
		fEventData->fTrackId->push_back( trackid );
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
	if(volList == G4String("NULL")){
		fOptPhSenDetVolNames.clear();
		return;
	}
	
	G4PhysicalVolumeStore *pPhysVolStore = G4PhysicalVolumeStore::GetInstance();
	
	G4int nVols = pPhysVolStore->size();
	
	if(nVols <= 0) return;
	
	stringstream hStream;
	hStream.str(volList);
	G4String hVolumeName;
	
	
	// store all the volume names
	set<G4String> candidatevolnames;
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
				fOptPhSenDetVolNames.insert(hName);
				fOptPhSenDetVolPtrs.insert(pPhysVolStore->at(iVol));
			}
		}
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
	
	if(fOptPhSenDetVolPtrsMap.size()==0){
		BuildPysVolDict();
		if(fOptPhSenDetVolPtrsMap.size()==0) return -2;
	}
	
	if(fOptPhSenDetVolPtrsMap.find((G4VPhysicalVolume*)aVolume)==fOptPhSenDetVolPtrsMap.end()){
		return -3;//This should not happen at this stage as the volume is not found in the list of all volumes
	}
	
	return fOptPhSenDetVolPtrsMap[(G4VPhysicalVolume*)aVolume];
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
	fOptPhSenDetVolPtrsMap.clear();
	
	std::string dictstr("");
	
	int iVol = 0;
	if( fSave<AnalysisManagerOptPh::kAll ){
		if(fOptPhSenDetVolPtrs.size()==0) return std::string("");
		
		vector<G4String> volnames;
		
		std::set<G4VPhysicalVolume*>::iterator it;
		for(it=fOptPhSenDetVolPtrs.begin(); it!=fOptPhSenDetVolPtrs.end(); ++it){
			fOptPhSenDetVolPtrsMap[(*it)] = iVol;
			fOptPhPhysVolsMap[iVol] = (*it)->GetName();
			iVol++;
		}
	}else{//Make a dictionary and maps with all the physical volumes
		G4PhysicalVolumeStore *pPhysVolStore = G4PhysicalVolumeStore::GetInstance();
		
		G4int nVols = pPhysVolStore->size();
		
		for(G4int iVol=0; iVol<nVols; iVol++){
			fOptPhSenDetVolPtrsMap[pPhysVolStore->at(iVol)] = iVol;
			fOptPhPhysVolsMap[iVol] = pPhysVolStore->at(iVol)->GetName();
		}
	}
	
	if(fOptPhPhysVolsMap.size()!=0){
		//Make the json dict
		json obj;
		
		stringstream ss_tmp; ss_tmp.str("");
		std::map<int, G4String>::iterator it;
		for(it=fOptPhPhysVolsMap.begin(); it!=fOptPhPhysVolsMap.end(); it++){
			ss_tmp << it->first;
			//obj[ std::stoi(it->first).c_str() ] = it->second;
			obj[ ss_tmp.str().c_str() ] = it->second;
			ss_tmp.str("");
		}
		
		dictstr = obj.dump();
	}
	
	
	
	return dictstr;
}






