#ifndef __EVENT_DATA_OPT_PH_HH__
#define __EVENT_DATA_OPT_PH_HH__

#include "TROOT.h"
#include "TRint.h"

#include <string>
#include <vector>

using std::string;
using std::vector;


class EventDataOptPh
{
public:
	
	EventDataOptPh();
	virtual ~EventDataOptPh();
	
	
	void Reset();
	
	
	
public:
	//Use root type to be sure to match the datatype of the branch
	
	Int_t fEventId; // Event ID
	
	Int_t fPrimaryVolumeIndex;
	Int_t fPrimaryVolumeCopyNum;
	Int_t fPrimaryVolumeId;
	
	vector<Int_t> *fPrimary_Id; //Id of the primary particle of the chain (used only for the hit mode)
	
	// Position of the primary particle
	Double_t fPrimary_Xpos;
	Double_t fPrimary_Ypos;
	Double_t fPrimary_Zpos;
	
	//Energy of the primary particle
	vector<Double_t> *fPrimEn;
	
	// Momentum of the primary particle
	vector<Double_t> *fPrimary_Xmom;
	vector<Double_t> *fPrimary_Ymom;
	vector<Double_t> *fPrimary_Zmom;
	
	// Polarization of the primary particle
	vector<Double_t> *fPrimary_Xpol;
	vector<Double_t> *fPrimary_Ypol;
	vector<Double_t> *fPrimary_Zpol;
	
	// Number of all the hits for the event
	// Quantities for hits and for making LUTs 
	Int_t fNbTotHits;
	
	vector<Int_t> *fVolIndex;
	vector<Int_t> *fHitVolCopyNum;
	vector<Int_t> *fHitVolId; //Id of the touchable where the step accurred (absorption volume if LUT table are needed)
	
	vector<Double_t> *fTime; //Time of the step
	
	vector<Int_t> *fTrackId; //Id of the track
	vector<Int_t> *fPartGener; //Generation of the particle
	//Extended hit quantities for LUTs: these quantities are general both for "extended info" data level and for more stepping details. The only difference is that in "hit mode" the stuff is saved only when there is an absorption (useful for LUTs making) while in "steping mode" they are saved always, for any kind of process
	
	vector<Double_t> *fXpos, *fYpos, *fZpos; //Coordinates where the step occurred
	vector<Double_t> *fXmom, *fYmom, *fZmom; //Momentum direction at the step
	vector<Double_t> *fXpol, *fYpol, *fZpol; //Polarisation of the particle at the step
	
	
	//Variables present only in the "step mode"
	vector<Int_t> *fParentId; //Id of the parent particle (-1 for primary particles)
	vector<Int_t> *fCreatProc;
	vector<Int_t> *fFirstParentId;
	vector<Int_t> *fDepProc;
};

#endif

