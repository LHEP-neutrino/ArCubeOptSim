#ifndef __TRAC_ACTION__
#define __TRAC_ACTION__

#include "G4UserTrackingAction.hh"
#include "AnalysisManager.hh"

class TrackAct: public G4UserTrackingAction
{
	AnalysisManagerOptPh *fAnalysisManager;
	
public:
	TrackAct(AnalysisManagerOptPh *pAnalysisManager):fAnalysisManager(pAnalysisManager){;};
	
	virtual ~TrackAct(){;};
	
	inline void PreUserTrackingAction(const G4Track* pTrack){
		if(fAnalysisManager) fAnalysisManager->PreUserTrackingAction(pTrack);
	};
	
};


#endif