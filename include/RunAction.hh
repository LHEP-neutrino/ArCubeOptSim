#ifndef __RUN_ACTION_OPT_PH_HH__
#define __RUN_ACTION_OPT_PH_HH__

#include <G4UserRunAction.hh>

class G4Run;

class AnalysisManagerOptPh;

class RunActionOptPh: public G4UserRunAction
{
public:
	RunActionOptPh(AnalysisManagerOptPh *pAnalysisManager=NULL);
	~RunActionOptPh();

public:
	void BeginOfRunAction(const G4Run *pRun);
	void EndOfRunAction(const G4Run *pRun);
	
	

private:
	AnalysisManagerOptPh *fAnalysisManager;
};

#endif

