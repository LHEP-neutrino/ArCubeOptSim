#include "AnalysisManager.hh"
#include "RunAction.hh"





RunActionOptPh::RunActionOptPh(AnalysisManagerOptPh *pAnalysisManager)
{
	fAnalysisManager = pAnalysisManager;
}

RunActionOptPh::~RunActionOptPh()
{

}

void RunActionOptPh::BeginOfRunAction(const G4Run *pRun)
{
	if(fAnalysisManager) fAnalysisManager->BeginOfRun(pRun);
}

void RunActionOptPh::EndOfRunAction(const G4Run *pRun)
{
	if(fAnalysisManager) fAnalysisManager->EndOfRun(pRun);
}

