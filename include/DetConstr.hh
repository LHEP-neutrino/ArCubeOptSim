#ifndef __DET_CONSTR_HH__
#define __DET_CONSTR_HH__

#include "OptPropManager.hh"

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

#include "G4OpticalSurface.hh"
#include "G4GDMLParser.hh"
#include "G4Material.hh"

#include <vector>
#include <map>

using std::vector;
using std::map;
using std::string;


class G4String;
class G4Colour;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4SubtractionSolid;
class DetectorMessenger;
//class OptPropManager;


class DetConstrOptPh: public G4VUserDetectorConstruction
{
public:
	
	typedef std::map<std::string, std::vector<G4VPhysicalVolume*> > PVmap;
	
	enum verbosity{
		kSilent,
		kInfo,
		kDetails,
		kDebug
	};
	
	DetConstrOptPh(G4String gdmlfilename);
	virtual ~DetConstrOptPh();
	
	G4VPhysicalVolume* Construct();
	
	
	const G4GDMLParser* GetGdmlParser() const {return fGDMLParser;};
	const G4VPhysicalVolume* GetWorldVolume() const {return fWorld;};
	
	
	inline void SetVerbosity(DetConstrOptPh::verbosity verb){
		fVerbose=verb;
		fOptPropManager->SetVerbosity( (OptPropManager::verbosity)verb );
	};
	inline DetConstrOptPh::verbosity GetVerbosity(){return fVerbose;};
	
	inline void SetTpbThickness(G4double thick){fTpbThick = thick;};
	inline G4double GetTpbThickness(){return fTpbThick;}
	
	const std::vector<G4VPhysicalVolume* >* GetPvList(G4String pvname) const;
	
	// functions used moslty by the DetectorMessenger
	//void SetCheckOverlap(G4bool dCheckOverlap){pCheckOverlap = dCheckOverlap;};
	
	//static bool GetGeometryParameter(G4String szParameter, G4double& outval);
	
	inline const PVmap* GetVolsMap() const {return (const PVmap*)(&fPVolsMap);}
	
	void PrintVolumeCoordinates(const G4String& VolName);
	void PrintVolumeInfo(const G4String& VolName);
	
	void PrintListOfPhysVols();
	void PrintListOfLogVols();
	
protected:
	//These methods are only used at the startup as default (and to debug)
	//They should go away when the full user interface for optical setting
	//Here also all the optical surfaces are defined
	virtual void BuildTPBlayer();
	virtual void BuildDefaultOptSurf();
	virtual void BuildDefaultLogSurfaces();
	virtual void SetDefaultOptProperties();
	virtual G4Material* FindMaterial(G4String matname);
	
	
	
private:
	
	//Using a std::set in the map in order to avoid having more entries with same volumes pointers (like for a std::vector). This fills the fPVolsMap object.
	void ScanVols(G4VPhysicalVolume* mvol, std::map<G4String, std::set<G4VPhysicalVolume*> > *map=NULL);
	
	DetectorMessenger *fDetectorMessenger;
	
	OptPropManager *fOptPropManager;
	
	G4GDMLParser *fGDMLParser;
	G4VPhysicalVolume *fWorld;
	
	DetConstrOptPh::verbosity fVerbose;
	
	PVmap fPVolsMap;
	std::map<std::string, int > fPVolsRecour; //How many time it is found in the tree
	
	const G4SurfacePropertyTable *fOptSurfTab;
		
	G4double fTpbThick;
};





#endif
