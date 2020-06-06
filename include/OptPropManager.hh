#ifndef OPT_PROP_MANAGER_HH
#define OPT_PROP_MANAGER_HH


#include "G4OpBoundaryProcess.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "nlohmann/json.hpp"

#include <map>
#include <set>
#include <vector>
#include <fstream>
#include <string>



class G4String;


// for convenience
using json = nlohmann::json;

//Definition of the functins that take a "processing unit" from the main json file structure
class DetConstrOptPh;
class OptPropManager;
typedef void (OptPropManager::*json_proc_memfunc)(const json keyval);


class OptPropManager
{
public:
	enum verbosity{
		kSilent,
		kInfo,
		kDetails,
		kDebug
	};
	
	
	static OptPropManager* GetInstance();
	
	void SetDetConstr(DetConstrOptPh* detConst){fDetConstr = detConst;};
	
	void ProcessJsonFile(const G4String& jsonfilename);
	
	void SetVerbosity(OptPropManager::verbosity verb){fVerbose=verb;};
	OptPropManager::verbosity GetVerbosity(){return fVerbose;};
	
	//loads the refractive index from a 2 column ascii file and gives it to a material if this exists
	//The first column should be the energy photon (not the wavelenght) and the second column the refraction index
	void SetMaterialRindex(const G4String& materialname, const G4String& asciifilename);
	
	//loads the absorption lenght from a 2 column ascii file and gives it to a material if this exists
	//The first column should be the energy photon (not the wavelenght) and the second column the lenght (in m)
	void SetMaterialAbsLenght(const G4String& materialname, const G4String& asciifilename);
	
	//loads the Rayleigh scattering lenght from a 2 column ascii file and gives it to a material if this exists
	//The first column should be the energy photon (not the wavelenght) and the second column the lenght (in m)
	void SetMaterialRayleighLenght(const G4String& materialname, const G4String& asciifilename);

	//loads the photon detection efficiency from a 2 column ascii file and gives it to a material if this exists
	//The first column should be the energy photon (not the wavelenght) and the second column the efficiency
	void SetMaterialEfficiency(const G4String& materialname, const G4String& asciifilename);

	//loads the WLS absorption lenght from a 2 column ascii file and gives it to a material if this exists
	//The first column should be the energy photon (not the wavelenght) and the second column the lenght (in m)
	void SetMaterialWLSAbsLenght(const G4String& materialname, const G4String& asciifilename);
	
	//loads the WLS emission from a 2 column ascii file and gives it to a material if this exists
	//The first column should be the energy photon (not the wavelenght) and the second column the lenght (in m)
	void SetMaterialWLSEmission(const G4String& materialname, const G4String& asciifilename);
	
	//loads the WLS delay from a 2 column ascii file and gives it to a material if this exists
	//The first column should be the energy photon (not the wavelenght) and the second column the lenght (in m)
	void SetMaterialWLSDelay(const G4String& materialname, const G4String& asciifilename);
	
	
	void SetMaterialRindex(const G4String& materialname, const G4int Nentries, const G4double* wavelenghts, const G4double* rindeces);
	void SetMaterialRindex(const G4String& materialname, const std::vector<G4double>& wavelenghts, const std::vector<G4double>& rindeces);
	
	
	void SetMaterialAbsLenght(const G4String& materialname, const G4int Nentries, const G4double* wavelenghts, const G4double* abslenghts);
	void SetMaterialAbsLenght(const G4String& materialname, const std::vector<G4double>& wavelenghts, const std::vector<G4double>& abslenghts);
	
	
	void SetMaterialRayleighLenght(const G4String& materialname, const G4int Nentries, const G4double* wavelenghts, const G4double* rayleighlenghts);
	void SetMaterialRayleighLenght(const G4String& materialname, const std::vector<G4double>& wavelenghts, const std::vector<G4double>& rayleighlenghts);
	
	
	void SetMaterialEfficiency(const G4String& materialname, const G4int Nentries, const G4double* wavelenghts, const G4double* efficiencies);
	void SetMaterialEfficiency(const G4String& materialname, const std::vector<G4double>& wavelenghts, const std::vector<G4double>& efficiencies);
	
	
	void SetMaterialWLSAbsLenght(const G4String& materialname, const G4int Nentries, const G4double* wavelenghts, const G4double* wlsabslenghts);
	void SetMaterialWLSAbsLenght(const G4String& materialname, const std::vector<G4double>& wavelenghts, const std::vector<G4double>& wlsabslenghts);
	
	
	void SetMaterialWLSEmission(const G4String& materialname, const G4int Nentries, const G4double* wavelenghts, const G4double* wlsemissions);
	void SetMaterialWLSEmission(const G4String& materialname, const std::vector<G4double>& wavelenghts, const std::vector<G4double>& wlsemissions);
	
	
	void SetMaterialWLSDelay(const G4String& materialname, const G4double* wlsdelay);
	void SetMaterialWLSDelay(const G4String& materialname, const std::vector<G4double>& wlsdelay);
	
	
	//To use these methods the "logical surface" should already have an optical surface assigned otherwise it doesn't apply the settings
	void SetSurfModel(const G4String& logsurfname, const G4String& model);
	void SetSurfType(const G4String& logsurfname, const G4String& type);
	void SetSurfFinish(const G4String& logsurfname, const G4String& finish);
	void SetSurfSigmaAlpha(const G4String& logsurfname, const G4double& s_a);
	void SetSurfReflectivity(const G4String& logsurfname, const G4int Nentries, const G4double* wavelenghts, const G4double* reflectivities);
	void SetSurfReflectivity(const G4String& logsurfname, const std::vector<G4double>& wavelenghts, const std::vector<G4double>& reflectivities);
	void SetSurfPropFromFile(const G4String& logsurfname, const G4String& filename, const G4String& propertyname);
	
	void BuildOpticalSurface(const G4String& optsurfname, const G4String& model, const G4String& type, const G4String& finish );
	
	//Defaults to glisur, dielectric_dielectric, polished
	void BuildOpticalSurface(const G4String& optsurfname)
	{
		BuildOpticalSurface(optsurfname, "glisur", "dielectric_dielectric", "polished");
	};
	
	
	//Builds the logical surface given the names of the PVs.
	//If no name for the optical surface is given it will build the logical surface without that. --> Do not do it as can generate big issues not immediately recognisable!!!
	int BuildLogicalBorderSurface(const G4String& logsurfname, const G4String& physvol1, const G4String& physvol2, const G4String& optsurfname=G4String(""));
	
	
	void SetOpticalSurface(const G4String& logsurfname, const G4String& optsurfname);
	
	std::set< G4LogicalSurface* >* FindLogSurf(const G4String& logsurfname);
	
	G4OpticalSurface* FindOptSurf(const G4String& optsurfname);
	
	
private:
	OptPropManager();
	~OptPropManager(){;};
	
	
	OptPropManager::verbosity fVerbose;
	
	static OptPropManager* gThis;
	
	std::map<G4String, G4SurfaceType> OptSurfTypeMap;
	std::map<G4String, G4OpticalSurfaceModel> OptSurfModelMap;
	std::map<G4String, G4OpticalSurfaceFinish> OptSurfFinishMap;
	
	std::map<G4String, std::set< G4LogicalSurface*> > fLogSurfMap;
	
	DetConstrOptPh* fDetConstr;
	
	//The following functions process the info of the "setmatprop" key word at the upper level of the jsonfile structure
	
	//Modifies the optical properties of an already existing material
	void setmatprop(const json keyval);
	
	//Modifies the properties of an already existing G4OpticalSurface
	void setoptsurf(const json keyval);
	
	//Modifies the properties of an already existing G4LogicalBorderSurface. If necessary accesses to its G4OpticalSurface to apply the properties
	void setlogbordersurf(const json keyval);
	
	//Modifies the properties of an already existing G4LogicalSkinSurface. If necessary accesses to its G4OpticalSurface to apply the properties
	void setlogskinsurf(const json keyval);
	
	//Builds a new G4OpticalSurface if it doesn't already exist (check by name). If it exists doesn't do anything
	void buildoptsurf(const json keyval);
	
	//Builds a new G4LogicalBorderSurface if it doesn't already exist (check by name). If it exists doesn't do anything. In the process an already instanced G4OpticalSurface must be given to the surface, otherwise this object won't be created.
	void buildlogbordersurf(const json keyval);
	
	//Builds a new G4LogicalSkinSurface if it doesn't already exist (check by name). If it exists doesn't do anything. In the process an already instanced G4OpticalSurface must be given to the surface, otherwise this object won't be created.
	void buildlogskinsurf(const json keyval);
	
	
	//Private versions of functions to build/set logical surfaces, should be only used with homogeneous lists of logical surfaces (same name)
	void SetSurfModel(const std::set<G4LogicalSurface* >* logsurflist, const G4String& model);
	void SetSurfType(const std::set<G4LogicalSurface* >* logsurflist, const G4String& type);
	void SetSurfFinish(const std::set<G4LogicalSurface* >* logsurflist, const G4String& finish);
	void SetSurfSigmaAlpha(const std::set<G4LogicalSurface* >* logsurflist, const G4double& s_a);
	void SetSurfPropFromFile(const std::set<G4LogicalSurface* >* logsurflist, const G4String& filename, const G4String& propertyname);
	void SetSurfConstProp(const std::set<G4LogicalSurface* >* logsurflist, const G4double& value, const G4String& constpropname);
	
	int BuildLogicalBorderSurface(const G4String& logsurfname, const std::vector<G4VPhysicalVolume*>& vol1_vec, const std::vector<G4VPhysicalVolume*>& vol2_vec, const G4OpticalSurface* optsurf );
	
	
	//Map from the key (a string obj) to the corresponding function
	std::map<G4String, json_proc_memfunc> json_proc_tab;
	
	
	//Service functions used by the stuff here above
	const std::vector<G4VPhysicalVolume* >* FindPhysVol(const G4String& physvolname) const;
	
	
	
	
	G4Material* FindMaterial(const G4String& materialname);
	
	void ReadValuesFromFile(const G4String& filename, std::vector<G4double>& ph_en, std::vector<G4double>& vals);
	
	void RegisterLogSurf(G4LogicalSurface* logsurf);
	void DeRegisterLogSurf(G4LogicalSurface* logsurf);
	
};

#endif
