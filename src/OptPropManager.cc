#ifndef OPT_PROP_MANAGER_CC
#define OPT_PROP_MANAGER_CC

#include "OptPropManager.hh"
#include "DetConstr.hh"

#include "G4String.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"

#include <numeric>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cassert>




OptPropManager* OptPropManager::gThis = NULL;


//std::map<G4String, G4OpticalSurfaceModel>* OptPropManager::OptSurfModelMap = NULL;

//std::map<G4String, G4SurfaceType>* OptPropManager::OptSurfTypeMap = NULL;

//std::map<G4String, G4OpticalSurfaceFinish>* OptPropManager::OptSurfFinishMap = NULL;



OptPropManager::OptPropManager()
{
	fVerbose = OptPropManager::kSilent;
	
	fDetConstr = NULL;
	
	//Map of the optical surfaces models
	(OptPropManager::OptSurfModelMap)["glisur"] = glisur;
	(OptPropManager::OptSurfModelMap)["unified"] = unified;
	(OptPropManager::OptSurfModelMap)["LUT"] = LUT;
	(OptPropManager::OptSurfModelMap)["DAVIS"] = DAVIS;
	(OptPropManager::OptSurfModelMap)["dichroic"] = dichroic;
	
	//Map of the optical surfaces types
	(OptPropManager::OptSurfTypeMap)["dielectric_metal"] = dielectric_metal;
	(OptPropManager::OptSurfTypeMap)["dielectric_dielectric"] = dielectric_dielectric;
	(OptPropManager::OptSurfTypeMap)["dielectric_LUT"] = dielectric_LUT;
	(OptPropManager::OptSurfTypeMap)["dielectric_LUTDAVIS"] = dielectric_LUTDAVIS;
	(OptPropManager::OptSurfTypeMap)["dielectric_dichroic"] = dielectric_dichroic;
	(OptPropManager::OptSurfTypeMap)["firsov"] = firsov;
	(OptPropManager::OptSurfTypeMap)["x_ray"] = x_ray;
	
	//Map of the optical surfaces finishes
	(OptPropManager::OptSurfFinishMap)["polished"] = polished;
	(OptPropManager::OptSurfFinishMap)["polishedfrontpainted"] = polishedfrontpainted;
	(OptPropManager::OptSurfFinishMap)["polishedbackpainted"] = polishedbackpainted;
	(OptPropManager::OptSurfFinishMap)["ground"] = ground;
	(OptPropManager::OptSurfFinishMap)["groundfrontpainted"] = groundfrontpainted;
	(OptPropManager::OptSurfFinishMap)["groundbackpainted"] = groundbackpainted;
	
	(OptPropManager::OptSurfFinishMap)["polishedlumirrorair"] = polishedlumirrorair;
	(OptPropManager::OptSurfFinishMap)["polishedlumirrorglue"] = polishedlumirrorglue;
	(OptPropManager::OptSurfFinishMap)["polishedair"] = polishedair;
	(OptPropManager::OptSurfFinishMap)["polishedteflonair"] = polishedteflonair;
	(OptPropManager::OptSurfFinishMap)["polishedtioair"] = polishedtioair;
	(OptPropManager::OptSurfFinishMap)["polishedtyvekair"] = polishedtyvekair;
	(OptPropManager::OptSurfFinishMap)["polishedvm2000air"] = polishedvm2000air;
	(OptPropManager::OptSurfFinishMap)["polishedvm2000glue"] = polishedvm2000glue;
	
	(OptPropManager::OptSurfFinishMap)["etchedlumirrorair"] = etchedlumirrorair;
	(OptPropManager::OptSurfFinishMap)["etchedlumirrorglue"] = etchedlumirrorglue;
	(OptPropManager::OptSurfFinishMap)["etchedair"] = etchedair;
	(OptPropManager::OptSurfFinishMap)["etchedteflonair"] = etchedteflonair;
	(OptPropManager::OptSurfFinishMap)["etchedtioair"] = etchedtioair;
	(OptPropManager::OptSurfFinishMap)["etchedtyvekair"] = etchedtyvekair;
	(OptPropManager::OptSurfFinishMap)["etchedvm2000air"] = etchedvm2000air;
	(OptPropManager::OptSurfFinishMap)["etchedvm2000glue"] = etchedvm2000glue;
	
	(OptPropManager::OptSurfFinishMap)["groundlumirrorair"] = groundlumirrorair;
	(OptPropManager::OptSurfFinishMap)["groundlumirrorglue"] = groundlumirrorglue;
	(OptPropManager::OptSurfFinishMap)["groundair"] = groundair;
	(OptPropManager::OptSurfFinishMap)["groundteflonair"] = groundteflonair;
	(OptPropManager::OptSurfFinishMap)["groundtioair"] = groundtioair;
	(OptPropManager::OptSurfFinishMap)["groundtyvekair"] = groundtyvekair;
	(OptPropManager::OptSurfFinishMap)["groundvm2000air"] = groundvm2000air;
	(OptPropManager::OptSurfFinishMap)["groundvm2000glue"] = groundvm2000glue;
	
	(OptPropManager::OptSurfFinishMap)["Rough_LUT"] = Rough_LUT;
	(OptPropManager::OptSurfFinishMap)["RoughTeflon_LUT"] = RoughTeflon_LUT;
	(OptPropManager::OptSurfFinishMap)["RoughESR_LUT"] = RoughESR_LUT;
	(OptPropManager::OptSurfFinishMap)["RoughESRGrease_LUT"] = RoughESRGrease_LUT;
	
	(OptPropManager::OptSurfFinishMap)["Polished_LUT"] = Polished_LUT;
	(OptPropManager::OptSurfFinishMap)["PolishedTeflon_LUT"] = PolishedTeflon_LUT;
	(OptPropManager::OptSurfFinishMap)["PolishedESR_LUT"] = PolishedESR_LUT;
	(OptPropManager::OptSurfFinishMap)["PolishedESRGrease_LUT"] = PolishedESRGrease_LUT;
	
	(OptPropManager::OptSurfFinishMap)["Detector_LUT"] = Detector_LUT;
	
	
	json_proc_tab["setmatprop"] = &OptPropManager::setmatprop;
	json_proc_tab["setoptsurf"] = &OptPropManager::setoptsurf;
	json_proc_tab["setlogbordersurf"] = &OptPropManager::setlogbordersurf;
	json_proc_tab["setlogskinsurf"] = &OptPropManager::setlogskinsurf;
	json_proc_tab["buildoptsurf"] = &OptPropManager::buildoptsurf;
	json_proc_tab["buildlogbordersurf"] = &OptPropManager::buildlogbordersurf;
	json_proc_tab["buildlogskinsurf"] = &OptPropManager::buildlogskinsurf;
	
}



OptPropManager* OptPropManager::GetInstance()
{
	if(!gThis) gThis = new OptPropManager;
	return gThis;
}



void OptPropManager::ProcessJsonFile(const G4String& jsonfilename)
{
  if(fVerbose>=OptPropManager::kDebug){
    std::cout << "Debug --> OptPropManager::ProcessJsonFile(...): processing <" << jsonfilename << ">" << std::endl;
  }
	json jsonObj;
	std::ifstream infile( jsonfilename.c_str() );
	if(!infile) return;
	
	infile >> jsonObj;
	
	int iCmd=0;
	
	if( jsonObj.contains("commands") ){
		
		if(jsonObj.at("commands").is_array()){
			
			std::string cmdname = "";
			json jsonCmds = jsonObj.at("commands");
			
			
			for (json::iterator it = jsonObj.at("commands").begin(); it != jsonObj.at("commands").end(); it++){

				bool validcmd = true;
				cmdname = "";
				
				if(!it.value().is_object()){
					std::cout << "\nERROR --> The json array element num. " << iCmd << " is not of json object type!" << std::endl;
					validcmd = false;
				}
				
				if( validcmd && (!it.value().contains("name")) ){
					std::cout << "\nERROR --> The json array element num. " << iCmd << " does not contain the <name> key!" << std::endl;
					validcmd = false;
				}
				
				
				if(validcmd && (!it.value().at("name").is_string()) ){
					std::cout << "\nERROR --> The json array element num. " << iCmd << " has the <name> value that is not of string type!" << std::endl;
					validcmd = false;
				}
				
				
				if(validcmd){
					cmdname = it.value().at("name").get<std::string>();
					if(fVerbose>=OptPropManager::kDebug){
						std::cout << "Debug --> OptPropManager::ProcessJsonFile(...): Found command key: <" << cmdname << ">" << std::endl;
					}
				}
				
				if(validcmd && (!it.value().contains("obj")) ){
					std::cout << "\nERROR --> The json array element num. " << iCmd << ", with name <" << cmdname << "> has not the <obj> key!" << std::endl;
					validcmd = false;
				}
				
				if(validcmd && (!it.value().at("obj").is_object()) ){
					std::cout << "\nERROR --> The json array element num. " << iCmd << ", with name <" << cmdname << "> has <obj> key not corresponding to a json object value type!" << std::endl;
					validcmd = false;
				}
				
				
				if(validcmd){
					
					if( json_proc_tab.find(cmdname) != json_proc_tab.end() ){
						json_proc_memfunc fncptr = json_proc_tab.at( cmdname );
						if(fVerbose>=OptPropManager::kDebug){
							std::cout << "Debug --> OptPropManager::ProcessJsonFile(...): Processing command: <" << cmdname << ">" << std::endl;
						}
						(this->*fncptr)( it.value().at("obj") );
					}else{
						std::cout << "\nERROR --> OptPropManager::ProcessJsonFile(...): The command <" << cmdname << "> does not correspond to any callback function in the dictionary.\n" << std::endl;
					}
					
				}
				
				iCmd++;
			}
			
		}else{
			std::cout << "\nERROR --> OptPropManager::ProcessJsonFile(...): The key <commands> must be a json array!\n" << std::endl;
		}
	}else{
		std::cout << "\nERROR --> OptPropManager::ProcessJsonFile(...): The key <commands> must be present in the json file. Could not find it!\n" << std::endl;
	}
	
	
	if(fVerbose>=OptPropManager::kDetails){
		std::cout << "Detail --> OptPropManager::ProcessJsonFile(...): Found " << iCmd << " commands to process." << std::endl;
	}
	
}



std::set<G4LogicalSurface*>* OptPropManager::FindLogSurf(const G4String& logsurfname)
{
	std::map<G4String, std::set<G4LogicalSurface*> >::iterator it = fLogSurfMap.find(logsurfname);
	if(it == fLogSurfMap.end()) return  NULL;
	
	return &(it->second);
}



G4OpticalSurface* OptPropManager::FindOptSurf(const G4String& optsurfname)
{
	G4SurfacePropertyTable* surftab = (G4SurfacePropertyTable*)G4OpticalSurface::GetSurfacePropertyTable();
	size_t nSurf = G4OpticalSurface::GetNumberOfSurfaceProperties();
		
	if( (!surftab) || (nSurf<=0) ) return NULL;
	
	for(size_t iSurf=0; iSurf<nSurf; iSurf++){
		if(surftab->at(iSurf)->GetName() == optsurfname){
			return dynamic_cast<G4OpticalSurface*> (surftab->at(iSurf));
		}
	}
	
	return NULL;
}


G4Material* OptPropManager::FindMaterial(const G4String& materialname)
{
	G4MaterialTable *pMatTable = G4Material::GetMaterialTable();
	size_t nMat = G4Material::GetNumberOfMaterials();
	
	if( (!pMatTable) || (nMat<=0) ) return NULL;
	
	for(size_t iMat=0; iMat<nMat; iMat++){
		if( (pMatTable->at(iMat)->GetName()) == materialname ){
			return pMatTable->at(iMat);
		}
	}
	
	return NULL;
}


const std::vector<G4VPhysicalVolume* >* OptPropManager::FindPhysVol(const G4String& physvolname) const
{
	G4PhysicalVolumeStore *pPhysVolStore = G4PhysicalVolumeStore::GetInstance();
	
	if(!pPhysVolStore) return NULL; //This is a big problem as the static method above returns a static class member!!!
	
	if(!fDetConstr) return NULL;
	
	return fDetConstr->GetPvList(physvolname);
}


void OptPropManager::ReadValuesFromFile(const G4String& filename, std::vector<G4double>& ph_en, std::vector<G4double>& vals)
{
	ph_en.resize(0);
	vals.resize(0);
	
	std::ifstream infile(filename.c_str());
	
	if(!infile){
		std::cout << "\nERROR --> OptPropManager::ReadValuesFromFile(...): Cannot find or open in read mode the file <" << filename << ">\n" << std::endl;
		return;
	}
	
	std::string str;
	std::stringstream ss_tmp;
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::ReadValuesFromFile(...): Reading values from file <" << filename << ">:" << std::endl;
	}
	int iLine =0;
	G4double ph_en_d;
	G4double val_d;
	while(getline(infile,str)){
		iLine++;
		if(fVerbose>=OptPropManager::kDebug){
			std::cout << "   Line: " << iLine;
		}
		
		ss_tmp.clear(); ss_tmp.str("");
		ss_tmp << str;
		
		ss_tmp >> str;
		if(ss_tmp){
			ph_en_d = std::stod(str);
			if(fVerbose>=OptPropManager::kDebug) std::cout << "; en: " << ph_en_d;
		}else{
			if(fVerbose>=OptPropManager::kDebug) std::cout << "; corrupted!" << std::endl;
			ph_en.resize(0);
			vals.resize(0);
			return;
		}
		
		ss_tmp >> str;
		if(ss_tmp){//There is only one value while the file format is defined with 2 columns
			val_d = std::stod(str);
			if(fVerbose>=OptPropManager::kDebug) std::cout << "; value: " << val_d << std::endl;
		}else{
			if(fVerbose>=OptPropManager::kDebug) std::cout << "; corrupted!" << std::endl;
			ph_en.resize(0);
			vals.resize(0);
			return;
		}
		
		
		/*
		ph_en_d = std::stod(str);
		
		
		ss_tmp >> str;
		val_d = std::stod(str);
		
		
		if(fVerbose>=OptPropManager::kDebug){
			std::cout << "; en: " << ph_en_d << "; value: " << val_d << std::endl;
		}
		*/
		
		ph_en.push_back(ph_en_d);
		vals.push_back(val_d);
	}
	
	if(fVerbose>=OptPropManager::kDebug) std::cout << std::endl;
}



void OptPropManager::setmatprop(const json keyval)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::setmatprop(...): start of routine." << std::endl;
	}
	
	//Requirements are that the material must exist, otherwise return with an error
	if( !keyval.contains("matname") ){
		std::cout << "\nERROR --> OptPropManager::setmatprop(...): the json object doesn't contain the \"matname\" key. Cannot set material optical properties without knowing the material name!" << std::endl;
		return;
	}
	
	if( !keyval.at("matname").is_string() ){
		std::cout << "\nERROR --> OptPropManager::setmatprop(...): The \"matname\" keymust be a string!" << std::endl;
		return;
	}
	
	
	G4Material* mat = FindMaterial(keyval.at("matname").get<std::string>());
	
	if(!mat){
		std::cout << "\nERROR --> OptPropManager::setmatprop(...): Cannot find the material \""<< keyval.at("matname") <<"\" in table of the instanced materials!" << std::endl;
		return;
	}
	
	//The materials object have only the key "propfile", containing a list of keys that link a property to a specific text file
	
	if( keyval.contains("propfile") && (keyval.at("propfile").is_object()) ){
		json propObj = keyval.at("propfile");
	
		std::vector<G4double> en_vec(0), val_vec(0);
		if(fVerbose>=OptPropManager::kDebug){
			std::cout << "Debug --> OptPropManager::setmatprop(...): starting iterator over properties." << std::endl;
		}
		
		for (json::iterator it = propObj.begin(); it != propObj.end(); ++it){
			if(!it.value().is_string()){
				std::cout << "\nERROR --> OptPropManager::setmatprop(...): The field corresponding to the property " << it.key() << " is not a string!" << std::endl;
				continue;
			}
			
			en_vec.clear(); en_vec.resize(0);
			val_vec.clear(); val_vec.resize(0);
			
			
			if(fVerbose>=OptPropManager::kDetails){
				std::cout << "Detail --> OptPropManager::setmatprop(...): Setting property " << it.key() << " for \"" << mat->GetName() << "\" from file \"" << it.value().get<std::string>() << "\"" << std::endl;
			}
			
			ReadValuesFromFile( it.value().get<std::string>(), en_vec, val_vec );
			
			if( (en_vec.size()==0) || (val_vec.size()==0) || (en_vec.size()!=val_vec.size()) ){
				std::cout << "\nERROR --> OptPropManager::setmatprop(...): Wrong or zero length vectors after reading property file <" << it.value().get<std::string>() << ">. Cannot set material property <" << it.key() << "> for material <" << mat->GetName() << ">\n" << std::endl;
				continue;
			}
			
			
			
			G4MaterialPropertiesTable* propTab = mat->GetMaterialPropertiesTable();
			
			if(!propTab){
				propTab = new G4MaterialPropertiesTable();
				mat->SetMaterialPropertiesTable(propTab);
			}
			
			if(propTab->GetProperty( it.key().c_str() )){
				propTab->RemoveProperty( it.key().c_str() );
			}
			propTab->AddProperty( it.key().c_str() ,(G4double*)&en_vec.at(0), (G4double*)&val_vec.at(0), en_vec.size() );
		}
	}
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::setmatprop(...): end of routine." << std::endl;
	}
}


void OptPropManager::setoptsurf(const json keyval)
{

	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::setoptsurf(...): start of routine." << std::endl;
	}
	
	//Requirements are that the optical surface exist, otherwise return with an error
	if( !keyval.contains("surfname") ){
		std::cout << "\nERROR --> OptPropManager::setoptsurf(...): the json object doesn't contain the <surfname> key. Cannot change settings of an optical surface without knowing its name!" << std::endl;
		return;
	}
	
	if( !keyval.at("surfname").is_string() ){
		std::cout << "\nERROR --> OptPropManager::setoptsurf(...): The json value of the <surfname> key is not of json string type!" << std::endl;
		return;
	}
	
	std::string surfname = keyval.at("surfname").get<std::string>();
	
	G4OpticalSurface* optsurf = FindOptSurf( surfname );
	
	if(!optsurf){
		std::cout << "\nERROR --> OptPropManager::setoptsurf(...): Cannot find the optical surface <"<< surfname <<"> in table of the instanced optical surfaces! You have to build a new one with the json command <buildoptsurf>.\n" << std::endl;
		return;
	}
	
	
	G4MaterialPropertiesTable* propTab = optsurf->GetMaterialPropertiesTable();
	if(!propTab){
		if(fVerbose>=OptPropManager::kDetails){
			std::cout << "Detail --> OptPropManager::setoptsurf(...): Optical surface <" << surfname << "> does not have a properties table (null pointer). Instancing and assigning a new one." << std::endl;
		}
		propTab = new G4MaterialPropertiesTable();
		optsurf->SetMaterialPropertiesTable(propTab);
	}
	
	
	if(keyval.contains("model")){
		if( keyval.at("model").is_string() && (OptSurfModelMap.find(keyval.at("model").get<std::string>())!=OptSurfModelMap.end()) ){
			if(fVerbose>=OptPropManager::kDetails){
				std::cout << "Detail --> OptPropManager::setoptsurf(...): Applying model <" << keyval.at("model").get<std::string>() << "> to the <" << surfname << "> optical surface." << std::endl;
			}
			optsurf->SetModel( OptSurfModelMap.at(keyval.at("model").get<std::string>()) );
		}
	}
	
	if(keyval.contains("type")){
		if( keyval.at("type").is_string() && (OptSurfTypeMap.find(keyval.at("type").get<std::string>())!=OptSurfTypeMap.end()) ){
			if(fVerbose>=OptPropManager::kDetails){
				std::cout << "Detail --> OptPropManager::setoptsurf(...): Applying surface type <" << keyval.at("type").get<std::string>() << "> to the <" << surfname << "> optical surface." << std::endl;
			}
			optsurf->SetType( OptSurfTypeMap.at(keyval.at("type").get<std::string>()) );
		}
	}
	
	if(keyval.contains("finish")){
		if( keyval.at("finish").is_string() && (OptSurfFinishMap.find(keyval.at("finish").get<std::string>())!=OptSurfFinishMap.end()) ){
			if(fVerbose>=OptPropManager::kDetails){
				std::cout << "Detail --> OptPropManager::setoptsurf(...): Applying surface finish <" << keyval.at("finish").get<std::string>() << "> to the <" << surfname << "> optical surface." << std::endl;
			}
			optsurf->SetFinish( OptSurfFinishMap.at(keyval.at("finish").get<std::string>()) );
		}
	}
	
	if(keyval.contains("sigma_alpha")){
		if( keyval.at("sigma_alpha").is_number_float() ){
			if(fVerbose>=OptPropManager::kDetails){
				std::cout << "Detail --> OptPropManager::setoptsurf(...): Applying value <sigma_alpha> = " << keyval.at("sigma_alpha").get<double>() << " to the <" << surfname << "> optical surface." << std::endl;
			}
			optsurf->SetSigmaAlpha( keyval.at("sigma_alpha").get<G4double>() );
		}
	}
	
	
	if( keyval.contains("propfile") && (keyval.at("propfile").is_object()) ){
		
		if(fVerbose>=OptPropManager::kDebug){
			std::cout << "Debug --> OptPropManager::setoptsurf(...): entering into <propfile> section for the <" << surfname << "> optical surface." << std::endl;
		}
		
		json propObj = keyval.at("propfile");
		
		std::vector<G4double> en_vec(0), val_vec(0);
		
		for (json::iterator it = propObj.begin(); it != propObj.end(); ++it){
			if(!it.value().is_string()){
				std::cout << "\nERROR --> OptPropManager::setoptsurf(...): The field corresponding to the property " << it.key() << " is not a json string type!" << std::endl;
				continue;
			}
			
			en_vec.clear(); en_vec.resize(0);
			val_vec.clear(); val_vec.resize(0);
			
			if(fVerbose>=OptPropManager::kDetails){
				std::cout << "Detail --> OptPropManager::buildoptsurf(...): Setting property  <" << it.key() << "> from file <" << it.value().get<std::string>() << "> for the optical surface <" << surfname << ">." << std::endl;
			}
			ReadValuesFromFile( it.value().get<std::string>(), en_vec, val_vec );
			
			if( (en_vec.size()==0) || (val_vec.size()==0) || (en_vec.size()!=val_vec.size()) ){
				std::cout << "\nERROR --> OptPropManager::buildoptsurf(...): Wrong dimensions of vectors after the readout of the file <" << it.value().get<std::string>() << ">. Cannot set the property <" << it.key() << "> for the <" << surfname << "> optical surface!\n" << std::endl;
				continue;
			}
			
			if(propTab->GetProperty( it.key().c_str() )){
				if(fVerbose>=OptPropManager::kDebug){
					std::cout << "Debug --> OptPropManager::buildoptsurf(...): Property <" << it.key() << "> already present in the table for the <" << surfname << "> optical surface. Removing and re-adding to it." << std::endl;
				}
				propTab->RemoveProperty( it.key().c_str() );
			}
			propTab->AddProperty( it.key().c_str() ,(G4double*)&en_vec.at(0), (G4double*)&val_vec.at(0), en_vec.size() );
		}
		
	}
	
	
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::setoptsurf(...): end of routine." << std::endl;
	}
}


void OptPropManager::setlogbordersurf(const json keyval)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::setlogbordersurf(...): start of routine." << std::endl;
	}
	
	
	//Requirements are that the logical border surface exists, otherwise return with an error
	if( !keyval.contains("surfname") ){
		std::cout << "\nERROR --> OptPropManager::setlogbordersurf(...): the json object doesn't contain the <surfname> key. Cannot change settings of a logical border surface without knowing its name!\n" << std::endl;
		return;
	}
	
	if( !keyval.at("surfname").is_string() ){
		std::cout << "\nERROR --> OptPropManager::setlogbordersurf(...): The <surfname> key must be a string!\n" << std::endl;
		return;
	}
	
	G4String logsurfname = keyval.at("surfname").get<std::string>();
	
	std::set<G4LogicalSurface* > *logsurflist = FindLogSurf(logsurfname);
	
	if(!logsurflist){
		std::cout << "\nERROR --> OptPropManager::setlogbordersurf(...): Cannot find the logical surface <"<< logsurfname <<"> in the table of the instanced logical border surfaces!\n" << std::endl;
		return;
	}
	if(fVerbose>=OptPropManager::kDetails){
		std::cout << "Detail --> OptPropManager::setlogbordersurf(...): Found "<< logsurflist->size() << " logical border surfaces named <"<< logsurfname << ">" << std::endl;
	}
	
	G4OpticalSurface *optsurf = NULL;
	std::string optsurfname = "";
	if( keyval.contains("optsurf") ){
		if( keyval.at("optsurf").is_string() ){
			
			optsurfname = keyval.at("optsurf").get<std::string>();
			optsurf = FindOptSurf( optsurfname );
			if(!optsurf->GetMaterialPropertiesTable()){
				optsurf->SetMaterialPropertiesTable( new G4MaterialPropertiesTable() );
			}
			
			if(optsurf){
				if(fVerbose>=OptPropManager::kDetails){
					std::cout << "Detail --> OptPropManager::setlogbordersurf(...): Assigning optical surface <"<< optsurfname << "> to the logical border surfaces named <"<< logsurfname << ">" << std::endl;
				}
				std::set<G4LogicalSurface* >::iterator iT;
				for(iT=logsurflist->begin(); iT!=logsurflist->end(); ++iT){
					(*iT)->SetSurfaceProperty(optsurf);
				}
				
			}else{
				std::cout << "\nWARNING --> OptPropManager::setlogbordersurf(...): Could not find the optical surface <" << keyval.at("optsurf").get<std::string>() << ">. The optical surface of the <" << logsurfname << "> logical surface(s) will not be changed." << std::endl;
			}
			
		}else{
			std::cout << "\nERROR --> OptPropManager::setlogbordersurf(...): The <optsurf> key must be a string!\n" << std::endl;
		}
	}
	
	
	if(keyval.contains("model")){
		if( keyval.at("model").is_string() && (OptSurfModelMap.find(keyval.at("model").get<std::string>())!=OptSurfModelMap.end()) ){
			if(optsurf){
				if(fVerbose>=OptPropManager::kDetails){
					std::cout << "Detail --> OptPropManager::setlogbordersurf(...): Setting logical border surfaces <" << logsurfname << ">: applying model <" << keyval.at("model").get<std::string>() << "> to the optical surface <"<< optsurfname << ">" << std::endl;
				}
				optsurf->SetModel( OptSurfModelMap.at(keyval.at("model").get<std::string>()) );
			}else{
				if(fVerbose>=OptPropManager::kDetails){
					std::cout << "Detail --> OptPropManager::setlogbordersurf(...): Setting logical border surfaces <" << logsurfname << ">: using iterative routine for applying model <" << keyval.at("model").get<std::string>() << ">" << std::endl;
				}
				SetSurfModel(logsurflist, keyval.at("model").get<std::string>() );
			}
		}
	}
	
	if(keyval.contains("type")){
		if( keyval.at("type").is_string() && (OptSurfTypeMap.find(keyval.at("type").get<std::string>())!=OptSurfTypeMap.end()) ){
			if(optsurf){
				if(fVerbose>=OptPropManager::kDetails){
					std::cout << "Detail --> OptPropManager::setlogbordersurf(...): Setting logical border surfaces <" << logsurfname << ">: applying type <" << keyval.at("type").get<std::string>() << "> to the optical surface <"<< optsurfname << ">" << std::endl;
				}
				optsurf->SetType( OptSurfTypeMap.at(keyval.at("type").get<std::string>()) );
			}else{
				if(fVerbose>=OptPropManager::kDetails){
					std::cout << "Detail --> OptPropManager::setlogbordersurf(...): Setting logical border surfaces <" << logsurfname << ">: using iterative routine for applying type <" << keyval.at("type").get<std::string>() << ">" << std::endl;
				}
				SetSurfType( logsurflist, keyval.at("type").get<std::string>() );
			}
		}
	}
	
	if(keyval.contains("finish")){
		if( keyval.at("finish").is_string() && (OptSurfFinishMap.find(keyval.at("finish").get<std::string>())!=OptSurfFinishMap.end()) ){
			
			if(optsurf){
				if(fVerbose>=OptPropManager::kDetails){
					std::cout << "Detail --> OptPropManager::setlogbordersurf(...): Setting logical border surfaces <" << logsurfname << ">: applying finish <" << keyval.at("finish").get<std::string>() << "> to the optical surface <"<< optsurfname << ">" << std::endl;
				}
				optsurf->SetFinish( OptSurfFinishMap.at(keyval.at("finish").get<std::string>()) );
			}else{
				if(fVerbose>=OptPropManager::kDetails){
					std::cout << "Detail --> OptPropManager::setlogbordersurf(...): Setting logical border surfaces <" << logsurfname << ">: using iterative routine for applying finish <" << keyval.at("finish").get<std::string>() << ">" << std::endl;
				}
				SetSurfFinish( logsurflist, keyval.at("finish").get<std::string>() );
			}
		}
	}
	
	if(keyval.contains("sigma_alpha")){
		if( keyval.at("sigma_alpha").is_number_float() ){
			if(optsurf){
				if(fVerbose>=OptPropManager::kDetails){
					std::cout << "Detail --> OptPropManager::setlogbordersurf(...): Setting logical border surfaces <" << logsurfname << ">: applying sigma_alpha value of <" << keyval.at("sigma_alpha").get<double>() << "> to the optical surface <"<< optsurfname << ">" << std::endl;
				}
				optsurf->SetSigmaAlpha( OptSurfFinishMap.at(keyval.at("sigma_alpha").get<double>()) );
			}else{
				if(fVerbose>=OptPropManager::kDetails){
					std::cout << "Detail --> OptPropManager::setlogbordersurf(...): Setting logical border surfaces <" << logsurfname << ">: using iterative routine for applying sigma_alpha value of <" << keyval.at("sigma_alpha").get<double>() << ">" << std::endl;
				}
				SetSurfSigmaAlpha( logsurflist, keyval.at("sigma_alpha").get<double>() );
			}
		}
	}
	
	
	if( keyval.contains("propfile") && (keyval.at("propfile").is_object()) ){
		
		if(fVerbose>=OptPropManager::kDebug){
			std::cout << "Debug --> OptPropManager::setlogbordersurf(...): entering into <propfile> section for the <" << logsurfname << "> logical border surfaces." << std::endl;
		}
		
		json propObj = keyval.at("propfile");
		
		for (json::iterator it = propObj.begin(); it != propObj.end(); ++it){
			if(!it.value().is_string()){
				std::cout << "\nERROR --> OptPropManager::setlogbordersurf(...): The field corresponding to the property <" << it.key() << "> is not a json string type!" << std::endl;
				continue;
			}
			
			std::string propertyname = it.key();
			std::string propertyfilename = it.value().get<std::string>();
			
			if(fVerbose>=OptPropManager::kDetails){
				std::cout << "Detail --> OptPropManager::setlogbordersurf(...): Setting property <" << propertyname << "> for logical surface(s) named <" << logsurfname << "> from file <" << propertyfilename << ">" << std::endl;
			}
			
			if(optsurf){
				std::vector<G4double> en_vec(0), val_vec(0);
				
				ReadValuesFromFile( propertyfilename, en_vec, val_vec );
				
				if( (en_vec.size()==0) || (val_vec.size()==0) || (en_vec.size()!=val_vec.size()) ){
					std::cout << "\nERROR --> OptPropManager::setlogbordersurf(...): Wrong size of the vectors read from file <" << propertyfilename << ">. The property <" << propertyname << "> will not be set for the logical surfaces named <" << logsurfname << ">!\n" << std::endl;
					return;
				}
				
				G4MaterialPropertiesTable* propTab = optsurf->GetMaterialPropertiesTable();
				if(!propTab){
					propTab = new G4MaterialPropertiesTable();
					optsurf->SetMaterialPropertiesTable(propTab);
				}
				
				if(propTab->GetProperty( propertyname.c_str() )){
					if(fVerbose>=OptPropManager::kDebug){
						std::cout << "Debug --> OptPropManager::setlogbordersurf(...): Property <" << propertyname << "> already present for in the list of properties of the <" << optsurf->GetName() << "> optical surface. Removing from and re-adding to it." << std::endl;
					}
					propTab->RemoveProperty( propertyname.c_str() );
				}
				propTab->AddProperty( propertyname.c_str(), (G4double*)&en_vec.at(0), (G4double*)&val_vec.at(0), en_vec.size() );
			}else{
				
				if(fVerbose>=OptPropManager::kDetails){
					std::cout << "Detail --> OptPropManager::setlogbordersurf(...): Setting logical border surfaces <" << logsurfname << ">: using iterative routine for applying property <" << propertyname << "> from file <" << propertyfilename << ">" << std::endl;
				}
				
				SetSurfPropFromFile(logsurflist, propertyfilename, propertyname );
			}
		}
	}
	
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::setlogbordersurf(...): end of routine." << std::endl;
	}
}


void OptPropManager::setlogskinsurf(const json keyval)
{
	//Not implemented yet
}


void OptPropManager::buildoptsurf(const json keyval)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::buildoptsurf(...): start of routine." << std::endl;
	}
	
	
	//Requirements are that the optical surface does not aready exist, otherwise return with an error
	if( !keyval.contains("surfname") ){
		std::cout << "\nERROR --> OptPropManager::buildoptsurf(...): the json object doesn't contain the \"surfname\" key. Cannot change settings of an optical surface without knowing its name!" << std::endl;
		return;
	}
	
	if( !keyval.at("surfname").is_string() ){
		std::cout << "\nERROR --> OptPropManager::buildoptsurf(...): The \"surfname\" key must be a string!" << std::endl;
		return;
	}
	
	std::string surfname = keyval.at("surfname").get<std::string>();
	
	G4OpticalSurface* optsurf = FindOptSurf(surfname);
	
	if(optsurf){
		std::cout << "\nERROR --> OptPropManager::buildoptsurf(...): The optical surface <" << surfname << "> already exists! Use the json command <setoptsurface> for changing the properties of this surface.\n" << std::endl;
		return;
	}
	
	
	if(fVerbose>=OptPropManager::kDetails){
		std::cout << "Detail --> OptPropManager::buildoptsurf(...): Building new optical surface with name <" << surfname << ">" << std::endl;
	}
	optsurf = new G4OpticalSurface( surfname) ;
	optsurf->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());
	G4MaterialPropertiesTable* propTab = optsurf->GetMaterialPropertiesTable();
	
	
	if(keyval.contains("model")){
		if( keyval.at("model").is_string() && (OptSurfModelMap.find(keyval.at("model").get<std::string>())!=OptSurfModelMap.end()) ){
			if(fVerbose>=OptPropManager::kDetails){
				std::cout << "Detail --> OptPropManager::buildoptsurf(...): Applying model <" << keyval.at("model").get<std::string>() << "> to the <" << surfname << "> optical surface." << std::endl;
			}
			optsurf->SetModel( OptSurfModelMap.at(keyval.at("model").get<std::string>()) );
		}
	}
	
	if(keyval.contains("type")){
		if( keyval.at("type").is_string() && (OptSurfTypeMap.find(keyval.at("type").get<std::string>())!=OptSurfTypeMap.end()) ){
			if(fVerbose>=OptPropManager::kDetails){
				std::cout << "Detail --> OptPropManager::buildoptsurf(...): Applying surface type <" << keyval.at("type").get<std::string>() << "> to the <" << surfname << "> optical surface." << std::endl;
			}
			optsurf->SetType( OptSurfTypeMap.at(keyval.at("type").get<std::string>()) );
		}
	}
	
	if(keyval.contains("finish")){
		if( keyval.at("finish").is_string() && (OptSurfFinishMap.find(keyval.at("finish").get<std::string>())!=OptSurfFinishMap.end()) ){
			if(fVerbose>=OptPropManager::kDetails){
				std::cout << "Detail --> OptPropManager::buildoptsurf(...): Applying surface finish <" << keyval.at("finish").get<std::string>() << "> to the <" << surfname << "> optical surface." << std::endl;
			}
			optsurf->SetFinish( OptSurfFinishMap.at(keyval.at("finish").get<std::string>()) );
		}
	}
	
	if(keyval.contains("sigma_alpha")){
		if( keyval.at("sigma_alpha").is_number_float() ){
			if(fVerbose>=OptPropManager::kDetails){
				std::cout << "Detail --> OptPropManager::buildoptsurf(...): Applying value <sigma_alpha> = " << keyval.at("sigma_alpha").get<double>() << " to the <" << surfname << "> optical surface." << std::endl;
			}
			optsurf->SetSigmaAlpha( keyval.at("sigma_alpha").get<double>() );
		}
	}
	
	
	if( keyval.contains("propfile") && (keyval.at("propfile").is_object()) ){
		
		if(fVerbose>=OptPropManager::kDebug){
			std::cout << "Debug --> OptPropManager::buildoptsurf(...): entering into <propfile> section for the <" << surfname << "> optical surface." << std::endl;
		}
		
		json propObj = keyval.at("propfile");
		
		std::vector<G4double> en_vec(0), val_vec(0);
		
		
		for (json::iterator it = propObj.begin(); it != propObj.end(); ++it){
			
			if(!it.value().is_string()){
				std::cout << "\nERROR --> OptPropManager::buildoptsurf(...): The field corresponding to the property " << it.key() << " is not a json string type!\n" << std::endl;
				continue;
			}
			
			en_vec.clear(); en_vec.resize(0);
			val_vec.clear(); val_vec.resize(0);
			
			
			if(fVerbose>=OptPropManager::kDetails){
				std::cout << "Detail --> OptPropManager::buildoptsurf(...): Setting property  <" << it.key() << "> from file <" << it.value().get<std::string>() << "> for the optical surface <" << surfname << ">." << std::endl;
			}
			ReadValuesFromFile( it.value().get<std::string>(), en_vec, val_vec );
			
			if( (en_vec.size()==0) || (val_vec.size()==0) || (en_vec.size()!=val_vec.size()) ){
				std::cout << "\nERROR --> OptPropManager::buildoptsurf(...): Wrong dimensions of vectors after the readout of the file <" << it.value().get<std::string>() << ">. Cannot set the property <" << it.key() << "> for the <" << surfname << "> optical surface!\n" << std::endl;
				continue;
			}
			
			
			if(propTab->GetProperty( it.key().c_str() )){
				if(fVerbose>=OptPropManager::kDebug){
					std::cout << "Debug --> OptPropManager::buildoptsurf(...): Property <" << it.key() << "> already present in the table for the <" << surfname << "> optical surface. Removing and re-adding to it." << std::endl;
				}
				propTab->RemoveProperty( it.key().c_str() );
			}
			propTab->AddProperty( it.key().c_str() ,(G4double*)(&en_vec.at(0)), (G4double*)(&val_vec.at(0)), en_vec.size() );
		}
		
	}
	
	
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::buildoptsurf(...): end of routine." << std::endl;
	}
}


void OptPropManager::buildlogbordersurf(const json keyval)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::buildlogbordersurf(...): start of routine." << std::endl;
	}
	
	
	//Requirements are that the logical border surface does not exist, otherwise return with an error
	if( !keyval.contains("surfname") ){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...): the json object doesn't contain the <surfname> key. Cannot change settings of alogical border surface without knowing its name!\n" << std::endl;
		return;
	}
	
	if( !keyval.at("surfname").is_string() ){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...): The <surfname> key must be a json string type!\n" << std::endl;
		return;
	}
	
	std::string logsurfname = keyval.at("surfname").get<std::string>();
	
	if( !keyval.contains("optsurf") ){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...): The <optsurf> key is mandatory! The logical border surfaces <" << logsurfname << "> will not be built!\n" << std::endl;
		return;
	}
	
	if( !keyval.at("optsurf").is_string() ){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...): The <optsurf> key must be a string!\n" << std::endl;
		return;
	}
	
	std::string optsurfname = keyval.at("optsurf").get<std::string>();
	
	
	G4OpticalSurface *optsurf = FindOptSurf( optsurfname );
	
	if(!optsurf){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...): Could not find the optical surface <" << optsurfname << ">. The optical surface must be built before than the logical surface that uses it. The logical border surfaces <" << logsurfname << "> will not be built.\n" << std::endl;
	}
	
	
	if( !keyval.contains("vol1") ){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...): The <vol1> key is mandatory to build a new logical border surface! The logical border surfaces <" << logsurfname << "> will not be built!\n" << std::endl;
	}
	
	if( !keyval.at("vol1").is_string() ){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...): The <vol1> key must be a json string type! The logical border surfaces <" << logsurfname << "> will not be built!\n" << std::endl;
		return;
	}
	
	
	if( !keyval.contains("vol2") ){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...): The <vol2> key is mandatory to build a new logical border surface! The logical border surfaces <" << logsurfname << "> will not be built!\n" << std::endl;
		return;
	}
	
	if( !keyval.at("vol2").is_string() ){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...): The <vol2> key must be a json string type! The logical border surfaces <" << logsurfname << "> will not be built!\n" << std::endl;
		return;
	}
	
	
	if( keyval.at("vol1").get<std::string>() == keyval.at("vol2").get<std::string>() ){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...): The <vol1> and <vol2> are the same. A logical border surface can be built only between 2 different physical volumes!\n" << std::endl;
		return;
	}
	
	
	const std::vector<G4VPhysicalVolume* > *vol_vec1 = FindPhysVol(keyval.at("vol1").get<std::string>());
	
	if(!vol_vec1){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...):Could not find the <vol1> physical volume with name <" << keyval.at("vol1").get<std::string>() << ">. The logical border surfaces <" << logsurfname << "> will not be built!\n" << std::endl;
		return;
	}
	
	if( (!vol_vec1->size()) ){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...):The <vol1> physical volume with name <" << keyval.at("vol1").get<std::string>() << "> corresponds to an empty list! This is an unexpected behaviour. The logical border surfaces <" << logsurfname << "> will not be built.\n" << std::endl;
		return;
	}
	
	
	const std::vector<G4VPhysicalVolume* > *vol_vec2 = FindPhysVol(keyval.at("vol2").get<std::string>());
	
	if(!vol_vec2){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...):Could not find the <vol2> physical volume with name <" << keyval.at("vol2").get<std::string>() << ">. The logical border surfaces <" << logsurfname << "> will not be built!\n" << std::endl;
		return;
	}
	
	if( !(vol_vec2->size()) ){
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...): The <vol2> physical volume with name <" << keyval.at("vol2").get<std::string>() << "> corresponds to an empty list! This is an unexpected behaviour. The logical border surfaces <" << logsurfname << "> will not be built.\n" << std::endl;
		return;
	}
	
	
	
	int nLBS = BuildLogicalBorderSurface(logsurfname, (*vol_vec1), (*vol_vec2), optsurf);
	if(nLBS<=0){
		
		std::cout << "\nERROR --> OptPropManager::buildlogbordersurf(...): The logical border surfaces <" << keyval.at("surfname").get<std::string>() << "> could not be built!\n" << std::endl;
		
	}else{
		if(fVerbose>=OptPropManager::kDetails){
			std::cout << "Detail --> OptPropManager::buildlogbordersurf(...): built " << nLBS << " logical border surfaces named <" << logsurfname << "> between vol1 <" << keyval.at("vol1").get<std::string>() << "> and vol2 <" << keyval.at("vol2").get<std::string>() << "> by using the optical surface <" << optsurfname << ">." << std::endl;
		}
	}
	
	
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::buildlogbordersurf(...): end of routine." << std::endl;
	}
}


void OptPropManager::buildlogskinsurf(const json keyval)
{
	//Not implemented yet
}


void OptPropManager::SetSurfReflectivity(const G4String& logsurfname, const G4int Nentries, const G4double* photonenergies, const G4double* reflectivities)
{
	const G4LogicalBorderSurfaceTable* surftab = G4LogicalBorderSurface::GetSurfaceTable();
	
	if(surftab){
		G4LogicalSurface* Surface = NULL;
		
		for(size_t iSurf=0; iSurf<surftab->size(); iSurf++){
			G4String name = surftab->at(iSurf)->GetName();
			if(name == logsurfname){
				Surface = surftab->at(iSurf);
			}
		}
		
		if(Surface){
			G4OpticalSurface* OpticalSurface = dynamic_cast <G4OpticalSurface*> (Surface->GetSurfaceProperty());
			if(OpticalSurface){
				G4MaterialPropertiesTable* propTab = OpticalSurface->GetMaterialPropertiesTable();
				
				if(!propTab){
					propTab = new G4MaterialPropertiesTable();
					OpticalSurface->SetMaterialPropertiesTable(propTab);
				}
				
				if(propTab->GetProperty("REFLECTIVITY")){
					propTab->RemoveProperty("REFLECTIVITY");
				}
				
				propTab->AddProperty("REFLECTIVITY",(G4double*)photonenergies,(G4double*)reflectivities,Nentries);
			}else{
				std::cout << "WARNING --> OptPropManager::SetSurfReflectivity(...): The G4LogicalBorder surface \"" << logsurfname << "\" returned null pointer to its G4OpticalSurface (maybe not yet set). The surface roughness cannot be set." << std::endl;
			}
		}else{
			std::cout << "WARNING --> OptPropManager::SetSurfReflectivity(...): The G4LogicalBorder surface \"" << logsurfname << "\" has not yet been created. The surface roughness cannot be set." << std::endl;
		}
	}
}

void OptPropManager::SetSurfReflectivity(const G4String& logsurfname, const std::vector<G4double>& photonenergies, const std::vector<G4double>& reflectivities)
{
  if(photonenergies.size()==reflectivities.size()) OptPropManager::SetSurfReflectivity(logsurfname, photonenergies.size(), &photonenergies.at(0), &reflectivities.at(0));
}



void OptPropManager::SetMaterialRindex(const G4String& materialname, const G4int Nentries, const G4double* photonenergies, const G4double* rindexes)
{
	G4MaterialTable *pMatTable = G4Material::GetMaterialTable();
	size_t nMat = G4Material::GetNumberOfMaterials();
	
	if(pMatTable && (nMat>0)){
		for(size_t iMat=0; iMat<nMat;iMat++){
			G4Material *mat = pMatTable->at(iMat);
			
			if( mat->GetName() == materialname ){
				
				G4MaterialPropertiesTable* propTab = mat->GetMaterialPropertiesTable();
				
				
				if(!propTab){
					propTab = new G4MaterialPropertiesTable();
					mat->SetMaterialPropertiesTable(propTab);
				}
				
				
				if(propTab->GetProperty("RINDEX")){
					propTab->RemoveProperty("RINDEX");
				}
				propTab->AddProperty("RINDEX",(G4double*)photonenergies,(G4double*)rindexes,Nentries);
			}
		}
	}
}

void OptPropManager::SetMaterialRindex(const G4String& materialname, const std::vector<G4double>& photonenergies, const std::vector<G4double>& rindexes)
{
	if(photonenergies.size()==rindexes.size()) OptPropManager::SetMaterialRindex(materialname, photonenergies.size(), &photonenergies.at(0), &rindexes.at(0));
}



void OptPropManager::SetMaterialAbsLenght(const G4String& materialname, const G4int Nentries, const G4double* photonenergies, const G4double* abslenghts)
{
	G4MaterialTable *pMatTable = G4Material::GetMaterialTable();
	size_t nMat = G4Material::GetNumberOfMaterials();
	
	if(pMatTable && (nMat>0)){
		for(size_t iMat=0; iMat<nMat;iMat++){
			G4Material *mat = pMatTable->at(iMat);
			if( mat->GetName() == materialname ){
				
				G4MaterialPropertiesTable* propTab = mat->GetMaterialPropertiesTable();
				
				if(!propTab){
					propTab = new G4MaterialPropertiesTable();
					mat->SetMaterialPropertiesTable(propTab);
				}
				
				if(propTab->GetProperty("ABSLENGTH")){
					propTab->RemoveProperty("ABSLENGTH");
				}
				propTab->AddProperty("ABSLENGTH",(G4double*)photonenergies,(G4double*)abslenghts,Nentries);
			}
		}
	}
}

void OptPropManager::SetMaterialAbsLenght(const G4String& materialname, const std::vector<G4double>& photonenergies, const std::vector<G4double>& abslenghts)
{
	if(photonenergies.size()==abslenghts.size()) OptPropManager::SetMaterialAbsLenght(materialname, photonenergies.size(), &photonenergies.at(0), &abslenghts.at(0));
}



void OptPropManager::SetMaterialRayleighLenght(const G4String& materialname, const G4int Nentries, const G4double* photonenergies, const G4double* rayleighlenghts)
{
	G4MaterialTable *pMatTable = G4Material::GetMaterialTable();
	size_t nMat = G4Material::GetNumberOfMaterials();
	
	if(pMatTable && (nMat>0)){
		for(size_t iMat=0; iMat<nMat;iMat++){
			G4Material *mat = pMatTable->at(iMat);
			
			if( mat->GetName() == materialname ){
				
				G4MaterialPropertiesTable* propTab = mat->GetMaterialPropertiesTable();
				
				if(!propTab){
					propTab = new G4MaterialPropertiesTable();
					mat->SetMaterialPropertiesTable(propTab);
				}
				
				if(propTab->GetProperty("RAYLEIGH")){
					propTab->RemoveProperty("RAYLEIGH");
				}
				propTab->AddProperty("RAYLEIGH",(G4double*)photonenergies,(G4double*)rayleighlenghts,Nentries);
			}
		}
	}
}

void OptPropManager::SetMaterialRayleighLenght(const G4String& materialname, const std::vector<G4double>& photonenergies, const std::vector<G4double>& rayleighlenghts)
{
	if(photonenergies.size()==rayleighlenghts.size()) OptPropManager::SetMaterialRayleighLenght(materialname, photonenergies.size(), &photonenergies.at(0), &rayleighlenghts.at(0));
}



void OptPropManager::SetMaterialEfficiency(const G4String& materialname, const G4int Nentries, const G4double* photonenergies, const G4double* efficiencies)
{
	G4MaterialTable *pMatTable = G4Material::GetMaterialTable();
	size_t nMat = G4Material::GetNumberOfMaterials();
	
	if(pMatTable && (nMat>0)){
		for(size_t iMat=0; iMat<nMat;iMat++){
			G4Material *mat = pMatTable->at(iMat);
			if( mat->GetName() == materialname ){
				
				G4MaterialPropertiesTable* propTab = mat->GetMaterialPropertiesTable();
				
				if(!propTab){
					propTab = new G4MaterialPropertiesTable();
					mat->SetMaterialPropertiesTable(propTab);
				}
				
				if(propTab->GetProperty("EFFICIENCY")){
					propTab->RemoveProperty("EFFICIENCY");
				}
				propTab->AddProperty("EFFICIENCY",(G4double*)photonenergies,(G4double*)efficiencies,Nentries);
			}
		}
	}
}

void OptPropManager::SetMaterialEfficiency(const G4String& materialname, const std::vector<G4double>& photonenergies, const std::vector<G4double>& efficiencies)
{
	if(photonenergies.size()==efficiencies.size()) OptPropManager::SetMaterialEfficiency(materialname, photonenergies.size(), &photonenergies.at(0), &efficiencies.at(0));
}



void OptPropManager::SetMaterialWLSAbsLenght(const G4String& materialname, const G4int Nentries, const G4double* photonenergies, const G4double* wlsabslenghts)
{
	G4MaterialTable *pMatTable = G4Material::GetMaterialTable();
	size_t nMat = G4Material::GetNumberOfMaterials();
	
	if(pMatTable && (nMat>0)){
		for(size_t iMat=0; iMat<nMat;iMat++){
			G4Material *mat = pMatTable->at(iMat);
			if( mat->GetName() == materialname ){
				
				G4MaterialPropertiesTable* propTab = mat->GetMaterialPropertiesTable();
				
				if(!propTab){
					propTab = new G4MaterialPropertiesTable();
					mat->SetMaterialPropertiesTable(propTab);
				}
				
				if(propTab->GetProperty("WLSABSLENGTH")){
					propTab->RemoveProperty("WLSABSLENGTH");
				}
				propTab->AddProperty("WLSABSLENGTH",(G4double*)photonenergies,(G4double*)wlsabslenghts,Nentries);
			}
		}
	}
}

void OptPropManager::SetMaterialWLSAbsLenght(const G4String& materialname, const std::vector<G4double>& photonenergies, const std::vector<G4double>& wlsabslenghts)
{
	if(photonenergies.size()==wlsabslenghts.size()) OptPropManager::SetMaterialWLSAbsLenght(materialname, photonenergies.size(), &photonenergies.at(0), &wlsabslenghts.at(0));
}



void OptPropManager::SetMaterialWLSEmission(const G4String& materialname, const G4int Nentries, const G4double* photonenergies, const G4double* wlsemissions)
{
	G4MaterialTable *pMatTable = G4Material::GetMaterialTable();
	size_t nMat = G4Material::GetNumberOfMaterials();
	
	if(pMatTable && (nMat>0)){
		for(size_t iMat=0; iMat<nMat;iMat++){
			G4Material *mat = pMatTable->at(iMat);
			if( mat->GetName() == materialname ){
				
				G4MaterialPropertiesTable* propTab = mat->GetMaterialPropertiesTable();
				
				if(!propTab){
					propTab = new G4MaterialPropertiesTable();
					mat->SetMaterialPropertiesTable(propTab);
				}
				
				if(propTab->GetProperty("WLSCOMPONENT")){
					propTab->RemoveProperty("WLSCOMPONENT");
				}
				propTab->AddProperty("WLSCOMPONENT",(G4double*)photonenergies,(G4double*)wlsemissions,Nentries);
			}
		}
	}
}

void OptPropManager::SetMaterialWLSEmission(const G4String& materialname, const std::vector<G4double>& photonenergies, const std::vector<G4double>& wlsemissions)
{
	if(photonenergies.size()==wlsemissions.size()) OptPropManager::SetMaterialWLSEmission(materialname, photonenergies.size(), &photonenergies.at(0), &wlsemissions.at(0));
}



void OptPropManager::SetMaterialWLSDelay(const G4String& materialname, const G4double* delay)
{
	G4MaterialTable *pMatTable = G4Material::GetMaterialTable();
	size_t nMat = G4Material::GetNumberOfMaterials();
	
	if(pMatTable && (nMat>0)){
		for(size_t iMat=0; iMat<nMat;iMat++){
			G4Material *mat = pMatTable->at(iMat);
			if( mat->GetName() == materialname ){
				
				G4MaterialPropertiesTable* propTab = mat->GetMaterialPropertiesTable();
				
				if(!propTab){
					propTab = new G4MaterialPropertiesTable();
					mat->SetMaterialPropertiesTable(propTab);
				}
				
				if(propTab->GetProperty("WLSTIMECONSTANT")){
					propTab->RemoveProperty("WLSTIMECONSTANT");
				}
				propTab->AddConstProperty("WLSTIMECONSTANT",(G4double)*delay);
			}
		}
	}
}



void OptPropManager::SetMaterialWLSDelay(const G4String& materialname, const std::vector<G4double>& delay)
{
	OptPropManager::SetMaterialWLSDelay(materialname, &delay.at(0));
}



void OptPropManager::SetSurfModel(const G4String& logsurfname, const G4String& model)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfModel(const G4String& logsurfname, ...): Entering the function." << std::endl;
	}
	
	std::set<G4LogicalSurface* > *logsurflist = FindLogSurf(logsurfname);
	
	if(!logsurflist){
		std::cout << "\nERROR --> OptPropManager::SetSurfModel(...): Cannot find the logical surface(s) named <" << logsurfname << "> in the table of the instanced logical surfaces!\n" << std::endl;
		return;
	}
	
	if(!(logsurflist->size())){
		std::cout << "\nERROR --> OptPropManager::SetSurfModel(...): The list of logical surface(s) named <" << logsurfname << "> is empty!\n" << std::endl;
		return;
	}
	
	
	if(OptSurfModelMap.find(model)==OptSurfModelMap.end()){
		std::cout << "\nERROR --> OptPropManager::SetSurfModel(...): The model <" << model << "> is not in the list of known models and can not be applied to logical surface(s) <" << logsurfname << ">.\n" << std::endl;
		return;
	}
	
	SetSurfModel(logsurflist, model);
	
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfModel(const G4String& logsurfname, ...): Exiting the function." << std::endl;
	}
}

void OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >* logsurflist, const G4String& model)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >*, ...): Entering the function." << std::endl;
	}
	
	if(!logsurflist){
		std::cout << "\nERROR --> OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >* , ...): Null pointer to the list of logical surfaces!\n" << std::endl;
		return;
	}
	
	if(!(logsurflist->size())){
		std::cout << "\nERROR --> OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >* , ...): The list of logical surfaces is empty!\n" << std::endl;
		return;
	}
	
	
	if(OptSurfModelMap.find(model)==OptSurfModelMap.end()){
		std::cout << "\nERROR --> OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >* , ...): The model <" << model << "> is not in the list of known models and can not be applied to the list of logical surfaces in input.\n" << std::endl;
		return;
	}
	
	
	G4OpticalSurface *lastoptsurf = NULL;
	std::string firstsurfname, lastsurfname;
	std::set<G4LogicalSurface* >::iterator iT;
	int nSurfs = logsurflist->size();
	int nAppl = 0;
	
	int iSurf = 0;
	for(iT=logsurflist->begin(); iT!=logsurflist->end(); ++iT, iSurf++){
		if(iSurf==0){
			firstsurfname = (*iT)->GetName();
			lastsurfname = firstsurfname;
			lastoptsurf = dynamic_cast<G4OpticalSurface*>((*iT)->GetSurfaceProperty());
			if(lastoptsurf){
				if(fVerbose>=OptPropManager::kDebug){
					std::cout << "Debug --> OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >* , ...): Applying model <" << model << "> to logical surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ")." << std::endl;
				}
				lastoptsurf->SetModel( OptSurfModelMap.at(model) );
				nAppl++;
			}else{
				std::cout << "\nERROR --> OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") has no optical surface assigned. Cannot assign model to this logical surface instance. Add to it an optical surface before!\n" << std::endl;
			}
		}else{
			if(!(*iT)->GetSurfaceProperty()){
				std::cout << "\nERROR --> OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") has no optical surface assigned. Cannot assign model to this logical surface instance. Add to it an optical surface  before!\n" << std::endl;
			}else{
				std::string logsurfname = (*iT)->GetName();
				if(logsurfname != lastsurfname){
					std::cout << "WARNING --> OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >* , ...): Change of name from <" << lastsurfname << "> to <" << logsurfname << "> for logical surface (list element " << iSurf << " at " << (*iT) << "). This function is supposed to apply the model to a set of logical surfaces with the same name. Unexpected behaviours of the simulation might occur!" << std::endl;
					lastsurfname = logsurfname;
				}
				if( lastoptsurf != ((*iT)->GetSurfaceProperty()) ){
					std::cout << "WARNING --> OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << lastsurfname << "> (list element " << iSurf << " at " << (*iT) << ") has an assigned optical surface different (fro pointer) to the one processed before. Applying anyway the model <" << model << "> to this optical surface!" << std::endl;
					lastoptsurf = dynamic_cast<G4OpticalSurface*>((*iT)->GetSurfaceProperty());
					lastoptsurf->SetModel( OptSurfModelMap.at(model) );
				}else{
					if(fVerbose>=OptPropManager::kDebug){
						std::cout << "Debug --> OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >* , ...): Surface model <" << model << "> to logical surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") was already applied at the previous iteration. NO NEED to actually apply it for this logical surface instance." << std::endl;
					}
				}
				nAppl++;
			}
		}
	}
	
	if(fVerbose>=OptPropManager::kDetails){
		std::cout << "Details --> OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >*, ...): Applied model <" << model << "> to " << nAppl << " out of " << nSurfs << " logical surfaces named <" << firstsurfname << "> (if list was homogeneous)." << std::endl;
	}
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfModel(const std::set<G4LogicalSurface* >*, ...): Exiting the function." << std::endl;
	}
}



void OptPropManager::SetSurfType(const G4String& logsurfname, const G4String& type)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfType(const G4String&, ...): Entering the function." << std::endl;
	}
	
	std::set<G4LogicalSurface* > *logsurflist = FindLogSurf(logsurfname);
	
	if(!logsurflist){
		std::cout << "\nERROR --> OptPropManager::SetSurfType(...): Cannot find the logical surface(s) named <" << logsurfname << "> in the table of the instanced logical surfaces!\n" << std::endl;
		return;
	}
	
	if(!(logsurflist->size())){
		std::cout << "\nERROR --> OptPropManager::SetSurfType(...): The list of logical surface(s) named <" << logsurfname << "> is empty!\n" << std::endl;
		return;
	}
	
	
	if(OptSurfTypeMap.find(type)==OptSurfTypeMap.end()){
		std::cout << "\nERROR --> OptPropManager::SetSurfType(...): The type <" << type << "> is not in the list of known surface types and can not be applied to the logical surface(s) <" << logsurfname << ">.\n" << std::endl;
		return;
	}
	
	SetSurfType(logsurflist, type);
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfType(const G4String&, ...): Exiting the function." << std::endl;
	}
}

void OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >* logsurflist, const G4String& type)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >*, ...): Entering the function." << std::endl;
	}
	
	if(!logsurflist){
		std::cout << "\nERROR --> OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >* , ...): Null pointer to the list of logical surfaces!\n" << std::endl;
		return;
	}
	
	if(!(logsurflist->size())){
		std::cout << "\nERROR --> OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >* , ...): The list of logical surfaces is empty!\n" << std::endl;
		return;
	}
	
	
	if(OptSurfTypeMap.find(type)==OptSurfTypeMap.end()){
		std::cout << "\nERROR --> OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >* , ...): The type <" << type << "> is not in the list of known types and can not be applied to the logical surfaces.\n" << std::endl;
		return;
	}
	
	
	G4OpticalSurface *lastoptsurf = NULL;
	std::string firstsurfname, lastsurfname;
	std::set<G4LogicalSurface* >::iterator iT;
	int nSurfs = logsurflist->size();
	int nAppl = 0;
	
	int iSurf = 0;
	for(iT=logsurflist->begin(); iT!=logsurflist->end(); ++iT, iSurf++){
		if(iSurf==0){
			firstsurfname = (*iT)->GetName();
			lastsurfname = firstsurfname;
			lastoptsurf = dynamic_cast<G4OpticalSurface*>((*iT)->GetSurfaceProperty());
			if(lastoptsurf){
				if(fVerbose>=OptPropManager::kDebug){
					std::cout << "Debug --> OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >* , ...): Applying type <" << type << "> to logical surfaces <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ")." << std::endl;
				}
				lastoptsurf->SetType( OptSurfTypeMap.at(type) );
				nAppl++;
			}else{
				std::cout << "\nERROR --> OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") has no optical surface assigned. Cannot assign type to this logical surface instance. Add to it an optical surface before!\n" << std::endl;
			}
		}else{
			if(!(*iT)->GetSurfaceProperty()){
				std::cout << "\nERROR --> OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") has no optical surface assigned. Cannot assign type to this logical surface instance. Add to it an optical surface before!\n" << std::endl;
			}else{
				std::string logsurfname = (*iT)->GetName();
				if(logsurfname != lastsurfname){
					std::cout << "WARNING --> OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >* , ...): Change of name from <" << lastsurfname << "> to <" << logsurfname << "> for logical surface (list element " << iSurf << " at " << (*iT) << "). This function is supposed to apply the type to a set of logical surfaces with the same name. Unexpected behaviours of the simulation might occur!" << std::endl;
					lastsurfname = logsurfname;
				}
				if( lastoptsurf != ((*iT)->GetSurfaceProperty()) ){
					std::cout << "WARNING --> OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << lastsurfname << "> (list element " << iSurf << " at " << (*iT) << ") has an assigned optical surface different (from pointer) to the one processed before. Applying anyway the type <" << type << "> to this optical surface!" << std::endl;
					lastoptsurf = dynamic_cast<G4OpticalSurface*>((*iT)->GetSurfaceProperty());
					lastoptsurf->SetType( OptSurfTypeMap.at(type) );
				}else{
					if(fVerbose>=OptPropManager::kDebug){
						std::cout << "Debug --> OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >* , ...): Surface type <" << type << "> to logical surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") was already applied at the previous iteration. NO NEED to actually apply it for this logical surface instance." << std::endl;
					}
				}
				nAppl++;
			}
		}
	}
	
	if(fVerbose>=OptPropManager::kDetails){
		std::cout << "Details --> OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >*, ...): Applied type <" << type << "> to " << nAppl << " out of " << nSurfs << " logical surfaces named <" << firstsurfname << "> (if list was homogeneous)." << std::endl;
	}
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfType(const std::set<G4LogicalSurface* >*, ...): Exiting the function." << std::endl;
	}
}



void OptPropManager::SetSurfFinish(const G4String& logsurfname, const G4String& finish)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfFinish(const G4String&, ...): Entering the function." << std::endl;
	}
	
	std::set<G4LogicalSurface* > *logsurflist = FindLogSurf(logsurfname);
	
	if(!logsurflist){
		std::cout << "\nERROR --> OptPropManager::SetSurfFinish(const G4String&, ...): Cannot find the logical surface(s) named <" << logsurfname << "> in the table of the instanced logical surfaces!\n" << std::endl;
		return;
	}
	
	if(!(logsurflist->size())){
		std::cout << "\nERROR --> OptPropManager::SetSurfFinish(const G4String&, ...): The list of logical surface(s) named <" << logsurfname << "> is empty!\n" << std::endl;
		return;
	}
	
	
	if(OptSurfFinishMap.find(finish)==OptSurfFinishMap.end()){
		std::cout << "\nERROR --> OptPropManager::SetSurfFinish(const G4String&, ...): The finish <" << finish << "> is not in the list of known surface finishes and can not be applied to the logical surface(s) <" << logsurfname << ">.\n" << std::endl;
		return;
	}
	
	
	SetSurfFinish(logsurflist, finish);
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfFinish(const G4String&, ...): Exiting the function." << std::endl;
	}
}

void OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >* logsurflist, const G4String& finish)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >*, ...): Entering the function." << std::endl;
	}
	
	if(!logsurflist){
		std::cout << "\nERROR --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >* , ...): Null pointer to the list of logical surfaces!\n" << std::endl;
		return;
	}
	
	if(!(logsurflist->size())){
		std::cout << "\nERROR --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >* , ...): The list of logical surfaces is empty!\n" << std::endl;
		return;
	}
	
	
	if(OptSurfFinishMap.find(finish)==OptSurfFinishMap.end()){
		std::cout << "\nERROR --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >* , ...): The finish <" << finish << "> is not in the list of known types and can not be applied to the logical surfaces.\n" << std::endl;
		return;
	}
	
	
	G4OpticalSurface *lastoptsurf = NULL;
	std::string firstsurfname, lastsurfname;
	std::set<G4LogicalSurface* >::iterator iT;
	int nSurfs = logsurflist->size();
	int nAppl = 0;
	
	int iSurf = 0;
	for(iT=logsurflist->begin(); iT!=logsurflist->end(); ++iT, iSurf++){
		if(iSurf==0){
			firstsurfname = (*iT)->GetName();
			lastsurfname = firstsurfname;
			lastoptsurf = dynamic_cast<G4OpticalSurface*>((*iT)->GetSurfaceProperty());
			if(lastoptsurf){
				if(fVerbose>=OptPropManager::kDebug){
					std::cout << "Debug --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >* , ...): Applying finish <" << finish << "> to logical surfaces <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ")." << std::endl;
				}
				lastoptsurf->SetFinish( OptSurfFinishMap.at(finish) );
				nAppl++;
			}else{
				std::cout << "\nERROR --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") has no optical surface assigned. Cannot assign finish to this logical surface instance. Add to it an optical surface before!\n" << std::endl;
			}
		}else{
			if(!(*iT)->GetSurfaceProperty()){
				std::cout << "\nERROR --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") has no optical surface assigned. Cannot assign type to this logical surface instance. Add to it an optical surface before!\n" << std::endl;
			}else{
				std::string logsurfname = (*iT)->GetName();
				if(logsurfname != lastsurfname){
					std::cout << "WARNING --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >* , ...): Change of name from <" << lastsurfname << "> to <" << logsurfname << "> for logical surface (list element " << iSurf << " at " << (*iT) << "). This function is supposed to apply the type to a set of logical surfaces with the same name. Unexpected behaviours of the simulation might occur!" << std::endl;
					lastsurfname = logsurfname;
				}
				if( lastoptsurf != ((*iT)->GetSurfaceProperty()) ){
					std::cout << "WARNING --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << lastsurfname << "> (list element " << iSurf << " at " << (*iT) << ") has an assigned optical surface different (from pointer) to the one processed before. Applying anyway the finish <" << finish << "> to this optical surface!" << std::endl;
					lastoptsurf = dynamic_cast<G4OpticalSurface*>((*iT)->GetSurfaceProperty());
					lastoptsurf->SetFinish( OptSurfFinishMap.at(finish) );
				}else{
					if(fVerbose>=OptPropManager::kDebug){
						std::cout << "Debug --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >* , ...): Surface finish <" << finish << "> to logical surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") was already applied at the previous iteration. NO NEED to actually apply it for this logical surface instance." << std::endl;
					}
				}
				nAppl++;
			}
		}
	}
	
	if(fVerbose>=OptPropManager::kDetails){
		std::cout << "Details --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >*, ...): Applied finish <" << finish << "> to " << nAppl << " out of " << nSurfs << " logical surfaces named <" << firstsurfname << "> (if list was homogeneous)." << std::endl;
	}
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >*, ...): Exiting the function." << std::endl;
	}
}



void OptPropManager::SetSurfSigmaAlpha(const G4String& logsurfname, const G4double& s_a)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfSigmaAlpha(const G4String&, ...): Entering the function." << std::endl;
	}
	
	std::set<G4LogicalSurface* > *logsurflist = FindLogSurf(logsurfname);
	
	
	if(!logsurflist){
		std::cout << "\nERROR --> OptPropManager::SetSurfSigmaAlpha(const G4String&, ...): Cannot find the logical surface(s) named <" << logsurfname << "> in the table of the instanced logical surfaces!\n" << std::endl;
		return;
	}
	
	if(!(logsurflist->size())){
		std::cout << "\nERROR --> OptPropManager::SetSurfSigmaAlpha(const G4String&, ...): The list of logical surface(s) named <" << logsurfname << "> is empty!\n" << std::endl;
		return;
	}
	
	
	SetSurfSigmaAlpha(logsurflist, s_a);
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfSigmaAlpha(const G4String&, ...): Exiting the function." << std::endl;
	}
}


void OptPropManager::SetSurfSigmaAlpha(const std::set<G4LogicalSurface* >* logsurflist, const G4double& s_a)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfSigmaAlpha(const std::set<G4LogicalSurface* >*, ...): Entering the function." << std::endl;
	}
	
	if(!logsurflist){
		std::cout << "\nERROR --> OptPropManager::SetSurfSigmaAlpha(const std::set<G4LogicalSurface* >* , ...): Null pointer to the list of logical surfaces!\n" << std::endl;
		return;
	}
	
	if(!(logsurflist->size())){
		std::cout << "\nERROR --> OptPropManager::SetSurfSigmaAlpha(const std::set<G4LogicalSurface* >* , ...): The list of logical surfaces is empty!\n" << std::endl;
		return;
	}
	
	
	G4OpticalSurface *lastoptsurf = NULL;
	std::string firstsurfname, lastsurfname;
	std::set<G4LogicalSurface* >::iterator iT;
	int nSurfs = logsurflist->size();
	int nAppl = 0;
	
	int iSurf = 0;
	
	for(iT=logsurflist->begin(); iT!=logsurflist->end(); ++iT, iSurf++){
		if(iSurf==0){
			firstsurfname = (*iT)->GetName();
			lastsurfname = firstsurfname;
			lastoptsurf = dynamic_cast<G4OpticalSurface*>((*iT)->GetSurfaceProperty());
			if(lastoptsurf){
				if(fVerbose>=OptPropManager::kDebug){
					std::cout << "Debug --> OptPropManager::SetSurfSigmaAlpha(const std::set<G4LogicalSurface* >* , ...): Applying sigma_alpha value of <" << s_a << "> to logical surfaces <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ")." << std::endl;
				}
				lastoptsurf->SetSigmaAlpha(s_a);
				nAppl++;
			}else{
				std::cout << "\nERROR --> OptPropManager::SetSurfSigmaAlpha(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") has no optical surface assigned. Cannot assign sigma_alpha value to this logical surface instance. Add to it an optical surface before!\n" << std::endl;
			}
		}else{
			if(!(*iT)->GetSurfaceProperty()){
				std::cout << "\nERROR --> OptPropManager::SetSurfSigmaAlpha(const std::set<G4LogicalSurface* >* , ...): logical surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") has no optical surface assigned. Cannot assign sigma_alpha value to this logical surface instance. Add to it an optical surface before!\n" << std::endl;
			}else{
				std::string logsurfname = (*iT)->GetName();
				if(logsurfname != lastsurfname){
					std::cout << "WARNING --> OptPropManager::SetSurfSigmaAlpha(const std::set<G4LogicalSurface* >* , ...): Change of name from <" << lastsurfname << "> to <" << logsurfname << "> for logical surface (list element " << iSurf << " at " << (*iT) << "). This function is supposed to apply the sigma_alpha value to a set of logical surfaces with the same name. Unexpected behaviours of the simulation might occur!" << std::endl;
					lastsurfname = logsurfname;
				}
				if( lastoptsurf != ((*iT)->GetSurfaceProperty()) ){
					std::cout << "WARNING --> OptPropManager::SetSurfSigmaAlpha(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << lastsurfname << "> (list element " << iSurf << " at " << (*iT) << ") has an assigned optical surface different (from pointer) to the one processed before. Applying anyway the sigma_alpha value of <" << s_a << "> to this optical surface!" << std::endl;
					lastoptsurf = dynamic_cast<G4OpticalSurface*>((*iT)->GetSurfaceProperty());
					lastoptsurf->SetSigmaAlpha(s_a);
				}else{
					if(fVerbose>=OptPropManager::kDebug){
						std::cout << "Debug --> OptPropManager::SetSurfSigmaAlpha(const std::set<G4LogicalSurface* >* , ...): Surface sigma_alpha value of <" << s_a << "> to logical surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") was already applied at the previous iteration. NO NEED to actually apply it for this logical surface instance." << std::endl;
					}
				}
				nAppl++;
			}
		}
	}
	
	
	if(fVerbose>=OptPropManager::kDetails){
		std::cout << "Details --> OptPropManager::SetSurfFinish(const std::set<G4LogicalSurface* >*, ...): Applied sigma_alpha value of <" << s_a << "> to " << nAppl << " out of " << nSurfs << " logical surfaces named <" << firstsurfname << "> (if list was homogeneous)." << std::endl;
	}
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfSigmaAlpha(const std::set<G4LogicalSurface* >*, ...): Exiting the function." << std::endl;
	}
}



void OptPropManager::SetSurfPropFromFile(const G4String& logsurfname, const G4String& filename, const G4String& propertyname)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfPropFromFile(const G4String& logsurfname, ...): Entering the function." << std::endl;
	}
	
	
	std::set<G4LogicalSurface* > *logsurflist = FindLogSurf(logsurfname);
	
	if(!logsurflist){
		std::cout << "\nERROR --> OptPropManager::SetSurfPropFromFile(const G4String& logsurfname, ...): Cannot find the logical surfaces named <"<< logsurfname <<"> in the table of the instanced logical surfaces!\n" << std::endl;
		return;
	}
	
	if(!(logsurflist->size())){
		std::cout << "\nERROR --> OptPropManager::SetSurfPropFromFile(const G4String& logsurfname, ...): The list of logical surface(s) named <" << logsurfname << "> is empty!\n" << std::endl;
		return;
	}
	
	SetSurfPropFromFile( logsurflist, filename, propertyname );
	
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfPropFromFile(const G4String& logsurfname, ...): Exiting the function." << std::endl;
	}
}



void OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >* logsurflist, const G4String& filename, const G4String& propertyname)
{
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >*, ...): Entering the function." << std::endl;
	}
	
	if(!logsurflist){
		std::cout << "\nERROR --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >* , ...): Null pointer to the list of logical surfaces!\n" << std::endl;
		return;
	}
	
	if(!(logsurflist->size())){
		std::cout << "\nERROR --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >* , ...): The list of logical surfaces is empty!\n" << std::endl;
		return;
	}
	
	
	std::vector<G4double> en_vec(0), val_vec(0);
	
	ReadValuesFromFile( filename, en_vec, val_vec );
	
	
	if( (en_vec.size()==0) || (val_vec.size()==0) || (en_vec.size()!=val_vec.size()) ){
		std::cout << "\nERROR --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >* , ...): Wrong size of the vectors read from file <"<< filename <<">. The property <" << propertyname << "> will not be set for the logical surfaces in the input list!\n" << std::endl;
		return;
	}
	
	
	G4OpticalSurface *lastoptsurf = NULL;
	std::string firstsurfname, lastsurfname;
	std::set<G4LogicalSurface* >::iterator iT;
	int nSurfs = logsurflist->size();
	int nAppl = 0;
	int iSurf = 0;
	
	for(iT=logsurflist->begin(); iT!=logsurflist->end(); ++iT, iSurf++){
		if(iSurf==0){
			firstsurfname = (*iT)->GetName();
			lastsurfname = firstsurfname;
			lastoptsurf = dynamic_cast<G4OpticalSurface*>((*iT)->GetSurfaceProperty());
			if(lastoptsurf){
				if(fVerbose>=OptPropManager::kDebug){
					std::cout << "Debug --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >* , ...): Applying property <" << propertyname << "> to logical surfaces <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ")." << std::endl;
				}
				
				G4MaterialPropertiesTable* propTab = lastoptsurf->GetMaterialPropertiesTable();if(!propTab){
					propTab = new G4MaterialPropertiesTable();
					lastoptsurf->SetMaterialPropertiesTable(propTab);
				}
				
				if(propTab->GetProperty( propertyname.c_str() )){
					if(fVerbose>=OptPropManager::kDebug){
						std::cout << "Debug --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >* , ...): Property <" << propertyname << "> already present for logical surface <" << lastsurfname << "> (list element " << iSurf << " at " << (*iT) << "). Removing from and re-adding to it." << std::endl;
					}
					propTab->RemoveProperty( propertyname.c_str() );
				}
				propTab->AddProperty( propertyname.c_str(), (G4double*)&en_vec.at(0), (G4double*)&val_vec.at(0), en_vec.size() );
				nAppl++;
			}else{
				std::cout << "\nERROR --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") has no optical surface assigned. Cannot set property <" << propertyname << "> logical surface instance. Add to it an optical surface before!\n" << std::endl;
			}
		}else{
			if(!(*iT)->GetSurfaceProperty()){
				std::cout << "\nERROR --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >* , ...): logical surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") has no optical surface assigned. Cannot apply the property <" << propertyname << "> to logical surface instance. Add to it an optical surface before!\n" << std::endl;
			}else{
				std::string logsurfname = (*iT)->GetName();
				if(logsurfname != lastsurfname){
					std::cout << "WARNING --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >* , ...): Change of name from <" << lastsurfname << "> to <" << logsurfname << "> for logical surface (list element " << iSurf << " at " << (*iT) << "). This function is supposed to apply the property <" << propertyname << "> to a set of logical surfaces with the same name. Unexpected behaviours of the simulation might occur!" << std::endl;
					lastsurfname = logsurfname;
				}
				if( lastoptsurf != ((*iT)->GetSurfaceProperty()) ){
					std::cout << "WARNING --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >* , ...): logical border surface <" << lastsurfname << "> (list element " << iSurf << " at " << (*iT) << ") has an assigned optical surface different (from pointer) to the one processed before. Applying anyway the property <" << propertyname << "> to this optical surface!" << std::endl;
					lastoptsurf = dynamic_cast<G4OpticalSurface*>((*iT)->GetSurfaceProperty());
					
					G4MaterialPropertiesTable* propTab = lastoptsurf->GetMaterialPropertiesTable();if(!propTab){
						propTab = new G4MaterialPropertiesTable();
						lastoptsurf->SetMaterialPropertiesTable(propTab);
					}
					
					if(propTab->GetProperty( propertyname.c_str() )){
						if(fVerbose>=OptPropManager::kDebug){
							std::cout << "Debug --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >* , ...): Property <" << propertyname << "> already present for logical surface <" << lastsurfname << "> (list element " << iSurf << " at " << (*iT) << "). Removing from and re-adding to it." << std::endl;
						}
						propTab->RemoveProperty( propertyname.c_str() );
					}
					propTab->AddProperty( propertyname.c_str(), (G4double*)&en_vec.at(0), (G4double*)&val_vec.at(0), en_vec.size() );
				}else{
					if(fVerbose>=OptPropManager::kDebug){
						std::cout << "Debug --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >* , ...): Surface property <" << propertyname << "> to logical surface <" << (*iT)->GetName() << "> (list element " << iSurf << " at " << (*iT) << ") was already applied at the previous iteration. NO NEED to actually apply it for this logical surface instance." << std::endl;
					}
				}
				nAppl++;
			}
		}
	}
	
	
	
	if(fVerbose>=OptPropManager::kDetails){
		std::cout << "Details --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >*, ...): Applied property <" << propertyname << "> to " << nAppl << " out of " << nSurfs << " logical surfaces named <" << firstsurfname << "> (if list was homogeneous)." << std::endl;
	}
	
	
	
	if(fVerbose>=OptPropManager::kDebug){
		std::cout << "Debug --> OptPropManager::SetSurfPropFromFile(const std::set<G4LogicalSurface* >*, ...): Exiting the function." << std::endl;
	}
}


void OptPropManager::BuildOpticalSurface(const G4String& optsurfname, const G4String& model, const G4String& type, const G4String& finish )
{
	G4SurfacePropertyTable *surftab = (G4SurfacePropertyTable*)G4OpticalSurface::GetSurfacePropertyTable();
	size_t nOptSurf = G4OpticalSurface::GetNumberOfSurfaceProperties();
	
	if( (!surftab) ){
		//This is a big problem as it is a static class member!!!
		std::cout << "\nERROR --> OptPropManager::BuildOpticalSurface(...): Null pointer for the static table of the optical surfaces!\n" << std::endl;
		return; 
	}
	
	if(OptSurfModelMap.find(model)==OptSurfModelMap.end()){
		std::cout << "\nERROR --> OptPropManager::BuildOpticalSurface(...): The model <" << model << "> is not a known model for optical surfaces! The surface <" << optsurfname << "> will not be built!\n" << std::endl;
		return;
	}
	
	if(OptSurfTypeMap.find(type)==OptSurfTypeMap.end()){
		std::cout << "\nERROR --> OptPropManager::BuildOpticalSurface(...): The type <" << type << "> is not a known type for optical surfaces! The surface <" << optsurfname << "> will not be built!\n" << std::endl;
		return;
	}
	
	if(OptSurfFinishMap.find(finish)==OptSurfFinishMap.end()){
		std::cout << "\nERROR --> OptPropManager::BuildOpticalSurface(...): The finish <" << finish << "> is not a known finish for optical surfaces! The surface <" << optsurfname << "> will not be built!\n" << std::endl;
		return;
	}
	
	//Check if it alredy exists. If it does do not proceed to the creation of this surface as this might overwrite another optical surface used by some logical surface
	G4OpticalSurface *optsurf = FindOptSurf(optsurfname);
	if( optsurf ){
		std::cout << "\nERROR --> OptPropManager::BuildOpticalSurface(...): The surface <" << optsurfname << "> already exists! Try to change the settings to this optical surface instead of building a new one.\n" << std::endl;
		return;
	}
	
	optsurf = new G4OpticalSurface(optsurfname, OptSurfModelMap[model], OptSurfFinishMap[finish], OptSurfTypeMap[type]);
	optsurf->SetMaterialPropertiesTable( new G4MaterialPropertiesTable() );
}



int OptPropManager::BuildLogicalBorderSurface(const G4String& logsurfname, const G4String& physvol1, const G4String& physvol2, const G4String& optsurfname )
{
	//All the procedure should avoid that two logical surfaces with the same name can be created by mistake
	
	G4PhysicalVolumeStore *pPhysVolStore = G4PhysicalVolumeStore::GetInstance();
	
	if(!pPhysVolStore) return -1;
	
	
	if(!fDetConstr) return -2;
	
	DetConstrOptPh::PVmap *pv_map = (DetConstrOptPh::PVmap*)fDetConstr->GetVolsMap();
	
	if(!pv_map) return -3;
	
	std::vector<G4VPhysicalVolume*> vol1_vec, vol2_vec;
	
	
	
	if( pv_map->find(physvol1)==pv_map->end() ){
		std::cout << "\nWARNING --> OptPropManager::BuildLogicalBorderSurface(...): The physical volume <" << physvol1 << "> does not exists in the list of volumes. The logical surface <" << logsurfname << "> can not be created.\n" << std::endl;
		return -4;
	}else{
		//const std::string vname = physvol1.data();
		vol1_vec = (*pv_map)[ physvol1 ];
		//vol1_vec = (*pv_map)[ vname ];
	}
	
	if( pv_map->find(physvol2)==pv_map->end() ){
		std::cout << "\nWARNING --> OptPropManager::BuildLogicalBorderSurface(...): The physical volume <" << physvol2 << "> does not exists in the list of volumes. The logical surface <" << logsurfname << "> can not be created.\n" << std::endl;
		return -5;
	}else{
		//const std::string vname = physvol2.data();
		vol2_vec = (*pv_map)[ physvol2 ];
		//vol2_vec = (*pv_map)[ vname ];
	}
	
	
	G4SurfacePropertyTable *optsurftab = (G4SurfacePropertyTable*)G4OpticalSurface::GetSurfacePropertyTable();
	if( (!optsurftab) ) return -6; //This is a big problem as it is a static class member
	
	
	//Check if the optical surface already exists
	G4OpticalSurface *optsurf = NULL;
	size_t nOptSurf = G4OpticalSurface::GetNumberOfSurfaceProperties();
	if(optsurfname != G4String("")){
		for(size_t iSurf=0; iSurf<nOptSurf; iSurf++){
			if( (optsurftab->at(iSurf)->GetName())==optsurfname ) optsurf = dynamic_cast<G4OpticalSurface*>(optsurftab->at(iSurf));
		}
		
		if(!optsurf){
			std::cout << "\nWARNING --> OptPropManager::BuildLogicalBorderSurface(...): The optical surface <" << optsurfname << "> could not be found in the table of all optical surfaces. A dummy surface with the given name will be created for building the logical surface <" << logsurfname << ">.\n" << std::endl;
		}
		
		optsurf = new G4OpticalSurface(optsurfname);
		optsurf->SetMaterialPropertiesTable( new G4MaterialPropertiesTable() );
	}
	
	return BuildLogicalBorderSurface(logsurfname, vol1_vec, vol2_vec, optsurf);
	
}


int OptPropManager::BuildLogicalBorderSurface(const G4String& logsurfname, const std::vector<G4VPhysicalVolume*>& vol1_vec, const std::vector<G4VPhysicalVolume*>& vol2_vec, const G4OpticalSurface* optsurf )
{
	//All the procedure should avoid that two logical surfaces with the same name can be created by mistake
	
	G4int iSurf=0;
	
	size_t nVols1 = vol1_vec.size();
	if(nVols1==0){
		std::cout << "\nERROR --> OptPropManager::BuildLogicalBorderSurface(...): The vector of physical volumes <vol1_vec> is empty!!! The logical surface <" << logsurfname << "> will not be built.\n" << std::endl;
		return -7;
	}
	
	size_t nVols2 = vol2_vec.size();
	if(nVols2==0){
		std::cout << "\nERROR --> OptPropManager::BuildLogicalBorderSurface(...): The vector of physical volumes <vol2_vec>" << " is empty!!! The logical surface <" << logsurfname << "> will not be built.\n" << std::endl;
		return -8;
	}
	
	if(!optsurf){
		std::cout << "\nERROR --> OptPropManager::BuildLogicalBorderSurface(...): Passed a null pointer as optical surface. The logical border surfaces <" << logsurfname << "> will not be built.\n" << std::endl;
		return -9;
	}
	
	
	for(size_t iVol1 = 0; iVol1<nVols1; iVol1++){
		for(size_t iVol2 = 0; iVol2<nVols2; iVol2++){
			
			G4VPhysicalVolume *vol1 = vol1_vec.at(iVol1);
			G4VPhysicalVolume *vol2 = vol2_vec.at(iVol2);
			
			G4LogicalBorderSurface *logsurf = G4LogicalBorderSurface::GetSurface(vol1, vol2);
			
			if(logsurf){
				
				std::cout << "\nWARNING --> OptPropManager::BuildLogicalBorderSurface(...): The logical surface between volume <" << vol1->GetName() << "> at address" << vol1 << " and volume <" << vol2->GetName() << "> at address" << vol2 << " already exists with name <" << logsurf->GetName() << ">. Renaming it as <" << logsurfname << ">" << std::endl;
				if(optsurf){
					std::cout << " and assigning the surface <" << optsurf->GetName() << "> to it.\n" << std::endl;
				}else{
					std::cout << ".\n" << std::endl;
				}
				
				DeRegisterLogSurf(logsurf);
				
				logsurf->SetName(logsurfname);
				logsurf->SetSurfaceProperty( (G4OpticalSurface*)optsurf );
				
				RegisterLogSurf(logsurf);
				
				iSurf++;
				
			}else{
				
				if(fVerbose >= OptPropManager::kInfo){
					std::cout << "Info --> OptPropManager::BuildLogicalBorderSurface(...): building logical surface <" << logsurfname << ">." << std::endl;
				}
				RegisterLogSurf( new G4LogicalBorderSurface(logsurfname, vol1, vol2, (G4OpticalSurface*)optsurf) );
				
				iSurf++;
			}
		}//Exiting from internal loop (jVol) of phys volumes
	}//Exiting from external loop (iVol) of phys volumes
	
	return iSurf;
}



void OptPropManager::SetOpticalSurface(const G4String& logsurfname, const G4String& optsurfname)
{
	G4OpticalSurface *optsurf = NULL;
	
	G4SurfacePropertyTable *optsurftab = (G4SurfacePropertyTable*)G4OpticalSurface::GetSurfacePropertyTable();
	
	if( !optsurftab ) return; //This is a big problem as they are static class members
	
	
	//size_t nLogSurfs = G4LogicalBorderSurface::GetNumberOfBorderSurfaces();
	size_t nOptSurfs = G4OpticalSurface::GetNumberOfSurfaceProperties();
	
	if( nOptSurfs==0 ) return;
	
	
	for(size_t iSurf=0; iSurf<nOptSurfs; iSurf++){
		if( (optsurftab->at(iSurf)->GetName())==optsurfname ) optsurf = static_cast<G4OpticalSurface*>(optsurftab->at(iSurf));
	}
	
	if(!optsurf) return;
	
	
	std::set< G4LogicalSurface* >* logsurflist = FindLogSurf(logsurfname);
	
	if(!logsurflist) return;
	
	if(logsurflist->size()) return;
	
	std::set< G4LogicalSurface* >::iterator iT;
	for(iT=logsurflist->begin(); iT!=logsurflist->end(); ++iT){
		(*iT)->SetSurfaceProperty(optsurf);
	}
}



void OptPropManager::RegisterLogSurf(G4LogicalSurface* logsurf)
{
	if(!logsurf) return;
	
	if(fLogSurfMap.find(logsurf->GetName())==fLogSurfMap.end()){
		std::set< G4LogicalSurface* > ls_set;
		ls_set.insert(logsurf);
		fLogSurfMap[logsurf->GetName()] = ls_set;
	}else{
		(fLogSurfMap[logsurf->GetName()]).insert(logsurf);
	}
}


void OptPropManager::DeRegisterLogSurf(G4LogicalSurface* logsurf)
{
	if(!logsurf) return;
	
	if(fLogSurfMap.find(logsurf->GetName())!=fLogSurfMap.end()){
		(fLogSurfMap[logsurf->GetName()]).erase(logsurf);
	}
}



#endif
