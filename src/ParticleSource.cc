#include "TGraph.h"
#include "TF1.h"
#include "Math/Integrator.h"

#include "ParticleSource.hh"
#include "ParticleSourceMessenger.hh"
#include "PrimGenAction.hh"

#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <sstream>
#include <cmath>

using std::stringstream;



ParticleSourceOptPh::ParticleSourceOptPh(PrimaryGeneratorActionOptPh *primGenAct, G4int nPrim, G4int verb):
fPrimGenAct(primGenAct),
fPrimNb(nPrim),
fEnergySpectLoaded(false),
fInvIntegralSpectrum(nullptr)
{
	//fPrimGenAct = primGenAct;
	
	//fPrimNb = nPrim;
	
	fParticleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
	
	if(!fParticleDefinition){
		G4Exception("ParticleSourceOptPh::SetParticleDef(...)","Event0100", FatalException,"Cannot find the \"opticalphoton\" definition in the G4ParticleTable.");
	}
	
	fVerbosityLevel = verb;
	
	SetInitialValues();
	
	fMessenger = new ParticleSourceOptPhMessenger(this);
	fNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
	
	G4cout << "\nParticle source initialized." << G4endl;
}


ParticleSourceOptPh::~ParticleSourceOptPh()
{
	delete fMessenger;
	if(fInvIntegralSpectrum) delete fInvIntegralSpectrum;
}


void ParticleSourceOptPh::SetInitialValues()
{
	fPartName = fParticleDefinition->GetParticleName();
	fCharge = fParticleDefinition->GetPDGCharge();
	fMass = fParticleDefinition->GetPDGMass();
	
	fMomentumDirection = G4ParticleMomentum(0., 0., 1.);
	fKinEnergy = 9.69*eV; //128 nm optical photons
	fTotEnergy = fMass + fKinEnergy;
	fMom = std::sqrt( fTotEnergy*fTotEnergy - fMass*fMass );
	
	
	fVol = NULL;
	fPosition = G4ThreeVector(0.,0.,0.);
	fTime = 0.;
	fPolarization = G4ThreeVector(1.,0.,0.);
	
	
	fSourcePosType = SourceType::kPoint;
	fShape = SourceShape::kNone;
	fCenterCoords=G4ThreeVector(0.,0.,0.);
	
	fHalfx = 0.;
	fHalfy = 0.;
	fHalfz = 0.;
	fRadius = 0.;
	fConfine = false;
	fAngDistType = AngDistType::kIso;
	fEnergyDisType = EnergyDistType::kMono;
	
	fMaxConfineLoop = 100000;
	
	fVolumeNames.clear();
	
	fVolPrim = NULL;
	fEnPrim.resize(fPrimNb,fTotEnergy);
	fPosPrim = G4ThreeVector(0.,0.,0.);
	fMomPrim.resize(fPrimNb);
	fPolPrim.resize(fPrimNb);
	
	fPhysVolPrim.resize(fPrimNb);
}


//THIS IS THE CENTRAL FUNCTION OF THIS CLASS
void ParticleSourceOptPh::GeneratePrimaryVertex(G4Event* evt)
{
	if(!fParticleDefinition)
	{
		G4cerr << "\nERROR --> ParticleSourceOptPh::GeneratePrimaryVertex: No particle is defined!" << G4endl;
		return;
	}
	
	// Position
	G4bool srcconf = false;
	
	
	if(fSourcePosType == SourceType::kPoint){
		GeneratePointSource();
	}else if(fSourcePosType == SourceType::kVolume){
		
		if(fConfine == true){
			
			G4int LoopCount = 0;
			do{
				GeneratePointsInVolume();
				LoopCount++;
				
				if(LoopCount == fMaxConfineLoop)
				{
					G4cerr << "\nParticleSourceOptPh::GeneratePrimaryVertex(...) --> ERROR:" << G4endl;
					G4cerr << "*************************************" << G4endl;
					G4cerr << "LoopCount = " << fMaxConfineLoop << G4endl;
					G4cerr << "Either the source distribution >> confinement" << G4endl;
					G4cerr << "or any confining volume may not overlap with" << G4endl;
					G4cerr << "the source distribution or any confining volumes" << G4endl;
					G4cerr << "may not exist" << G4endl;
					G4cerr << "If you have set confine then this will be ignored" << G4endl;
					G4cerr << "for this event." << G4endl;
					G4cerr << "*************************************" << G4endl;
					break; //Avoids an infinite loop
				}
				
			}while(!IsSourceConfined());
			
		}else{
			GeneratePointsInVolume();
		}
	}
	
	fPosPrim = fPosition;
	fVolPrim = fVol;
	
	
	// create a new vertex
	G4PrimaryVertex *vertex = new G4PrimaryVertex(fPosition, fTime);
	
	G4String PhysVolName = fNavigator->LocateGlobalPointAndSetup(fPosition)->GetName();
	
	
	if(fVerbosityLevel > 1)
		G4cout << "\nDetail --> ParticleSourceOptPh::GeneratePrimaryVertex: Creating primaries and assigning to vertex." << G4endl;
	
	
	if(fVerbosityLevel > 0){
		G4cout << "\nInfo --> ParticleSourceOptPh::GeneratePrimaryVertex: " << G4endl;
		G4cout << " Number of prymaries per event: " << fPrimNb << G4endl;
		G4cout << "                 Particle name: " << fParticleDefinition->GetParticleName() << G4endl;
		if(fEnergyDisType==EnergyDistType::kMono){
			G4cout << "                  Energy(Mono): " << fTotEnergy/eV << " eV" << G4endl;
		}
		G4cout << "                          Mass: " << fMass/MeV << " MeV" << G4endl;
		G4cout << "                      Position: " << fPosition << G4endl;
		if(fAngDistType == AngDistType::kDirection){
			G4cout << "                     Direction: " << fMomentumDirection << G4endl;
			G4cout << "                  Polarization: " << fPolarization << G4endl;
		}else if(fAngDistType == AngDistType::kIso){
			G4cout << "                     Direction: isotropic" << G4endl;
			G4cout << "                  Polarization: random" << G4endl;
		}
	}
	
	
	if(fVerbosityLevel > 1) G4cout << "\nDetail --> ParticleSourceOptPh::GeneratePrimaryVertex: " << G4endl;
	
	for(G4int iPart = 0; iPart < fPrimNb; iPart++)
	{
		if(fAngDistType == AngDistType::kIso){
			GenerateIsotropic();
			if(fVerbosityLevel > 1){
				G4cout << "                  Direction["<<iPart<<"]: " << fMomentumDirection.unit() << G4endl;
				G4cout << "               Polarization["<<iPart<<"]: " << fPolarization << G4endl;
			}
		}
		
		if(fEnergyDisType==EnergyDistType::kSpectrum){
			
			if( (!fInvIntegralSpectrum) || (!fEnergySpectLoaded) ){
				G4Exception("ParticleSourceOptPh::GeneratePrimaryVertex","Event0102",FatalException,"Energy spectrum not loaded.");
			}
			
			SetKinEnergy(SampleEnergyFromSpectrum());
			
			if(fVerbosityLevel > 1){
					G4cout << "                     Energy["<<iPart<<"]: " << fTotEnergy/eV << " eV" << G4endl;
			}
		}
		
		G4PrimaryParticle *particle = new G4PrimaryParticle(fParticleDefinition);
		particle->SetMass(fMass);
		particle->SetMomentumDirection(fMomentumDirection.unit());
		particle->SetTotalEnergy(fTotEnergy);
		particle->SetCharge(fCharge);
		particle->SetPolarization(fPolarization.unit());
		vertex->SetPrimary(particle);
		
		fMomPrim.at(iPart) = particle->GetMomentum();
		fPolPrim.at(iPart) = particle->GetPolarization();
	}
	evt->AddPrimaryVertex(vertex);
	if(fVerbosityLevel > 0) G4cout << "\nInfo --> ParticleSourceOptPh::GeneratePrimaryVertex: Primary vertex generated with " << fMomPrim.size() << " particles." << G4endl;
}


void ParticleSourceOptPh::SetPrimNb(G4int nprim)
{
	fPrimNb = nprim;
	
	if(fPrimGenAct) fPrimGenAct->SetPrimNb(fPrimNb);
	
	fEnPrim.resize(fPrimNb,fTotEnergy);
	fMomPrim.resize(fPrimNb);
	fPolPrim.resize(fPrimNb);
	fPhysVolPrim.resize(fPrimNb);
	
	return;
}


void ParticleSourceOptPh::SetParticleDef(G4ParticleDefinition *aParticleDefinition)
{
	if(!aParticleDefinition){
		G4Exception("ParticleSourceOptPh::SetParticleDef","Event0101",FatalException,"Null pointer is given.");
	}
	fParticleDefinition = aParticleDefinition;
	fCharge = fParticleDefinition->GetPDGCharge();
	fMass = fParticleDefinition->GetPDGMass();
	fPartName = fParticleDefinition->GetParticleName();
}


void ParticleSourceOptPh::SetKinEnergy(G4double dKinEnergy)
{
	fKinEnergy = dKinEnergy;
	if(!fParticleDefinition){
		//Assuming zero mass
		fMass = 0;
		fTotEnergy = fKinEnergy;
		fMom = fTotEnergy;
		return;
	}
	fTotEnergy = fMass + fKinEnergy;
	fMom = std::sqrt( fTotEnergy*fTotEnergy - fMass*fMass );
}


void ParticleSourceOptPh::SetMomentum(G4double aMomentum)
{
	fMom = aMomentum;
	if(!fParticleDefinition){
		//Assuming zero mass
		fMass = 0;
		fTotEnergy = fMom;
		fKinEnergy = fTotEnergy;
	}
	
	fTotEnergy = std::sqrt( fMom*fMom + fMass*fMass );
	fKinEnergy = fTotEnergy - fMass;
}


void ParticleSourceOptPh::SetMomentum(G4ParticleMomentum aMomentum)
{
	fMomentumDirection = aMomentum.unit();
	if(fEnergyDisType==EnergyDistType::kMono) SetMomentum(aMomentum.mag());
}


void ParticleSourceOptPh::PrintParticle(){ 
	G4cout << "\nSelected particle: " << fParticleDefinition->GetParticleName() << G4endl;
	//std::cout << "\nSelected particle: " << fParticleDefinition->GetParticleName() << std::endl;
}


void ParticleSourceOptPh::PrintDirection()
{
	if(fAngDistType==AngDistType::kIso){
		G4cout << "\nParticle direction: isotropic" << G4endl;
	}else{
		G4cout << "\nParticle direction: " << fMomentumDirection.unit() << G4endl;
	}
}


void ParticleSourceOptPh::PrintPolar()
{
	if(fAngDistType==AngDistType::kIso){
		G4cout << "\nParticle polarization: random" << G4endl;
	}else{
		G4cout << "\nParticle polarization: " << fPolarization.unit() << G4endl;
	}
}


void ParticleSourceOptPh::ConfineSourceToVolume(G4String hVolumeList)
{
	stringstream hStream;
	hStream.str(hVolumeList);
	G4String hVolumeName;
	
	// store all the volume names
	while(!hStream.eof())
	{
		hStream >> hVolumeName;
		fVolumeNames.insert(hVolumeName);
	}

	// checks if the selected volumes exist and store all volumes that match
	G4PhysicalVolumeStore *PVStore = G4PhysicalVolumeStore::GetInstance();
	G4bool bFoundAll = true;

	set<G4String> hActualVolumeNames;
	for(set<G4String>::iterator pIt = fVolumeNames.begin(); pIt != fVolumeNames.end(); pIt++){
		G4String hRequiredVolumeName = *pIt;
		G4bool bMatch = false;

		if(bMatch = (hRequiredVolumeName.last('*') != std::string::npos))
			hRequiredVolumeName = hRequiredVolumeName.strip(G4String::trailing, '*');

		G4bool bFoundOne = false;
		for(G4int iIndex = 0; iIndex < (G4int) PVStore->size(); iIndex++)
		{
			G4String hName = (*PVStore)[iIndex]->GetName();

			if((bMatch && (hName.substr(0, hRequiredVolumeName.size())) == hRequiredVolumeName) || hName == hRequiredVolumeName)
			{
				hActualVolumeNames.insert(hName);
				bFoundOne = true;
			}
		}

		bFoundAll = bFoundAll && bFoundOne;
	}

	if(bFoundAll)
	{
		fVolumeNames = hActualVolumeNames;
		fConfine = true;

		if(fVerbosityLevel >= 1)
			G4cout << "Source confined to volumes: " << hVolumeList << G4endl;

		if(fVerbosityLevel >= 2)
		{
			G4cout << "Volume list: " << G4endl;

			for(set<G4String>::iterator pIt = fVolumeNames.begin(); pIt != fVolumeNames.end(); pIt++)
				G4cout << *pIt << G4endl;
		}
	}
	else if(fVolumeNames.empty())
		fConfine = false;
	else
	{
		G4cout << "ParticleSourceOptPh::ConfineSourceToVolume(...) --> ERROR: One or more volumes do not exist! " << G4endl;
		G4cout << " Ignoring confine condition" << G4endl;
		fVolumeNames.clear();
		fConfine = false;
	}
}


void ParticleSourceOptPh::GeneratePointSource()
{
	// Generates Points given the point source.
	fPosition = fCenterCoords;
	G4ThreeVector nullvect(0., 0., 0.);
	G4ThreeVector *ptr = &nullvect;
	fVol = fNavigator->LocateGlobalPointAndSetup(fPosition, ptr, true);
}


void ParticleSourceOptPh::GeneratePointsInVolume()
{
	if(fShape == SourceShape::kNone){
		G4Exception("ParticleSourceOptPh::GeneratePointsInVolume","Event0103",FatalException,"Volume shape not defined. Define it with the proper UI command. Allowed values are [\"Spere\", \"Cylinder\", \"Box\"]!");
	}
	
	
	G4double x = 0., y = 0., z = 0.;

	if(fShape == SourceShape::kSphere)
	{
		x = fRadius * 2.;
		y = fRadius * 2.;
		z = fRadius * 2.;
		
		do{
			x = G4UniformRand();
			y = G4UniformRand();
			z = G4UniformRand();

			x = (2*x-1)*fRadius;
			y = (2*y-1)*fRadius;
			z = (2*z-1)*fRadius;
		}while(((x * x) + (y * y) + (z * z)) > (fRadius * fRadius));
		
	}

	else if(fShape == SourceShape::kCylinder)
	{
		x = fRadius * 2.;
		y = fRadius * 2.;
		
		z = G4UniformRand();
		z = (2*z - 1)*fHalfz;
		
		
		do{
			x = G4UniformRand();
			y = G4UniformRand();
			
			x = (2*x-1)*fRadius;
			y = (2*y-1)*fRadius;
		}while(((x * x) + (y * y)) > (fRadius * fRadius));
	}

	else if(fShape == SourceShape::kBox)
	{
		x = 2*(G4UniformRand()-0.5)*fHalfx;
		y = 2*(G4UniformRand()-0.5)*fHalfy;
		z = 2*(G4UniformRand()-0.5)*fHalfz;
	}
		
	
	fPosition = fCenterCoords + G4ThreeVector(x,y,z);
	G4ThreeVector nullvect(0., 0., 0.);
	G4ThreeVector *ptr = &nullvect;
	fVol = fNavigator->LocateGlobalPointAndSetup(fPosition, ptr, true);
}


G4bool ParticleSourceOptPh::IsSourceConfined()
{
	// Method to check point is within the volume specified
	if(fConfine == false)
		G4cerr << "\nParticleSourceOptPh::IsSourceConfined() --> ERROR: Confine is false" << G4endl;
	
	G4ThreeVector nullvect(0., 0., 0.);
	G4ThreeVector *ptr = &nullvect;

	// Check fParticlePosition is within a volume in our list
	G4VPhysicalVolume *theVolume;

	theVolume = fNavigator->LocateGlobalPointAndSetup(fPosition, ptr, true);
	G4String theVolName = theVolume->GetName();

	set<G4String>::iterator pIt;
	if((pIt = fVolumeNames.find(theVolName)) != fVolumeNames.end())
	{
		if(fVerbosityLevel >= 1)
			G4cout << "Particle is in volume " << *pIt << G4endl;
		return (true);
	}
	else
		return (false);
}


void ParticleSourceOptPh::GenerateIsotropic()
{
	G4double px, py, pz;

	G4double sintheta, sinphi, costheta, cosphi;

	costheta = 1 - 2*G4UniformRand();
	sintheta = std::sqrt(1. - std::pow(costheta,2));

	G4double phi = CLHEP::twopi * G4UniformRand();
	sinphi = std::sin(phi);
	cosphi = std::cos(phi);

	px = sintheta * cosphi;
	py = sintheta * sinphi;
	pz = costheta;

	fMomentumDirection= G4ThreeVector(px,py,pz).unit();
	
	if( (fPartName=="gamma") || (fPartName=="opticalphoton") ){
		
		//Generate the random polarization as a 2pi around the momentum direction
		G4double alpha = CLHEP::twopi * G4UniformRand();
	
		G4double sinalpha = std::sin(alpha);
		G4double cosalpha = std::cos(alpha);
	
		//This is always orthogonal to the momentum direction (try to make the scalar product)
		px = cosalpha*costheta*cosphi - sinalpha*sinphi;
		py = cosalpha*costheta*sinphi + sinalpha*cosphi;
		pz = -cosalpha*sintheta;
	
		fPolarization= G4ThreeVector(px,py,pz).unit();
		
	}
	
	// m_hParticleMomentumDirection now holds unit momentum vector.
	if(fVerbosityLevel > 2) G4cout << "Debug --> ParticleSourceOptPh::GenerateIsotropic: Generating isotropic vector: " << fMomentumDirection << G4endl;
}




void ParticleSourceOptPh::SetEnergySpectrum(G4String filename){
	
	fEnergySpectLoaded = false;
	
	std::ifstream infile(filename.c_str());
	
	if(!infile){
		std::cerr << "\nERROR --> ParticleSourceOptPh::SetEnergySpectrum: Cannot find or open in read mode the file <" << filename << "> with the primary spectrum.\n" << std::endl;
		return;
	}
	
	std::string str;
	std::stringstream ss_tmp;
	
	if(fVerbosityLevel>1){
		std::cout << "Info --> Reading optical photons primary spectrum from file <" << filename << ">:" << std::endl;
	}
	
	TGraph *gr = new TGraph();
	
	int iLine =0;
	G4double ph_en_d;
	G4double val_d;
	while(getline(infile,str)){
		
		ss_tmp.clear(); ss_tmp.str("");
		ss_tmp << str;
		
		ss_tmp >> str;
		if(ss_tmp){
			ph_en_d = std::stod(str);
			
		}else{
			std::cerr << "\nERROR --> ParticleSourceOptPh::SetEnergySpectrum: Line " << (iLine-1) << " of input file <" << filename << "> is corrupted. The spectrum of primary optical photons will not be loaded from this file." << std::endl;
			delete gr;
			return;
		}
		
		ss_tmp >> str;
		if(ss_tmp){//There is only one value while the file format is defined with 2 columns
			val_d = std::stod(str);
		}else{
			std::cerr << "\nERROR --> ParticleSourceOptPh::SetEnergySpectrum: Line " << (iLine-1) << " of input file <" << filename << "> is corrupted. The spectrum of primary optical photons will not be loaded from this file." << std::endl;
			delete gr;
			return;
		}
		
		gr->SetPoint(iLine, ph_en_d*eV, val_d);
		
		iLine++;
	}
	
	gr->Sort(); //Usually it is redundant.
	
	Int_t nPts = gr->GetN();
	
	ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
	
	//Make a TF1 to integrate the spectral distribution of the primary optical photons
	TF1 *func = new TF1("func", [&](Double_t *x, Double_t *y){return gr->Eval(x[0]);}, gr->GetX()[0], gr->GetX()[nPts-1], 0);
	
	
	//Make the already inverted spectrum by swapping the x (the energy) and the y (the integrated spectral density)
	if(fInvIntegralSpectrum) delete fInvIntegralSpectrum;
	fInvIntegralSpectrum = new TGraph();
	
	for(int iPt=0; iPt<nPts; iPt++){
		fInvIntegralSpectrum->SetPoint(iPt, func->Integral( gr->GetX()[0], gr->GetX()[iPt] ), gr->GetX()[iPt] );
	}
	
	//Normalise to 1 the integral
	Double_t maxIntegral = fInvIntegralSpectrum->GetX()[nPts-1];
	
	for(int iPt=0; iPt<nPts; iPt++){
		fInvIntegralSpectrum->GetX()[iPt] /= maxIntegral;
	}
	
	fEnergySpectLoaded = true;
	
	delete func;
	delete gr;
}


double ParticleSourceOptPh::SampleEnergyFromSpectrum(){
	
	G4double en = fInvIntegralSpectrum->Eval(G4UniformRand());
	//if(fVerbosityLevel>1){
	//	G4cout << "Detail --> Sampled energy: " << en/eV << " eV" << G4endl;
	//}
	return en;
	
}

