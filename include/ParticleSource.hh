#ifndef __PARTICLE_SOURCE_OPT_PH_HH__
#define __PARTICLE_SOURCE_OPT_PH_HH__

#include "TGraph.h"

#include "globals.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4Navigator.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"


#include <set>
#include <vector>

using std::set;
using std::vector;

class PrimaryGeneratorActionOptPh;
class ParticleSourceOptPhMessenger;
class G4VPhysicalVolume;



class ParticleSourceOptPh: public G4VPrimaryGenerator
{
public:
	
	enum class EnergyDistType{
		kMono,
		kSpectrum
	};
	
	enum class AngDistType{
		kDirection,
		kIso
	};
	
	enum class SourceType{
		kPoint,
		kVolume
	};
	
	enum class SourceShape{
		kNone,
		kSphere,
		kCylinder,
		kBox
	};
	
	
	ParticleSourceOptPh(PrimaryGeneratorActionOptPh *primGenAct, G4int nPrim=1, G4int verb=0);
	virtual ~ParticleSourceOptPh();
	
	
	//THIS IS THE CENTRAL FUNCTION OF ALL THE CLASS
	void GeneratePrimaryVertex(G4Event *pEvent);
	
	//This method is eventually used to switch to geantino (for debugging) otherwise the standard "opticalphoton" particle is defined already in the constructor
	void SetParticleDef(G4ParticleDefinition* aParticleDefinition);
	void SetKinEnergy(G4double KinEnergy);
	inline void SetDirection(G4ThreeVector aDirection){fMomentumDirection = aDirection.unit();};
	void SetMomentum(G4double aMomentum);
	void SetMomentum(G4ParticleMomentum aMomentum);
	
	inline void SetMaxConfineLoop(G4int _max){fMaxConfineLoop = _max;};
	inline void SetPosDisType(SourceType posType) { fSourcePosType = posType; }
	inline void SetPosDisShape(SourceShape distShape) { fShape = distShape; }
	inline void SetCenterCoords(G4ThreeVector hCenterCoords) { fCenterCoords = hCenterCoords; }

	inline void SetHalfX(G4double dHalfx) { fHalfx = dHalfx; }
	inline void SetHalfY(G4double dHalfy) { fHalfy = dHalfy; }
	inline void SetHalfZ(G4double dHalfz) { fHalfz = dHalfz; }
	inline void SetRadius(G4double dRadius) { fRadius = dRadius; }
	
	void SetEnergySpectrum(G4String filename);
	
	inline void SetVerbosity(G4int iVerbosityLevel) { fVerbosityLevel = iVerbosityLevel; };
	
	inline void SetAngDistType(AngDistType angDistType) { fAngDistType = angDistType; }
	
	inline void SetPhotonPolar(G4ThreeVector hPol) { fPolarization = hPol.unit(); }
	inline void SetEnergyDisType(EnergyDistType enerDistType){ fEnergyDisType = enerDistType; };
	
	void SetPrimNb(G4int nprim);
	inline G4int GetPrimNb()const{return fPrimNb;};
	
	const G4VPhysicalVolume* GetPrimVol(){return fVolPrim;};
	const G4ThreeVector& GetPrimPos()const{return fPosPrim;};
	const vector<G4ParticleMomentum>& GetPrimMom()const{return fMomPrim;};
	const vector<G4ThreeVector>& GetPrimPol()const{return fPolPrim;};
	
	
	
	const G4String &GetParticleType()const{ return fParticleDefinition->GetParticleName(); };
	G4double GetTotalEnergy()const{ return fTotEnergy; }
	G4double GetKinEnergy()const{ return fKinEnergy; }
	G4double GetMomentum()const{ return fMom; }
	const G4ThreeVector &GetDirection()const{ return fMomentumDirection; }
	const G4ThreeVector &GetPosition()const{ return fPosition; }
	const G4ThreeVector &GetPolarization()const{ return fPolarization; }
	G4String GetPosDisType() const; //Inline
	G4String GetShapeType() const; //Inline
	G4String GetAngDistrType() const; //Inline
	G4String GetEnergyDistType() const; //Inline
	const G4bool IsEnergySpectrumLoaded(){return fEnergySpectLoaded;};
	
	
	void PrintParticle();
	void PrintDirection();
	void PrintPolar();

	void GeneratePointSource();
	void GeneratePointsInVolume();
	G4bool IsSourceConfined();
	void ConfineSourceToVolume(G4String);

	void GenerateIsotropic();
	
	double SampleEnergyFromSpectrum();

protected:
	virtual void SetInitialValues();
	
private:
	PrimaryGeneratorActionOptPh *fPrimGenAct;
	ParticleSourceOptPhMessenger *fMessenger;
	G4Navigator *fNavigator;

protected:
	G4int fPrimNb;
	
	G4bool fEnergySpectLoaded; //This is a flag to load the energy spectrum from a txt file
	TGraph *fInvIntegralSpectrum; //Inverted spectral integral (x and y swapped)
	
	G4ParticleDefinition *fParticleDefinition;
	
	G4String fPartName;
	G4double fMass, fCharge;
	G4double fTotEnergy, fKinEnergy, fMom;
	G4ParticleMomentum fMomentumDirection;
	G4ThreeVector fPolarization;
	
	
	SourceType fSourcePosType;
	SourceShape fShape;
	G4ThreeVector fCenterCoords;
	G4double fHalfx;
	G4double fHalfy;
	G4double fHalfz;
	G4double fRadius;
	G4bool fConfine;
	set<G4String> fVolumeNames;
	AngDistType fAngDistType;
	EnergyDistType fEnergyDisType;
	G4int fMaxConfineLoop;
	
	//Stuff to be used internally for generating the primaries
	G4VPhysicalVolume *fVol;
	G4ThreeVector fPosition;
	G4double fTime;
	
	//Stuff that is given to the outside world from the class interface (e.g. made for the analysis manager)
	vector<G4double> fEnPrim;
	G4VPhysicalVolume *fVolPrim; //This should never ever be deleted inside this class, only reassigned when needed
	G4ThreeVector fPosPrim;
	vector<G4ThreeVector> fPolPrim;
	vector<G4ParticleMomentum> fMomPrim;
	vector<G4String> fPhysVolPrim;
	
	G4int fVerbosityLevel;
};


//Definitions of inline functions

inline G4String ParticleSourceOptPh::GetPosDisType() const
{
	if(fSourcePosType == SourceType::kPoint) return G4String("Point");
	if(fSourcePosType == SourceType::kVolume) return G4String("Volume");
};

inline G4String ParticleSourceOptPh::GetShapeType() const
{
	if(fShape == SourceShape::kNone) return G4String("None");
	if(fShape == SourceShape::kSphere) return G4String("Sphere");
	if(fShape == SourceShape::kCylinder) return G4String("Cylinder");
	if(fShape == SourceShape::kBox) return G4String("Box");
};

inline G4String ParticleSourceOptPh::GetAngDistrType() const
{
	if(fAngDistType == AngDistType::kDirection) return G4String("Direction");
	if(fAngDistType == AngDistType::kIso) return G4String("Iso");
};

inline G4String ParticleSourceOptPh::GetEnergyDistType() const
{
	if(fEnergyDisType == EnergyDistType::kMono) return G4String("Mono");
	if(fEnergyDisType == EnergyDistType::kSpectrum) return G4String("Spectrum");
};

#endif

