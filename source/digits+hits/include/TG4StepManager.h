#ifndef TG4_STEP_MANAGER_H
#define TG4_STEP_MANAGER_H

//------------------------------------------------
// The Geant4 Virtual Monte Carlo package
// Copyright (C) 2007 - 2015 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file TG4StepManager.h
/// \brief Definition of the TG4StepManager class 
///
/// \author I. Hrivnacova; IPN, Orsay

#include <Rtypes.h>

#include "TG4StepStatus.h"

#include <G4Step.hh>
#include <G4GFlashSpot.hh>
#include <G4TransportationManager.hh>
#include <G4SteppingManager.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>

#include <TString.h>
#include <TArrayI.h>
#include <TMCProcess.h>

// for inline functions
#include "TG4G3Units.h"
#include "TLorentzVector.h"
#include "TG4SteppingAction.h"
#include "TG4PhysicsManager.h"
#include "TG4TrackInformation.h"
#include "TG4ParticlesManager.h"
#include "TG4TrackManager.h"

class TG4Limits;
class TG4TrackManager;
class TG4SteppingAction;

class G4Track;
class G4SteppingManager;
class G4VPhysicalVolume;
class G4VTouchable;

class TLorentzVector;

/// \ingroup digits_hits
/// \brief Geant4 implementation of the TVirtualMC interface methods                    
/// for access to Geant4 at step level.
///
/// The public methods that do not implement TVirtualMC methods
/// are commented as G4 specific
///
/// \author I. Hrivnacova; IPN, Orsay

class TG4StepManager
{
  public:
    TG4StepManager(const TString& userGeometry);
    virtual ~TG4StepManager();

    // static access method
    static TG4StepManager* Instance();
        
    // methods
    void LateInitialize();
    void StopTrack();
    void StopEvent();
    void StopRun();
    
    // set methods
    void SetStep(G4Step* step, TG4StepStatus status);    // G4 specific
    void SetStep(G4Track* track, TG4StepStatus status);  // G4 specific
    void SetStep(G4GFlashSpot* gflashSpot, TG4StepStatus status);  // G4 specific
    void SetSteppingManager(G4SteppingManager* manager); // G4 specific
    void SetMaxStep(Double_t step);
    void SetMaxStepBack();                               // G4 specific
    void SetMaxNStep(Int_t maxNofSteps); 
    void SetCollectTracks(Bool_t collectTracks);
    void ForceDecayTime(Float_t pdg);
    
    // get methods
    G4Track* GetTrack() const;                            // G4 specific
    G4Step*  GetStep() const;                             // G4 specific
    TG4StepStatus GetStepStatus() const;                  // G4 specific
    TG4Limits*    GetLimitsModifiedOnFly() const;         // G4 specific
    Bool_t   IsCollectTracks() const;
        
        // tracking volume(s) 
    G4VPhysicalVolume* GetCurrentPhysicalVolume() const;  // G4 specific
    TG4Limits* GetCurrentLimits() const;  // G4 specific
    Int_t CurrentVolID(Int_t& copyNo) const;
    Int_t CurrentVolOffID(Int_t off, Int_t& copyNo) const;
    const char* CurrentVolName() const;
    const char* CurrentVolOffName(Int_t off) const;
    const char* CurrentVolPath();
    Bool_t CurrentBoundaryNormal(
                    Double_t &x, Double_t &y, Double_t &z) const;
    Int_t  CurrentMaterial(Float_t &a, Float_t &z, Float_t &dens, 
                    Float_t &radl, Float_t &absl) const;
    Int_t CurrentMedium() const;
    void Gmtod(Double_t* xm, Double_t* xd, Int_t iflag);
    void Gmtod(Float_t* xm, Float_t* xd, Int_t iflag);
    void Gdtom(Double_t* xd, Double_t* xm, Int_t iflag);
    void Gdtom(Float_t* xd, Float_t* xm, Int_t iflag);
    Double_t MaxStep() const;
    Int_t GetMaxNStep() const;

        // tracking particle 
        // dynamic properties
    void TrackPosition(TLorentzVector& position) const;
    void TrackPosition(Double_t& x, Double_t& y, Double_t& z) const;
    void TrackPosition(Float_t& x, Float_t& y, Float_t& z) const;
    void TrackMomentum(TLorentzVector& momentum) const;
    void TrackMomentum(Double_t& px, Double_t& py, Double_t&pz, 
                       Double_t& etot) const;
    void TrackMomentum(Float_t& px, Float_t& py, Float_t&pz,
                       Float_t& etot) const;
    Double_t TrackStep() const;  
    Double_t TrackLength() const;   
    Double_t TrackTime() const;  
    Double_t Edep() const;
        // static properties
    Int_t TrackPid() const;
    Double_t TrackCharge() const;
    Double_t TrackMass() const;
    Double_t Etot() const;

        // track status
    Bool_t IsTrackInside() const;
    Bool_t IsTrackEntering() const;
    Bool_t IsTrackExiting() const;
    Bool_t IsTrackOut() const;
    Bool_t IsTrackDisappeared() const;
    Bool_t IsTrackStop() const;
    Bool_t IsTrackAlive() const;
    Bool_t IsNewTrack() const;

        // secondaries
    Int_t NSecondaries() const;
    void GetSecondary(Int_t index, Int_t& particleId,
                      TLorentzVector& position, TLorentzVector& momentum);      
    TMCProcess ProdProcess(Int_t isec) const; 
    Int_t StepProcesses(TArrayI &proc) const;

  private:
    /// Not implemented
    TG4StepManager(const TG4StepManager& right);
    /// Not implemented
    TG4StepManager& operator=(const TG4StepManager& right);

    // methods
    void CheckTrack() const;
    void CheckStep(const G4String& method) const;
    void CheckGflashSpot(const G4String& method) const;
    void CheckSteppingManager() const;
    void SetTLorentzVector(G4ThreeVector xyz, G4double t, 
                           TLorentzVector& lv) const;    
    const G4VTouchable* GetCurrentTouchable() const; 
    G4VPhysicalVolume*  GetCurrentOffPhysicalVolume(
                           G4int off, G4bool warn = false) const;

    Double_t EdepBoundaryAndKilled() const;

    // static data members
    static G4ThreadLocal TG4StepManager*  fgInstance;   ///< this instance

    
    //
    // data members
    
    /// current track
    G4Track*            fTrack;
    
    /// current step
    G4Step*             fStep;

    /// current Gflash spot
    G4GFlashSpot*       fGflashSpot;
    
    /// \brief step status 
    /// \details that decides whether track properties will be returned from 
    /// PreStepPoint or PostStepPoint
    TG4StepStatus       fStepStatus; 
    
    /// \brief limits which step limit was modified during tracking
    TG4Limits*          fLimitsModifiedOnFly;
    
    /// G4SteppingManager  
    G4SteppingManager*  fSteppingManager;

    /// buffer for current volume name or path
    mutable G4String    fNameBuffer;
    
    /// volume copy number offset
    G4int               fCopyNoOffset;

    /// division copy number offset
    G4int               fDivisionCopyNoOffset;

    /// Cached pointer to thread-local track manager
    TG4TrackManager*    fTrackManager;
};

// inline methods

inline TG4StepManager* TG4StepManager::Instance() { 
  /// Return this instance.
  return fgInstance;
}

inline void TG4StepManager::SetStep(G4Step* step, TG4StepStatus status) { 
  /// Set current step and step status. 
  fTrack = step->GetTrack(); fStep = step; fStepStatus = status; fGflashSpot = 0;
}

inline void TG4StepManager::SetStep(G4Track* track, TG4StepStatus status) { 
  /// Set current track and step status. 
  fTrack = track; fStep = 0; fStepStatus = status;  fGflashSpot = 0;
}

inline void TG4StepManager::SetStep(G4GFlashSpot* gflashSpot, TG4StepStatus status) {
  /// Set current track and step status.
  fTrack = const_cast<G4Track*>(gflashSpot->GetOriginatorTrack()->GetPrimaryTrack());
  fStep = 0; fStepStatus = status;  fGflashSpot = gflashSpot;
}

inline void TG4StepManager::SetSteppingManager(G4SteppingManager* manager) { 
  /// Set G4 stepping manger. 
  fSteppingManager = manager; 

  /// Set navigator !!!
  //G4cout << "SetNavigator:"
  //       << G4TransportationManager::GetTransportationManager()
  //		   ->GetNavigatorForTracking() << G4endl;
  fSteppingManager
    ->SetNavigator(G4TransportationManager::GetTransportationManager()
		   ->GetNavigatorForTracking());
}

inline G4Track* TG4StepManager::GetTrack() const { 
  /// Return current track manger. 
  return fTrack; 
}

inline G4Step* TG4StepManager::GetStep() const { 
  /// Return current step. 
  return fStep; 
}

inline TG4StepStatus TG4StepManager::GetStepStatus() const { 
  /// Return current step status. 
  return fStepStatus; 
}

inline TG4Limits* TG4StepManager::GetLimitsModifiedOnFly() const {
  /// Return limits that has been modified on fly
  return fLimitsModifiedOnFly;
}

//_____________________________________________________________________________
inline Int_t TG4StepManager::GetMaxNStep() const
{
/// Return the maximum number of steps.

  return TG4SteppingAction::Instance()->GetMaxNofSteps();
}

//_____________________________________________________________________________
inline void TG4StepManager::TrackPosition(TLorentzVector& position) const
{
/// Fill the current particle position in the world reference frame
/// and the global time since the event in which the track belongs is created.
/// (position in the PostStepPoint).

#ifdef MCDEBUG
  CheckTrack();
#endif

  G4ThreeVector positionVector;
  if ( fStepStatus == kGflashSpot ) {
    positionVector = fGflashSpot->GetEnergySpot()->GetPosition();
  } else {
    // get position
    // check if this is == to PostStepPoint position !!
    positionVector = fTrack->GetPosition();
  }
  positionVector *= TG4G3Units::InvLength();

  // global time
  G4double time = fTrack->GetGlobalTime();
  time *= TG4G3Units::InvTime();

  SetTLorentzVector(positionVector, time, position);
}

//_____________________________________________________________________________
inline void TG4StepManager::TrackPosition(Double_t& x, Double_t& y, Double_t& z) const
{
/// Fill the current particle position in the world reference frame
/// (position in the PostStepPoint).


#ifdef MCDEBUG
  CheckTrack();
#endif

  G4ThreeVector positionVector;
  if ( fStepStatus == kGflashSpot ) {
    positionVector = fGflashSpot->GetEnergySpot()->GetPosition();
  } else {
  // get position
  // check if this is == to PostStepPoint position !!
    positionVector = fTrack->GetPosition();
  }
  positionVector *= TG4G3Units::InvLength();

  x = positionVector.x();
  y = positionVector.y();
  z = positionVector.z();
}

//_____________________________________________________________________________
inline void TG4StepManager::TrackPosition(Float_t& x, Float_t& y, Float_t& z) const
{
/// Fill the current particle position in the world reference frame
/// (position in the PostStepPoint) as float.

  Double_t dx, dy, dz;
  TrackPosition(dx, dy, dz);

  x = static_cast<float>(dx);
  y = static_cast<float>(dy);
  z = static_cast<float>(dz);
}

//_____________________________________________________________________________
inline void TG4StepManager::TrackMomentum(TLorentzVector& momentum) const
{
/// Fill the current particle momentum (px, py, pz, Etot)
/// Not updated in Gflash fast simulation.

#ifdef MCDEBUG
  CheckTrack();
#endif

  G4ThreeVector momentumVector = fTrack->GetMomentum();
  momentumVector *= TG4G3Units::InvEnergy();

  G4double energy = fTrack->GetDynamicParticle()->GetTotalEnergy();
  energy *= TG4G3Units::InvEnergy();

  SetTLorentzVector(momentumVector, energy, momentum);
}

//_____________________________________________________________________________
inline void TG4StepManager::TrackMomentum(Double_t& px, Double_t& py, Double_t&pz,
                                   Double_t& etot) const
{
/// Fill the current particle momentum
/// Not updated in Gflash fast simulation.

#ifdef MCDEBUG
  CheckTrack();
#endif

  G4ThreeVector momentumVector = fTrack->GetMomentum();
  momentumVector *= TG4G3Units::InvEnergy();

  px = momentumVector.x();
  py = momentumVector.y();
  pz = momentumVector.z();

  etot = fTrack->GetDynamicParticle()->GetTotalEnergy();
  etot *= TG4G3Units::InvEnergy();
}

//_____________________________________________________________________________
inline void TG4StepManager::TrackMomentum(Float_t& px, Float_t& py, Float_t&pz,
                                   Float_t& etot) const
{
/// Fill the current particle momentum as float.
/// Not updated in Gflash fast simulation.

  Double_t dpx, dpy, dpz, detot;
  TrackMomentum(dpx, dpy, dpz, detot);

  px = static_cast<float>(dpx);
  py = static_cast<float>(dpy);
  pz = static_cast<float>(dpz);
  etot = static_cast<float>(detot);
}

//_____________________________________________________________________________
inline Double_t TG4StepManager::TrackStep() const
{
/// Return the current step length.
/// Not updated in Gflash fast simulation.

  if ( fStepStatus == kNormalStep ) {
#ifdef MCDEBUG
    CheckStep("TrackStep");
#endif
    return fStep->GetStepLength()*TG4G3Units::InvLength();
  }
  else
    return 0;
}

//_____________________________________________________________________________
inline Double_t TG4StepManager::TrackLength() const
{
/// Return the length of the current track from its origin.
/// Not updated in Gflash fast simulation.

#ifdef MCDEBUG
  CheckTrack();
#endif

  return fTrack->GetTrackLength()*TG4G3Units::InvLength();
}

//_____________________________________________________________________________
inline Double_t TG4StepManager::TrackTime() const
{
/// Return the global track time = time since the event in which
/// the track belongs is created.                                           \n
/// Note that in Geant4: there is also defined proper time as
/// the proper time of the dynamical particle of the current track.
/// Not updated in Gflash fast simulation.

#ifdef MCDEBUG
  CheckTrack();
#endif

  return fTrack->GetGlobalTime()*TG4G3Units::InvTime();
}

//_____________________________________________________________________________
inline Double_t TG4StepManager::Edep() const
{
/// Return the total energy deposit in this step.

  if ( fStepStatus == kNormalStep ) {

#ifdef MCDEBUG
    CheckStep("Edep");
#endif

    return fStep->GetTotalEnergyDeposit()*TG4G3Units::InvEnergy();
  }

  if ( fStepStatus == kBoundary &&
       fTrack->GetTrackStatus() == fStopAndKill ) {
	// in this case dispatch to a specialized (non-inlined) function
	return EdepBoundaryAndKilled();
  }

  if ( fStepStatus == kGflashSpot ) {

#ifdef MCDEBUG
    CheckGflashSpot("Edep");
#endif

    return fGflashSpot->GetEnergySpot()->GetEnergy()*TG4G3Units::InvEnergy();
  }

  return 0;
}

//_____________________________________________________________________________
inline Int_t TG4StepManager::TrackPid() const
{
/// Return the current particle PDG encoding.

#ifdef MCDEBUG
  CheckTrack();
#endif

  G4ParticleDefinition* particle
    = fTrack->GetDynamicParticle()->GetDefinition();

  // Ask TG4ParticlesManager to get PDG encoding
  // (in order to get PDG from extended TDatabasePDG
  // in case the standard PDG code is not defined)
  G4int pdgEncoding
    = TG4ParticlesManager::Instance()->GetPDGEncoding(particle);

  // Make difference between optical photon from Cerenkov and
  // feedback photon generated by user
  if ( pdgEncoding == 50000050 ) {
    TG4TrackInformation* trackInformation
      = fTrackManager->GetTrackInformation(fTrack);
    if ( trackInformation && trackInformation->GetPDGEncoding() )
      pdgEncoding = trackInformation->GetPDGEncoding();
  }

  return pdgEncoding;
}

//_____________________________________________________________________________
inline Double_t TG4StepManager::TrackCharge() const
{
/// Return the current particle charge.

#ifdef MCDEBUG
  CheckTrack();
#endif

  return fTrack->GetDynamicParticle()->GetDefinition()
           ->GetPDGCharge()*TG4G3Units::InvCharge();
}

//_____________________________________________________________________________
inline Double_t TG4StepManager::TrackMass() const
{
/// Return the current particle mass at rest.

#ifdef MCDEBUG
  CheckTrack();
#endif

  return fTrack->GetDynamicParticle()->GetDefinition()
           ->GetPDGMass()*TG4G3Units::InvMass();
}

//_____________________________________________________________________________
inline Double_t TG4StepManager::Etot() const
{
/// Return the total energy of the current particle.

#ifdef MCDEBUG
  CheckTrack();
#endif

  return fTrack->GetDynamicParticle()->GetTotalEnergy()*TG4G3Units::InvEnergy();
}

// TO DO: revise these with added kGflashSpot status

//_____________________________________________________________________________
inline Bool_t TG4StepManager::IsTrackInside() const
{
/// Return true if the particle does not cross a geometrical boundary
/// and is not in the vertex.

  if ( fStepStatus == kNormalStep  && ! ( IsTrackExiting() ) ) {
    // track is always inside during a normal step
    return true;
  }

  return false;
}

//_____________________________________________________________________________
inline Bool_t TG4StepManager::IsTrackEntering() const
{
/// Return true if the particle crosses a geometrical boundary
/// or is in the vertex.

  if ( fStepStatus != kNormalStep ) {
    // track is entering during a vertex or boundary step
    return true;
  }

  return false;
}

//_____________________________________________________________________________
inline Bool_t TG4StepManager::IsTrackExiting() const
{
/// Return true if the particle crosses a geometrical boundary.

  if (fStepStatus == kNormalStep) {

#ifdef MCDEBUG
    CheckStep("IsTrackExiting");
#endif

    if (fStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary)
       return true;
  }

  return false;
}

//_____________________________________________________________________________
inline Bool_t TG4StepManager::IsTrackOut() const
{
/// Return true if the particle crosses the world boundary
/// at the post-step point.

  if ( fStepStatus == kVertex ) return false;

#ifdef MCDEBUG
  CheckStep("IsTrackOut");
#endif

  if ( fStep->GetPostStepPoint()->GetStepStatus() == fWorldBoundary )
    return true;
  else
    return false;
}

//_____________________________________________________________________________
inline Bool_t TG4StepManager::IsTrackStop() const
{
/// Return true if the particle has stopped
/// or has been killed, suspended or postponed to the next event.
///
/// Possible track status from G4:
///   - fAlive,              // Continue the tracking
///   - fStopButAlive,       // Invoke active rest physics processes and
///                          // and kill the current track afterward
///   - fStopAndKill,        // Kill the current track
///   - fKillTrackAndSecondaries, // Kill the current track and also associated
///                          // secondaries.
///   - fSuspend,            // Suspend the current track
///   - fPostponeToNextEvent // Postpones the tracking of thecurrent track
///                          // to the next event.

#ifdef MCDEBUG
  CheckTrack();
#endif

  // check
  G4TrackStatus status
     = fTrack->GetTrackStatus();
  if ( ( status == fStopAndKill ) ||
       ( status == fKillTrackAndSecondaries ) ||
       ( status == fSuspend ) ||
       ( status == fPostponeToNextEvent ) ) {
    return true;
  }
  else
    return false;
}

//_____________________________________________________________________________
inline Bool_t TG4StepManager::IsTrackDisappeared() const
{
/// Return true if particle has disappeared
/// (due to any physical process)
/// or has been killed or postponed to next event.

#ifdef MCDEBUG
  CheckTrack();
#endif

  // check
  G4TrackStatus status
     = fTrack->GetTrackStatus();
  if ( ( status == fStopAndKill ) ||
       ( status == fKillTrackAndSecondaries ) ||
       ( status == fPostponeToNextEvent ) ) {
    return true;
  }
  else
    return false;
}

//_____________________________________________________________________________
inline Bool_t TG4StepManager::IsTrackAlive() const
{
/// Return true if particle continues tracking.

#ifdef MCDEBUG
  CheckTrack();
#endif

  G4TrackStatus status
     = fTrack->GetTrackStatus();
  if ( (status == fAlive) ||
       (status == fStopButAlive) )
    return true;
  else
    return false;
}

//_____________________________________________________________________________
inline Bool_t TG4StepManager::IsNewTrack() const
{
/// Return true when the track performs the first step.

  if ( fStepStatus == kVertex )
    return true;
  else
    return false;
}

//_____________________________________________________________________________
inline Int_t TG4StepManager::NSecondaries() const
{
/// Return the number of secondary particles generated
/// in the current step.

  if ( fStepStatus == kVertex || fStepStatus == kGflashSpot ) return 0;

#ifdef MCDEBUG
  CheckSteppingManager();
#endif

  G4int nofSecondaries = 0;
  nofSecondaries += fSteppingManager->GetfN2ndariesAtRestDoIt();
  nofSecondaries += fSteppingManager->GetfN2ndariesAlongStepDoIt();
  nofSecondaries += fSteppingManager->GetfN2ndariesPostStepDoIt();

  return nofSecondaries;
}

//_____________________________________________________________________________
inline void TG4StepManager::SetTLorentzVector(G4ThreeVector xyz, G4double t,
                                              TLorentzVector& lv) const
{
  /// Fill TLorentzVector with G4ThreeVector and G4double.

  lv[0] = xyz.x();
  lv[1] = xyz.y();
  lv[2] = xyz.z();
  lv[3] = t;
}

#endif //TG4_STEP_MANAGER_H

