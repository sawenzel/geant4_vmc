//------------------------------------------------
// The Geant4 Virtual Monte Carlo package
// Copyright (C) 2007 - 2015 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file TG4StepManager.cxx
/// \brief Implementation of the TG4StepManager class 
///
/// \author I. Hrivnacova; IPN, Orsay

#include "TG4StepManager.h"
#include "TG4SteppingAction.h"
#include "TG4GeometryServices.h"
#include "TG4SDServices.h"
#include "TG4ParticlesManager.h"
#include "TG4PhysicsManager.h"
#include "TG4TrackManager.h"
#include "TG4TrackInformation.h"
#include "TG4Limits.h"
#include "TG4Globals.h"
#include "TG4G3Units.h"

#include <G4SteppingManager.hh>
#include <G4UserLimits.hh>
#include <G4UImanager.hh>
#include <G4AffineTransform.hh>
#include <G4TransportationManager.hh>
#include <G4Navigator.hh>
#include <G4VProcess.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessVector.hh>
#include <G4OpticalPhoton.hh>
#include <G4VTouchable.hh>

#include <TLorentzVector.h>
#include <TMath.h>

G4ThreadLocal TG4StepManager* TG4StepManager::fgInstance = 0;

//_____________________________________________________________________________
TG4StepManager::TG4StepManager(const TString& userGeometry) 
  : fTrack(0),
    fStep(0),
    fGflashSpot(0),
    fStepStatus(kNormalStep),
    fLimitsModifiedOnFly(0),
    fSteppingManager(0),
    fNameBuffer(),
    fCopyNoOffset(0),
    fDivisionCopyNoOffset(0),
    fTrackManager(0)
{
/// Standard constructor
/// \param userGeometry  User selection of geometry definition and navigation 

  //G4cout << "TG4StepManager::TG4StepManager " << this << G4endl;
  
  if ( fgInstance ) {
    TG4Globals::Exception(
      "TG4StepManager", "TG4StepManager", 
      "Cannot create two instances of singleton.");
  }
      
  fgInstance = this;  
  
  /// Set offset for passing copyNo to 1;
  /// as G3toG4 decrement copyNo passed by user by 1
  if ( userGeometry == "VMCtoGeant4") fCopyNoOffset = 1;

  /// Set offset for passing copyNo to 1;
  /// to be equivalent to Root geometrical model
  /// (Root starts numbering from 1, while Geant4 from 0)
  if ( userGeometry == "RootToGeant4"  || userGeometry == "Geant4") 
    fDivisionCopyNoOffset = 1;
}

//_____________________________________________________________________________
TG4StepManager::~TG4StepManager() 
{
/// Destructor

  fgInstance = 0;
}

//
// private methods
//

//_____________________________________________________________________________
void TG4StepManager::CheckTrack() const
{
/// Give exception in case the track is not defined.

  if ( ! fTrack ) 
    TG4Globals::Exception(
      "TG4StepManager", "CheckTrack", "Track is not defined.");
}     


//_____________________________________________________________________________
void TG4StepManager::CheckStep(const G4String& method) const
{
/// Give exception in case the step is not defined.

  if ( ! fStep ) {
    TG4Globals::Exception(
      "TG4StepManager", method, "Step is not defined.");
  }
}     

//_____________________________________________________________________________
void TG4StepManager::CheckGflashSpot(const G4String& method) const
{
/// Give exception in case the step is not defined.

  if ( ! fGflashSpot ) {
    TG4Globals::Exception(
      "TG4StepManager", method, "Gflash spot is not defined.");
  }
}

//_____________________________________________________________________________
void TG4StepManager::CheckSteppingManager() const
{
/// Give exception in case the step is not defined.

  if ( ! fSteppingManager ) 
    TG4Globals::Exception(
      "TG4StepManager", "CheckSteppingManager", 
      "Stepping manager is not defined.");
}     


//_____________________________________________________________________________
const G4VTouchable* TG4StepManager::GetCurrentTouchable() const 
{
/// Return the current touchable.
 
#ifdef MCDEBUG
    CheckTrack();
#endif    

  if ( fStepStatus == kGflashSpot ) {
    G4ReferenceCountedHandle<G4VTouchable> touchableHandle
      = fGflashSpot->GetTouchableHandle();
    return touchableHandle();
  }
  else if ( fStepStatus != kBoundary )
    return fTrack->GetTouchable();
  else
    return fTrack->GetNextTouchable();
}  

//_____________________________________________________________________________
G4VPhysicalVolume* 
TG4StepManager::GetCurrentOffPhysicalVolume(G4int off, G4bool warn) const 
{
/// Return the physical volume of the off-th mother
/// of the current volume.
 
  // Get current touchable
  //
  const G4VTouchable* touchable = GetCurrentTouchable(); 

  // Check touchable depth
  //
  if ( touchable->GetHistoryDepth() < off ) {
    if ( warn ) {
      TString text = "level=";
      text += off;
      TG4Globals::Warning(
        "TG4StepManager", "GetCurrentOffPhysicalVolume", 
        "Volume " + TString(touchable->GetVolume()->GetName()) +
        " has not defined mother in " + text + ".");  
    }    
    return 0;
  }  

  return touchable->GetVolume(off);
}     

//
// public methods
//

//_____________________________________________________________________________
void TG4StepManager::LateInitialize()
{
  fTrackManager = TG4TrackManager::Instance();
}

//_____________________________________________________________________________
void TG4StepManager::StopTrack()
{
/// Stop the current track and skips to the next.                          \n
/// Possible "stop" track status from G4:
///  - fStopButAlive       // Invoke active rest physics processes and
///                        // and kill the current track afterward
///  - fStopAndKill        // Kill the current track (used)
///  - fKillTrackAndSecondaries  // Kill the current track and also associated
///                       // secondaries.

  if ( fTrack ) {
    fTrack->SetTrackStatus(fStopAndKill);
    // fTrack->SetTrackStatus(fStopButAlive);
    // fTrack->SetTrackStatus(fKillTrackAndSecondaries);
  }
  else {
    TG4Globals::Warning(
      "TG4StepManager", "StopTrack()",
      "There is no current track to be stopped.");
  }    
}

//_____________________________________________________________________________
void TG4StepManager::StopEvent()
{
/// Abort the current event processing.

  if ( fTrack ) {
    fTrack->SetTrackStatus(fKillTrackAndSecondaries);
            //StopTrack();   // cannot be used as it keeps secondaries
  }
              
  G4UImanager::GetUIpointer()->ApplyCommand("/event/abort");
}

//_____________________________________________________________________________
void TG4StepManager::StopRun()
{
/// Abort the current run processing.

  TG4SDServices::Instance()->SetIsStopRun(true);

  StopEvent();
  G4UImanager::GetUIpointer()->ApplyCommand("/run/abort");
}

//_____________________________________________________________________________
void TG4StepManager::SetMaxStep(Double_t step)
{
/// Set the maximum step allowed in the current logical volume;
/// the value is restored after exiting from the current tracking
/// medium

  TG4Limits* userLimits = GetCurrentLimits();

  if ( ! userLimits )  return;

  //G4cout << "TG4StepManager::SetMaxStep  in " 
  //       << GetCurrentPhysicalVolume()->GetLogicalVolume()->GetName() << "  "
  //       << userLimits->GetName() << G4endl;

  // set max step
  userLimits->SetCurrentMaxAllowedStep(step*TG4G3Units::Length());
  fLimitsModifiedOnFly = userLimits;
}

//_____________________________________________________________________________
void TG4StepManager::SetMaxStepBack()
{
/// Restore back the maximum step after exiting from the tracking
/// medium where it has been changed of fly

  if ( ! fLimitsModifiedOnFly ) {
    TG4Globals::Warning(
      "TG4StepManager", "SetMaxStepBack", 
      "No limits modified on fly found.");
    return;  
  }

  // set max step
  fLimitsModifiedOnFly->SetMaxAllowedStepBack();
  fLimitsModifiedOnFly = 0;
}

//_____________________________________________________________________________
void TG4StepManager::SetMaxNStep(Int_t maxNofSteps)
{
/// Set the maximum number of steps.

  TG4SteppingAction::Instance()->SetMaxNofSteps(TMath::Abs(maxNofSteps));
}

//_____________________________________________________________________________
void TG4StepManager::SetCollectTracks(Bool_t collectTracks)
{
/// (In)Activate collecting TGeo tracks 

  TG4SteppingAction::Instance()->SetCollectTracks(collectTracks);
}  

//_____________________________________________________________________________
void TG4StepManager::ForceDecayTime(Float_t time)
{
/// Force decay time.                                                        \n
/// Not yet implemented.

#ifdef MCDEBUG
  CheckTrack();
#endif

  G4ParticleDefinition* particle
    = fTrack->GetDynamicParticle()->GetDefinition();
  
  // Store the original particle lifetime in track information
  // (as it has to be set back after track is finished)
  TG4TrackInformation* trackInformation
    =  fTrackManager->GetTrackInformation(fTrack);
  trackInformation->SetPDGLifetime(particle->GetPDGLifeTime()); 
    
  // Set new lifetime value
  particle->SetPDGLifeTime(time*TG4G3Units::Time());
}

//_____________________________________________________________________________
Bool_t  TG4StepManager::IsCollectTracks() const
{
/// Return the info if collecting tracks is activated

  return TG4SteppingAction::Instance()->GetCollectTracks();
}  

//_____________________________________________________________________________
G4VPhysicalVolume* TG4StepManager::GetCurrentPhysicalVolume() const 
{
/// Return the current physical volume.
/// According to fStepStatus the volume from track vertex,
/// pre step point or post step point is returned.

#ifdef MCDEBUG
  CheckTrack();
#endif 
 
  if ( fStepStatus == kGflashSpot )
    return fGflashSpot->GetTouchableHandle()->GetVolume();
  else if ( fStepStatus != kBoundary )
    return fTrack->GetVolume();
  else
    return fTrack->GetNextVolume();
}

//_____________________________________________________________________________
TG4Limits* TG4StepManager::GetCurrentLimits() const 
{
/// Return the current limits.

#ifdef MCDEBUG
  TG4Limits* userLimits 
     = TG4GeometryServices::Instance()
         ->GetLimits(GetCurrentPhysicalVolume()->GetLogicalVolume()->GetUserLimits()); 
#else  
  TG4Limits* userLimits 
    = (TG4Limits*) GetCurrentPhysicalVolume()->GetLogicalVolume()->GetUserLimits();
#endif    

  if ( ! userLimits ) {
    TG4Globals::Warning(
      "TG4StepManager", "Get current limits", "User limits not defined.");
    return 0;  
  }
  
  return userLimits;  
}

//_____________________________________________________________________________
Int_t TG4StepManager::CurrentVolID(Int_t& copyNo) const
{
/// Return the current sensitive detector ID
/// and fill the copy number of the current physical volume 

  G4VPhysicalVolume* physVolume = GetCurrentPhysicalVolume();
  if ( ! physVolume ) {
    TG4Globals::Exception(
      "TG4StepManager", "CurrentVolID", "No current physical volume found");
    return 0;  
  }
  copyNo = physVolume->GetCopyNo() + fCopyNoOffset;
  
  if ( physVolume->IsParameterised() ||
       physVolume->IsReplicated() )  copyNo += fDivisionCopyNoOffset;

  // sensitive detector ID
  return TG4SDServices::Instance()->GetVolumeID(physVolume->GetLogicalVolume());
} 

//_____________________________________________________________________________
Int_t TG4StepManager::CurrentVolOffID(Int_t off, Int_t&  copyNo) const
{ 
/// Return the  the sensitive detector ID of the off-th mother of the current  
/// volume and  fill the copy number of its physical volume 

  if ( off == 0 ) return CurrentVolID(copyNo);
#ifdef MCDEBUG
   G4VPhysicalVolume* mother = GetCurrentOffPhysicalVolume(off, true); 
#else
   G4VPhysicalVolume* mother = GetCurrentOffPhysicalVolume(off); 
#endif   

  if ( mother ) {
    copyNo = mother->GetCopyNo() + fCopyNoOffset;

    if ( mother->IsParameterised() ||
         mother->IsReplicated() )  copyNo += fDivisionCopyNoOffset;

    // sensitive detector ID
    return TG4SDServices::Instance()->GetVolumeID(mother->GetLogicalVolume());
  }
  else {
    copyNo = 0;
    return 0;
  }  
}

//_____________________________________________________________________________
const char* TG4StepManager::CurrentVolName() const
{
/// Return the current physical volume name.

  fNameBuffer
    = TG4GeometryServices::Instance()->UserVolumeName(
        GetCurrentPhysicalVolume()->GetLogicalVolume()->GetName());

  return fNameBuffer.data();
}

//_____________________________________________________________________________
const char* TG4StepManager::CurrentVolOffName(Int_t off) const
{ 
/// Return the off-th mother's physical volume name.

  if ( off == 0) return CurrentVolName();

  G4VPhysicalVolume* mother = GetCurrentOffPhysicalVolume(off); 

  if ( mother ) {
    fNameBuffer
      = TG4GeometryServices::Instance()->UserVolumeName(
          mother->GetLogicalVolume()->GetName());
  }             
  else {
    fNameBuffer = "";
  }
  return fNameBuffer.data();
}

//_____________________________________________________________________________
const char* TG4StepManager::CurrentVolPath()
{ 
/// Return the current volume path.

  TG4GeometryServices* geometryServices = TG4GeometryServices::Instance();

  // Get current touchable
  const G4VTouchable* touchable = GetCurrentTouchable();
  
  // Check touchable depth
  //
  G4int depth = touchable->GetHistoryDepth();
  
  // Compose the path
  //
  fNameBuffer = "";
  for ( G4int i=0; i<depth; i++ ) {
    G4VPhysicalVolume* physVolume = touchable->GetHistory()->GetVolume(i);
    fNameBuffer += "/";
    fNameBuffer
      += geometryServices->UserVolumeName(physVolume->GetName());
    fNameBuffer += "_";
    TG4Globals::AppendNumberToString(fNameBuffer, physVolume->GetCopyNo());
  }     

  // Add current volume to the path
  G4VPhysicalVolume* curPhysVolume = GetCurrentPhysicalVolume(); 
  fNameBuffer += "/";
  fNameBuffer += geometryServices->UserVolumeName(curPhysVolume->GetName());
  fNameBuffer += "_";
  TG4Globals::AppendNumberToString(fNameBuffer, curPhysVolume->GetCopyNo());

  return fNameBuffer.data();
}

//_____________________________________________________________________________
Bool_t TG4StepManager::CurrentBoundaryNormal(
                              Double_t& x, Double_t& y, Double_t& z) const
{
/// Return the he normal vector of the surface of the last volume exited

   G4Navigator* theNavigator 
     = G4TransportationManager::GetTransportationManager()
       ->GetNavigatorForTracking();

   G4bool valid;
   G4ThreeVector theLocalNormal 
     = theNavigator->GetLocalExitNormal(&valid);
   if ( ! valid ) return false;
   
   G4ThreeVector theGlobalNormal 
     = theNavigator->GetLocalToGlobalTransform().TransformAxis(theLocalNormal);
     
   x = theGlobalNormal.x();
   y = theGlobalNormal.y();
   z = theGlobalNormal.z();
   
   return true;
}    

//_____________________________________________________________________________
Int_t TG4StepManager::CurrentMaterial(Float_t &a, Float_t &z, Float_t &dens, 
                          Float_t &radl, Float_t &absl) const
{
/// Get parameters of the current material material during transport.
/// Return the number of elements in the mixture
/// \param a     The atomic mass in au
/// \param z     The atomic number
/// \param dens  The density in g/cm3
/// \param radl  The radiation length in cm
/// \param absl  The absorption length in cm

  G4VPhysicalVolume* physVolume = GetCurrentPhysicalVolume(); 
    
  G4Material* material 
    = physVolume->GetLogicalVolume()->GetMaterial();

  G4int nofElements = material->GetNumberOfElements();
  TG4GeometryServices* geometryServices = TG4GeometryServices::Instance();
  a = geometryServices->GetEffA(material);
  z = geometryServices->GetEffZ(material);
      
  // density 
  dens = material->GetDensity();
  dens /= TG4G3Units::MassDensity();      
      
  // radiation length
  radl = material->GetRadlen();
  radl /= TG4G3Units::Length();
      
  absl = 0.;  // this parameter is not defined in Geant4
  return nofElements;
}

//_____________________________________________________________________________
Int_t TG4StepManager::CurrentMedium() const
{   
/// Return the medium ID 

  return TG4SDServices::Instance()->GetMediumID(
            GetCurrentPhysicalVolume()->GetLogicalVolume());
}

//_____________________________________________________________________________
void TG4StepManager::Gmtod(Float_t* xm, Float_t* xd, Int_t iflag) 
{ 
/// Transform a position from the world reference frame
/// to the current volume reference frame.
/// \param xm    Known coordinates in the world reference system
/// \param xd    Computed coordinates in the daughter reference system
/// \param iflag The option: 
///              - IFLAG=1  convert coordinates,                                 \n
///              - IFLAG=2  convert direction cosinus
///
 

  G4double* dxm = TG4GeometryServices::Instance()->CreateG4doubleArray(xm, 3);
  G4double* dxd = TG4GeometryServices::Instance()->CreateG4doubleArray(xd, 3);

  Gmtod(dxm, dxd, iflag);

  for ( G4int i=0; i<3; i++ ) {
    xm[i] = dxm[i]; 
    xd[i] = dxd[i];
  }   

  delete [] dxm;
  delete [] dxd;
} 
 
//_____________________________________________________________________________
void TG4StepManager::Gmtod(Double_t* xm, Double_t* xd, Int_t iflag) 
{ 
/// Transform a position from the world reference frame
/// to the current volume reference frame.
/// \param xm    Known coordinates in the world reference system
/// \param xd    Computed coordinates in the daughter reference system
/// \param iflag The option: 
///              - IFLAG=1  convert coordinates,                                 \n
///              - IFLAG=2  convert direction cosinus
///

#ifdef MCDEBUG
  if ( iflag != 1 && iflag != 2 ) {
      TString text = "iflag=";
      text += iflag;
      TG4Globals::Warning(
        "TG4StepManager", "Gmtod", text + " is different from 1..2.");
      return;        
  }        
#endif

  const G4AffineTransform& affineTransform
    = GetCurrentTouchable()->GetHistory()->GetTopTransform();

  G4ThreeVector theGlobalPoint(xm[0]* TG4G3Units::Length(),
                               xm[1]* TG4G3Units::Length(),                       
                               xm[2]* TG4G3Units::Length()); 
  G4ThreeVector theLocalPoint;
  if ( iflag == 1 ) 
    theLocalPoint = affineTransform.TransformPoint(theGlobalPoint);
  else {
    // if ( iflag == 2)
    theLocalPoint = affineTransform.TransformAxis(theGlobalPoint);
  }  

  xd[0] = theLocalPoint.x()/TG4G3Units::Length();
  xd[1] = theLocalPoint.y()/TG4G3Units::Length();
  xd[2] = theLocalPoint.z()/TG4G3Units::Length();
} 
 
//_____________________________________________________________________________
void TG4StepManager::Gdtom(Float_t* xd, Float_t* xm, Int_t iflag) 
{ 
/// Transform a position from the current volume reference frame
/// to the world reference frame.
/// \param xd    Known coordinates in the daughter reference system
/// \param xm    Computed coordinates in the world reference system
/// \param iflag The option: 
///              - IFLAG=1  convert coordinates,                                 \n
///              - IFLAG=2  convert direction cosinus


  G4double* dxd = TG4GeometryServices::Instance()->CreateG4doubleArray(xd, 3);
  G4double* dxm = TG4GeometryServices::Instance()->CreateG4doubleArray(xm, 3);

  Gdtom(dxd, dxm, iflag);

  for ( G4int i=0; i<3; i++ ) {
    xd[i] = dxd[i];
    xm[i] = dxm[i]; 
  }   

  delete [] dxd;
  delete [] dxm;
} 
 
//_____________________________________________________________________________
void TG4StepManager::Gdtom(Double_t* xd, Double_t* xm, Int_t iflag) 
{ 
/// Transform a position from the current volume reference frame
/// to the world reference frame.
/// \param xd    Known coordinates in the daughter reference system
/// \param xm    Computed coordinates in the world reference system
/// \param iflag The option: 
///              - IFLAG=1  convert coordinates,                                 \n
///              - IFLAG=2  convert direction cosinus

#ifdef MCDEBUG
  if ( iflag != 1 && iflag != 2 ) {
      TString text = "iflag=";
      text += iflag;
      TG4Globals::Warning(
        "TG4StepManager", "Gmtod", text + " is different from 1..2.");
      return;        
  }        
#endif

  const G4AffineTransform& affineTransform
    = GetCurrentTouchable()->GetHistory()->GetTopTransform().Inverse();

  G4ThreeVector theLocalPoint(xd[0]*TG4G3Units::Length(),
                              xd[1]*TG4G3Units::Length(),
                              xd[2]*TG4G3Units::Length()); 
  G4ThreeVector theGlobalPoint;
  if ( iflag == 1 )
    theGlobalPoint = affineTransform.TransformPoint(theLocalPoint);
  else {
    // if( iflag == 2)
    theGlobalPoint = affineTransform.TransformAxis(theLocalPoint);
  }    

  xm[0] = theGlobalPoint.x()/TG4G3Units::Length();
  xm[1] = theGlobalPoint.y()/TG4G3Units::Length();
  xm[2] = theGlobalPoint.z()/TG4G3Units::Length();
} 
 
//_____________________________________________________________________________
Double_t TG4StepManager::MaxStep() const
{   
/// Return the maximum step allowed in the current logical volume
/// by user limits.

  G4LogicalVolume* curLogVolume
    = GetCurrentPhysicalVolume()->GetLogicalVolume();

  // check this
  G4UserLimits* userLimits 
    = curLogVolume->GetUserLimits();

  G4double maxStep;
  if ( userLimits == 0 ) { 
    TG4Globals::Warning(
      "TG4StepManager", "MaxStep",
      "User Limits are not defined for the current logical volume " + 
       TString(curLogVolume->GetName()) + ".");
    return FLT_MAX;
  }
  else { 
    const G4Track& trackRef = *(fTrack);
    maxStep = userLimits->GetMaxAllowedStep(trackRef); 
    maxStep /= TG4G3Units::Length(); 
    return maxStep;
  }  
}


//_____________________________________________________________________________
void TG4StepManager::GetSecondary(Int_t index, Int_t& particleId, 
                          TLorentzVector& position, TLorentzVector& momentum)
{
/// Fill the parameters of the generated secondary particle 
/// !! Check if indexing of secondaries is same !!
/// \param index      The secondary particle index
/// \param particleId The PDG encoding
/// \param position   The position 
/// \param momentum   The momentum 

#ifdef MCDEBUG
  CheckSteppingManager();
#endif

  G4int nofSecondaries = NSecondaries();
  if ( ! nofSecondaries ) return;

  const G4TrackVector* secondaryTracks = fSteppingManager->GetSecondary();
#ifdef MCDEBUG
  if ( ! secondaryTracks ) {
    TG4Globals::Exception(
      "TG4StepManager", "GetSecondary", "Secondary tracks vector is empty");
  }
  
  if ( index >= nofSecondaries ) {
    TG4Globals::Exception(
      "TG4StepManager", "GetSecondary", "Wrong secondary track index.");
  }
#endif
  
  // the index of the first secondary of this step
  G4int startIndex 
    = secondaryTracks->size() - nofSecondaries;
         // (the secondaryTracks vector contains secondaries 
         // produced by the track at previous steps, too)
  G4Track* track 
    = (*secondaryTracks)[startIndex + index]; 
   
  // particle encoding
  particleId 
    = track->GetDynamicParticle()->GetDefinition()->GetPDGEncoding();
 
  // position & time
  G4ThreeVector positionVector = track->GetPosition();
  positionVector *= 1./(TG4G3Units::Length());
  G4double time = track->GetGlobalTime();
  time /= TG4G3Units::Time();
  SetTLorentzVector(positionVector, time, position);

  // momentum & energy
  G4ThreeVector momentumVector = track->GetMomentum();        
  momentumVector *= 1./(TG4G3Units::Energy());   
  G4double energy = track->GetDynamicParticle()->GetTotalEnergy();
  energy /= TG4G3Units::Energy();
  SetTLorentzVector(momentumVector, energy, momentum);
}

//_____________________________________________________________________________
TMCProcess TG4StepManager::ProdProcess(Int_t isec) const
{
/// Return the VMC code of the process that has produced the secondary particle
/// specified by its index

  G4int nofSecondaries = NSecondaries();
  if ( fStepStatus == kVertex || ! nofSecondaries ) return kPNoProcess;

#ifdef MCDEBUG
  CheckStep("ProdProcess");
#endif

  const G4TrackVector* secondaryTracks = fSteppingManager->GetSecondary();
 
#ifdef MCDEBUG
  // should never happen
  if ( ! secondaryTracks ) {
    TG4Globals::Exception(
      "TG4StepManager", "ProdProcess", "Secondary tracks vector is empty.");

    return kPNoProcess;  
  }    

  if ( isec >= nofSecondaries ) {
    TG4Globals::Exception(
      "TG4StepManager", "ProdProcess", "Wrong secondary track index.");

    return kPNoProcess;  
  }
#endif

  // the index of the first secondary of this step
  G4int startIndex 
    = secondaryTracks->size() - nofSecondaries;
         // the secondaryTracks vector contains secondaries 
         // produced by the track at previous steps, too

  // the secondary track with specified isec index
  G4Track* track = (*secondaryTracks)[startIndex + isec]; 
   
  const G4VProcess* kpProcess = track->GetCreatorProcess(); 
  
  TMCProcess mcProcess 
   = TG4PhysicsManager::Instance()->GetMCProcess(kpProcess);
  
  // distinguish kPDeltaRay from kPEnergyLoss  
  if (mcProcess == kPEnergyLoss) mcProcess = kPDeltaRay;
  
  return mcProcess;
}

//_____________________________________________________________________________
Int_t TG4StepManager::StepProcesses(TArrayI& processes) const
{
/// Fill the array of processes that were active in the current step.
/// The array is filled with the process VMC codes (TMCProcess).
/// Return the number of active processes    
/// (TBD: Distinguish between kPDeltaRay and kPEnergyLoss)

 if ( fStepStatus == kVertex || fStepStatus == kBoundary || fStepStatus == kGflashSpot ) {
   G4int nofProcesses = 1;
   processes.Set(nofProcesses);
   processes[0] = kPNull;
   return nofProcesses;
 }  
   
#ifdef MCDEBUG
  CheckSteppingManager();
  CheckStep("StepProcesses");
#endif

  // along step processes
  G4ProcessVector* processVector 
    = fStep->GetTrack()->GetDefinition()->GetProcessManager()
        ->GetAlongStepProcessVector();
  G4int nofAlongStep = processVector->entries();
  
  // process defined step
  const G4VProcess* kpLastProcess 
    = fStep->GetPostStepPoint()->GetProcessDefinedStep();

  // set array size
  processes.Set(nofAlongStep+2);
     // maximum number of processes:
     // nofAlongStep (along step) - 1 (transportations) + 1 (post step process)
     // + possibly 2 (additional processes if OpBoundary )
     // => nofAlongStep + 2
 
  // fill array with (nofAlongStep-1) along step processes 
  TG4PhysicsManager* physicsManager = TG4PhysicsManager::Instance();
  G4int counter = 0;  
  for ( G4int i=0; i<nofAlongStep; i++ ) {
    G4VProcess* g4Process = (*processVector)[i];    
    // do not fill transportation along step process
    if ( g4Process && g4Process->GetProcessName() != "Transportation" )
      processes[counter++] = physicsManager->GetMCProcess(g4Process);
  }
    
  // fill array with optical photon information
  if ( fStep->GetTrack()->GetDefinition() == G4OpticalPhoton::Definition() &&
       kpLastProcess->GetProcessName() == "Transportation" &&
       physicsManager->IsOpBoundaryProcess() ) {
       
     // add light scattering anbd reflection/absorption as additional processes
     processes[counter++] = kPLightScattering;
     processes[counter++] = physicsManager->GetOpBoundaryStatus();
  }        

  // fill array with last process
  processes[counter++] = physicsManager->GetMCProcess(kpLastProcess);

  return counter;  
}

Double_t TG4StepManager::EdepBoundaryAndKilled() const {
  G4VProcess* proc = fSteppingManager->GetfCurrentProcess();
  TG4PhysicsManager* physicsManager = TG4PhysicsManager::Instance();
  if (proc && physicsManager->GetMCProcess(proc) == kPLightScattering &&
      physicsManager->GetOpBoundaryStatus() == kPLightDetection) {
    return fTrack->GetTotalEnergy() * TG4G3Units::InvEnergy();
  }
}
