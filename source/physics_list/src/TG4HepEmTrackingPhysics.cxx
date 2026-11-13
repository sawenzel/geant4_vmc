//------------------------------------------------
// The Geant4 Virtual Monte Carlo package
// Copyright (C) 2007 - 2025 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file TG4HepEmTrackingPhysics.cxx
/// \brief Implementation of the TG4HepEmTrackingPhysics class
///
/// \author I. Hrivnacova; IPN, Orsay

#ifdef USE_G4HEPEM

#include "TG4HepEmTrackingPhysics.h"

#include "G4HepEmTrackingManager.hh"
#include "G4HepEmConfig.hh"

#include "G4EmParameters.hh"
#include "G4Version.hh"

#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4ParticleDefinition.hh"
#include "G4Positron.hh"
#include "G4BuilderType.hh"

#include "G4EmBuilder.hh"
#include "G4hMultipleScattering.hh"
#include "G4NuclearStopping.hh"
#include "G4ionIonisation.hh"
#include "G4EmModelActivator.hh"

//_____________________________________________________________________________
TG4HepEmTrackingPhysics::TG4HepEmTrackingPhysics(const G4String& name)
   :  G4VPhysicsConstructor(name)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();

  param->SetMscRangeFactor(0.04);
  SetPhysicsType(bElectromagnetic);

  verboseLevel = 1;
}

//_____________________________________________________________________________
TG4HepEmTrackingPhysics::~TG4HepEmTrackingPhysics()
{}

//_____________________________________________________________________________
void TG4HepEmTrackingPhysics::ConstructProcess()
{
  if(verboseLevel > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }

  // G4HepEm
  // Register custom tracking manager for e-/e+ and gammas.
  G4cout << "Constructing G4HepEmTrackingManager with verbose level " << verboseLevel << G4endl;
  auto* trackingManager = new G4HepEmTrackingManager(verboseLevel);
  // Configuration of G4HepEm
  // Several paramaters can be configured per detector region. These are:
  //  MSC parameters, continuous energy loss step limit function parameters,
  //  MSC minimal/default step limit, Woodcock tracking of photons, energy loss
  //  fluctuation, multiple steps in the combined MSC with Transportation
  // Here we activate only one: Woodcock tracking in the calorimeter region (Woodcock_Region)
  G4HepEmConfig* config = trackingManager->GetConfig();
  config->SetWoodcockTrackingRegion("Woodcock_Region");

  // several msc steps per transportation step only if requested, as in Geant4
  // (/process/em/transportationWithMsc)
#if G4VERSION_NUMBER >= 1110
  config->SetMultipleStepsInMSCWithTransportation(
    G4EmParameters::Instance()->TransportationWithMsc() ==
    G4TransportationWithMscType::fMultipleSteps);
#else
  config->SetMultipleStepsInMSCWithTransportation(false);
#endif
  G4Electron::Definition()->SetTrackingManager(trackingManager);
  G4Positron::Definition()->SetTrackingManager(trackingManager);
  G4Gamma::Definition()->SetTrackingManager(trackingManager);

  if (G4Threading::IsMasterThread() && verboseLevel > 0) {
    G4EmParameters::Instance()->Dump();
  }
  // end G4HepEm

  /// Added by us: EM processes for other particles (from G4EmStandardPhysics)

  G4EmBuilder::PrepareEMPhysics();

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4EmParameters* param = G4EmParameters::Instance();

  // processes used by several particles
  G4hMultipleScattering* hmsc = new G4hMultipleScattering("ionmsc");

  // nuclear stopping is enabled if th eenergy limit above zero
  G4double nielEnergyLimit = param->MaxNIELEnergy();
  G4NuclearStopping* pnuc = nullptr;
  if(nielEnergyLimit > 0.0) {
    pnuc = new G4NuclearStopping();
    pnuc->SetMaxKinEnergy(nielEnergyLimit);
  }

  // high energy limit for e+- scattering models and bremsstrahlung
  G4double highEnergyLimit = param->MscEnergyLimit();

  // Other particles
  // (From G4EmStandardPhysics)

  // generic ion
  auto particle = G4GenericIon::GenericIon();
  G4ionIonisation* ionIoni = new G4ionIonisation();
  ph->RegisterProcess(hmsc, particle);
  ph->RegisterProcess(ionIoni, particle);
  if(nullptr != pnuc) { ph->RegisterProcess(pnuc, particle); }

  // muons, hadrons ions
  G4EmBuilder::ConstructCharged(hmsc, pnuc);

  // extra configuration
  G4EmModelActivator mact(GetPhysicsName());

  if (G4Threading::IsMasterThread() && verboseLevel > 0) {
    G4EmParameters::Instance()->Dump();
  }

  G4cout << "### G4HepEmTrackingPhysics constructed." << G4endl;
}

#endif