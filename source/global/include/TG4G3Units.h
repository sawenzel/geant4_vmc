#ifndef TG4_G3_UNITS_H
#define TG4_G3_UNITS_H

//------------------------------------------------
// The Geant4 Virtual Monte Carlo package
// Copyright (C) 2007 - 2014 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file TG4G3Units.h
/// \brief Definition of the TG4G3Units class 
///
/// \author I. Hrivnacova; IPN, Orsay

#include <globals.hh>

/// \ingroup global
/// \brief G3 default units of physical quantities; 
///
/// All physical quantities returned by VMC are expressed in these units.
///
/// \author I. Hrivnacova; IPN, Orsay

class TG4G3Units
{
  public:

    // static get methods
    static G4double Length(); 
    static G4double Angle(); 
    static G4double Time(); 
    static G4double Charge(); 
    static G4double Energy(); 
    static G4double Mass(); 
    static G4double MassDensity(); 
    static G4double AtomicWeight();     
    static G4double Field();

    static G4double InvLength();
    static G4double InvAngle();
    static G4double InvTime();
    static G4double InvCharge();
    static G4double InvEnergy();
    static G4double InvMass();
    static G4double InvMassDensity();
    static G4double InvAtomicWeight();
    static G4double InvField();

   private:
    TG4G3Units();      

    // static data members  
    static const G4double fgkLength;       ///< G3 length unit 
    static const G4double fgkAngle;        ///< G3 angle unit 
    static const G4double fgkTime;         ///< G3 time unit 
    static const G4double fgkCharge;       ///< G3 charge unit  
    static const G4double fgkEnergy;       ///< G3 energy unit  
    static const G4double fgkMass;         ///< G3 mass unit
    static const G4double fgkMassDensity;  ///< G3 mass density unit 
    static const G4double fgkAtomicWeight; ///< G3 atomic weight unit  
    static const G4double fgkField;        ///< G3 magnetic field unit

    // the inverses to avoid division at runtime
    static const G4double fgkInvLength;        ///< G3 length unit
    static const G4double fgkInvAngle;         ///< G3 angle unit
    static const G4double fgkInvTime;          ///< G3 time unit
    static const G4double fgkInvCharge;        ///< G3 charge unit
    static const G4double fgkInvEnergy;        ///< G3 energy unit
    static const G4double fgkInvMass;          ///< G3 mass unit
    static const G4double fgkInvMassDensity;   ///< G3 mass density unit
    static const G4double fgkInvAtomicWeight;  ///< G3 atomic weight unit
    static const G4double fgkInvField;         ///< G3 magnetic field unit
};     

// inline methods

inline G4double TG4G3Units::Length() { 
  /// Return G3 length unit
  return fgkLength; 
}

inline G4double TG4G3Units::Angle() { 
  /// Return G3 angle unit
  return fgkAngle; 
}

inline G4double TG4G3Units::Time() { 
  /// Return G3 time unit
  return fgkTime; 
}

inline G4double TG4G3Units::Charge() { 
  /// Return G3 chanrge unit
  return fgkCharge; 
}

inline G4double TG4G3Units::Energy() { 
  /// Return G3 energy unit
  return fgkEnergy; 
}

inline G4double TG4G3Units::Mass() { 
  /// Return G3 mass unit
  return fgkMass; 
}

inline G4double TG4G3Units::MassDensity()  
{ 
  /// Return G3 mass density unit
  return fgkMassDensity; 
}

inline G4double TG4G3Units::AtomicWeight() { 
  /// Return G3 atomic weight unit
  return fgkAtomicWeight; 
}

inline G4double TG4G3Units::Field() 
{ 
  /// Return G3 magnetic field unit
  return fgkField; 
}

inline G4double TG4G3Units::InvLength() {
  /// Return G3 length unit
  return fgkInvLength;
}

inline G4double TG4G3Units::InvAngle() {
  /// Return G3 angle unit
  return fgkInvAngle;
}

inline G4double TG4G3Units::InvTime() {
  /// Return G3 time unit
  return fgkInvTime;
}

inline G4double TG4G3Units::InvCharge() {
  /// Return G3 chanrge unit
  return fgkInvCharge;
}

inline G4double TG4G3Units::InvEnergy() {
  /// Return G3 energy unit
  return fgkInvEnergy;
}

inline G4double TG4G3Units::InvMass() {
  /// Return G3 mass unit
  return fgkInvMass;
}

inline G4double TG4G3Units::InvMassDensity()
{
  /// Return G3 mass density unit
  return fgkInvMassDensity;
}

inline G4double TG4G3Units::InvAtomicWeight() {
  /// Return G3 atomic weight unit
  return fgkInvAtomicWeight;
}

inline G4double TG4G3Units::InvField()
{
  /// Return G3 magnetic field unit
  return fgkInvField;
}

#endif //TG4_G3_UNITS_H
