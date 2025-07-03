//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B4/B4a/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B4a::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include <cmath>

using namespace B4a;

namespace B4a
{

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    fParticleSource(new G4GeneralParticleSource())
{
  // Define Eu-152 ion
  G4int Z = 63;
  G4int A = 152;
  G4double E = 0.*keV;

  auto ion = G4IonTable::GetIonTable()->GetIon(Z, A, E);
  fParticleSource->SetParticleDefinition(ion);
  fParticleSource->SetParticleCharge(0.0);
  fParticleSource->GetCurrentSource()->GetIonSource()->SetEnableDecay(true);

  // Set position at origin
  fParticleSource->SetParticlePosition(G4ThreeVector(0., 0., 0.));
  fParticleSource->SetParticleTime(0.);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleSource;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Generate random isotropic direction
  G4double u = G4UniformRand();
  G4double theta = std::acos(1. - 2. * u);
  G4double phi = 2. * CLHEP::pi * G4UniformRand();

  G4double dirX = std::sin(theta) * std::cos(phi);
  G4double dirY = std::sin(theta) * std::sin(phi);
  G4double dirZ = std::cos(theta);

  fParticleSource->SetParticleMomentumDirection(G4ThreeVector(dirX, dirY, dirZ));
  fParticleSource->GeneratePrimaryVertex(anEvent);
}

}
