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

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

namespace B4a
{

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleSource(new G4GeneralParticleSource)
{

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleSource;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Random direction
  G4double u = G4UniformRand();
  G4double theta = std::acos(1. - 2. * u);
  G4double phi = 2. * CLHEP::pi * G4UniformRand();

  G4double dirX = std::sin(theta) * std::cos(phi);
  G4double dirY = std::sin(theta) * std::sin(phi);
  G4double dirZ = std::cos(theta);

  G4ThreeVector direction(dirX, dirY, dirZ);
  fParticleSource->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(direction);

  fParticleSource->GeneratePrimaryVertex(anEvent);
}

}
