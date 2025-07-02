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
/// \file B4/B4a/src/DetectorConstruction.cc
/// \brief Implementation of the B4::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Box.hh"

namespace B4
{

DetectorConstruction::DetectorConstruction()
{}

DetectorConstruction::~DetectorConstruction()
{}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  auto nist = G4NistManager::Instance();

  // Define materials
  fAbsorberMaterial = nist->FindOrBuildMaterial("G4_Al");
  fGapMaterial = nist->FindOrBuildMaterial("G4_AIR");

  // World
  auto worldSize = 20 * cm;
  auto solidWorld = new G4Box("World", worldSize, worldSize, worldSize);
  auto logicWorld = new G4LogicalVolume(solidWorld, fGapMaterial, "World");
  auto physWorld = new G4PVPlacement(nullptr, {}, logicWorld, "World", nullptr, false, 0);

  // Cylindrical detector
  auto height = 2 * 2.54 * cm;  // 2 inches in cm
  auto radius = 1 * 2.54 * cm;  // 1 inch in cm

  auto solidDetector = new G4Tubs("Absorber", 0., radius, height / 2., 0., 360. * deg);
  auto logicDetector = new G4LogicalVolume(solidDetector, fAbsorberMaterial, "Absorber");
  fAbsorberPV = new G4PVPlacement(nullptr, {}, logicDetector, "Absorber", logicWorld, false, 0);

  fGapPV = nullptr;

  return physWorld;
}

}  // namespace B4