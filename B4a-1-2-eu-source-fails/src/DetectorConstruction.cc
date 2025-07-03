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
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

namespace B4
{

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fAbsorberPV(nullptr), fGapPV(nullptr)
{}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();

  // Materials
  G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* cylinderMat = nist->FindOrBuildMaterial("G4_Si");

  // Geometry parameters
  G4double worldSize = 20*cm; // World box
  G4double cylRadius = 2.54*cm; // 2 inches diameter
  G4double cylHeight = 2.54*cm; // 2 inches height

  // World volume
  G4Box* solidWorld = new G4Box("World", worldSize/2, worldSize/2, worldSize/2);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(
    nullptr, {}, logicWorld, "World", nullptr, false, 0);

  // Cylindrical detector
  G4Tubs* solidCyl = new G4Tubs("Cyl", 0., cylRadius, cylHeight/2, 0., 360.*deg);
  G4LogicalVolume* logicCyl = new G4LogicalVolume(solidCyl, cylinderMat, "Cyl");
  fAbsorberPV = new G4PVPlacement(
    nullptr, {}, logicCyl, "Cyl", logicWorld, false, 0);

  return physWorld;
}

const G4VPhysicalVolume* DetectorConstruction::GetAbsorberPV() const
{
  return fAbsorberPV;
}

const G4VPhysicalVolume* DetectorConstruction::GetGapPV() const
{
  return nullptr; // No gap
}

}