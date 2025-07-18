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
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Materials defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Ge");
  nistManager->FindOrBuildMaterial("G4_Al");
  nistManager->FindOrBuildMaterial("G4_Li");
  nistManager->FindOrBuildMaterial("G4_B");

  // Vacuum
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4double crystalLength = 71.1 * mm; // A
  G4double crystalDiameter = 56.8 * mm; // B
  G4double crystalHoleDepth = 57.7 * mm; // C
  G4double crystalHoleDiameter = 20.6 * mm; // D
  G4double endcapWellDepth = 51.3 * mm; // E
  G4double endcapWellDiameter = 15.5 * mm; // F
  G4double endcapToCrystalGap = 9 * mm; // G

  G4double endcapTopThickness = 2 * mm; // 1
  G4double endcapSideThickness = 1.6 * mm; // 2
  G4double endcapWellSideThickness = 0.5 * mm; // 3
  G4double endcapWellBottomThickness = 2 * mm; // 4
  G4double mountCupThickness = 1.6 * mm; // 5
  G4double crystalDeadLayerTopSide = 0.7 * mm; // 6
  G4double crystalDeadLayerHole = 0.3 * um; // 7

  // World size (add margins)
  G4double worldSizeXY = 1.5 * (crystalDiameter / 2 + endcapSideThickness + mountCupThickness);
  G4double worldSizeZ = 1.5 * (crystalLength + endcapTopThickness + mountCupThickness);

  // Get materials
  auto vacuum = G4Material::GetMaterial("Galactic");
  auto aluminium = G4Material::GetMaterial("G4_Al");
  auto germanium = G4Material::GetMaterial("G4_Ge");
  auto geLiDeadLayer = G4Material::GetMaterial("G4_Li");
  auto geBDeadLayer = G4Material::GetMaterial("G4_B");

  if ( ! vacuum || ! aluminium || ! germanium || ! geLiDeadLayer || !geBDeadLayer) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 vacuum,           // its material
                 "World");         // its name

  auto worldPV = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                         // at (0,0,0)
    worldLV,                                 // its logical volume
    "World",                                 // its name
    nullptr,                                 // its mother  volume
    false,                                   // no boolean operation
    0,                                       // copy number
    fCheckOverlaps);                         // checking overlaps

  //
  // Mount Cup (outermost Aluminium casing)
  //
  G4double mountCupOuterRadius = crystalDiameter / 2 + endcapSideThickness + mountCupThickness;
  G4double mountCupHeight = crystalLength + endcapTopThickness + mountCupThickness + endcapToCrystalGap;
  G4double mountCupInnerRadius = crystalDiameter / 2 + endcapSideThickness;

  auto mountCupSolid = new G4Tubs("MountCup",   // its name
      mountCupInnerRadius,                      // its inner radius
      mountCupOuterRadius,                      // its outer radius
      mountCupHeight / 2,                       // hz
      0. * deg,                                 // its start angle
      360. * deg);                              // its spanning angle

  auto mountCupLV = new G4LogicalVolume(mountCupSolid,  // its solid
      aluminium,                                        // its material
      "MountCup");                                      // its name

  new G4PVPlacement(nullptr,        // no rotation
      G4ThreeVector(),              // at (0,0,0)
      mountCupLV,                   // its logical volume
      "MountCup",                   // its name
      worldLV,                      // its mother  volume
      false,                        // no boolean operation
      0,                            // copy number
      fCheckOverlaps);              // checking overlaps

  //
  // Endcap (Aluminium)
  //
  // The endcap is a cylinder with a well (hole)
  G4double endcapOuterRadius = crystalDiameter / 2 + endcapSideThickness;
  G4double endcapTotalHeight = endcapTopThickness + endcapWellDepth + endcapWellBottomThickness;

  auto endcapBaseSolid = new G4Tubs("EndcapBase",
      0., endcapOuterRadius,
      endcapTotalHeight / 2,
      0. * deg, 360. * deg);

  G4double endcapWellInnerRadius = endcapWellDiameter / 2;
  G4double endcapWellSubtractedHeight = endcapWellDepth;

  auto endcapWellSubtractor = new G4Tubs("EndcapWellSubtractor",
      0., endcapWellInnerRadius,
      endcapWellSubtractedHeight / 2,
      0. * deg, 360. * deg);

  // Position the well subtractor relative to the base of the endcap
  // The well is at the 'bottom' of the endcap when considering its full extent
  G4ThreeVector wellPosition = G4ThreeVector(0., 0., -endcapTotalHeight / 2 + endcapWellSubtractedHeight / 2 + endcapWellBottomThickness);

  auto endcapSolid = new G4SubtractionSolid("Endcap",   // pName
      endcapBaseSolid,                                  // pSolidA
      endcapWellSubtractor,                             // pSolidB
      nullptr,                                          // rotMatrix
      wellPosition);                                    // transVector

  auto endcapLV = new G4LogicalVolume(endcapSolid,
      aluminium,
      "Endcap");

  // Position the endcap inside the mount cup. It sits at the 'top' of the detector assembly.
  G4double endcapZPositionInMountCup = mountCupHeight / 2 - endcapTotalHeight / 2 - mountCupThickness; // Top of mount cup, minus its thickness for the lip, minus half endcap height

  new G4PVPlacement(nullptr,
      G4ThreeVector(0., 0., endcapZPositionInMountCup),
      endcapLV,
      "Endcap",
      mountCupLV,
      false,
      0,
      fCheckOverlaps);

  //
  // Crystal (Germanium) - with an inner hole
  //
  G4double crystalOuterRadius = crystalDiameter / 2;
  G4double crystalInnerRadius = crystalHoleDiameter / 2;
  G4double crystalTotalLength = crystalLength;
  G4double crystalHoleTotalDepth = crystalHoleDepth;

  auto crystalFullSolid = new G4Tubs("CrystalFull",
      0., crystalOuterRadius,
      crystalTotalLength / 2,
      0. * deg, 360. * deg);

  auto crystalHoleSubtractor = new G4Tubs("CrystalHoleSubtractor",
      0., crystalInnerRadius,
      crystalHoleTotalDepth / 2,
      0. * deg, 360. * deg);

  // Position the hole subtractor relative to the crystal
  // The hole starts from one end of the crystal
  G4ThreeVector holePositionInCrystal = G4ThreeVector(0., 0., crystalTotalLength / 2 - crystalHoleTotalDepth / 2);

  auto crystalSolid = new G4SubtractionSolid("Crystal",
      crystalFullSolid,
      crystalHoleSubtractor,
      nullptr,
      holePositionInCrystal);

  auto crystalLV = new G4LogicalVolume(crystalSolid,
      germanium,
      "Crystal");

  // Position the crystal below the endcap, with the specified gap.
  G4double crystalZPositionInMountCup = endcapZPositionInMountCup - endcapTotalHeight / 2 - endcapToCrystalGap - crystalTotalLength / 2;

  new G4PVPlacement(nullptr,
      G4ThreeVector(0., 0., crystalZPositionInMountCup),
      crystalLV,
      "Crystal",
      mountCupLV,
      false,
      0,
      fCheckOverlaps);

  //
  // Crystal Dead Layer - Top and Side (Ge/Li)
  // This layer covers the top surface and the outer cylindrical surface of the crystal.
  //
  // Side dead layer (cylindrical shell)
  G4double crystalActiveOuterRadius = crystalOuterRadius - crystalDeadLayerTopSide;
  auto crystalSideDeadLayerSolid = new G4Tubs("CrystalSideDeadLayer",
      crystalActiveOuterRadius, crystalOuterRadius,
      crystalTotalLength / 2,
      0. * deg, 360. * deg);

  auto crystalSideDeadLayerLV = new G4LogicalVolume(crystalSideDeadLayerSolid,
      geLiDeadLayer,
      "CrystalSideDeadLayer");

  new G4PVPlacement(nullptr,
      G4ThreeVector(0., 0., 0.), // Placed within the crystal volume, at the same Z
      crystalSideDeadLayerLV,
      "CrystalSideDeadLayer",
      crystalLV,
      false,
      0,
      fCheckOverlaps);

  // Top dead layer (disc)
  // It's a disc at the very top of the crystal, with thickness crystalDeadLayerTopSide
  // and outer radius crystalOuterRadius, and inner radius crystalInnerRadius (because of the hole).
  auto crystalTopDeadLayerSolid = new G4Tubs("CrystalTopDeadLayer",
      crystalInnerRadius, crystalOuterRadius,
      crystalDeadLayerTopSide / 2,
      0. * deg, 360. * deg);

  auto crystalTopDeadLayerLV = new G4LogicalVolume(crystalTopDeadLayerSolid,
      geLiDeadLayer,
      "CrystalTopDeadLayer");

  new G4PVPlacement(nullptr,
      G4ThreeVector(0., 0., crystalTotalLength / 2 - crystalDeadLayerTopSide / 2), // At the top of the crystal
      crystalTopDeadLayerLV,
      "CrystalTopDeadLayer",
      crystalLV,
      false,
      0,
      fCheckOverlaps);


  //
  // Crystal Dead Layer - Hole Inner and Bottom (Ge/B)
  // This layer covers the inner cylindrical surface and the bottom of the crystal hole.
  //
  // Hole inner dead layer (cylindrical shell inside the hole)
  G4double crystalHoleActiveInnerRadius = crystalInnerRadius - crystalDeadLayerHole;
  auto crystalHoleInnerDeadLayerSolid = new G4Tubs("CrystalHoleInnerDeadLayer",
      crystalHoleActiveInnerRadius, crystalInnerRadius,
      crystalHoleTotalDepth / 2,
      0. * deg, 360. * deg);

  auto crystalHoleInnerDeadLayerLV = new G4LogicalVolume(crystalHoleInnerDeadLayerSolid,
      geBDeadLayer,
      "CrystalHoleInnerDeadLayer");

  new G4PVPlacement(nullptr,
      G4ThreeVector(0., 0., holePositionInCrystal.z()), // Same position as hole
      crystalHoleInnerDeadLayerLV,
      "CrystalHoleInnerDeadLayer",
      crystalLV,
      false,
      0,
      fCheckOverlaps);


  // Hole bottom dead layer (disc at the bottom of the hole)
  auto crystalHoleBottomDeadLayerSolid = new G4Tubs("CrystalHoleBottomDeadLayer",
      0., crystalInnerRadius,
      crystalDeadLayerHole / 2,
      0. * deg, 360. * deg);

  auto crystalHoleBottomDeadLayerLV = new G4LogicalVolume(crystalHoleBottomDeadLayerSolid,
      geBDeadLayer,
      "CrystalHoleBottomDeadLayer");

  new G4PVPlacement(nullptr,
      G4ThreeVector(0., 0., holePositionInCrystal.z() - crystalHoleTotalDepth / 2 + crystalDeadLayerHole / 2), // At the bottom of the hole
      crystalHoleBottomDeadLayerLV,
      "CrystalHoleBottomDeadLayer",
      crystalLV,
      false,
      0,
      fCheckOverlaps);

  //
  // print parameters
  //
  G4cout
      << G4endl
      << "------------------------------------------------------------" << G4endl
      << "---> The detector is [UPDATE HERE]" << G4endl;

  //
  // Visualization attributes
  //
  G4VisAttributes* aluminiumVisAtt = new G4VisAttributes(G4Colour(0.6, 0.6, 0.6)); // Grey
  mountCupLV->SetVisAttributes(aluminiumVisAtt);
  endcapLV->SetVisAttributes(aluminiumVisAtt);
  // No separate vis attributes for endcapWellLV since it's now part of endcapSolid

  G4VisAttributes* germaniumVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // Red
  crystalLV->SetVisAttributes(germaniumVisAtt);

  G4VisAttributes* deadLayerVisAtt = new G4VisAttributes(G4Colour(0.2, 0.8, 0.2)); // Green
  crystalSideDeadLayerLV->SetVisAttributes(deadLayerVisAtt);
  crystalTopDeadLayerLV->SetVisAttributes(deadLayerVisAtt);
  crystalHoleInnerDeadLayerLV->SetVisAttributes(deadLayerVisAtt);
  crystalHoleBottomDeadLayerLV->SetVisAttributes(deadLayerVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

