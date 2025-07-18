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

  // Calculate derived dimensions
  G4double crystalRadius = crystalDiameter / 2.0;
  G4double crystalHoleRadius = crystalHoleDiameter / 2.0;
  G4double endcapWellRadius = endcapWellDiameter / 2.0;

  // World size (add margins)
  G4double worldSizeXY = 2.0 * (crystalRadius + endcapSideThickness + endcapToCrystalGap + 10. * mm);
  G4double worldSizeZ = 2.0 * (crystalLength / 2.0 + endcapTopThickness + endcapToCrystalGap + 10. * mm);

  // Get materials
  auto vacuum = G4Material::GetMaterial("Galactic");
  auto crystalMaterial = G4Material::GetMaterial("G4_Ge");
  auto endcapMaterial = G4Material::GetMaterial("G4_Al");
  auto mountCupMaterial = G4Material::GetMaterial("G4_Al");
  auto crystalDeadLayerMaterialTopSide = G4Material::GetMaterial("G4_Li");
  auto crystalDeadLayerMaterialHole = G4Material::GetMaterial("G4_B");

  if ( ! vacuum || ! crystalMaterial || ! endcapMaterial || ! mountCupMaterial || ! crystalDeadLayerMaterialTopSide || ! crystalDeadLayerMaterialHole) {
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

  // Well Detector Assembly (Mother Volume for the entire detector)
  // This will be a logical volume representing the whole detector system (crystal, endcap, vacuum)
  // Define its size based on the outermost dimensions.
  G4double detectorOuterRadius = crystalRadius + endcapSideThickness + endcapToCrystalGap;
  G4double detectorTotalLength = crystalLength + endcapTopThickness + endcapToCrystalGap;

  auto detectorAssemblyS = new G4Tubs("DetectorAssembly",           // its name
                                      0.,                           // its inner radius
                                      detectorOuterRadius,          // its outer radius
                                      detectorTotalLength / 2.0,    // hz
                                      0. * deg,                     // its start angle
                                      360. * deg);                  // its spanning angle

  auto detectorAssemblyLV = new G4LogicalVolume(detectorAssemblyS, vacuum, "DetectorAssembly");

  new G4PVPlacement(nullptr,
      G4ThreeVector(),
      detectorAssemblyLV,
      "DetectorAssembly",
      worldLV,
      false,
      0,
      fCheckOverlaps);

  // Endcap (Outer container)
  // Modeled as a hollow cylinder with a central well.
  // We'll use G4Tubs and subtract to create the well.
  G4double endcapOuterRadius = crystalRadius + endcapSideThickness + endcapToCrystalGap;
  G4double endcapInnerRadius = crystalRadius + endcapToCrystalGap;
  G4double endcapTotalLength = crystalLength + endcapToCrystalGap + endcapTopThickness;
  G4double endcapBodyLength = endcapTotalLength - endcapTopThickness;

  // Solid for the main cylindrical part of the endcap
  auto endcapBodyS = new G4Tubs("EndcapBody",
      endcapInnerRadius, // Inner radius is where the crystal gap starts
      endcapOuterRadius,
      endcapBodyLength / 2.0,
      0. * deg, 360. * deg);

  // Solid for the top plate of the endcap
  auto endcapTopS = new G4Tubs("EndcapTop",
      0.,
      endcapOuterRadius,
      endcapTopThickness / 2.0,
      0. * deg, 360. * deg);

  // Solid for the endcap well - this will be subtracted from the top part
  G4double endcapWellOuterRadius = endcapWellRadius + endcapWellSideThickness;
  G4double endcapWellTotalDepth = endcapWellDepth + endcapWellBottomThickness;

  auto endcapWellHoleS = new G4Tubs("EndcapWellHole",
      0.,
      endcapWellOuterRadius,
      endcapWellTotalDepth / 2.0,
      0. * deg, 360. * deg);

  // Combine endcap components.
  // For simplicity, we'll create the endcap as a main hollow cylinder
  // and then position the top with a hole.

  // Main endcap cylinder (sides and bottom where crystal rests)
  auto endcapS = new G4Tubs("Endcap",
      crystalRadius + endcapToCrystalGap,
      endcapOuterRadius,
      endcapTotalLength / 2.0, // length of the whole endcap structure
      0. * deg, 360. * deg);

  auto endcapLV = new G4LogicalVolume(endcapS,
      endcapMaterial,
      "EndcapLV");

  // Position the endcap within the detector assembly.
  // The bottom of the endcap will be at -detectorTotalLength/2.0 in the assembly,
  // so its center is offset.
  new G4PVPlacement(nullptr,
      G4ThreeVector(0., 0., 0.), // Centered within the assembly for now
      endcapLV,
      "EndcapPV",
      detectorAssemblyLV,
      false,
      0,
      fCheckOverlaps);

  // Crystal
  // Modeled as a hollow cylinder.
  G4double activeCrystalOuterRadius = crystalRadius - crystalDeadLayerTopSide;
  G4double activeCrystalInnerRadius = crystalHoleRadius + crystalDeadLayerHole;
  G4double activeCrystalLength = crystalLength - crystalDeadLayerTopSide; // Assuming top dead layer affects length
  G4double activeCrystalHoleDepth = crystalHoleDepth - crystalDeadLayerHole; // Assuming bottom dead layer affects hole depth

  auto crystalS = new G4Tubs("Crystal",
      activeCrystalInnerRadius, // Inner radius of the active crystal
      activeCrystalOuterRadius, // Outer radius of the active crystal
      activeCrystalLength / 2.0, // Half length
      0. * deg, 360. * deg);

  auto crystalLV = new G4LogicalVolume(crystalS,
      crystalMaterial,
      "CrystalLV");

  // Position the crystal within the endcap, considering the gap and endcap thickness.
  // The crystal's center needs to be calculated relative to the endcap.
  // We assume the bottom of the crystal aligns with the bottom of the endcap's inner cavity.
  G4double crystalZPos = -(endcapTotalLength / 2.0) + crystalLength / 2.0 + endcapWellBottomThickness; // Simplified initial positioning

  new G4PVPlacement(nullptr,
      G4ThreeVector(0., 0., 0.), // Position relative to detectorAssemblyLV
      "CrystalPV",
      detectorAssemblyLV,
      false,
      0,
      fCheckOverlaps);

  // Crystal Dead Layers
  // These will be thin shells.
  // Outer/Top Dead Layer (6)
  auto crystalDeadLayerOuterS = new G4Tubs("CrystalDeadLayerOuter",
      activeCrystalOuterRadius,
      crystalRadius,
      crystalLength / 2.0, // Full crystal length
      0. * deg, 360. * deg);
  auto crystalDeadLayerOuterLV = new G4LogicalVolume(crystalDeadLayerOuterS,
      crystalDeadLayerMaterialTopSide,
      "CrystalDeadLayerOuterLV");
  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), crystalDeadLayerOuterLV,
      "CrystalDeadLayerOuterPV", detectorAssemblyLV, false, 0, fCheckOverlaps);

  // Hole Dead Layer (7)
  auto crystalDeadLayerHoleS = new G4Tubs("CrystalDeadLayerHole",
      crystalHoleRadius,
      activeCrystalInnerRadius,
      crystalHoleDepth / 2.0, // Full hole depth
      0. * deg, 360. * deg);
  auto crystalDeadLayerHoleLV = new G4LogicalVolume(crystalDeadLayerHoleS,
      crystalDeadLayerMaterialHole,
      "CrystalDeadLayerHoleLV");
  // Position it inside the crystal hole
  G4double crystalHoleZOffset = crystalLength / 2.0 - crystalHoleDepth / 2.0; // Assuming hole is from the top
  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., crystalHoleZOffset), crystalDeadLayerHoleLV,
      "CrystalDeadLayerHolePV", detectorAssemblyLV, false, 0, fCheckOverlaps);


  // Mount Cup (5)
  // A cylindrical ring at the bottom of the crystal.
  G4double mountCupInnerRadius = crystalRadius; // Assumed to be same as crystal outer radius
  G4double mountCupOuterRadius = mountCupInnerRadius + mountCupThickness;
  G4double mountCupLength = mountCupThickness; // Assuming it's a short ring

  auto mountCupS = new G4Tubs("MountCup",
      mountCupInnerRadius,
      mountCupOuterRadius,
      mountCupLength / 2.0,
      0. * deg, 360. * deg);

  auto mountCupLV = new G4LogicalVolume(mountCupS,
      mountCupMaterial,
      "MountCupLV");

  // Position the mount cup at the bottom of the crystal
  G4double mountCupZPos = -(crystalLength / 2.0) - (mountCupLength / 2.0); // Adjust this based on exact position
  new G4PVPlacement(nullptr,
      G4ThreeVector(0., 0., mountCupZPos),
      mountCupLV,
      "MountCupPV",
      detectorAssemblyLV, // Place it in the detector assembly
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
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VisAttributes* crystalVisAtt = new G4VisAttributes(G4Colour::Blue());
  crystalVisAtt->SetForceSolid(true);
  crystalLV->SetVisAttributes(crystalVisAtt);

  G4VisAttributes* endcapVisAtt = new G4VisAttributes(G4Colour::Grey());
  endcapVisAtt->SetForceSolid(true);
  endcapLV->SetVisAttributes(endcapVisAtt);

  G4VisAttributes* deadLayerVisAttOuter = new G4VisAttributes(G4Colour::Red());
  deadLayerVisAttOuter->SetForceSolid(true);
  crystalDeadLayerOuterLV->SetVisAttributes(deadLayerVisAttOuter);

  G4VisAttributes* deadLayerVisAttHole = new G4VisAttributes(G4Colour::Magenta());
  deadLayerVisAttHole->SetForceSolid(true);
  crystalDeadLayerHoleLV->SetVisAttributes(deadLayerVisAttHole);

  G4VisAttributes* mountCupVisAtt = new G4VisAttributes(G4Colour::Green());
  mountCupVisAtt->SetForceSolid(true);
  mountCupLV->SetVisAttributes(mountCupVisAtt);

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

