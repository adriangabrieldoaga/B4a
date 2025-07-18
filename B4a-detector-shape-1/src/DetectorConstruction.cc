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
  G4double A = 71.1 * mm;
  G4double B = 56.8 * mm;
  G4double C = 57.7 * mm;
  G4double D = 20.6 * mm;
  G4double E = 51.3 * mm;
  G4double F = 15.5 * mm;
  G4double G = 9.0 * mm;

  G4double t1 = 2.0 * mm;     // endcap top thickness
  G4double t2 = 1.6 * mm;     // endcap side thickness
  G4double t3 = 0.5 * mm;     // endcap well side thickness
  G4double t4 = 2.0 * mm;     // endcap well bottom thickness
  G4double t5 = 1.6 * mm;     // mount cup thickness
  G4double t6 = 0.7 * mm;     // crystal top/side dead layer (Ge/Li)
  G4double t7 = 0.0003 * mm;  // crystal hole inner/bottom dead layer (Ge/B)

  // World size (add margins)
  G4double worldSizeX = B + 20 * mm;
  G4double worldSizeY = B + 20 * mm;
  G4double worldSizeZ = A + t1 + t2 + G + t4 + 20 * mm;

  // Get materials
  auto vacuum = G4Material::GetMaterial("Galactic");
  auto germanium = G4Material::GetMaterial("G4_Ge");
  auto aluminium = G4Material::GetMaterial("G4_Al");

  // Dead layers use the same material for now
  auto deadLayer_Li = germanium;
  auto deadLayer_B = germanium;

  if ( ! vacuum || ! aluminium || ! germanium) {
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
                 worldSizeX/2, worldSizeY/2, worldSizeZ/2); // its size

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 vacuum,  // its material
                 "World");         // its name

  auto worldPV = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                         // at (0,0,0)
    worldLV,                                 // its logical volume
    "World",                                 // its name
    nullptr,                                 // its mother  volume
    false,                                   // no boolean operation
    0,                                       // copy number
    fCheckOverlaps);                         // checking overlaps

  // Crystal outer cylinder
  auto crystalOuterS = new G4Tubs("CrystalOuter", // its name
                                  0,              // its inner radius
                                  B / 2,          // its outer radius
                                  A / 2,          // hz
                                  0,              // its start angle
                                  360 * deg);     // its spanning angle

  auto crystalOuterLV = new G4LogicalVolume(crystalOuterS,   // its solid
                                            germanium,       // its material
                                            "CrystalOuter"); // its name

  new G4PVPlacement(nullptr,            // no rotation
                    G4ThreeVector(),    // at (0,0,0)
                    crystalOuterLV,     // its logical volume
                    "CrystalOuter",     // its name
                    worldLV,            // its mother  volume
                    false,              // no boolean operation
                    0,                  // copy number
                    fCheckOverlaps);    // checking overlaps

  // Crystal inner hole
  auto holeS = new G4Tubs("CrystalHole", 0, D / 2, C / 2, 0, 360 * deg);
  auto holeLV = new G4LogicalVolume(holeS, vacuum, "CrystalHole");
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (A - C) / 2), holeLV, "CrystalHole", crystalOuterLV, false, 0, fCheckOverlaps);

  // Dead layer - top and side (Li)
  auto deadOuterS = new G4Tubs("DeadLayerOuter", 0, (B / 2), A / 2, 0, 360 * deg);
  auto deadInnerS = new G4Tubs("DeadLayerInner", 0, (B / 2 - t6), A / 2 - t6, 0, 360 * deg);
  auto deadLayerShell = new G4SubtractionSolid("DeadLayerSideTop", // pName
                                               deadOuterS,         // pSolidA
                                               deadInnerS);        // pSolidB
  auto deadLayerLV = new G4LogicalVolume(deadLayerShell, deadLayer_Li, "DeadLayerLi");
  new G4PVPlacement(nullptr, G4ThreeVector(), deadLayerLV, "DeadLayerLi", crystalOuterLV, false, 0, fCheckOverlaps);

  // Dead layer - inner hole (B)
  auto deadHoleOuterS = new G4Tubs("DeadHoleOuter", 0, D / 2, C / 2, 0, 360 * deg);
  auto deadHoleInnerS = new G4Tubs("DeadHoleInner", 0, (D / 2 - t7), C / 2 - t7, 0, 360 * deg);
  auto deadHoleShell = new G4SubtractionSolid("DeadLayerHole", deadHoleOuterS, deadHoleInnerS);
  auto deadHoleLV = new G4LogicalVolume(deadHoleShell, deadLayer_B, "DeadLayerB");
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (A - C) / 2), deadHoleLV, "DeadLayerB", crystalOuterLV, false, 0, fCheckOverlaps);

  // Endcap well (outer aluminium cup)
  G4double wellOuterR = (F / 2) + t3;
  G4double wellOuterH = E + t4;
  auto endcapWellOuter = new G4Tubs("EndcapWellOuter", 0, wellOuterR, wellOuterH / 2, 0, 360 * deg);
  auto endcapWellInner = new G4Tubs("EndcapWellInner", 0, F / 2, E / 2, 0, 360 * deg);
  auto endcapWellShell = new G4SubtractionSolid("EndcapWell",                       // pName
                                                endcapWellOuter,                    // pSolidA
                                                endcapWellInner,                    // pSolidB
                                                nullptr,                            // rotMatrix
                                                G4ThreeVector(0, 0, (t4) / 2));     // transVector
  auto endcapWellLV = new G4LogicalVolume(endcapWellShell, aluminium, "EndcapWell");

  G4double zWell = A / 2 + G + E / 2;
  // new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zWell), endcapWellLV, "EndcapWell", worldLV, false, 0, fCheckOverlaps);
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), endcapWellLV, "EndcapWell", worldLV, false, 0, fCheckOverlaps);

  // Endcap top and side shell (around crystal)
  auto endcapOuter = new G4Tubs("EndcapOuter", 0, (B / 2 + t2), (A + t1) / 2, 0, 360 * deg);
  auto endcapInner = new G4Tubs("EndcapInner", 0, B / 2, A / 2, 0, 360 * deg);
  auto endcapShell = new G4SubtractionSolid("Endcap", endcapOuter, endcapInner);
  auto endcapLV = new G4LogicalVolume(endcapShell, aluminium, "Endcap");

  G4double zCap = t1 / 2;
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zCap), endcapLV, "Endcap", worldLV, false, 0, fCheckOverlaps);

  // Mount cup (bottom support)
  auto mountCupOuter = new G4Tubs("MountCupOuter", 0, B / 2, t5 / 2, 0, 360 * deg);
  auto mountCupLV = new G4LogicalVolume(mountCupOuter, aluminium, "MountCup");
  G4double zMount = -(A + t5) / 2;
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zMount), mountCupLV, "MountCup", worldLV, false, 0, fCheckOverlaps);

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
  crystalOuterLV->SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
  deadLayerLV->SetVisAttributes(new G4VisAttributes(G4Colour::Blue()));
  deadHoleLV->SetVisAttributes(new G4VisAttributes(G4Colour::Brown()));
  endcapLV->SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
  endcapWellLV->SetVisAttributes(new G4VisAttributes(G4Colour::White()));
  mountCupLV->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));

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

