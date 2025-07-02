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
/// \file B4/B4a/include/ActionInitialization.hh
#ifndef B4A_ACTIONINITIALIZATION_H
#define B4A_ACTIONINITIALIZATION_H

#include "G4VUserActionInitialization.hh"

namespace B4 {
  class DetectorConstruction;
}

namespace B4a {

class ActionInitialization : public G4VUserActionInitialization
{
public:
  explicit ActionInitialization(const B4::DetectorConstruction* det);
  virtual ~ActionInitialization() = default;

  virtual void Build() const override;
  virtual void BuildForMaster() const override;

private:
  const B4::DetectorConstruction* fDetConstruction;
};

}  // namespace B4a

#endif  // B4A_ACTIONINITIALIZATION_H