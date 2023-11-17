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
// $Id: BNSimRunAction.hh 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file BNSimRunAction.hh
/// \brief Definition of the BNSimRunAction class

#ifndef BNSimRunAction_h
#define BNSimRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"
// #include "ListAccumulable.hh"

#include <vector>
#include <list>
#include <deque>

using namespace std;

class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class BNSimRunAction : public G4UserRunAction
{
public:
	BNSimRunAction();
	virtual ~BNSimRunAction();

	virtual G4Run* GenerateRun();
	virtual void BeginOfRunAction(const G4Run*);
	virtual void   EndOfRunAction(const G4Run*);

	// inline void PushBackEnergyDepositOfEvent(G4double edep) { energyDeposit.push_back(edep); };
	inline void AddNeutronCount() { fNeutronCount += 1; }
	inline void AddNeutronNumber() { fNeutronNumber += 1; }

private:
	//G4Accumulable<G4double> fEdep;
	//G4Accumulable<G4double> fEdep2;

	// static list<G4double> energyDeposit;
	// deque<G4double> EnergyDeposit;
	//static G4int neutronCount;
	//static G4int neutronNumber;
	G4Accumulable<G4int> fNeutronCount;
	G4Accumulable<G4int> fNeutronNumber;
	G4String runCondition;
	G4String runConditionWithoutThickness;
};

#endif
