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
// $Id: B1Run.cc 66536 2012-12-19 14:32:36Z ihrivnac $
//
/// \file B1Run.cc
/// \brief Implementation of the B1Run class

#include "BNSimRun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimRun::BNSimRun()
: G4Run(),
  energyDeposit(0),
  neutronEnergy(0)
  /*fEdep(0.), 
  fEdep2(0.)*/
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimRun::~BNSimRun()
{} 
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BNSimRun::Merge(const G4Run* run)
{
	// G4cout << "---------------------------Merge-----------------------" << G4endl;
	// getchar();

	const BNSimRun* localRun = static_cast<const BNSimRun*>(run);

	std::list<G4double> localEnergyDeposit = localRun->energyDeposit;
	 //G4cout << "size of eneryDeposit before merge: " << energyDeposit.size() << G4endl;
	 //G4cout << "size of localEneryDeposit: " << localEnergyDeposit.size() << G4endl;
	 //for (std::list<G4double>::iterator itr = energyDeposit.begin(); itr != energyDeposit.end(); ++itr)
	 //{
	 //  if (*itr != 0)
	 //  {
	 //	  G4cout << *itr << G4endl;
	 //  }
	 //}
	energyDeposit.merge(localEnergyDeposit);
	// G4cout << "size of eneryDeposit after merge: " << energyDeposit.size() << G4endl;
	// getchar();
	// for (std::list<G4double>::iterator itr = energyDeposit.begin(); itr != energyDeposit.end(); ++itr)
	// {
	//     if (*itr != 0)
	//     {
	//   	  G4cout << *itr << G4endl;
	//     }
	// }	

	std::list<G4double> localNeutronEnergy = localRun->neutronEnergy;
	neutronEnergy.merge(localNeutronEnergy);

	/*fEdep  += localRun->fEdep;
	fEdep2 += localRun->fEdep2;*/

	G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BNSimRun::PushBackEnergyDeposit(G4double edep)
{
	energyDeposit.push_back(edep);
	//for (std::list<G4double>::iterator itr = energyDeposit.begin(); itr != energyDeposit.end(); ++itr)
	//{
	//	G4cout << *itr << G4endl;
	//}
}

void BNSimRun::PushbackNeutronEnergy(G4double nenergy)
{
	neutronEnergy.push_back(nenergy);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void BNSimRun::AddEdep (G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


