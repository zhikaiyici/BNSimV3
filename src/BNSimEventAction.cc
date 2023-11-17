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
// $Id: BNSimEventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file BNSimEventAction.cc
/// \brief Implementation of the BNSimEventAction class

#include "BNSimAnalysis.hh"
#include "BNSimRun.hh"
#include "BNSimEventAction.hh"
#include "BNSimRunAction.hh"
#include "BNSimSteppingAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimEventAction::BNSimEventAction(BNSimRunAction* runAction)
	: G4UserEventAction(),
	fRunAction(runAction),
	fEdep(0.),
	neutronCount(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimEventAction::~BNSimEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BNSimEventAction::BeginOfEventAction(const G4Event* /*event*/)
{
	fEdep = 0.;

	//G4cout << "----------------------BeginOfEventAction------------------------------" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BNSimEventAction::EndOfEventAction(const G4Event* event)
{
	/*// get analysis manager
	auto analysisManager = G4AnalysisManager::Instance();

	if (fEdep != 0)
	{
		// fill histograms
		// analysisManager->FillH1(0, fEdep);

		// fill ntuple
		analysisManager->FillNtupleDColumn(0, fEdep);
		analysisManager->AddNtupleRow();
	}*/

	G4String detectorMaterial = UserDataInput::GetDetectorMaterial();
	if (detectorMaterial != "AIR")
	{
		// accumulate statistics in B1Run
		BNSimRun* fBNSimRun
			= static_cast<BNSimRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
		if (fEdep != 0)
		{
			fBNSimRun->PushBackEnergyDeposit(fEdep);
		}
		/*if (fEdep != 0)
		{
			fRunAction->PushBackEnergyDeposit(fEdep);
		}*/
	}

	G4int eventID = event->GetEventID();
	G4int numberofevents = UserDataInput::GetNumberOfEvents();
	// G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
	if (eventID == 0 || ((eventID + 1) % (numberofevents / 10) == 0))
	{
		if (eventID == 0)
		{
			G4int numofevents = UserDataInput::GetNumberOfEvents();
			G4double doubNumOfEvents = numofevents;
			G4double dtctrx = UserDataInput::GetDectorDimensionX() / mm;
			G4double dtctry = UserDataInput::GetDectorDimensionY() / mm;
			G4double dtctrz = UserDataInput::GetDectorDimensionZ() / um;
			G4cout << G4endl << " Initialization completed. " << doubNumOfEvents << " event(s) will be simulated."
				   << G4endl
				   << " The dimension of detector is " << dtctrx << " mm x " << dtctry << " mm x " << dtctrz << " um"
				   << G4endl
				   << G4endl;
		}
		G4int per =(int) ((1. * eventID + 1) / (numberofevents * 0.01));
		//G4cout << " numberofevents: "<< numberofevents << " eventID: "<< eventID <<G4endl;
		G4long seconds = time(NULL); // 格林威治时间
		seconds = seconds + 8 * 3600; // 北京时间
		G4int secondnow = seconds % 60;
		G4int minutes = (seconds - secondnow) / 60;
		G4int minutenow = minutes % 60;
		G4int hours = (minutes - minutenow) / 60;
		G4int hournow = hours % 24;
		G4cout << " Time now: " << setw(2) << hournow << ":" << setw(2) << minutenow << ":" << setw(2) << secondnow
			<< ". " << setw(3) << per << "% of simulation completed."
			<< G4endl;
		//getchar();
	}
	//G4cout << "-----------------------EndOfEventAction------------------------------" << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
