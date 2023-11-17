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
// $Id: BNSimRunAction.cc 99560 2016-09-27 07:03:29Z gcosmo $
//
/// \file BNSimRunAction.cc
/// \brief Implementation of the BNSimRunAction class

#include "BNSimRun.hh"
#include "BNSimRunAction.hh"
#include "BNSimPrimaryGeneratorAction.hh"
#include "BNSimDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "BNSimAnalysis.hh"
#include "BNSimSteppingAction.hh"
#include "BNSimEventAction.hh"
#include "UserDataInput.hh"

#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
/*#include <direct.h>
#include <io.h>*/

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//G4Accumulable<G4int> BNSimRunAction::fNeutronCount = 0;
//G4int BNSimRunAction::neutronCount = 0;
//G4int BNSimRunAction::neutronNumber = 0;
// list<G4double> BNSimRunAction::energyDeposit = *(new list<G4double>);

BNSimRunAction::BNSimRunAction()
	: G4UserRunAction(),
	fNeutronCount(0),
	fNeutronNumber(0)

{
	G4int nofEvents = UserDataInput::GetNumberOfEvents();
	G4double dtctrz = UserDataInput::GetDectorDimensionZ() / um;
	G4bool moderatorStatus = UserDataInput::GetModeratorStatus();
	G4String detectorMaterial = UserDataInput::GetDetectorMaterial();
	G4double gammaPercentage = UserDataInput::GetGammaPercentage() * 100;
	G4double numofevent = nofEvents;
	ostringstream ostrsThickness, ostrsEventNumber, ostrsGammaPercentage;
	ostrsThickness << dtctrz;
	ostrsEventNumber << setprecision(1) << numofevent;
	ostrsGammaPercentage << gammaPercentage;
	G4String strThickness = ostrsThickness.str();
	G4String strEventNumber = ostrsEventNumber.str();
	G4String strGammaPercentage = ostrsGammaPercentage.str();
	G4String strModeratorStatus;
	if (moderatorStatus != TRUE)
	{
		strModeratorStatus = "OFF";
	}
	else
	{
		G4double moderatorThickness = UserDataInput::GetModeratorThickness() / cm;
		ostringstream ostrsModeratorThickness; ostrsModeratorThickness << moderatorThickness;
		G4String strModeratorThickness = ostrsModeratorThickness.str();
		G4String strModeratorPosition = UserDataInput::GetModeratorPosition();
		strModeratorStatus = "ON" + strModeratorThickness + "cm_" + strModeratorPosition;
	}
	runCondition ="(" + strThickness + "um_" + strEventNumber
		+ "_moder" + strModeratorStatus + "_gamma" + strGammaPercentage + "%_" + detectorMaterial + ")";
	runConditionWithoutThickness = "(" + strEventNumber
		+ "_moder" + strModeratorStatus + "_gamma" + strGammaPercentage + "%_" + detectorMaterial + ")";

	// set printing event number per each event
	// G4RunManager::GetRunManager()->SetPrintProgress(nofEvents / 10);

	/*// Create analysis manager
	// The choice of analysis technology is done via selectin of a namespace
	// in B4Analysis.hh
	auto analysisManager = G4AnalysisManager::Instance();
	G4cout << "Using " << analysisManager->GetType() << G4endl;

	// Create directories 
	//analysisManager->SetHistoDirectoryName("histograms");
	//analysisManager->SetNtupleDirectoryName("ntuple");
	analysisManager->SetVerboseLevel(1);
	analysisManager->SetNtupleMerging(true);
	// Note: merging ntuples is available only with Root output

	// Book histograms, ntuple
	//

	// Creating histograms
	// analysisManager->CreateH1("Edep", "Edep in BN", 2000, 0., 20 * MeV);

	// Creating ntuple
	//
	analysisManager->CreateNtuple("Edep", "Edep in BN");
	analysisManager->CreateNtupleDColumn("Edep in BN");
	analysisManager->FinishNtuple();*/

	// Register accumulable to the accumulable manager
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->RegisterAccumulable(fNeutronCount);
	accumulableManager->RegisterAccumulable(fNeutronNumber);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimRunAction::~BNSimRunAction()
{
	// delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* BNSimRunAction::GenerateRun()
{
	return new BNSimRun;
}

void BNSimRunAction::BeginOfRunAction(const G4Run* /*run*/)
{
	// inform the runManager to save random number seed
	// G4RunManager::GetRunManager()->SetRandomNumberStore(false);

	/*// Get analysis manager
	auto analysisManager = G4AnalysisManager::Instance();

	// Open an output file
	//
	G4String edepFileName =
		"output/edep" + runCondition;
	analysisManager->OpenFile(edepFileName);*/

	// reset accumulables to their initial values
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->Reset();

	// energyDeposit.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BNSimRunAction::EndOfRunAction(const G4Run* run)
{
	G4int nofEvents = run->GetNumberOfEvent();
	if (nofEvents == 0) return;

	G4String detectorMaterial = UserDataInput::GetDetectorMaterial();
	if (detectorMaterial != "AIR")
	{
		// Merge accumulables
		G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
		accumulableManager->Merge();

		//G4int neutronCount = fNeutronCount.GetValue();

		/*// 判断输出目标文件夹是否存在，不存在则创建
		if (_access("output", 4) != 0)
		{
			_mkdir("output");
			//G4cout << "已创建output" << G4endl;
		}*/

		// Print
		//  
		if (IsMaster()) {
			/*// Get analysis manager
			auto analysisManager = G4AnalysisManager::Instance();

			// save histograms & ntuple
			//
			analysisManager->Write();
			analysisManager->CloseFile();*/

			const BNSimRun* fBNSimRun = static_cast<const BNSimRun*>(run);
			list<G4double> energyDeposit = fBNSimRun->GetEnergyDeposit();

			ofstream edepFile, captureEfficiencyFile;
			G4String edepFileName = "output/edep" + runCondition + ".txt";
			edepFile.open(edepFileName, ios::out);
			for (list<G4double>::iterator itr = energyDeposit.begin(); itr != energyDeposit.end(); ++itr)
			{
				if (*itr != 0)
				{
					edepFile << *itr << G4endl;
				}
			}
			edepFile.close();

			G4int neutronCount = fNeutronCount.GetValue();
			G4double gammaPercentage = UserDataInput::GetGammaPercentage();
			G4int neutronNumber = nofEvents * (1 - gammaPercentage); // fNeutronNumber.GetValue();

			G4double neutronCaptureEfficiency = 100. * neutronCount / neutronNumber;
			G4double dtctrz = UserDataInput::GetDectorDimensionZ() / um;

			G4String captureEfficiencyFileName =
				"output/effofneucap" + runConditionWithoutThickness + ".txt";
			captureEfficiencyFile.open(captureEfficiencyFileName, ios_base::app);
			captureEfficiencyFile << dtctrz << setw(10) << neutronCaptureEfficiency << "%" << G4endl;
			captureEfficiencyFile.close();

			G4cout
				<< G4endl
				<< "--------------------End of Global Run-----------------------";
			G4cout
				<< G4endl
				<< " The count of captured neutron is " << neutronCount
				<< G4endl
				<< " The efficiency of neutron capture is " << neutronCaptureEfficiency << "%"
				<< G4endl
				<< "------------------------------------------------------------"
				<< G4endl
				<< G4endl;
		}
		else {
			G4cout
				<< G4endl
				<< "--------------------End of Local Run------------------------"
				<< G4endl;
		}
	}
	else
	{
		const BNSimRun* fBNSimRun = static_cast<const BNSimRun*>(run);
		list<G4double> neutronEnergy = fBNSimRun->GetNeutronEnergy();

		ofstream neutronEnergyFile;
		G4String neutronEnergyFileName = "output/kn" + runCondition + ".txt";
		neutronEnergyFile.open(neutronEnergyFileName, ios::out);
		for (list<G4double>::iterator itr = neutronEnergy.begin(); itr != neutronEnergy.end(); ++itr)
		{
			if (*itr != 0)
			{
				neutronEnergyFile << *itr << G4endl;
			}
		}
		neutronEnergyFile.close();
	}
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
