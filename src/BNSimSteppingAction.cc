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
// $Id: BNSimSteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file BNSimSteppingAction.cc
/// \brief Implementation of the BNSimSteppingAction class

#include "BNSimRun.hh"
#include "BNSimSteppingAction.hh"
#include "BNSimDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

#include "UserDataInput.hh"

#include <algorithm>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

/*
vector<G4double> BNSimSteppingAction::EnergyDepositOf2nd;// = *(new vector<G4double>);
vector<G4int> BNSimSteppingAction::trackIDs;// = *(new vector<G4int>);
*/

BNSimSteppingAction::BNSimSteppingAction(BNSimEventAction* eventAction, BNSimRunAction* runAction)
	: G4UserSteppingAction(),
	fRunAction(runAction),
	fEventAction(eventAction),
	fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimSteppingAction::~BNSimSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BNSimSteppingAction::UserSteppingAction(const G4Step* step)
{
	//G4cout << "TrackStatus: " << trackStatus << G4endl;
	//G4cout << "trackID: " << trackID << G4endl;
	//G4cout << "parentID: " << parentID << G4endl;

	//G4String strPrtclName = theTrack->GetDefinition()->GetParticleName();
	//G4cout << "particalname: " << strPrtclName << G4endl;

	//G4float stepLength = step->GetStepLength();		//cm
	//G4int stepNumber = theTrack->GetCurrentStepNumber();
	//G4StepPoint* prestepPoint = step->GetPreStepPoint();
	//G4float PreKineticEnergy = prestepPoint->GetKineticEnergy(); 	//MeV
	//if (/*parentID != 0 */stepNumber == 1)
	//{
	//	//G4cout << "stepNumber: " << stepNumber << G4endl;
	//	if ((PreKineticEnergy > 0.04 * MeV && PreKineticEnergy < 0.2 * MeV)/* || (PreKineticEnergy > 0.04 * MeV && PreKineticEnergy < 0.05 * MeV)*/)
	//	{
	//		G4cout << "stepNumber: " << stepNumber << G4endl;
	//		G4cout << "PreKineticEnergy: " << PreKineticEnergy / keV << G4endl;
	//		G4String strPrtclName = theTrack->GetDefinition()->GetParticleName();
	//		G4cout << "particalname: " << strPrtclName << G4endl;
	//		G4int trackID = theTrack->GetTrackID();
	//		G4cout << "trackID: " << trackID << G4endl;
	//		G4cout << "parentID: " << parentID << G4endl;
	//		G4String modelName = theTrack->GetCreatorModelName();
	//		G4cout << "modelName: " << modelName << G4endl;
	//		getchar();
	//	}
	//}
	//G4StepPoint* poststepPoint = step->GetPostStepPoint();
	//G4StepStatus poststepStatus = poststepPoint->GetStepStatus();
	//G4float PostKineticEnergy = poststepPoint->GetKineticEnergy(); 	//MeV
	//G4cout << "PostKineticEnergy: " << PostKineticEnergy / keV << G4endl;
	//G4cout << "stepLength: " << stepLength /um << G4endl;
	//G4ThreeVector prepoint = prestepPoint->GetPosition();
	//G4ThreeVector postpoint = poststepPoint->GetPosition();
	//G4cout << "prepoint: " << prepoint / um << G4endl;
	//G4cout << "postpoint: " << postpoint / um << G4endl << G4endl;

	////通过慢化剂后的中子能量
	//G4StepPoint* poststepPoint = step->GetPostStepPoint();
	//G4StepStatus poststepStatus = poststepPoint->GetStepStatus();
	//if (poststepStatus == fGeomBoundary)
	//{
	//    G4TouchableHandle touch = step->GetPreStepPoint()->GetTouchableHandle();
	//    G4String PhyVolName = touch->GetVolume()->GetName();
	//    G4String strPrtclName = theTrack->GetDefinition()->GetParticleName();
	//    if (PhyVolName == "Moderator" && strPrtclName == "neutron")
	//    {
	//        G4ThreeVector postpoint = poststepPoint->GetPosition();
	//        //G4cout << "postpoint: " << postpoint / cm << G4endl << G4endl;
	//        G4double mdrtrz = UserDataInput::GetModeratorThickness();
	//        G4double srctdtctrdstnc = UserDataInput::GetDistanceOfSD();
	//        G4double nz = postpoint.getZ();
	//        if (abs(nz) >= mdrtrz / 2.)
	//        {
	//            //G4cout << "postpoint: " << postpoint  << G4endl;
	//            //G4cout << "nzz: " << - mdrtrz / 2. << G4endl;
	//            //getchar();
	//            G4float PostKineticEnergy = poststepPoint->GetKineticEnergy(); 	//MeV
	//            if (PostKineticEnergy <= 1. * keV)
	//            {
	//                ofstream moderatedneutron;
	//                G4String MDRTEDNTRN = "moderatedneutron.txt";
	//                moderatedneutron.open(MDRTEDNTRN, ios::app);
	//                moderatedneutron << PostKineticEnergy / keV << G4endl;
	//                moderatedneutron.close();
	//            }
	//        }
	//    }
	//}

	////G4double dtctrx = UserDataInput::GetDectorDimensionX();
	////G4double dtctry = UserDataInput::GetDectorDimensionY();
	////G4double dtctrz = UserDataInput::GetDectorDimensionZ();
	////if (abs(postpoint.getZ()) >= dtctrz / 2. || abs(postpoint.getY()) >= dtctry / 2. || abs(postpoint.getX()) >= dtctrx / 2.)
	////    G4cout << "粒子射出探测器" << G4endl;
	////else
	////{
	////    G4cout << "粒子被吸收" << G4endl;
	////    getchar();
	////}

	//if (poststepStatus == fGeomBoundary)
	//{
	//    G4cout << "粒子射出探测器" << G4endl;
	//}
	//else if (poststepStatus != fGeomBoundary && trackStatus == fAlive)
	//{
	//    G4cout << "粒子继续碰撞" << G4endl;
	//}
	//else
	//{
	//    G4cout << "粒子被吸收" << G4endl;
	//    //getchar();
	//}

	if (!fScoringVolume) {
		const BNSimDetectorConstruction* detectorConstruction
			= static_cast<const BNSimDetectorConstruction*>
			(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		fScoringVolume = detectorConstruction->GetScoringVolume();
	}

	// get volume of the current step
	G4TouchableHandle touch = step->GetPreStepPoint()->GetTouchableHandle();
	G4LogicalVolume* volume = touch->GetVolume()->GetLogicalVolume();

	// check if we are in scoring volume
	if (volume != fScoringVolume) return;

	G4String logVolName = volume->GetName();
	if (logVolName != "EmptyDetector")
	{
		G4Track* theTrack = step->GetTrack();
		G4int parentID = theTrack->GetParentID();

		// 判断是否为中子
		G4String strPrtclName = theTrack->GetDefinition()->GetParticleName();
		G4TrackStatus trackStatus = theTrack->GetTrackStatus();
		if (strPrtclName == "neutron")
		{
			// fRunAction->AddNeutronNumber();
			if (trackStatus == fStopAndKill) //若中子在探测器内终止则中子俘获计数+1
			{
				//getchar();
				fRunAction->AddNeutronCount();
			}
		}

		//if (strPrtclName == "alpha" || strPrtclName == "Li7")
		//{
		//	// if (strPrtclType ==)
		//	// collect energy deposited in this step
		//	G4double edepStep = step->GetTotalEnergyDeposit() / keV;
		//
		//	//统计每个event的能量沉积，即中子反应产生的所有次级粒子的能量沉积的总和
		//	fEventAction->AddEdep(edepStep);
		//
		//	G4cout << "particalname: " << strPrtclName << G4endl;
		//	G4int stepNumber = theTrack->GetCurrentStepNumber();
		//	G4cout << "stepNumber: " << stepNumber << G4endl;
		//	G4cout << "edepStep: " << edepStep << G4endl << G4endl;
		//	getchar();
		//}

		// if (strPrtclType ==)
		// collect energy deposited in this step
		G4double edepStep = step->GetTotalEnergyDeposit() / keV;

		//BNSimRun* fBNSimRun
		//	= static_cast<BNSimRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
		//G4double tt = step->GetPreStepPoint()->GetPosition().getX();
		//fBNSimRun->PushBackEnergyDeposit(tt);

		//统计每个event的能量沉积，即中子反应产生的所有次级粒子的能量沉积的总和
		fEventAction->AddEdep(edepStep);

		//G4cout << "edepStep: " << edepStep << G4endl << G4endl;
		//while (edepStep != 0.)
		//{
		//    getchar();
		//    break;
		//}
	}
	else
	{
		G4StepPoint* poststepPoint = step->GetPostStepPoint();
		G4StepStatus poststepStatus = poststepPoint->GetStepStatus();
		if (poststepStatus == fGeomBoundary)
		{
			G4Track* theTrack = step->GetTrack();
			G4String strPrtclName = theTrack->GetDefinition()->GetParticleName();
			if (strPrtclName == "neutron")
			{
				G4double postKineticEnergy = poststepPoint->GetKineticEnergy(); 	//MeV
				BNSimRun* fBNSimRun
					= static_cast<BNSimRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
				fBNSimRun->PushbackNeutronEnergy(postKineticEnergy);

				//G4cout << "strPrtclName:" << strPrtclName << G4endl;
				//G4cout << "postKineticEnergy:" << postKineticEnergy << G4endl;
				//getchar();
			}
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
