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
// $Id: BNSimDetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file BNSimDetectorConstruction.cc
/// \brief Implementation of the BNSimDetectorConstruction class

#include "BNSimDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"
#include "G4BooleanSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimDetectorConstruction::BNSimDetectorConstruction()
	: G4VUserDetectorConstruction(),
	fScoringVolume(0), checkOverlaps(true), userData(0), 
	dtctrx(userData->GetDectorDimensionX()),
	dtctry(userData->GetDectorDimensionY()),
	dtctrz(userData->GetDectorDimensionZ()),
	sourceDetectorDistance(userData->GetDistanceOfSD()),
	moderatorThickness(userData->GetModeratorThickness())
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BNSimDetectorConstruction::~BNSimDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* BNSimDetectorConstruction::Construct()
{
	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();

	//     
	// World
	//

	G4double world_sizeX = 1.1 * (dtctrx + 2. * moderatorThickness);
	G4double world_sizeY = 1.1 * (dtctry + 2. * moderatorThickness);
	G4double world_sizeZ = 1.1 * (sourceDetectorDistance + dtctrz + 2. * moderatorThickness);
	G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

	G4Box* solidWorld =
		new G4Box("World",                       //its name
			0.5 * world_sizeX, 0.5 * world_sizeY, 0.5 * world_sizeZ);     //its size

	G4LogicalVolume* logicWorld =
		new G4LogicalVolume(solidWorld,          //its solid
			world_mat,           //its material
			"World");            //its name

	G4VPhysicalVolume* physWorld =
		new G4PVPlacement(0,                     //no rotation
			G4ThreeVector(),       //at (0,0,0)
			logicWorld,            //its logical volume
			"World",               //its name
			0,                     //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking 

	G4bool moderatorStatus = userData->GetModeratorStatus();
	if (moderatorStatus != FALSE)
		ConstructModerator(logicWorld);

	G4LogicalVolume* logicDetectorContainer = ConstructDetectorContainer(logicWorld);

	G4String detectorMaterial = userData->GetDetectorMaterial();
	if (detectorMaterial != "AIR")
		ConstructDetector(logicDetectorContainer);
	else
		ConstructEmptyDetector(logicDetectorContainer);

	G4cout << *(G4Material::GetMaterialTable()) << G4endl;

	return physWorld;
}

void BNSimDetectorConstruction::ConstructModerator(G4LogicalVolume* logicWorld)
{
	//
	// 慢化剂HDPE
	//

	G4double z, a, density;
	G4String name, symbol;
	G4int numberofelements, numberofatoms;

	// 定义HDPE
	G4Element* elC = new G4Element(name = "Carbon",
		symbol = "C",
		z = 6.0, a = 12.0107 * g / mole);
	G4Element* elH = new G4Element(name = "Hydrogen",
		symbol = "H",
		z = 1.0, a = 1.00784 * g / mole);

	G4Material* HDPE = new G4Material(name = "HDPE",
		density = 0.95 * g / cm3,
		numberofelements = 2);
	HDPE->AddElement(elH, numberofatoms = 2);
	HDPE->AddElement(elC, numberofatoms = 1);

	// G4double pi = 3.14159265358979323846;
	G4double rAnti = 0.5 * sqrt(dtctrx * dtctrx + dtctry * dtctry + dtctrz * dtctrz);
	G4double rOrb = moderatorThickness/* + rAnti*/;
	G4double zAnti = rOrb + 2 * rAnti;

	G4Orb* solidOrbModerator = new G4Orb("OrbModerator", rOrb);
	G4EllipticalTube* solidAntiModerator = new G4EllipticalTube("AntiModerator", rAnti, rAnti, 0.5 * zAnti);
	G4SubtractionSolid* solidModerator =
		new G4SubtractionSolid("Moderator", solidOrbModerator, solidAntiModerator, 
							   G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., 0.5 * rOrb)));

	G4LogicalVolume* logicModerator = new G4LogicalVolume(solidModerator, HDPE, "Moderator");
	logicModerator->SetVisAttributes(new G4VisAttributes(G4Color(1., 1., 1., 0.6)));

	G4RotationMatrix* rotModerator = new G4RotationMatrix;
	rotModerator->rotateX(pi / 2.);
	G4ThreeVector moderatorPosition;
	G4String strModeratorPosition = userData->GetModeratorPosition();
	if (strModeratorPosition == "Detector")
		moderatorPosition = G4ThreeVector(0., 0., -0.5 * (dtctrz + sourceDetectorDistance));
	else
		moderatorPosition = G4ThreeVector(0., 0., 0.5 * sourceDetectorDistance);

	G4VPhysicalVolume* physModerator =
		new G4PVPlacement(rotModerator, moderatorPosition, 
						  logicModerator, "Moderator", logicWorld, false, 0, checkOverlaps);
}

G4LogicalVolume* BNSimDetectorConstruction::ConstructDetectorContainer(G4LogicalVolume* logicWorld)
{
	G4NistManager* nist = G4NistManager::Instance();
	G4Material* containerMat = nist->FindOrBuildMaterial("G4_Fe");
	G4Material* antiContainerMat = nist->FindOrBuildMaterial("G4_AIR");

	G4double containerX = 1. * cm;
	G4double containerY = 1. * cm;
	G4double containerZ = 1. * cm;
	G4double containerThickness = 0.2 * cm;

	G4Box* solidDetectorContainerSolid = 
		new G4Box("DetectorContainerSolid", 0.5 * containerX, 0.5 * containerY, 0.5 * containerZ);
	G4Box* solidAntiDetectorContainer = 
		new G4Box("AntiDetectorContainer", 0.5 * (containerX - containerThickness), 0.5 * (containerX - containerThickness), 0.5 * (containerX - containerThickness));
	//G4SubtractionSolid* solidDetectorContainer =
	//	new G4SubtractionSolid("DetectorContainer", solidDetectorContainerSolid, solidAntiDetectorContainer);
	G4LogicalVolume* logicDetectorContainerSolid = 
		new G4LogicalVolume(solidDetectorContainerSolid, containerMat, "DetectorContainerSolid");
	G4LogicalVolume* logicAntiDetectorContainer =
		new G4LogicalVolume(solidAntiDetectorContainer, antiContainerMat, "AntiDetectorContainer");

	//logicDetectorContainer->SetVisAttributes(new G4VisAttributes(G4Color(300. / 255.0, 200. / 255.0, 100. / 255.0)));

	G4VPhysicalVolume* physDetectorContainerSolid =
		new G4PVPlacement(0, G4ThreeVector(0., 0., -0.5 * (sourceDetectorDistance + containerZ)),
						  logicDetectorContainerSolid, "DetectorContainerSolid", 
						  logicWorld, false, 0, checkOverlaps);
	G4VPhysicalVolume* physAntiDetectorContainer =
		new G4PVPlacement(0, G4ThreeVector(),
			logicAntiDetectorContainer, "AntiDetectorContainer",
			logicDetectorContainerSolid, false, 0, checkOverlaps);

	return logicAntiDetectorContainer;
}

void BNSimDetectorConstruction::ConstructDetector(G4LogicalVolume* logicDetectorContainer)
{
	G4double z, a, density;
	G4String name, symbol;
	G4int numberofelements, numberofatoms;

	//
	// 探测器几何
	//

	// 定义探测器材料hBN

	G4Element* elN = new G4Element(name = "Nitrogen",
		symbol = "N",
		z = 7.0, a = 14.0067 * g / mole);
	G4Element* elB;

	G4String detectorMaterial = userData->GetDetectorMaterial();
	if (detectorMaterial == "B10")
	{
		// B-10
		/* 构建指定丰度的某元素，只能通过相对丰度，由同位素定义，
		 * 不能通过指定a来定义，即不能使用下面的方法，
		 * 下面的方法貌似还是该元素的天然丰度。
		elB = new G4Element(name = "Boron10",
			symbol = "B10",
			z = 5.0, a = 10.012926862 * g / mole);*/

		G4int iz, in, numberofisotopes;
		G4Isotope* isoB10 = new G4Isotope(name = "B10",
			iz = 5, in = 10, a = 10.012937068613 * g / mole);
		elB = new G4Element(name = "Boron",
			symbol = "B",
			numberofisotopes = 1);
		elB->AddIsotope(isoB10, 100. * perCent);

		// h-B10N
		density = 2.213 * g / cm3;
	}
	else if (detectorMaterial == "B11")
	{
		// B-11
		/*elB = new G4Element(name = "Boron11",
			symbol = "B11",
			z = 5.0, a = 11.0093 * g / mole);*/

		G4int iz, in, numberofisotopes;
		G4Isotope* isoB11 = new G4Isotope(name = "B11",
			iz = 5, in = 11, a = 11.0093054826847 * g / mole);
		elB = new G4Element(name = "Boron",
			symbol = "B",
			numberofisotopes = 1);
		elB->AddIsotope(isoB11, 100. * perCent);

		// h-B11N
		density = 2.305 * g / cm3;
	}
	else
	{
		// 天然B
		elB = new G4Element(name = "Boron",
			symbol = "B",
			z = 5.0, a = 10.811 * g / mole);

		//G4int iz, in, numberofisotopes;
		//
		//// 数据来源
		//// a: http://www.nuclear.csdb.cn/texing.html
		//// 丰度：https://www.nndc.bnl.gov/nudat2/
		//
		//G4Isotope* isoB10 = new G4Isotope(name = "B10",
		//    iz = 5, in = 10, a = 10.012937068613 * g/mole);
		//G4Isotope* isoB11 = new G4Isotope(name = "B11",
		//    iz = 5, in = 11, a = 11.0093054826847 * g/mole);
		//G4Element* elB = new G4Element(name = "Boron",
		//    symbol = "B",
		//    numberofisotopes = 2);
		//elB->AddIsotope(isoB10, 19.9 * perCent);
		//elB->AddIsotope(isoB11, 80.1 * perCent);
		//
		//G4Isotope* isoN14 = new G4Isotope(name = "N14",
		//    iz = 7, in = 14, a = 14.0030739869773 * g/mole);
		//G4Isotope* isoN15 = new G4Isotope(name = "N15",
		//    iz = 7, in = 15, a = 15.0001088574001 * g/mole);
		//G4Element* elN = new G4Element(name = "Nitrogen",
		//    symbol = "N",
		//    numberofisotopes = 2);
		//elN->AddIsotope(isoN14, 99.636 * perCent);
		//elN->AddIsotope(isoN15, 0.364 * perCent);

		//G4Element* elB = nist->FindOrBuildElement(z = 5);
		//G4Element* elN = nist->FindOrBuildElement(z = 7);

		//h-BN
		density = 2.287 * g / cm3;
	}

	G4Material* hBN = new G4Material(name = "h-BN",
		density,
		numberofelements = 2);
	hBN->AddElement(elB, numberofatoms = 1);
	hBN->AddElement(elN, numberofatoms = 1);

	//G4double dtctrx = userData.GetDectorDimensionX();
	//G4double dtctry = userData.GetDectorDimensionY();
	//G4double dtctrz = userData.GetDectorDimensionZ();
	//G4double sourceDetectorDistance = userData.GetDistanceOfSD();

	G4Box* solidDetector = new G4Box("Detector", 0.5 * dtctrx, 0.5 * dtctry, 0.5 * dtctrz);

	G4LogicalVolume* logicDetector = new G4LogicalVolume(solidDetector, hBN, "Detector");

	logicDetector->SetVisAttributes(new G4VisAttributes(G4Color(300. / 255.0, 200. / 255.0, 100. / 255.0)));

	G4VPhysicalVolume* physDetector =
		new G4PVPlacement(0, G4ThreeVector(/*0., 0., -0.5 * (1. * cm + sourceDetectorDistance)/*母体为世界时的坐标*/), 
						  logicDetector, "Detector", logicDetectorContainer, false, 0, checkOverlaps);

	// Set Detector as scoring volume 
	fScoringVolume = logicDetector;
}

void BNSimDetectorConstruction::ConstructEmptyDetector(G4LogicalVolume* logicDetectorContainer)
{
	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();
	G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

	G4double density = 1.e-5 * g / cm3;
	G4double pressure = 2.e-2 * bar;
	G4double temperature = STP_Temperature; // from PhysicalConstants.h
	G4double fractionmass;
	G4int numberofcomponents;
	G4String name;
	G4Material* vacuum = 
		new G4Material(name = "Vacuum ", density, numberofcomponents = 1, kStateGas, temperature, pressure);
	vacuum->AddMaterial(air, fractionmass = 1.);

	G4Box* solidEmptyDetector = new G4Box("EmptyDetector", 0.5 * dtctrx, 0.5 * dtctry, 0.5 * dtctrz);

	G4LogicalVolume* logicEmptyDetector = new G4LogicalVolume(solidEmptyDetector, vacuum, "EmptyDetector");

	logicEmptyDetector->SetVisAttributes(new G4VisAttributes(G4Color(50. / 255.0, 200. / 255.0, 100. / 255.0)));

	G4VPhysicalVolume* physDetector =
		new G4PVPlacement(0, G4ThreeVector(/*0., 0., -0.5 * (dtctrz + sourceDetectorDistance)/*母体为世界时的坐标*/),
						  logicEmptyDetector, "EmptyDetector", logicDetectorContainer, false, 0, checkOverlaps);

	// Set Detector as scoring volume 
	fScoringVolume = logicEmptyDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
