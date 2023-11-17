
#include "UserDataInput.hh"

#include "G4SystemOfUnits.hh"

#include <fstream>
#include <iomanip>

using namespace std;

G4int UserDataInput::numberOfThreads = 4;
G4int UserDataInput::numberofevents = 1000;

G4bool UserDataInput::uiStatus = TRUE;
G4bool UserDataInput::moderatorStatus = FALSE;

G4String UserDataInput::source = "Point";
G4String UserDataInput::detectorMaterial = "Natural";
G4String UserDataInput::moderatorPosition = "Detector";
/*G4String UserDataInput::gammaSpectrumName = "gammaspectrum.spc";
G4String UserDataInput::neutronSpectrumName = "neutronspectrum.spc";*/

G4double UserDataInput::sourceDetectorDistance = 20.;
G4double UserDataInput::dtctrx = 2.;
G4double UserDataInput::dtctry = 2.;
G4double UserDataInput::dtctrz = 50.;
//G4double UserDataInput::StartThickness = 5.;
//G4double UserDataInput::ThicknessIncrement = 5.;
//G4double UserDataInput::EndThickness = 50.;
G4double UserDataInput::gammaPercentage = 0.5;
G4double UserDataInput::moderatorThickness = 0.9 * sourceDetectorDistance;

vector<G4double> UserDataInput::gammaEnergy = *(new vector<G4double>); //粒子能量MeV
vector<G4double> UserDataInput::gammaNormSpectrum = *(new vector<G4double>); //归一化能谱
vector<G4double> UserDataInput::neutronEnergy = *(new vector<G4double>); //粒子能量MeV
vector<G4double> UserDataInput::neutronNormSpectrum = *(new vector<G4double>); //归一化能谱

UserDataInput::UserDataInput()
{
}

UserDataInput::~UserDataInput()
{
}

void UserDataInput::ReadInputData()
{
	G4String prmtrFile = "InputSet.txt";
	ifstream prameterFile(prmtrFile, ios_base::in);		//输入文件
	if (prameterFile.good() != 1)
		G4cout << "FAILED to open file " << prmtrFile << G4endl;
	else
	{
		G4cout << "File " << prmtrFile << " has been opened SUCCESSFULLY." << G4endl;

		G4String input;
		G4String gammaSpectrumName; //gamma能谱文件名
		G4String neutronSpectrumName; //中子能谱文件名

		while (prameterFile >> input)			//读取参数
		{
			if (input == "Number_Of_Events:")
			{
				G4double doubNumberOfEvents;
				prameterFile >> doubNumberOfEvents;
				numberofevents = (G4int)doubNumberOfEvents;
				G4cout << doubNumberOfEvents << " events will be simulated." << G4endl;
			}
			else if (input == "Number_Of_Threads:")
			{
				prameterFile >> numberOfThreads;
			}
			else if (input == "Source_Type:")
			{
				prameterFile >> source;
				G4cout << "The source tpye is " << source << " source." << G4endl;
			}
			else if (input == "Distance_of_SD:")
			{
				prameterFile >> sourceDetectorDistance;
				G4cout << "The distance of SD is " << sourceDetectorDistance << " cm." << G4endl;
			}
			else if (input == "XYZ:")
			{
				prameterFile >> dtctrx >> dtctry >> dtctrz;
				G4cout << "The dimension of detector is " << dtctrx << " mm x " << dtctry << " mm x " << dtctrz << " um" << G4endl;
			}
			else if (input == "Material_Of_Detector:")
			{
				prameterFile >> detectorMaterial;
				G4cout << "The detector material is " << detectorMaterial << " h-BN." << G4endl;
			}
			else if (input == "Gamma_Spectrum_File_Name:")
			{
				prameterFile >> gammaSpectrumName;
			}
			else if (input == "Neutron_Spectrum_File_Name:")
			{
				prameterFile >> neutronSpectrumName;
			}
			else if (input == "Percentage_Of_Gamma:")
			{
				prameterFile >> gammaPercentage;
				G4cout << "The Percentage_Of_Gamma is " << gammaPercentage << G4endl;
			}
			else if (input == "HDPE_Moderator:")
			{
				G4String HDPE_Moderator;
				prameterFile >> HDPE_Moderator;
				if (HDPE_Moderator == "ON")
				{
					moderatorStatus = TRUE;
					G4cout << "The moderator is ON." << G4endl;
				}
				else
				{
					moderatorStatus = FALSE;
					G4cout << "The moderator is OFF." << G4endl;
				}
			}
			else if (input == "Thickness_Of_Moderator:")
			{
				prameterFile >> moderatorThickness;
				G4cout << "The moderator thickness is " << moderatorThickness << " cm." << G4endl;
			}
			else if (input == "Position_Of_Moderator:")
			{
				prameterFile >> moderatorPosition;
				G4cout << "The moderator position is " << moderatorPosition << G4endl; 
			}
			else if (input == "Status_Of_UI:")
			{
				G4String Status_Of_UI;
				prameterFile >> Status_Of_UI;
				if (Status_Of_UI == "ON")
				{
					uiStatus = TRUE;
				}
				else
				{
					uiStatus = FALSE;
				}
			}
			//else if (input == "Start_Thickness:")
			//{
			//	prameterFile >> StartThickness;
			//	G4cout << "" << StartThickness << " um" << G4endl;
			//}
			//else if (input == "Increment_Of_Thickness:")
			//{
			//	prameterFile >> ThicknessIncrement;
			//	G4cout << "" << ThicknessIncrement << " um" << G4endl;
			//}
			//else if (input == "End_Thickness:")
			//{
			//	prameterFile >> EndThickness;
			//	G4cout << "" << EndThickness << " um" << G4endl;
			//}
		}
		prameterFile.close();

		ReadSpectra(gammaSpectrumName, gammaEnergy, gammaNormSpectrum);
		ReadSpectra(neutronSpectrumName, neutronEnergy, neutronNormSpectrum);
	}
}


void UserDataInput::ReadSpectra(G4String spectrumName, vector<G4double>& energy, vector<G4double>& normSpectrum)
{
	ifstream spcFile;
	spectrumName = "spectra/" + spectrumName;

	spcFile.open(spectrumName, ios_base::in);

	if (spcFile.good() != 1)
		G4cout << "FAILED to open spectrum file " << spectrumName << G4endl;
	else
	{
		G4cout << "Spectrum file " << spectrumName << " has been opened SUCCESSFULLY." << G4endl;

		//读取能谱,存入energy和spectrum中
		unsigned int i = 0;
		G4double temp, sum = 0;

		vector<G4double> spectrum;

		while (spcFile >> temp)
		{
			energy.push_back(temp);
			spcFile >> temp;
			spectrum.push_back(temp);
			i++;
		}

		spcFile.close();

		//归一化能谱
		for (i = 0; i < energy.size(); ++i)
			sum += spectrum[i]; //求能谱总和
		for (i = 0; i < energy.size(); ++i)
			normSpectrum.push_back(spectrum[i] / sum);  //归一化谱数据 
	}
}
