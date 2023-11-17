#ifndef UserDataInput_h
#define UserDataInput_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

class UserDataInput
{
public:
	UserDataInput();
	~UserDataInput();

	static void ReadInputData(); //读取外部输入文件

	static inline G4int GetNumberOfThreads() { return numberOfThreads; }
	static inline G4int GetNumberOfEvents() { return numberofevents; }

	static inline G4bool GetStatusOfUI() { return uiStatus; }
	static inline G4bool GetModeratorStatus() { return moderatorStatus; }

	static inline G4String GetSourceType() { return source; }
	static inline G4String GetDetectorMaterial() { return detectorMaterial; }
	static inline G4String GetModeratorPosition() { return moderatorPosition; }

	static inline G4double GetDistanceOfSD() { return sourceDetectorDistance * cm; }
	static inline G4double GetDectorDimensionX() { return dtctrx * mm; }
	static inline G4double GetDectorDimensionY() { return dtctry * mm; }
	static inline G4double GetDectorDimensionZ() { return dtctrz * um; }
	static inline G4double GetGammaPercentage() { return gammaPercentage; }
	static inline G4double GetModeratorThickness() { return moderatorThickness * cm; }

	static inline vector<G4double> GetGammaEnergy() { return gammaEnergy; }
	static inline vector<G4double> GetGammaNormSpectrum() { return gammaNormSpectrum; }
	static inline vector<G4double> GetNeutronEnergy() { return neutronEnergy; }
	static inline vector<G4double> GetNeutronNormSpectrum() { return neutronNormSpectrum; }

private:

	static G4int numberOfThreads;
	static G4int numberofevents; //要模拟的离子数

	static G4bool uiStatus;
	static G4bool moderatorStatus; //中子慢化是否打开

	static G4String source; //源的类型
	static G4String detectorMaterial;
	static G4String moderatorPosition;

	static G4double sourceDetectorDistance; //源到探测器距离
	static G4double dtctrx, dtctry, dtctrz; //探测器的长宽高
	//static G4double StartThickness, ThicknessIncrement, EndThickness;
	static G4double gammaPercentage; //总谱中gamma占比
	static G4double moderatorThickness; //慢化剂厚度

	static vector<G4double> gammaEnergy; // gamma能量
	static vector<G4double> gammaNormSpectrum; //gamma归一化能谱
	static vector<G4double> neutronEnergy; //中子能量
	static vector<G4double> neutronNormSpectrum; //中子归一化能谱

	static void ReadSpectra(G4String spectrumName, vector<G4double>& energy, vector<G4double>& normSpectrum); //读取能谱文件
};

#endif // !UserDataInput_h
