#include <vector>
#include <map>
#include "TGraph.h"
#include "TROOT.h"
#include "TTree.h"
#include "TF1.h"
#include "TRandom.h"
#include "TH1D.h"

class SNF
{
 public:
  SNF();
  SNF(double);
  ~SNF();
  static const int kError = -999999;
  static const int kDetector = 4;
  static const int kReactor = 6;
  static const int kDaysInWeek = 7;
  static double EbinWidth(){return 0.25;}//width of energy bin in MeV for referece SNF spectrum
  static double ReactorOnPower(){return 0.9;}//reactor fully on above this
  static double MinOnPower(){return 0.01;}//reactor fully off below this
  static double SinSq2Theta12(){return 0.846;}
  static double SinSq2Theta13(){return 0.0841;}
  static double DeltaM21sq(){return 7.53e-5;}
  static double DeltaM31sq(){return 2.44e-3;}
  static double DeltaMeeSq(){return 2.5e-3;}
  static double IBDthresh(){return 1.80;}//threshold for IBD interaction in MeV
  vector<double>vReactorPower;
  vector<double>vRefuelFreq;
  vector<double>vEnergy;
  vector<double>vDB_Spectrum;//measured Daya Bay cross section spectrum
  vector<double>vDB_EbinCenter;//energy bins for measured cross section
  vector<double>vDB_EbinWidth;//energy bin width for measured cross section
  vector<double>vTime;
  vector<vector<double> >vAvgProbSurvival;//[det][core]
  vector<vector<double> >vRefuelTimes;//[core][n] week when refuel starts 
  vector<vector<double> >vSpectrum;//[time][energy]
  vector<vector<double> >vParameter;//[par_num][energy]
  vector<vector <double> >vDetPos;//[det][axis] axis=0: x; axis=1: y; axis=2: z
  vector<vector <double> >vCorePos;//[det][axis] axis=0: x; axis=1: y; axis=2: z
  vector<vector <double> >vSNFPos;//[det][axis] axis=0: x; axis=1: y; axis=2: z
  vector<TF1*>vSNFofT;
  TGraph *gr;
  TH1D *hSpec;
  double AbsDayaBaySpectrum(double);//measured Daya Bay cross section spectrum
  double AbsDetectorResponseSpectrum(int, double, double, int week);//absolute SNF flux between lE and hE (in MeV) and for week since beginning of P15A dataset at specified detector location convoluted with measured Daya Bay spectrum. 
  int AddRandomError();//add Gaussian error to the various elements of spectrum to determine spread for error calculation
  double Alpha(int);//scale factor accounts for SNF with burnup != 45GWd/tU
  double BaselineCore(int, int);//distance from detector to core
  double BaselineSNF(int, int);//distance from detector to SNF pool near core
  double Beta(int, double, double);//fractional contribution to neutrino flux spectrum from a given reactor as a function of energy (in GeV) and time (in days) 
  double Delta(int, int);//(core, week)scale factor to account for reactor off SNF
  int Configure(const char *);//Read parameters from config file
  void FindAvgProbSurvival(bool);//find average survival probabilities for each detector-reactor combination
  int FindRefuelTimes(int);//Find refueling times for each reactor
  TGraph* GetAbsDayaBaySpectrumGraph();//return TGraph of Daya Bay cross section spectrum. Units are 10^-42 cm^2/MeV/fission
  TH1D* GetAbsDayaBaySpectrumHisto();//return TH1D of Daya Bay cross section spectrum. Units are 10^-42 cm^2/MeV/fission
  TGraph* GetAbsSNFSpectrumGraph(int, int);//return TGraph of absolute spectrum for det at week
  TH1D* GetAbsSNFSpectrumHisto(int, int);//return TH1D of absolute spectrum for det at week
  double GetAlphaError();//Relative uncertainty in Alpha(): burnup per cycle 
  double GetAvgProbSurvival(int, int);//return average survival probabilities for each detector-reactor combination
  double GetBaselineError();//Absolute uncertainty in baseline
  double GetDBSpectrumIntegral();//return integrated DB spectrum 
  TString GetDetectorName(int);//Useful for labeling graphs
  TString GetReactorName(int);//Useful for labeling graphs
  TGraph* GetRelSNFSpectrumGraph(int, int);//return TGraph of fractional spectrum for det at week
  TH1D* GetRelSNFSpectrumHisto(int, int);//return TH1D of fractional spectrum for det at week
  double GetSpectrumError();//Relative uncertainty in calculated SNF spectrum
  void IncludeProbSurvival(bool);//Use to include survival probability in reactor weighting
  void InsertRefuelTime(int);//function to add refuel times before start of P15A
  double IntegrateDBfluxSpectrum();//total DB spectrum flux for normalization
  int IsolateSingleReactorContribution(int);//Calculate fractional contribution of single reactor 
  int LoadReactorPowerTree(const char *, const char *);//Root tree of reactor power
  int ParameterizeSNFvsT(bool);//parameterize decay versus time for each energy bin
  double ProbSurvival(double, double);//survival probability function of L in km and E in MeV
  double Ptilde(int, int);//if reactor power > MinOnPower Ptilde = reactor power;otherwise Ptilde = 1
  double RatedOutputPower(int);//reactor rated output power in GW
  double ReactorOffScaling();//Scale SNF by this when reactor  off
  double ReactorData(const char*, int);//Reactor data from tree by name of leaf for week
  double ReactorPower(int, int);//Reactor power as function of core# and week#
  double RelDetectorResponseSpectrum(int, double, double, int);//fractional SNF flux at E (in MeV) and week since beginning of P15A dataset at specified detector location convoluted with measured Daya Bay spectrum. 
  double RelDetectorResponseSpectrum(int, double, double, int, int);//fractional SNF flux at E (in MeV) and week*7+day days since beginning of P15A dataset at specified detector location convoluted with measured Daya Bay spectrum. 
  double RelSpectrumAtDet(int, double, int);//fractional(%) SNF flux at E (in MeV) and week since beginning of P15A dataset at specified detector location (not convoluted with IBD X-section)
  double RelSpectrumAtDet(int, double, int, int);//(det, E, week, day)fractional(%) SNF flux at E (in MeV) for week and day(since beginning of week 0-6)for P15A dataset at specified detector location (not convoluted with IBD X-section)
  double RelSpectrumSim(double, double);//Simulated fractional SNF spectrum of E (in MeV) and time (in days) 
  void SetAlpha(int, double);//set value returned by Alpha(): burnup per cycle for a reactor
  void SetAlphaError(double);//Uncertainty in Alpha(): burnup per cycle for a reactor  
  void SetBaselineError(double);//Absolute uncertainty in baseline
  void SetDayLevelPrecision(bool);//set whether or not to use day-level precision of refuel times. Default is week-level precision.
  void SetDeltaE(double);//set width dE of bin for sum approximation of convolution integral
  void SetReactorOffScaling(double);//Scale SNF by this when reactor  off
  void SetRelativeSpectrumError(double);//Relative uncertainty in reference SNF spectrum 
  void SetSpectrumError(double);//Relative uncertainty in reference SNF spectrum
  double Theta13();
private:
  TRandom *rand;
  TTree *RPTree;
  TString RName[kReactor];
  TString PName[kReactor];
  TString DName[kDetector];
  TString RProperName[kReactor];
  bool fAddRandomError;
  bool fAtDayLevelPrecision;
  bool fIsolateReactorContribution;
  bool fIncludeProbSurvival;
  int fSingleReactor;
  int fNEnergyBins;
  int fNTimeBins;
  int fNWeeksInP15A;
  double fRefuelTimeError;//relative error in time between refuels before P15A
  double fBurnupScaleFactor[kReactor];
  double fSpectrumError;//Relative uncertainty in reference SNF spectrum 
  double fBaselineError;//Absolute uncertainty in baseline
  double fAlphaError;//Relative uncertainty in Alpha(): burnup per cycle for a reactor
  double fDBSpectrumIntegral;
  double fDeltaE;//width dE of bin for sum approximation of convolution integral
  double fReactorOffScaling;//scale SNF fraction by this when reactor off 
};
