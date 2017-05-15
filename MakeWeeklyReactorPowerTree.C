#include <iostream>
#include <sstream>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>
#include "time.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TString.h"
#include "TDatime.h"
const int SIZE = 255, nReactor = 6, EPOCHtoROOT = 788918400;
int MakeWeeklyReactorPowerTree(const char *filename="ReactorPowerWeeklyAvg_P15A.dat")
{
  struct tm *startdate, *enddate;
  ifstream file(filename);
  if(!(file.good() && file.is_open())){
    cout<<"File not found. Exiting.\n";
    return -1;
  }
  int week, reactor;
  int syear, eyear, smonth, emonth, sdom, edom, sdow, edow;
  int s_epoch, e_epoch;
  double power[nReactor], u235[nReactor], u238[nReactor];
  double pu239[nReactor], pu241[nReactor];
  string line;
  TString name[nReactor] = {"DB1","DB2","LA1","LA2","LA3","LA4"};
  const int SIZE = 255;
  TFile *newfile = new TFile(Form("WeeklyReactorData.root"), "recreate");
  TTree *tr = new TTree("tree","tree");
  TDatime stime, etime; 
  TBranch *brSEpoch = tr->Branch("start_epoch_time", &s_epoch, 
			     "start_epoch_time/I");
  TBranch *brEEpoch = tr->Branch("end_epoch_time", &e_epoch, 
			     "end_epoch_time/I");
  TBranch *brSY = tr->Branch("start_year", &syear, "start_year/I");
  TBranch *brEY = tr->Branch("end_year", &eyear, "end_year/I");
  TBranch *brSM = tr->Branch("start_month", &smonth);
  TBranch *brEM = tr->Branch("end_month", &emonth);
  TBranch *brSDom = tr->Branch("start_day_of_month", &sdom, 
			       "start_day_of_month/I");
  TBranch *brEDom = tr->Branch("end_day_of_month", &edom,
			       "end_day_of_month/I");
  TBranch *brSDow = tr->Branch("start_day_of_week", &sdow,
			       "start_day_of_week/I");
  TBranch *brEDow = tr->Branch("end_day_of_week", &edow,
			       "end_day_of_week/I");
  TBranch *brWeek = tr->Branch("weeks_since_Dec_24_2011", &week,
			       "weeks_since_Dec_24_2011/I");

  TBranch *brPow[nReactor];
  TBranch *brU235[nReactor];
  TBranch *brU238[nReactor];
  TBranch *brPu239[nReactor];
  TBranch *brPu241[nReactor];

  for(int i=0;i<nReactor;++i){
    brPow[i] = tr->Branch(Form("%s_Power",name[i].Data()), &power[i],
				Form("%s_Power/D",name[i].Data()));
    brU235[i] = tr->Branch(Form("fracU235_%s",name[i].Data()), &u235[i],
				Form("fracU235_%s/D",name[i].Data()));
    brU238[i] = tr->Branch(Form("fracU238_%s",name[i].Data()), &u238[i],
				Form("fracU238_%s/D",name[i].Data()));
    brPu239[i] = tr->Branch(Form("fracPu239_%s",name[i].Data()), &pu239[i],
				Form("fracPu239_%s/D",name[i].Data()));
    brPu241[i] = tr->Branch(Form("fracPu241_%s",name[i].Data()), &pu241[i],
				Form("fracPu241_%s/D",name[i].Data()));
  }
  while(!file.eof()){
    string temp;
    char a[SIZE];
    for(int i=0;i<nReactor;++i){
      getline(file, line);
      stringstream ss(line);
      ss.getline(a,SIZE, ' ');
      week = atoi(a);
      ss.getline(a,SIZE, ' ');
      reactor = atoi(a);
      ss.getline(a,SIZE, ' ');
      ss.getline(a,SIZE, ' ');

      ss.getline(a,SIZE, ' ');
      s_epoch = atoi(a);
      time_t st(s_epoch);
      startdate = localtime(&st);
      syear = startdate->tm_year + 1900;
      smonth = startdate->tm_mon + 1;//months counted as Jan=0, Feb=1...
      sdom = startdate->tm_mday;
      sdow = startdate->tm_wday;
      ss.getline(a,SIZE, ' ');
      e_epoch = atoi(a);
      time_t et(e_epoch);
      enddate = localtime(&et);
      eyear = enddate->tm_year + 1900;
      emonth = enddate->tm_mon;
      edom = enddate->tm_mday;
      edow = enddate->tm_wday;

      ss.getline(a,SIZE, ' ');
      power[i] = atof(a);
      ss.getline(a,SIZE, ' ');
      u235[i] = atof(a);
      ss.getline(a,SIZE, ' ');
      u238[i] = atof(a);
      ss.getline(a,SIZE, ' ');
      pu239[i] = atof(a);
      ss.getline(a,SIZE, ' ');
      pu241[i] = atof(a);
      cout<<week<<" "<<syear<<" "<<smonth<<" "<<sdom<<" "<<sdow<<" ";
      cout<<eyear<<" "<<emonth<<" "<<edom<<" "<<edow<<" ";
      cout<<power[i]<<" "<<u235[i]<<" ";
      cout<<u238[i]<<" "<<pu239[i]<<" "<<pu241[i]<<endl;
    }
    tr->Fill();

    file.peek();
  }
  tr->Write("",TObject::kOverwrite);
  newfile->Close();
  delete newfile;
  return 0;
}
