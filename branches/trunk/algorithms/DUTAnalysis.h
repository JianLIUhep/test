#ifndef DUTAnalysis_H
#define DUTAnalysis_H 1

#include "Algorithm.h"
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

class DUTAnalysis : public Algorithm {
  
public:
  // Constructors and destructors
  DUTAnalysis(bool);
  ~DUTAnalysis(){}

  // Functions
  void initialise(Parameters*);
  StatusCode run(Clipboard*);
  void finalise();
  
  // Histograms
  TH1F* tracksVersusTime;
  TH1F* associatedTracksVersusTime;
  TH1F* residualsX;
  TH1F* residualsY;
  TH1F* hTrackCorrelationX;
  TH1F* hTrackCorrelationY;
  TH1F* hTrackCorrelationTime;
  TH1F* residualsTime;
  TH2F* clusterToTVersusTime;
  TH2F* residualsTimeVsTime;

  TH2F* hAssociatedTracksGlobalPosition;
  TH2F* hUnassociatedTracksGlobalPosition;

  TH1F* tracksVersusPowerOnTime;
  TH1F* associatedTracksVersusPowerOnTime;
  
  // Member variables
  int m_eventNumber;
  int m_nAlignmentClusters;
  long long int m_powerOnTime;
  long long int m_powerOffTime;
  long long int m_shutterOpenTime;
  long long int m_shutterCloseTime;
  bool m_digitalPowerPulsing;

};

#endif // DUTAnalysis_H
