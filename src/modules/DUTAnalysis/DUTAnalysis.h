#ifndef DUTAnalysis_H
#define DUTAnalysis_H 1

#include <iostream>
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "core/module/Module.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class DUTAnalysis : public Module {

    public:
        // Constructors and destructors
        DUTAnalysis(Configuration config, Detector* detector);

        ~DUTAnalysis() {}

        // Functions
        void initialise();
        StatusCode run(Clipboard* clipboard);
        void finalise();

    private:
        // Histograms
        TH1F* tracksVersusTime;
        TH1F* associatedTracksVersusTime;
        TH1F* residualsX;
        TH1F* residualsXfine;
        TH1F* residualsX1pix;
        TH1F* residualsX2pix;
        TH1F* residualsY;
        TH1F* clusterTotAssociated;
        TH1F* clusterSizeAssociated;
        TH1F* clusterSizeAssociated_X;
        TH1F* clusterSizeAssociated_Y;
        TH1F* hTrackCorrelationX;
        TH1F* hTrackCorrelationY;
        TH1F* hTrackCorrelationTime;
        TH1F* residualsTime;
        TH2F* clusterToTVersusTime;
        TH2F* residualsTimeVsTime;
        TH2F* residualsTimeVsSignal;

        TH1F* residualsXMCtruth;
        TH1F* telescopeResolution;

        TH2F* hAssociatedTracksGlobalPosition;
        TH2F* hUnassociatedTracksGlobalPosition;

        TH1F* tracksVersusPowerOnTime;
        TH1F* associatedTracksVersusPowerOnTime;

        // Member variables
        int m_eventNumber;
        int m_nAlignmentClusters;
        bool m_useMCtruth;
        long long int m_powerOnTime;
        long long int m_powerOffTime;
        double m_shutterOpenTime;
        double m_shutterCloseTime;
        bool m_digitalPowerPulsing;
        double chi2ndofCut;
    };
} // namespace corryvreckan
#endif // DUTAnalysis_H
