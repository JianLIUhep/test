/**
 * @file
 * @brief Definition of [AnalysisTiming] module
 * @copyright Copyright (c) 2018 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * *
 * Refer to the User's Manual for more details.
 */

#include <iostream>
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "core/module/Module.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     * @brief Module to do function
     *
     * More detailed explanation of module
     */
    class AnalysisTiming : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detector Pointer to the detector for this module instance
         */
        AnalysisTiming(Configuration config, std::shared_ptr<Detector> detector);
        void initialise();
        StatusCode run(std::shared_ptr<Clipboard> clipboard);
        void finalise();

    private:
        std::shared_ptr<Detector> m_detector;

        // timing correction functions:
        void correctClusterTimestamp(Cluster*, int mode);

        // 1D histograms:
        TH1F* hTrackCorrelationTime;
        TH1F* hTrackCorrelationTimeAssoc;
        TH1F* hTrackCorrelationTime_rowCorr;
        TH1F* hTrackCorrelationTime_rowAndTimeWalkCorr;
        TH1F* hTrackCorrelationTime_rowAndTimeWalkCorr_l25;
        TH1F* hTrackCorrelationTime_rowAndTimeWalkCorr_l40;
        TH1F* hTrackCorrelationTime_rowAndTimeWalkCorr_g40;
        TH1D* hTrackCorrelationTime_example;

        // 2D histograms:
        TH2F* hTrackCorrelationTimeVsCol; // control plot only
        TH2F* hTrackCorrelationTimeVsRow;
        TH2F* hTrackCorrelationTimeVsRow_1px;
        TH2F* hTrackCorrelationTimeVsRow_npx;
        TH2F* hTrackCorrelationTimeVsRow_rowCorr;
        TH2F* hTrackCorrelationTimeVsTot;
        TH2F* hTrackCorrelationTimeVsTot_1px;
        TH2F* hTrackCorrelationTimeVsTot_npx;
        TH1F* hClusterTimeMinusPixelTime;
        TH2F* hTrackCorrelationTimeVsTot_rowCorr;
        TH2F* hTrackCorrelationTimeVsTot_rowCorr_1px;
        TH2F* hTrackCorrelationTimeVsTot_rowCorr_npx;
        TH2F* hTrackCorrelationTimeVsRow_rowAndTimeWalkCorr;
        TH2F* hTrackCorrelationTimeVsTot_rowAndTimeWalkCorr;

        TH2F* hClusterSizeVsTot_Assoc;

        TH2F* hHitMapAssoc;
        TH2F* hHitMapAssoc_highTot;
        TH2F* hHitMapAssoc_inPixel;
        TH2F* hHitMapAssoc_inPixel_highTot;
        TH2F* hClusterMapAssoc;

        TH2F* hTotVsTime_low;
        TH2F* hTotVsTime_high;

        // Control Plots for "left tail":
        TH2F* hClusterMap_leftTail;
        TH1F* hTot_leftTail;
        TH1F* hPixelTS1_leftTail;
        TH1F* hPixelTS2_leftTail;
        TH1F* hPixelTS1bits_leftTail;
        TH1F* hPixelTS2bits_leftTail;
        TH1F* hClusterSize_leftTail;

        TH2F* hClusterMap_rightTail;
        TH1F* hTot_rightTail;
        TH1F* hPixelTS1_rightTail;
        TH1F* hPixelTS2_rightTail;
        TH1F* hPixelTS1bits_rightTail;
        TH1F* hPixelTS2bits_rightTail;
        TH1F* hClusterSize_rightTail;

        // TGraphErrors:
        TGraphErrors* gTimeCorrelationVsRow;
        TGraphErrors* gTimeCorrelationVsTot_rowCorr;
        TGraphErrors* gTimeCorrelationVsTot_rowCorr_1px;
        TGraphErrors* gTimeCorrelationVsTot_rowCorr_npx;

        TGraphErrors* gRowCorr;
        TGraphErrors* gTimeWalkCorr;

        // Member Variables:
        std::string m_DUT;
        double m_timingCut;
        double m_chi2ndofCut;
        double m_timeCutFrameEdge;
        double m_clusterTotCut;
        size_t m_clusterSizeCut;

        std::string m_correctionFile_row;
        std::string m_correctionFile_timewalk;
        bool m_calcCorrections;
        bool m_pointwise_correction_row;
        bool m_pointwise_correction_timewalk;
        int m_totBinExample;

        int total_tracks_uncut;
        int tracks_afterChi2Cut;
        int tracks_hasIntercept;
        int tracks_isWithinROI;
        int tracks_afterMasking;
        int total_tracks;
        int matched_tracks;
        int tracks_afterClusterTotCut;
        int tracks_afterClusterSizeCut;
    };

} // namespace corryvreckan
