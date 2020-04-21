/**
 * @file
 * @brief Definition of module AnalysisMaterialBudget
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile2D.h>
#include <iostream>
#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

namespace corryvreckan {
    /** @ingroup Modules
     * @brief Module to do function
     *
     * More detailed explanation of module
     */
    class AnalysisMaterialBudget : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detectors Vector of pointers to the detectors
         */
        AnalysisMaterialBudget(Configuration config, std::vector<std::shared_ptr<Detector>> detectors);

        /**
         * @brief [Initialise this module]
         */
        void initialise();

        /**
         * @brief [Run the function of this module]
         */
        StatusCode run(std::shared_ptr<Clipboard> clipboard);

        /**
         * @brief [Finalise module]
         */
        void finalise();

    private:
        int m_eventNumber;

        ROOT::Math::XYVector cell_size;
        ROOT::Math::XYVector image_size;

        double angle_cut;
        int min_cell_content;

        bool write_3D;

        int n_cells_x, n_cells_y;

        TH1F* entriesPerCell;

        TH1F* sixletKinkX;
        TH1F* sixletKinkY;
        TProfile2D* MBIpreview;
        TH2F* MBI;
        TH2F* meanAngles;
        TH3F* MBI3D;

        std::map<std::pair<uint8_t, uint8_t>, std::vector<double>> m_all_kinks;
        std::map<std::pair<uint8_t, uint8_t>, double> m_all_sum;
        std::map<std::pair<uint8_t, uint8_t>, int> m_all_entries;

        /**
         * @brief Method re-calculating the AAD of a given image cell
         * @param Cell ID in x
         * @param Cell ID in y
         */
        double getAAD(int cell_x, int cell_y);
    };

} // namespace corryvreckan
