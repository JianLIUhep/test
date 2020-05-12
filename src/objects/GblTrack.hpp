/**
 * @file
 * @brief Definition of GBL track object
 *
 * @copyright Copyright (c) 2019-2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef CORRYVRECKAN_GBLTRACK_H
#define CORRYVRECKAN_GBLTRACK_H 1

#include "Track.hpp"
namespace corryvreckan {
    /**
     * @ingroup Objects
     * @brief GblTrack object
     *
     * This class is a general broken line track which knows how to fit itself. It is dervied from Track
     */

    class GblTrack : public Track {

    public:
        void print(std::ostream& out) const override;

        /**
         * @brief The fitting routine
         */
        void fit() override;

        /**
         * @brief Get the track position for a certain z position
         * @param z positon
         * @return ROOT::Math::XYZPoint at z position
         */
        ROOT::Math::XYZPoint getIntercept(double z) const override;

        /**
         * @brief Get the track state at a detector
         * @param name of detector
         * @return ROOT::Math::XYZPoint state at detetcor layer
         */
        ROOT::Math::XYZPoint getState(std::string detectorID) const override;

        /**
         * @brief Get the track direction at a detector
         * @param name of detector
         * @return ROOT::Math::XYZPoint direction at detetcor layer
         */
        ROOT::Math::XYZVector getDirection(std::string detectorID) const override;

        /**
         * @brief Return kink of track at given detector
         * @param  detectorID Detector ID at which the kink should be evaluated
         * @return            Kink at given detector
         */
        ROOT::Math::XYPoint getKinkAt(std::string detectorID) const override;

        void setVolumeScatter(double length) override;

    private:
        /**
         * @brief Set seedcluster used for track fitting
         * @param Pointer to seedcluster of the GblTrack
         */
        void set_seed_cluster(const Cluster* cluster);

        /**
         * @brief Get seedcluster used for track fitting
         * @return Pointer to seedcluster of the GblTrack if set, nullptr otherwise
         */
        Cluster* get_seed_cluster() const;

        // Member variables
        TRef seed_cluster_{nullptr};
        double scattering_length_volume_{};
        bool use_volume_scatter_{};

        std::map<std::string, ROOT::Math::XYPoint> local_track_points_{};
        std::map<std::string, ROOT::Math::XYPoint> initital_residual{};
        // ROOT I/O class definition - update version number when you change this class!
        ClassDefOverride(GblTrack, 3);
    };
} // namespace corryvreckan

#endif // CORRYVRECKAN_GBLTRACK_H
