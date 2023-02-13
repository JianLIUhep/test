/** @file
 *  @brief Detector model class
 *  @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#ifndef CORRYVRECKAN_BENTPIXELDETECTOR_H
#define CORRYVRECKAN_BENTPIXELDETECTOR_H

#include <fstream>
#include <map>
#include <string>

#include "Math/DisplacementVector2D.h"
#include "Math/Point3Dfwd.h"
#include "Math/Transform3D.h"
#include "Math/Vector2D.h"
#include "Math/Vector3D.h"

#include "Detector.hpp"
#include "core/config/Configuration.hpp"
#include "core/utils/ROOT.h"
#include "core/utils/log.h"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

namespace corryvreckan {

    /**
     * @brief PixelDetector representation derived from Detector interface in the reconstruction chain
     *
     * Contains the PixelDetector with all its properties such as position and orientation, pitch, spatial resolution
     * etc.
     */
    class BentPixelDetector : public PixelDetector {
    public:
        /**
         * Delete default constructor
         */
        BentPixelDetector() = delete;

        /**
         * Default destructor
         */
        ~BentPixelDetector() = default;

        /**
         * @brief Constructs a detector in the geometry
         * @param config Configuration object describing the detector
         */
        BentPixelDetector(const Configuration& config);

        // Function to get global intercept with a track
        PositionVector3D<Cartesian3D<double>> getIntercept(const Track* track) const override;

        // Function to check if a track intercepts with a plane
        bool hasIntercept(const Track* track, double pixelTolerance = 0.) const override;

        /**
         * @brief Transform local coordinates of this detector into global coordinates
         * @param  local Local coordinates in the reference frame of this detector
         * @return       Global coordinates
         */
        XYZPoint localToGlobal(XYZPoint local) const override;

        /**
         * @brief Transform global coordinates into detector-local coordinates
         * @param  global Global coordinates
         * @return        Local coordinates in the reference frame of this detector
         */
        XYZPoint globalToLocal(XYZPoint global) const override;

        PositionVector3D<Cartesian3D<double>> getLocalPosition(double column, double row) const override;
    private:
        /**
         * @brief Different bent axis types
         */
        enum class BentAxis {
            COLUMN,
            ROW,
        } m_bent_axis;

        struct InterceptParameters {
        public:
            void setParams(double param1, double param2) {
                if(std::isnan(param1) || std::isnan(param2)) {
                    resetParams();
                    return;
                }
                if(param1 < param2) {
                    intercept_param1 = param1;
                    intercept_param2 = param2;
                } else if(param1 > param2) {
                    intercept_param1 = param2;
                    intercept_param2 = param1;
                } else {
                    intercept_param1 = param1;
                    intercept_param2 = param1;
                }
                valid_params = true;
            }
            void resetParams() {
                intercept_param1 = 0.;
                intercept_param2 = 0.;
                valid_params = false;
            }
            double getParam1() { return intercept_param1; }
            double getParam2() { return intercept_param2; }
            bool isValid() { return valid_params; }

        private:
            double intercept_param1 = 0.;
            double intercept_param2 = 0.;
            bool valid_params = false;
        };

        // Build axis, for devices which are not auxiliary
        // Different in Pixel/Strip Detector
        void build_axes(const Configuration& config) override;

        // Config position, orientation, mode of detector
        // Different in Pixel/Strip Detector
        void configure_pos_and_orientation(Configuration& config) const override;

        // Function to facilitate the incept calculation of sraight tracks with cylindrical sensor surface
        void
        get_intercept_parameters(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>>& state_track,
                                 const ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>>& direction_track,
                                 const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>>& state_cylinder,
                                 const ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>>& direction_cylinder,
                                 InterceptParameters& result) const;

        // bent geometry configuration parameters
        double m_radius;
        double m_rotate_by;
    };
} // namespace corryvreckan

#endif // CORRYVRECKAN_BENTPIXELDETECTOR_H
