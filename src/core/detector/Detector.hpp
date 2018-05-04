/** @file
 *  @brief Detector model class
 *  @copyright Copyright (c) 2017 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef CORRYVRECKAN_DETECTOR_H
#define CORRYVRECKAN_DETECTOR_H

#include <fstream>
#include <map>
#include <string>

#include <Math/DisplacementVector2D.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include "Math/Transform3D.h"
#include "Math/Vector3D.h"

#include "core/config/Configuration.hpp"
#include "core/utils/ROOT.h"
#include "core/utils/log.h"
#include "objects/Track.h"

using namespace ROOT::Math;

namespace corryvreckan {
    // Class containing just the detector parameters
    class Detector {
    public:
        // Constructors and desctructors
        Detector() = delete;
        Detector(const Configuration& config);
        ~Detector() {}

        // Functions to retrieve basic information
        std::string type() { return m_detectorType; }
        std::string name() { return m_detectorName; }

        Configuration getConfiguration();

        ROOT::Math::XYVector size() { return ROOT::Math::XYVector(m_pitch.X() * m_nPixelsX, m_pitch.Y() * m_nPixelsY); }
        ROOT::Math::XYVector pitch() { return m_pitch; }

        int nPixelsX() { return m_nPixelsX; }
        int nPixelsY() { return m_nPixelsY; }
        double timingOffset() { return m_timingOffset; }

        // Functions to set and retrieve basic translation parameters
        void displacementX(double x) { m_displacement.SetX(x); }
        void displacementY(double y) { m_displacement.SetY(y); }
        void displacementZ(double z) { m_displacement.SetZ(z); }
        void displacement(ROOT::Math::XYZPoint displacement) { m_displacement = displacement; }
        ROOT::Math::XYZPoint displacement() { return m_displacement; }

        // Functions to set and retrieve basic rotation parameters
        void rotationX(double rx) { m_orientation.SetX(rx); }
        void rotationY(double ry) { m_orientation.SetY(ry); }
        void rotationZ(double rz) { m_orientation.SetZ(rz); }
        ROOT::Math::XYZVector rotation() { return m_orientation; }
        void rotation(ROOT::Math::XYZVector rotation) { m_orientation = rotation; }

        PositionVector3D<Cartesian3D<double>> normal() { return m_normal; };

        // Functions to set and check channel masking
        void setMaskFile(std::string file);
        void processMaskFile();

        std::string maskFile() { return m_maskfile; }
        void maskChannel(int chX, int chY);
        bool masked(int chX, int chY);

        // Function to initialise transforms
        void initialise();

        // Function to update transforms (such as during alignment)
        void update();

        // Function to get global intercept with a track
        PositionVector3D<Cartesian3D<double>> getIntercept(Track* track);

        // Function to check if a track intercepts with a plane
        bool hasIntercept(Track* track, double pixelTolerance = 0.);

        // Function to check if a track goes through/near a masked pixel
        bool hitMasked(Track* track, int tolerance = 0.);

        // Functions to get row and column from local position
        double getRow(PositionVector3D<Cartesian3D<double>> localPosition);
        double getColumn(PositionVector3D<Cartesian3D<double>> localPosition);

        // Function to get local position from row and column
        PositionVector3D<Cartesian3D<double>> getLocalPosition(double row, double column);

        // Function to get in-pixel position
        double inPixelX(PositionVector3D<Cartesian3D<double>> localPosition);
        double inPixelY(PositionVector3D<Cartesian3D<double>> localPosition);

        ROOT::Math::XYZPoint localToGlobal(ROOT::Math::XYZPoint local) { return m_localToGlobal * local; };
        ROOT::Math::XYZPoint globalToLocal(ROOT::Math::XYZPoint global) { return m_globalToLocal * global; };

    private:
        // Member variables
        // Detector information
        std::string m_detectorType;
        std::string m_detectorName;
        ROOT::Math::XYVector m_pitch;
        int m_nPixelsX;
        int m_nPixelsY;
        double m_timingOffset;

        // Displacement and rotation in x,y,z
        ROOT::Math::XYZPoint m_displacement;
        ROOT::Math::XYZVector m_orientation;
        std::string m_orientation_mode;

        // Transforms from local to global and back
        Transform3D m_localToGlobal;
        Transform3D m_globalToLocal;

        // Normal to the detector surface and point on the surface
        PositionVector3D<Cartesian3D<double>> m_normal;
        PositionVector3D<Cartesian3D<double>> m_origin;

        // List of masked channels
        std::map<int, bool> m_masked;
        std::string m_maskfile;
        std::string m_maskfile_name;
    };
} // namespace corryvreckan

#endif // CORRYVRECKAN_DETECTOR_H