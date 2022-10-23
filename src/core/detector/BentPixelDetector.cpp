/** @file
 *  @brief Implementation of the detector model
 *  @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <fstream>
#include <map>
#include <string>

#include "Math/DisplacementVector2D.h"
#include "Math/PositionVector3D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/RotationZYX.h"

#include "BentPixelDetector.hpp"
#include "core/utils/log.h"
#include "exceptions.h"

using namespace ROOT::Math;
using namespace corryvreckan;

BentPixelDetector::BentPixelDetector(const Configuration& config) : PixelDetector(config) {
    build_axes(config);
}

void BentPixelDetector::build_axes(const Configuration& config) {

    m_flat_part = config.get<double>("flat_part", 0); // A bent sensor can be attached to the PCB and have a local flat part
    m_radius = config.get<double>("radius", 0);       // Radius of the bent sensor
    m_bent_axis = config.get<BentAxis>("bent_axis");  // Axis along which the detector is bent

    LOG(DEBUG) << "  Using bent coordinate system for detector '" << m_detectorName << "'.";
    if(m_bent_axis == BentAxis::COLUMN) {
        LOG(TRACE) << "\": COLUMN is bent axis, ";
    } else {
        LOG(TRACE) << "\": ROW is bent axis, ";
    }
    if(m_flat_part > 0 && m_bent_axis == BentAxis::COLUMN) {
        throw InvalidCombinationError(config, {"m_flat_part", "m_bent_axis"}, "Bending along column axis not implemented.");
    }
    if(m_radius == 0) {
        throw InvalidValueError(config, "m_radius", "A radius bigger than 0 must be given.");
    }
    LOG(TRACE) << "Bending radius of " << Units::display(m_radius, {{"mm"}}) << ", flat part of "
               << Units::display(m_flat_part, {{"mm"}});
}

void BentPixelDetector::configure_pos_and_orientation(Configuration& config) const {
    PixelDetector::configure_pos_and_orientation(config);
    config.set("radius", m_radius, {{"mm"}});
    config.set("flat_part", m_flat_part, {{"mm"}});
    config.set("bent_axis", m_bent_axis);
}

XYZPoint BentPixelDetector::localToGlobal(XYZPoint local) const {
    ROOT::Math::XYVector detector_dimensions = getSize();
    // Axis of the sensor that is bent
    if(m_bent_axis == BentAxis::COLUMN) {
        // origin of coordinate system in the middle of the flat chip (tangential point)
        double col_arc_length = local.X();
        ROOT::Math::XYZPoint local_transformed(
            m_radius *
                sin(col_arc_length / m_radius), // >0 for col_arc_length>0 and vice versa independent of sign(m_radius)
            local.y(),
            local.z() - m_radius * (1.0 - cos(col_arc_length / m_radius))); // local.z() + ... for r<0
        local_transformed = m_localToGlobal * local_transformed;
        LOG(TRACE) << "Transformed local point (" << local.x() << "|" << local.y() << "|" << local.z()
                   << ") to global point (" << local_transformed.x() << "|" << local_transformed.y() << "|"
                   << local_transformed.z() << "). Bent column";
        return local_transformed;
    } else {
        // origin of coordinate system in the middle of the flat chip
        double row_arc_length = detector_dimensions.Y() / 2 - local.Y(); // distance from top of the chip
        if(row_arc_length < m_flat_part) {                               // flat sensor
            ROOT::Math::XYZPoint local_transformed(local.x(), -row_arc_length + detector_dimensions.Y() / 2, local.z());
            local_transformed = m_localToGlobal * local_transformed;
            LOG(TRACE) << "Transformed local point (" << local.x() << "|" << local.y() << "|" << local.z()
                       << ") to global point (" << local_transformed.x() << "|" << local_transformed.y() << "|"
                       << local_transformed.z() << "). Bent row";
            return local_transformed;
        } else { // bent sensor
            ROOT::Math::XYZPoint local_transformed(
                local.x(),
                -m_radius * sin((row_arc_length - m_flat_part) / m_radius) - m_flat_part +
                    detector_dimensions.Y() /
                        2, // first term with the sin reflects the change in dimensions due to bending; together with the
                           // second term they describe the distance downward from the sensor edge; the last term translates
                           // local bent coordinates to local Corryvreckan coordinates
                m_radius * (cos((row_arc_length - m_flat_part) / m_radius) - 1.0) + local.z());
            local_transformed = m_localToGlobal * local_transformed;
            LOG(TRACE) << "Transformed local point (" << local.x() << "|" << local.y() << "|" << local.z()
                       << ") to global point (" << local_transformed.x() << "|" << local_transformed.y() << "|"
                       << local_transformed.z() << "). Bent row";
            return local_transformed;
        }
    }
}

XYZPoint BentPixelDetector::globalToLocal(XYZPoint global) const {
    ROOT::Math::XYVector detector_dimensions = getSize();
    // transform from global to bent local
    ROOT::Math::XYZPoint local_transformed = m_globalToLocal * global;
    // which axis of the sensor is bent
    if(m_bent_axis == BentAxis::COLUMN) {
        // origin of coordinate system in the middle of the chip (tangential point)
        double col_arc_length =
            m_radius * asin(local_transformed.x() / m_radius); // inverse trigonometric functions are needed to obtain an
                                                               // angle from the angle's trigonometric ratios
        ROOT::Math::XYZPoint local(col_arc_length,
                                   local_transformed.y(),
                                   m_radius * (cos(col_arc_length / m_radius) - 1)); // z is always 0 in local coordinates
        LOG(TRACE) << "Transformed global point (" << global.x() << "|" << global.y() << "|" << global.z()
                   << ") to local point (" << local.x() << "|" << local.y() << "|" << local.z() << "). Bent column";
        return local;
    } else {
        double row_arc_length =
            m_radius * asin((detector_dimensions.Y() / 2 - m_flat_part - local_transformed.Y()) / m_radius) +
            m_flat_part; // inverse of the def in localToGlobal
        // origin of coordinate in the middle of the flat chip
        if(row_arc_length < m_flat_part) { // flat part
            ROOT::Math::XYZPoint local(local_transformed.x(), local_transformed.y(), 0);
            LOG(TRACE) << "Transformed global point (" << global.x() << "|" << global.y() << "|" << global.z()
                       << ") to local point (" << local.x() << "|" << local.y() << "|" << local.z() << "). Bent row";
            return local;
        } else {                                              // bent part
            ROOT::Math::XYZPoint local(local_transformed.x(), // not modified
                                       detector_dimensions.Y() / 2 - row_arc_length,
                                       0);
            LOG(TRACE) << "Transformed global point (" << global.x() << "|" << global.y() << "|" << global.z()
                       << ") to local point (" << local.x() << "|" << local.y() << "|" << local.z() << "). Bent row";
            return local;
        }
    }
}

PositionVector3D<Cartesian3D<double>> BentPixelDetector::getIntercept(const Track* track) const {
    // FIXME: only works for straight line tracks
    if(track->getType() == "GblTrack") {
        return track->getState(getName());
    } else {
        // Get the distance from the plane to the track initial state
        double distance = (m_origin.X() - track->getState(m_detectorName).X()) * m_normal.X();
        distance += (m_origin.Y() - track->getState(m_detectorName).Y()) * m_normal.Y();
        distance += (m_origin.Z() - track->getState(m_detectorName).Z()) * m_normal.Z();
        distance /= (track->getDirection(m_detectorName).X() * m_normal.X() +
                     track->getDirection(m_detectorName).Y() * m_normal.Y() +
                     track->getDirection(m_detectorName).Z() * m_normal.Z());

        // Propagate the track
        PositionVector3D<Cartesian3D<double>> globalPlanarIntercept(
            track->getState(m_detectorName).X() + distance * track->getDirection(m_detectorName).X(),
            track->getState(m_detectorName).Y() + distance * track->getDirection(m_detectorName).Y(),
            track->getState(m_detectorName).Z() + distance * track->getDirection(m_detectorName).Z());

        // Get and transform track state and direction
        PositionVector3D<Cartesian3D<double>> state_track = track->getState(m_detectorName);
        DisplacementVector3D<Cartesian3D<double>> direction_track = track->getDirection(m_detectorName);

        // Bring track to local (transformed) coordinate system
        state_track = m_globalToLocal * state_track;
        direction_track = m_globalToLocal.Rotation() * direction_track;
        direction_track = direction_track.Unit();

        // From globalPlanarIntercept get intercept with bent surface of pixel detector
        InterceptParameters intercept_parameters;
        if(m_bent_axis == BentAxis::COLUMN) {
            // Define/initialise detector cylinder in local_transformed coordinates
            ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> state_cylinder(0, 0, -m_radius);
            ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> direction_cylinder(
                0, 1, 0); // cylinder axis along y (row)

            get_intercept_parameters(state_track, direction_track, state_cylinder, direction_cylinder, intercept_parameters);

        } else {
            double row_arc_length =
                m_radius * asin((this->getSize().Y() / 2 - m_flat_part - globalPlanarIntercept.Y()) / m_radius) +
                m_flat_part;
            if(row_arc_length < m_flat_part) { // flat part
                return globalPlanarIntercept;
            } else { // bent part
                ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> state_cylinder(
                    0, this->getSize().Y() / 2 - m_flat_part, -m_radius);
                ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> direction_cylinder(
                    1, 0, 0); // cylinder axis along x (column)

                get_intercept_parameters(
                    state_track, direction_track, state_cylinder, direction_cylinder, intercept_parameters);
            }
        }

        // Select solution according to bending direction
        PositionVector3D<Cartesian3D<double>> localBentIntercept;
        double inf = std::numeric_limits<double>::infinity();
        if(m_radius < 0) { // for negative radius select smaller solution
            if(!intercept_parameters.isValid()) {
                LOG(INFO) << "No intercept of track and cylinder sensor found. Return intercept outside of sensor area.";
                localBentIntercept.SetCoordinates(inf, inf, inf);
                return localBentIntercept;
            } else {
                localBentIntercept =
                    ROOT::Math::XYZPoint(state_track.x() + intercept_parameters.getParam1() * direction_track.x(),
                                         state_track.y() + intercept_parameters.getParam1() * direction_track.y(),
                                         state_track.z() + intercept_parameters.getParam1() * direction_track.z());
            }
        } else {
            if(!intercept_parameters.isValid()) {
                LOG(INFO) << "No intercept of track and cylinder sensor found. Return intercept outside of sensor area.";
                localBentIntercept.SetCoordinates(inf, inf, inf);
                return localBentIntercept;
            } else {
                localBentIntercept =
                    ROOT::Math::XYZPoint(state_track.x() + intercept_parameters.getParam2() * direction_track.x(),
                                         state_track.y() + intercept_parameters.getParam2() * direction_track.y(),
                                         state_track.z() + intercept_parameters.getParam2() * direction_track.z());
            }
        }

        return m_localToGlobal * localBentIntercept;
    }
}

void BentPixelDetector::get_intercept_parameters(const PositionVector3D<Cartesian3D<double>>& state_track,
                                                 const DisplacementVector3D<Cartesian3D<double>>& direction_track,
                                                 const PositionVector3D<Cartesian3D<double>>& state_cylinder,
                                                 const DisplacementVector3D<Cartesian3D<double>>& direction_cylinder,
                                                 InterceptParameters& result) const {
    // Get factors for quadratic equation: alpha z^2 + beta * z + gamma
    double alpha;
    double beta;
    double gamma;
    if(direction_cylinder.X() != 0 && direction_cylinder.Y() == 0 &&
       direction_cylinder.Z() == 0) { // cylinder along x (column)
        alpha = (direction_track.Y() * direction_track.Y()) + (direction_track.Z() * direction_track.Z());
        beta = 2 * (state_track.Y() * direction_track.Y() + state_track.Z() * direction_track.Z() -
                    direction_track.Y() * state_cylinder.Y() - direction_track.Z() * state_cylinder.Z());
        gamma =
            ((state_track.Y() * state_track.Y()) + (state_track.Z() * state_track.Z()) +
             (state_cylinder.Y() * state_cylinder.Y()) + (state_cylinder.Z() * state_cylinder.Z()) -
             2 * state_track.Y() * state_cylinder.Y() - 2 * state_track.Z() * state_cylinder.Z() - (m_radius * m_radius));
    } else if(direction_cylinder.X() == 0 && direction_cylinder.Y() != 0 &&
              direction_cylinder.Z() == 0) { // cylinder along y (row)
        alpha = (direction_track.X() * direction_track.X()) + (direction_track.Z() * direction_track.Z());
        beta = 2 * (state_track.X() * direction_track.X() + state_track.Z() * direction_track.Z() -
                    direction_track.X() * state_cylinder.X() - direction_track.Z() * state_cylinder.Z());
        gamma =
            ((state_track.X() * state_track.X()) + (state_track.Z() * state_track.Z()) +
             (state_cylinder.X() * state_cylinder.X()) + (state_cylinder.Z() * state_cylinder.Z()) -
             2 * state_track.X() * state_cylinder.X() - 2 * state_track.Z() * state_cylinder.Z() - (m_radius * m_radius));
    } else {
        throw InvalidSettingError(this,
                                  "direction_cylinder",
                                  "Unknown cylinder direction! Cylinder can only be bent along 'x' (column) or 'y' (row).");
    }

    // Check if the quadratic equation has a solution
    double discriminant = beta * beta - 4 * alpha * gamma;
    if(discriminant < 0 || alpha == 0) { // must have real solution
        return;
    }

    // Solve equation
    result.setParams((-beta + sqrt(discriminant)) / (2 * alpha), (-beta - sqrt(discriminant)) / (2 * alpha));
    return;
}

bool BentPixelDetector::hasIntercept(const Track* track, double pixelTolerance) const {

    // First, get the track intercept in global coordinates with the plane
    PositionVector3D<Cartesian3D<double>> globalIntercept = this->getIntercept(track);

    // Convert to local coordinates
    PositionVector3D<Cartesian3D<double>> localIntercept = this->m_globalToLocal * globalIntercept;

    // Get the row and column numbers
    double row = this->getRow(localIntercept);
    double column = this->getColumn(localIntercept);

    // Check if the row and column are outside of the chip
    // Chip reaches from -0.5 to nPixels-0.5
    bool intercept = true;
    if(row < pixelTolerance - 0.5 || row > (this->m_nPixels.Y() - pixelTolerance - 0.5) || column < pixelTolerance - 0.5 ||
       column > (this->m_nPixels.X() - pixelTolerance - 0.5)) {
        intercept = false;
    }

    // Check if the cylinder surface - track intercept parameters are valid
    std::vector<double> coords(3);
    globalIntercept.GetCoordinates(coords.begin(), coords.end());
    for(auto coord : coords) {
        if(coord == std::numeric_limits<double>::infinity()) {
            LOG(TRACE) << "Invalid cylinder surface - track intercept parameters!";
            intercept = false;
        }
    }

    return intercept;
}
