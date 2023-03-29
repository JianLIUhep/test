/** @file
 *  @brief Implementation of the detector model
 *  @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include <fstream>
#include <map>
#include <string>

#include "BentPixelDetector.hpp"
#include "core/utils/log.h"

using namespace ROOT::Math;
using namespace corryvreckan;

BentPixelDetector::BentPixelDetector(const Configuration& config) : PixelDetector(config) {
    build_axes(config);

}

void BentPixelDetector::build_axes(const Configuration& config) {
    m_rotate_by = config.get<double>("rotate_by", 0);
    m_radius = config.get<double>("radius");         // Radius of the bent sensor
    m_bent_axis = config.get<BentAxis>("bent_axis"); // Axis along which the detector is bent

    LOG(DEBUG) << "  Using bent coordinate system for detector '" << m_detectorName << "'.";
    if(m_bent_axis == BentAxis::COLUMN) {
        LOG(TRACE) << "\": COLUMN is bent axis, ";
    } else {
        LOG(TRACE) << "\": ROW is bent axis, ";
    }
    if(m_radius == 0) {
        throw InvalidValueError(config, "m_radius", "A radius bigger than 0 must be given.");
    }
    LOG(TRACE) << "Bending radius of " << Units::display(m_radius, {{"mm"}}) << ", rotated by "
               << Units::display(m_rotate_by, {{"deg"}});
}

void BentPixelDetector::configure_pos_and_orientation(Configuration& config) const {
    PixelDetector::configure_pos_and_orientation(config);
    config.set("radius", m_radius, {{"mm"}});
    config.set("m_rotate_by", m_rotate_by, {{"deg"}});
    config.set("bent_axis", m_bent_axis);
}

// XYZPoint BentPixelDetector::localToGlobal(XYZPoint local) const {
//     double locx, locy, locz;
//     // Axis of the sensor that is bent
//     if(m_bent_axis == BentAxis::COLUMN) {
//         if(m_rotate_by != 0) {
//             local.SetX(local.x() + static_cast<double>(Units::convert(m_rotate_by, "rad")) * m_radius);
//         }
//         locx = m_radius * sin(local.x() / m_radius);
//         locy = local.y();
//         locz = local.z() - m_radius * (1.0 - cos(local.x() / m_radius));
//     } else {
//         if(m_rotate_by != 0) {
//             local.SetY(local.y() + static_cast<double>(Units::convert(m_rotate_by, "rad")) * m_radius);
//         }
//         locx = local.x();
//         locy = m_radius * sin(local.y() / m_radius);
//         locz = local.z() - m_radius * (1.0 - cos(local.y() / m_radius));
//     }
//     ROOT::Math::XYZPoint local_transformed = toGlobal() * ROOT::Math::XYZPoint(locx, locy, locz);
//     LOG(TRACE) << "Transformed local point (" << local.x() << "|" << local.y() << "|" << local.z() << ") to global point ("
//                << local_transformed.x() << "|" << local_transformed.y() << "|" << local_transformed.z() << ").";
//     return local_transformed;
// }

// XYZPoint BentPixelDetector::globalToLocal(XYZPoint global) const {
//     double lx, ly, lz;
//     // Transform from global to bent local
//     ROOT::Math::XYZPoint local_transformed = toLocal() * global;
//     // Axis of the sensor that is bent
//     if(m_bent_axis == BentAxis::COLUMN) {
//         if(m_rotate_by != 0) {
//             local_transformed.SetX(local_transformed.x() +
//                                    static_cast<double>(Units::convert(m_rotate_by, "rad")) * m_radius);
//         }
//         lx = m_radius * asin(local_transformed.x() / m_radius);
//         ly = local_transformed.y();
//         lz = m_radius * (cos(asin(local_transformed.x() / m_radius) / m_radius) -
//                          1); // TODO: should z be 0 in local coords(i.e. add loc_tr.z())?
//     } else {
//         if(m_rotate_by != 0) {
//             local_transformed.SetY(local_transformed.y() +
//                                    static_cast<double>(Units::convert(m_rotate_by, "rad")) * m_radius);
//         }
//         lx = local_transformed.x();
//         ly = m_radius * asin(local_transformed.y() / m_radius);
//         lz = m_radius * (cos(asin(local_transformed.y() / m_radius) / m_radius) -
//                          1); // TODO: should z be 0 in local coords(i.e. add loc_tr.z())?
//     }
//     ROOT::Math::XYZPoint local = ROOT::Math::XYZPoint(lx, ly, lz);
//     LOG(TRACE) << "Transformed global point (" << global.x() << "|" << global.y() << "|" << global.z()
//                << ") to local point (" << local.x() << "|" << local.y() << "|" << local.z() << ").";
//     return local;
// }


// no need to implement as everything is in build_axes (See l78, l62 in pxdet.cpp)?
XYZVector BentPixelDetector::getSpatialResolutionXYZ(double column, double row) {

    double theta = 0.0;
    if(m_bent_axis == BentAxis::COLUMN) {
        theta = m_pitch.X() * (column - (m_nPixels.x() - 1) / 2.) / (m_radius);
    }
    else {
        theta = m_pitch.Y() * (row - (m_nPixels.y() - 1) / 2.) / (m_radius);
    }
    LOG(INFO) << "getSpatialResolutionXYZ: theta = " << theta;
    LOG(INFO) << "getSpatialResolutionXYZ: col, row = " << column << "," << row ;
    LOG(INFO) << "getSpatialResolutionXYZ: m_spatial_resolution = " << m_spatial_resolution;
    XYZVector sp_reso{m_spatial_resolution.x(), m_spatial_resolution.y(), m_spatial_resolution.x()};

    sp_reso.SetX(m_spatial_resolution.x() * cos (theta)); 
    sp_reso.SetZ(m_spatial_resolution.x() * sin (theta));
    // m_spatial_resolution.SetY(m_spatial_resolution.y() * sin (theta));
    LOG(INFO) << "getSpatialResolutionXYZ: sp_reso X  = " << sp_reso.x();
    LOG(INFO) << "getSpatialResolutionXYZ: sp_reso Y  = " << sp_reso.y();
    LOG(INFO) << "getSpatialResolutionXYZ: sp_reso Z  = " << sp_reso.z();
    return sp_reso;
}

XYVector BentPixelDetector::getSpatialResolution(double column, double row) {
    XYZVector res = this->getSpatialResolutionXYZ(column,row);
    LOG(INFO) << "getSpatialResolution: m_sp_res X  = " << res.x();
    LOG(INFO) << "getSpatialResolution: m_sp_res Y  = " << res.y();
    return XYVector(res.x(), res.y());
}

TMatrixD BentPixelDetector::getSpatialResolutionMatrixGlobal(double column, double row) {
    TMatrixD errorMatrix(3, 3);
    TMatrixD locToGlob(3, 3), globToLoc(3, 3);

    double theta = 0.0;
    if(m_bent_axis == BentAxis::COLUMN) {
        theta = m_pitch.X() * (column - (m_nPixels.x() - 1) / 2.) / (m_radius);
        LOG(INFO) << "getSpatialResolutionMatrixGlobal: rad = " << m_radius;
        LOG(INFO) << "getSpatialResolutionMatrixGlobal: col = " << column;
        // std::cout << "COLUMN = " << column << std::endl;
        LOG(INFO) << "getSpatialResolutionMatrixGlobal: theta  = " << theta;
        LOG(INFO) << "getSpatialResolutionMatrixGlobal: getSpatialResolution(column,row) = " << getSpatialResolution(column,row);
    }
    else {
        theta = m_pitch.Y() * (row - (m_nPixels.y() - 1) / 2.) / (m_radius);
    }
    errorMatrix(0, 0) = getSpatialResolutionXYZ(column,row).x() * getSpatialResolutionXYZ(column,row).x();
    errorMatrix(1, 1) = getSpatialResolutionXYZ(column,row).y() * getSpatialResolutionXYZ(column,row).y();
    // errorMatrix(2, 2) = getSpatialResolution(column,row).x() * getSpatialResolution(column,row).x() * sin (theta);
    errorMatrix(2, 2) = getSpatialResolutionXYZ(column,row).z() * getSpatialResolutionXYZ(column,row).z();

    LOG(WARNING) << "getSpatialResolutionMatrixGlobal: errorMatrix = ";
    errorMatrix.Print();
    LOG(INFO) << "getSpatialResolutionMatrixGlobal: l2g = " << alignment_->local2global();
    LOG(INFO) << "getSpatialResolutionMatrixGlobal: l2g rot = " << alignment_->local2global().Rotation();
    LOG(INFO) << "getSpatialResolutionMatrixGlobal: g2l = " << alignment_->global2local();
    LOG(INFO) << "getSpatialResolutionMatrixGlobal: g2l rot = " << alignment_->global2local().Rotation();
    alignment_->local2global().Rotation().GetRotationMatrix(locToGlob);
    alignment_->global2local().Rotation().GetRotationMatrix(globToLoc);

    m_spatial_resolution_matrix_global = locToGlob * errorMatrix * globToLoc;
    LOG(WARNING) << "getSpatialResolutionMatrixGlobal: m_spatial_resolution_matrix_global:";
    m_spatial_resolution_matrix_global.Print();

    return m_spatial_resolution_matrix_global;
}

PositionVector3D<Cartesian3D<double>> BentPixelDetector::getLocalPosition(double column, double row) const {
    ROOT::Math::RhoZPhiVector local;
    // (almost) almighty  cylinders - define "zero degree" at center of chip
    if(m_bent_axis == BentAxis::COLUMN) {
        // eversion possible by change or sign of the radius; but.. does it matter?
        local = ROOT::Math::RhoZPhiVector(m_radius, (m_pitch.Y() * (row - (m_nPixels.y() - 1) / 2.)), m_pitch.X() * (column - (m_nPixels.x() - 1) / 2.)/(-m_radius));
        LOG(INFO) << "getLocalPosition: col, row = " << column << "," << row ;
        LOG(INFO) << "getLocalPosition: local = " << local;
    } else {
        // not yet corrected; should just be a swap
        local = ROOT::Math::RhoZPhiVector(m_radius, (m_pitch.X() * (column - (m_nPixels.x() - 1) / 2.))/(-m_radius), m_pitch.Y() * (row - (m_nPixels.y() - 1) / 2.));
    }
    LOG(INFO) << "local cart = " << static_cast<PositionVector3D<Cartesian3D<double>>>(local);
    return static_cast<PositionVector3D<Cartesian3D<double>>>(local);
}

// PositionVector3D<Cartesian3D<double>> BentPixelDetector::getIntercept(const Track* track) const {
//     // Get and transform track state and direction
//     PositionVector3D<Cartesian3D<double>> state_track = track->getState(m_detectorName);
//     DisplacementVector3D<Cartesian3D<double>> direction_track = track->getDirection(m_detectorName);

//     // Bring track to local (transformed) coordinate system
//     state_track = toLocal() * state_track;
//     direction_track = toLocal().Rotation() * direction_track;
//     direction_track = direction_track.Unit();

//     // From globalPlanarIntercept get intercept with bent surface of pixel detector
//     InterceptParameters intercept_parameters;
//     ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> state_cylinder;
//     ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>> direction_cylinder;

//     if(m_bent_axis == BentAxis::COLUMN) {
//         state_cylinder = {0, 0, m_radius};
//         direction_cylinder = {0, 1, 0}; // cylinder axis along y (row)
//     } else {
//         state_cylinder = {0, this->getSize().Y() / 2, m_radius};
//         direction_cylinder = {1, 0, 0}; // cylinder axis along x (column)
//     }
//     LOG(INFO) << " st,dt,sc,dc = " << state_track << ", " << direction_track << ", " << state_cylinder << ", "<< direction_cylinder;
//     get_intercept_parameters(state_track, direction_track, state_cylinder, direction_cylinder, intercept_parameters);

//     // Select solution according to bending direction
//     PositionVector3D<Cartesian3D<double>> localBentIntercept;
//     double inf = std::numeric_limits<double>::infinity();

//     if(!intercept_parameters.isValid()) {
//         LOG(INFO) << "No intercept of track and cylinder sensor found. Return intercept outside of sensor area.";
//         localBentIntercept.SetCoordinates(inf, inf, inf);
//         return localBentIntercept;
//     }
//     auto param = (m_radius < 0) ? intercept_parameters.getParam1()
//                                 : intercept_parameters.getParam2(); // for negative radius select smaller solution
//     localBentIntercept = ROOT::Math::XYZPoint(state_track.x() + param * direction_track.x(),
//                                               state_track.y() + param * direction_track.y(),
//                                               state_track.z() + param * direction_track.z());
//     LOG(INFO) << " localBentIntercept = " << localBentIntercept;
//     return toGlobal() * localBentIntercept;
// }

// void BentPixelDetector::get_intercept_parameters(const PositionVector3D<Cartesian3D<double>>& state_track,
//                                                  const DisplacementVector3D<Cartesian3D<double>>& direction_track,
//                                                  const PositionVector3D<Cartesian3D<double>>& state_cylinder,
//                                                  const DisplacementVector3D<Cartesian3D<double>>& direction_cylinder,
//                                                  InterceptParameters& result) const {
//     // Get factors for quadratic equation: alpha z^2 + beta * z + gamma
//     double alpha;
//     double beta;
//     double gamma;

//     int component =
//         (direction_cylinder.X() != 0 && direction_cylinder.Y() == 0)
//             ? 0
//             : 1; // no other possibilities (see getIntercept's if); cylinder along x (column) : cylinder along y (row)
//     alpha = direction_track.Mag2() -
//             (component == 0 ? direction_track.X() * direction_track.X() : direction_track.Y() * direction_track.Y());
//     beta = 2 * (state_track - state_cylinder).Dot(direction_track) -
//            2 * (component == 0 ? direction_track.X() * (state_track.X() - state_cylinder.X())
//                                : direction_track.Y() * (state_track.Y() - state_cylinder.Y()));
//     gamma = state_track.Mag2() + state_cylinder.Mag2() -
//             (component == 0 ? state_track.X() * state_track.X() + state_cylinder.X() * state_cylinder.X()
//                             : state_track.Y() * state_track.Y() + state_cylinder.Y() * state_cylinder.Y()) -
//             2 * ((component == 0 ? state_track.Y() * state_cylinder.Y() : state_track.X() * state_cylinder.X()) +
//                  state_track.Z() * state_cylinder.Z()) -
//             (m_radius * m_radius);
//     // Check if the quadratic equation has a solution
//     double discriminant = beta * beta - 4 * alpha * gamma;
//     LOG(INFO) << " alpha,beta,gamma, disc = " << alpha << ", " << beta << ", " << gamma << ", " << discriminant;

//     if(discriminant < 0 || alpha == 0) { // must have real solution
//         return;
//     }
    
//     // Solve equation
//     result.setParams((-beta + sqrt(discriminant)) / (2 * alpha), (-beta - sqrt(discriminant)) / (2 * alpha));
//     return;
// }

// bool BentPixelDetector::hasIntercept(const Track* track, double pixelTolerance) const {

//     bool intercept = PixelDetector::hasIntercept(track, pixelTolerance);
//     if(!intercept)
//         return intercept;

//     // Get the global intercept
//     PositionVector3D<Cartesian3D<double>> globalIntercept = this->getIntercept(track);

//     // Check if the cylinder surface - track intercept parameters are valid
//     std::vector<double> coords(3);
//     globalIntercept.GetCoordinates(coords.begin(), coords.end());
//     for(auto coord : coords) {
//         if(coord == std::numeric_limits<double>::infinity()) {
//             LOG(TRACE) << "Invalid cylinder surface - track intercept parameters!";
//             intercept = false;
//         }
//     }

//     return intercept;
// }
