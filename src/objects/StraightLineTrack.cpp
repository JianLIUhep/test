/**
 * @file
 * @brief Implementation of StraightLine track object
 *
 * @copyright Copyright (c) 2017-2022 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "StraightLineTrack.hpp"
#include "Eigen/Dense"
#include "Track.hpp"
#include "core/detector/BentPixelDetector.hpp"
#include "core/utils/log.h"
#include "exceptions.h"

using namespace corryvreckan;

// PositionVector3D<Cartesian3D<double>> StraightLineTrack::getIntercept(const Track* track) {
//     return BentPixelDetector::getIntercept(track);
// }

ROOT::Math::XYPoint StraightLineTrack::distance(const Cluster* cluster) const {
    LOG(INFO) << "I am here: distance()!";
    if(get_plane(cluster->detectorID()) == nullptr) {
        throw MissingReferenceException(typeid(*this), typeid(Plane));
    }
    auto trackIntercept = get_plane(cluster->detectorID())->getToLocal() * getState(cluster->detectorID());
    LOG(INFO) << "slt: DETID = " << cluster->getDetectorID();
    LOG(INFO) << "slt: trackIntercept = " << trackIntercept;

    auto dist = cluster->local() - trackIntercept;
    LOG(INFO) << "slt: dist = " << dist;
    LOG(INFO) << "slt: dist^2 = " << ROOT::Math::XYPoint(dist.x(), dist.y());

    // Return the distance^2
    return ROOT::Math::XYPoint(dist.x(), dist.y());
}

ROOT::Math::XYPoint StraightLineTrack::getKinkAt(const std::string&) const {
    LOG(INFO) << "I am here: getKinkAt()!";
    return ROOT::Math::XYPoint(0, 0);
}

ROOT::Math::XYZPoint StraightLineTrack::get_m_state() const {
    return m_state;
}

ROOT::Math::XYZPoint StraightLineTrack::getState(const std::string& detectorID) const {
    LOG(TRACE) << "Requesting state at: " << detectorID;
    auto plane =
        std::find_if(planes_.begin(), planes_.end(), [&detectorID](Plane const& p) { return p.getName() == detectorID; });
    if(plane == planes_.end()) {
        throw MissingReferenceException(typeid(*this), typeid(Plane));
    }
    auto toGlobal = plane->getToGlobal();
    ROOT::Math::XYZVector planeU, planeV, planeN;
    toGlobal.Rotation().GetComponents(planeU, planeV, planeN);
    ROOT::Math::XYZPoint origin = toGlobal.Translation() * ROOT::Math::XYZPoint(0, 0, 0);

    ROOT::Math::XYZVector distance = m_state - origin;
    double pathLength;
    pathLength = -distance.Dot(planeN) / m_direction.Dot(planeN);

    // if (detectorID == "ALPIDE_3") {
    //     ROOT::Math::XYZPoint origin_cyl{0.,0.,18.11};
    //     ROOT::Math::XYZVector blabla = m_state - origin_cyl;
    //     LOG(INFO) << " blabla = " << blabla;
    //     pathLength = -distance.Dot(blabla) / m_direction.Dot(blabla);
    // }
    ROOT::Math::XYZPoint position = m_state + pathLength * m_direction;
    if(detectorID == "ALPIDE_3") {
        LOG(INFO) << " DETECTOR: " << detectorID;
        LOG(INFO) << " toglobal = " << toGlobal;
        LOG(INFO) << "planeU = " << planeU << ", "
                  << "planeV = " << planeV << ", "
                  << "planeN = " << planeN;
        LOG(INFO) << "origin = " << origin;
        LOG(INFO) << "m_state = " << m_state;
        LOG(INFO) << "distance = " << distance;
        LOG(INFO) << "pathLength = " << pathLength;
        LOG(INFO) << "position = " << position;
        LOG(INFO) << "m_direction = " << m_direction;
        LOG(INFO) << "-distance.Dot(planeN) = " << -distance.Dot(planeN);
        LOG(INFO) << "m_direction.Dot(planeN) = " << m_direction.Dot(planeN);
        LOG(INFO) << "pathLength * m_direction = " << pathLength * m_direction;
        LOG(INFO) << "m_direction.Dot(planeN.Unit()) = " << m_direction.Dot(planeN.Unit());
    }

    return position;
}

ROOT::Math::XYZVector StraightLineTrack::getDirection(const std::string&) const {
    LOG(INFO) << "I am here: getDirection()!";
    return m_direction;
}

ROOT::Math::XYZVector StraightLineTrack::getDirection(const double&) const {
    LOG(INFO) << "I am here: getDirection2()!";
    return m_direction;
}

void StraightLineTrack::calculateChi2() {
    LOG(INFO) << "I am here: calculateChi2()!";

    // Get the number of clusters
    // We do have a 2-dimensional offset(x_0,y_0) and slope (dx,dy). Each hit provides two measurements.
    // ndof_ = 2*num_planes - 4 = 2 * (num_planes -2)
    ndof_ = (track_clusters_.size() - 2) * 2;
    chi2_ = 0.;
    chi2ndof_ = 0.;

    // Loop over all clusters
    for(auto& cl : track_clusters_) {
        auto* cluster = cl.get();
        if(cluster == nullptr) {
            throw MissingReferenceException(typeid(*this), typeid(Cluster));
        }
        if(get_plane(cluster->detectorID()) == nullptr) {
            throw MissingReferenceException(typeid(*this), typeid(Plane));
        }

        // Get the distance and the error
        LOG(WARNING) << "chi2: detectorID = " << cluster->detectorID();
        auto intercept = get_plane(cluster->detectorID())->getToLocal() * getState(cluster->detectorID());
        if(cluster->detectorID() == "ALPIDE_3") {
            LOG(INFO) << "issue HERE!";
            // intercept = get_plane(cluster->detectorID())->getToLocal() * get_m_state();
            // intercept = cluster->getIntercept(*this);
            // Track* test = nullptr;
            // BentPixelDetector* test2 = nullptr;
            // test2->getIntercept(test);
            if(!bentPixelDetector_) {
                LOG(INFO) << "empty ptr";
            }
            intercept = bentPixelDetector_->getIntercept(this);
        }
        LOG(WARNING) << "chi2: intercept = " << intercept;
        LOG(INFO) << "chi2: get_m_state = " << get_m_state();
        LOG(INFO) << "chi2: getState(cluster->detectorID()) = " << getState(cluster->detectorID());
        auto dist = cluster->local() - intercept;
        LOG(WARNING) << "chi2: cluster->local() = " << cluster->local();
        LOG(INFO) << "chi2: dist = " << dist;
        double ex2 = cluster->errorX() * cluster->errorX();
        double ey2 = cluster->errorY() * cluster->errorY();
        chi2_ += ((dist.x() * dist.x() / ex2) + (dist.y() * dist.y() / ey2));
        LOG(INFO) << "chi2: chi2_ ++ = " << chi2_;
    }

    LOG(INFO) << "chi2: chi2_ = " << chi2_;
    // Store also the chi2/degrees of freedom
    chi2ndof_ = (ndof_ <= 0) ? -1 : (chi2_ / static_cast<double>(ndof_));
}

void StraightLineTrack::calculateResiduals() {
    LOG(INFO) << "I am here: calculateResiduals()!";
    for(const auto& c : track_clusters_) {
        auto* cluster = c.get();
        if(get_plane(cluster->detectorID()) == nullptr) {
            throw MissingReferenceException(typeid(*this), typeid(Plane));
        }
        residual_global_[cluster->detectorID()] = cluster->global() - getState(cluster->detectorID());
        residual_local_[cluster->detectorID()] =
            cluster->local() - get_plane(cluster->detectorID())->getToLocal() * getState(cluster->detectorID());

        if(cluster->detectorID() == "ALPIDE_3") {
            residual_global_[cluster->detectorID()] = cluster->global() - get_m_state();
            residual_local_[cluster->detectorID()] =
                cluster->local() - get_plane(cluster->detectorID())->getToLocal() * get_m_state();
        }
    }
}

double StraightLineTrack::operator()(const double* parameters) {
    LOG(INFO) << "I am here: operator()!";
    LOG(INFO) << "parameters: " << parameters;
    // Update the StraightLineTrack gradient and intercept
    this->m_direction.SetX(parameters[0]);
    this->m_state.SetX(parameters[1]);
    this->m_direction.SetY(parameters[2]);
    this->m_state.SetY(parameters[3]);

    // Calculate the chi2
    this->calculateChi2();

    // Return this to minuit
    return chi2_;
}

void StraightLineTrack::fit() {
    LOG(INFO) << "I am here: fit()!";
    isFitted_ = false;
    Eigen::Matrix4d mat(Eigen::Matrix4d::Zero());
    Eigen::Vector4d vec(Eigen::Vector4d::Zero());

    // Loop over all clusters and fill the matrices
    for(auto& cl : track_clusters_) {
        auto* cluster = cl.get();
        if(cluster == nullptr) {
            throw MissingReferenceException(typeid(*this), typeid(Cluster));
        }

        // Get the measurement and its errors
        double x = cluster->global().x();
        double y = cluster->global().y();
        double z = cluster->global().z();
        // LOG(WARNING) << " CLUSTER x,y,z = " << x << ", " << y << ", " << z;
        Eigen::Vector2d pos(x, y);
        auto errorMatrix = cluster->errorMatrixGlobal();
        Eigen::Matrix2d V;
        // LOG(INFO) << "eM = ";
        // errorMatrix.Print();
        if(errorMatrix(0, 0) < 1e-8) {
            errorMatrix(0, 0) = 1e-7;
            V << errorMatrix(0, 0), errorMatrix(0, 1), errorMatrix(1, 0), errorMatrix(1, 1);
        } else {
            V << errorMatrix(0, 0), errorMatrix(0, 1), errorMatrix(1, 0), errorMatrix(1, 1);
        }

        // if(cluster->getDetectorID()=="ALPIDE_3") V << errorMatrix(2, 2), errorMatrix(0, 1), errorMatrix(1, 0),
        // errorMatrix(1, 1);
        Eigen::Matrix<double, 2, 4> C;
        // C relates the position of a cluster in space to the track parameters
        // has two ros; so z is always fixed; do I need to change this?? otherwise it's always perpendicular
        C << 1., z, 0., 0., 0., 0., 1., z;
        // LOG(WARNING) << " C = " << C;
        // LOG(WARNING) << " V = " << V;

        // Fill the matrices
        if(fabs(V.determinant()) < std::numeric_limits<double>::epsilon()) {
            std::cout << "V = \n" << V << std::endl;
            std::cout << "epsilon = " << std::numeric_limits<double>::epsilon() << std::endl;
            std::cout << "determinant of V = " << V.determinant() << std::endl;
            throw TrackFitError(typeid(this), "Error matrix inversion in straight line fit failed");
        }
        vec += C.transpose() * V.inverse() * pos;
        mat += C.transpose() * V.inverse() * C;
    }
    // LOG(WARNING) << "mat f: " << mat;
    // LOG(WARNING) << " vec f: " << vec;

    // Check for singularities.
    if(fabs(mat.determinant()) < std::numeric_limits<double>::epsilon()) {
        throw TrackFitError(typeid(this), "Martix inversion in straight line fit failed");
    }
    // Get the StraightLineTrack parameters
    Eigen::Vector4d res = mat.inverse() * vec;
    // LOG(WARNING) << "res = " << res;
    // Set the StraightLineTrack parameters
    m_state.SetX(res(0));
    m_state.SetY(res(2));
    m_state.SetZ(0.);

    m_direction.SetX(res(1));
    m_direction.SetY(res(3));
    m_direction.SetZ(1.);

    // Calculate the chi2
    this->calculateChi2();
    this->calculateResiduals();
    isFitted_ = true;
}

ROOT::Math::XYZPoint StraightLineTrack::getIntercept(double z) const {
    LOG(INFO) << "I am here: getIntercept()!";
    LOG(INFO) << " slt: intercept = " << m_state + m_direction * z;
    return m_state + m_direction * z;
}

void StraightLineTrack::print(std::ostream& out) const {
    out << "StraightLineTrack " << this->m_state << ", " << this->m_direction << ", " << this->chi2_ << ", " << this->ndof_
        << ", " << this->chi2ndof_ << ", " << this->timestamp();
}
