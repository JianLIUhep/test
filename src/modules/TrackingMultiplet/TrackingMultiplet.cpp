/**
 * @file
 * @brief Implementation of module TrackingMultiplet
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "TrackingMultiplet.h"
#include <TDirectory.h>

using namespace corryvreckan;

TrackingMultiplet::TrackingMultiplet(Configuration config, std::vector<std::shared_ptr<Detector>> detectors)
    : Module(std::move(config), std::move(detectors)) {

    // timing cut, relative (x * time_resolution) or absolute:
    if(m_config.count({"time_cut_rel", "time_cut_abs"}) > 1) {
        throw InvalidCombinationError(
            m_config, {"time_cut_rel", "time_cut_abs"}, "Absolute and relative time cuts are mutually exclusive.");
    } else if(m_config.has("time_cut_abs")) {
        double time_cut_abs_ = m_config.get<double>("time_cut_abs");
        for(auto& detector : get_detectors()) {
            time_cuts_[detector] = time_cut_abs_;
        }
    } else {
        double time_cut_rel_ = m_config.get<double>("time_cut_rel", 3.0);
        for(auto& detector : get_detectors()) {
            time_cuts_[detector] = detector->getTimeResolution() * time_cut_rel_;
        }
    }

    // spatial cut, relative (x * spatial_resolution) or absolute:
    if(m_config.count({"spatial_cut_rel", "spatial_cut_abs"}) > 1) {
        throw InvalidCombinationError(
            m_config, {"spatial_cut_rel", "spatial_cut_abs"}, "Absolute and relative spatial cuts are mutually exclusive.");
    } else if(m_config.has("spatial_cut_abs")) {
        auto spatial_cut_abs_ = m_config.get<XYVector>("spatial_cut_abs");
        for(auto& detector : get_detectors()) {
            spatial_cuts_[detector] = spatial_cut_abs_;
        }
    } else {
        // default is 3.0 * spatial_resolution
        auto spatial_cut_rel_ = m_config.get<double>("spatial_cut_rel", 3.0);
        for(auto& detector : get_detectors()) {
            spatial_cuts_[detector] = detector->getSpatialResolution() * spatial_cut_rel_;
        }
    }

    // Read the scatterer position and the up- and downstream detectors
    // FIXME: Use DUT position
    auto dut_vector = get_duts();
    if(dut_vector.size() == 1) {
        m_config.setDefault<double>("scatterer_position", dut_vector.at(0)->displacement().Z());
    }

    scatterer_position_ = m_config.get<double>("scatterer_position");
    LOG(DEBUG) << "Set scatterer position: " << Units::display(scatterer_position_, {"mm", "m"});

    // Use detectors before and after scatterer as up- and downstream detectors
    std::vector<std::string> default_upstream_detectors, default_downstream_detectors;
    for(auto& detector : get_detectors()) {
        if(detector->isDUT() || detector->isAuxiliary()) {
            continue;
        }
        if(detector->displacement().Z() < scatterer_position_) {
            default_upstream_detectors.push_back(detector->getName());
        } else if(detector->displacement().Z() > scatterer_position_) {
            default_downstream_detectors.push_back(detector->getName());
        }
    }

    // Get the strings of detectors and translate it to shared_ptr<Detector>
    std::vector<std::string> upstream_detectors_str =
        m_config.getArray<std::string>("upstream_detectors", default_upstream_detectors);
    std::vector<std::string> downstream_detectors_str =
        m_config.getArray<std::string>("downstream_detectors", default_downstream_detectors);

    for(auto detectorID : upstream_detectors_str) {
        m_upstream_detectors.push_back(get_detector(detectorID));
    }
    for(auto detectorID : downstream_detectors_str) {
        m_downstream_detectors.push_back(get_detector(detectorID));
    }

    if(m_upstream_detectors.size() < 2) {
        throw InvalidValueError(m_config, "upstream_detectors", "At least two upstream detectors have to be provided.");
    }
    if(m_downstream_detectors.size() < 2) {
        throw InvalidValueError(m_config, "downstream_detectors", "At least two downstream detectors have to be provided.");
    }

    double max_z_upstream = std::numeric_limits<double>::min();
    for(auto detector : m_upstream_detectors) {
        if(detector->isDUT()) {
            LOG(WARNING) << "DUT listed as upstream detector. Update of configuration or geometry should be considered.";
        }
        if(detector->isAuxiliary()) {
            throw InvalidValueError(
                m_config, "upstream_detectors", "Auxiliary device listed as upstream detector. This is not supported.");
        }
        if(detector->displacement().Z() > max_z_upstream) {
            max_z_upstream = detector->displacement().Z();
        }
    }

    double min_z_downstream = std::numeric_limits<double>::max();
    for(auto detector : m_downstream_detectors) {
        if(detector->isDUT()) {
            LOG(WARNING) << "DUT listed as downstream detector. Update of configuration or geometry should be considered.";
        }
        if(detector->isAuxiliary()) {
            throw InvalidValueError(
                m_config, "downstream_detectors", "Auxiliary device listed as downstream detector. This is not supported.");
        }
        if(std::find(m_upstream_detectors.begin(), m_upstream_detectors.end(), detector) != m_upstream_detectors.end()) {
            throw InvalidCombinationError(m_config,
                                          {"upstream_detectors", "downstream_detectors"},
                                          "Detector " + detector->getName() +
                                              " is listed both as upstream and downstream detector.");
        }
        if(detector->displacement().Z() < min_z_downstream) {
            min_z_downstream = detector->displacement().Z();
        }
    }

    if(max_z_upstream > min_z_downstream) {
        throw InvalidCombinationError(m_config,
                                      {"upstream_detectors", "downstream_detectors"},
                                      "Last upstream detector is located behind first downstream detector.");
    }

    min_hits_upstream_ = m_config.get<size_t>("min_hits_upstream", m_upstream_detectors.size());
    min_hits_downstream_ = m_config.get<size_t>("min_hits_downstream", m_downstream_detectors.size());
    if(min_hits_upstream_ > m_upstream_detectors.size() || min_hits_upstream_ < 2) {
        throw InvalidValueError(
            m_config, "min_hits_upstream", "Number has to be 2 <= n <= " + to_string(m_upstream_detectors.size()));
    }
    if(min_hits_downstream_ > m_downstream_detectors.size() || min_hits_downstream_ < 2) {
        throw InvalidValueError(
            m_config, "min_hits_downstream", "Number has to be 2 <= n <= " + to_string(m_downstream_detectors.size()));
    }
    if(min_hits_upstream_ == 2) {
        LOG(WARNING) << "Number of required upstream hits equals 2. This leads to an underconstrained track fit.";
    }
    if(min_hits_downstream_ == 2) {
        LOG(WARNING) << "Number of required downstream hits equals 2. This leads to an underconstrained track fit.";
    }

    if(scatterer_position_ < max_z_upstream) {
        throw InvalidCombinationError(m_config,
                                      {"upstream_detectors", "scatterer_position"},
                                      "Scatterer position is located in front of last upstream detector.");
    }
    if(scatterer_position_ > min_z_downstream) {
        throw InvalidCombinationError(m_config,
                                      {"downstream_detectors", "scatterer_position"},
                                      "Scatterer position is located behind first downstream detector.");
    }

    scatterer_matching_cut_ = m_config.get<double>("scatterer_matching_cut");
    isolation_cut_ = m_config.get<double>("isolation_cut", scatterer_matching_cut_ * 2.);
}

void TrackingMultiplet::initialise() {

    std::string title = "Multiplet multiplicity;multiplets;events";
    multipletMultiplicity = new TH1F("multipletMultiplicity", title.c_str(), 40, 0, 40);

    title = "Track #chi^{2};#chi^{2};events";
    trackChi2 = new TH1F("trackChi2", title.c_str(), 150, 0, 150);

    title = "Track #chi^{2}/ndof;#chi^{2}/ndof;events";
    trackChi2ndof = new TH1F("trackChi2ndof", title.c_str(), 100, 0, 50);

    title = "Matching distance X at scatterer;distance x [mm];multiplet candidates";
    matchingDistanceAtScattererX = new TH1F("matchingDistanceAtScattererX", title.c_str(), 200, -10., 10.);
    title = "Matching distance Y at scatterer;distance y [mm];multiplet candidates";
    matchingDistanceAtScattererY = new TH1F("matchingDistanceAtScattererY", title.c_str(), 200, -10., 10.);

    title = "Multiplet offset X at scatterer;offset x [um];multiplets";
    multipletOffsetAtScattererX = new TH1F("multipletOffsetAtScattererX", title.c_str(), 200, -300., 300.);
    title = "Multiplet offset Y at scatterer;offset y [um];multiplets";
    multipletOffsetAtScattererY = new TH1F("multipletOffsetAtScattererY", title.c_str(), 200, -300., 300.);

    title = "Multiplet kink X at scatterer;kink x [mrad];multiplets";
    multipletKinkAtScattererX = new TH1F("multipletKinkAtScattererX", title.c_str(), 200, -20., 20.);
    title = "Multiplet kink Y at scatterer;kink y [mrad];multiplets";
    multipletKinkAtScattererY = new TH1F("multipletKinkAtScattererY", title.c_str(), 200, -20., 20.);

    for(auto stream : {upstream, downstream}) {
        std::string stream_name = stream == upstream ? "upstream" : "downstream";
        std::string stream_name_caps = stream == upstream ? "Upstream" : "Downstream";

        TDirectory* directory = getROOTDirectory();
        TDirectory* local_directory = directory->mkdir(stream_name.c_str());

        if(local_directory == nullptr) {
            throw RuntimeError("Cannot create or access local ROOT directory for module " + this->getUniqueName());
        }
        local_directory->cd();

        title = "";
        std::string hist_name = "";

        title = stream_name_caps + " tracklet multiplicity;" + stream_name + " tracklets;events";
        hist_name = stream_name + "Multiplicity";
        trackletMultiplicity[stream] = new TH1F(hist_name.c_str(), title.c_str(), 40, 0, 40);

        title = "Clusters per " + stream_name_caps + " tracklet;clusters;" + stream_name + " tracklets";
        hist_name = stream_name + "ClustersPerTracklet";
        clustersPerTracklet[stream] = new TH1F(hist_name.c_str(), title.c_str(), 10, 0, 10);

        title = stream_name_caps + " tracklet angle X;angle x [mrad];" + stream_name + " tracklets";
        hist_name = stream_name + "AngleX";
        trackletAngleX[stream] = new TH1F(hist_name.c_str(), title.c_str(), 250, -25., 25.);
        title = stream_name_caps + " tracklet angle Y;angle y [mrad];" + stream_name + " tracklets";
        hist_name = stream_name + "AngleY";
        trackletAngleY[stream] = new TH1F(hist_name.c_str(), title.c_str(), 250, -25., 25.);

        title = stream_name_caps + " tracklet X at scatterer;position x [mm];" + stream_name + " tracklets";
        hist_name = stream_name + "PositionAtScattererX";
        trackletPositionAtScattererX[stream] = new TH1F(hist_name.c_str(), title.c_str(), 200, -10., 10.);
        title = stream_name_caps + " tracklet Y at scatterer;position y [mm];" + stream_name + " tracklets";
        hist_name = stream_name_caps + "PositionAtScattererY";
        trackletPositionAtScattererY[stream] = new TH1F(hist_name.c_str(), title.c_str(), 200, -10., 10.);
    }

    // Loop over all up- and downstream planes
    std::vector<std::shared_ptr<Detector>> all_detectors;
    all_detectors.insert(all_detectors.end(), m_upstream_detectors.begin(), m_upstream_detectors.end());
    all_detectors.insert(all_detectors.end(), m_downstream_detectors.begin(), m_downstream_detectors.end());
    for(auto& detector : all_detectors) {
        std::string detectorID = detector->getName();

        TDirectory* directory = getROOTDirectory();
        TDirectory* local_directory = directory->mkdir(detectorID.c_str());

        if(local_directory == nullptr) {
            throw RuntimeError("Cannot create or access local ROOT directory for module " + this->getUniqueName());
        }
        local_directory->cd();

        title = detectorID + " Residual X;x_{track}-x [mm];events";
        residualsX[detectorID] = new TH1F("residualsX", title.c_str(), 500, -0.1, 0.1);
        title = detectorID + " Residual Y;y_{track}-y [mm];events";
        residualsY[detectorID] = new TH1F("residualsY", title.c_str(), 500, -0.1, 0.1);
    }
}

double TrackingMultiplet::calculate_average_timestamp(const Track* track) {
    double sum_weighted_time = 0;
    double sum_weights = 0;
    for(auto& cluster : track->clusters()) {
        double weight = 1 / (time_cuts_[get_detector(cluster->getDetectorID())]);
        double time_of_flight = static_cast<double>(Units::convert(cluster->global().z(), "mm") / (299.792458));
        sum_weights += weight;
        sum_weighted_time += (static_cast<double>(Units::convert(cluster->timestamp(), "ns")) - time_of_flight) * weight;
    }
    return (sum_weighted_time / sum_weights);
}

// Method containing the straight line tracklet finding for the arms of the multiplets
TrackVector TrackingMultiplet::find_multiplet_tracklets(const streams& stream,
                                                        std::map<std::shared_ptr<Detector>, KDTree*>& cluster_trees,
                                                        std::shared_ptr<Detector> reference_first,
                                                        std::shared_ptr<Detector> reference_last) {

    // Define upstream/downstream dependent variables
    size_t min_hits = stream == upstream ? min_hits_upstream_ : min_hits_downstream_;
    std::string stream_name = stream == upstream ? "upstream" : "downstream";

    // Choose reference detectors (first and last hit detector in the list)
    LOG(DEBUG) << "Start finding " + stream_name + " tracklets";

    TrackVector tracklets;

    // Tracklet finding
    for(auto& clusterFirst : cluster_trees[reference_first]->getAllClusters()) {
        for(auto& clusterLast : cluster_trees[reference_last]->getAllClusters()) {

            double time_cut = std::max(time_cuts_[reference_first], time_cuts_[reference_last]);
            if(std::fabs(clusterFirst->timestamp() - clusterLast->timestamp()) > time_cut) {
                LOG(DEBUG) << "Reference clusters not within time cuts.";
                continue;
            }

            auto trackletCandidate = new StraightLineTrack();
            trackletCandidate->addCluster(clusterFirst);
            trackletCandidate->addCluster(clusterLast);

            auto averageTimestamp = calculate_average_timestamp(trackletCandidate);
            trackletCandidate->setTimestamp(averageTimestamp);

            size_t detector_nr = 2;
            for(const auto& detector_tree : cluster_trees) {
                auto detector = detector_tree.first;
                if(detector == reference_first || detector == reference_last) {
                    continue;
                }

                detector_nr++;
                if(trackletCandidate->nClusters() + (cluster_trees.size() - detector_nr + 1) < min_hits) {
                    LOG(DEBUG) << "No chance to find a track - too few detectors left: " << trackletCandidate->nClusters()
                               << " + " << cluster_trees.size() << " - " << detector_nr << " < " << min_hits;
                    continue;
                }

                double timeCut =
                    std::max(std::min(time_cuts_[reference_first], time_cuts_[reference_last]), time_cuts_[detector]);
                LOG(DEBUG) << "Using timing cut of " << Units::display(timeCut, {"ns", "us", "s"});
                auto neighbours = detector_tree.second->getAllClustersInTimeWindow(trackletCandidate->timestamp(), timeCut);

                if(neighbours.empty()) {
                    LOG(DEBUG) << "No neighbours found within the correct time window.";
                    continue;
                }

                LOG(DEBUG) << "- found " << neighbours.size() << " neighbours within the correct time window";

                // Now let's see if there's a cluster matching in time and space.
                Cluster* closestCluster = nullptr;

                // Use spatial cut only as initial value (check if cluster is ellipse defined by cuts is done below):
                double closestClusterDistance = sqrt(spatial_cuts_[detector].x() * spatial_cuts_[detector].x() +
                                                     spatial_cuts_[detector].y() * spatial_cuts_[detector].y());

                // Now look for the spatially closest cluster on the next plane
                trackletCandidate->fit();

                double interceptX, interceptY;
                PositionVector3D<Cartesian3D<double>> interceptPoint = detector->getIntercept(trackletCandidate);
                interceptX = interceptPoint.X();
                interceptY = interceptPoint.Y();

                for(size_t ne = 0; ne < neighbours.size(); ne++) {
                    Cluster* newCluster = neighbours[ne];

                    // Calculate the distance to the previous plane's cluster/intercept
                    double distanceX = interceptX - newCluster->global().x();
                    double distanceY = interceptY - newCluster->global().y();
                    double distance = sqrt(distanceX * distanceX + distanceY * distanceY);

                    // Check if newCluster lies within ellipse defined by spatial cuts around intercept,
                    // following this example:
                    // https://www.geeksforgeeks.org/check-if-a-point-is-inside-outside-or-on-the-ellipse/
                    //
                    // ellipse defined by: x^2/a^2 + y^2/b^2 = 1: on ellipse,
                    //                                       > 1: outside,
                    //                                       < 1: inside
                    // Continue if outside of ellipse:

                    double norm = (distanceX * distanceX) / (spatial_cuts_[detector].x() * spatial_cuts_[detector].x()) +
                                  (distanceY * distanceY) / (spatial_cuts_[detector].y() * spatial_cuts_[detector].y());

                    if(norm > 1) {
                        LOG(DEBUG) << "Cluster outside the cuts. Normalized distance: " << norm;
                        continue;
                    }

                    // If this is the closest keep it for now
                    if(distance < closestClusterDistance) {
                        closestClusterDistance = distance;
                        closestCluster = newCluster;
                    }
                }

                if(closestCluster == nullptr) {
                    LOG(DEBUG) << "No cluster within spatial cut";
                    continue;
                }

                // Add the cluster to the tracklet
                trackletCandidate->addCluster(closestCluster);
                averageTimestamp = calculate_average_timestamp(trackletCandidate);
                trackletCandidate->setTimestamp(averageTimestamp);
                LOG(DEBUG) << "Added cluster to tracklet candidate";
            }

            if(trackletCandidate->nClusters() < min_hits) {
                LOG(DEBUG) << "Not enough clusters on the tracklet, found " << trackletCandidate->nClusters() << " but "
                           << min_hits << " required";
                delete trackletCandidate;
                continue;
            }

            LOG(DEBUG) << "Found good tracklet. Keeping this one.";
            trackletCandidate->fit();
            tracklets.push_back(trackletCandidate);
        }
    }

    // Check for isolation of tracklets
    std::vector<TrackVector::iterator> unisolatedTracklets;

    if(tracklets.size() > 1 && isolation_cut_ != 0) {
        for(TrackVector::iterator it0 = tracklets.begin(); it0 != tracklets.end(); ++it0) {
            auto positionAtScatterer = (*it0)->intercept(scatterer_position_);
            for(TrackVector::iterator it1 = it0 + 1; it1 != tracklets.end(); ++it1) {
                auto otherPositionAtScatterer = (*it1)->intercept(scatterer_position_);

                auto distance = otherPositionAtScatterer - positionAtScatterer;

                if(sqrt(distance.Mag2()) < isolation_cut_) {
                    LOG(DEBUG) << "Tracklet is not isolated. Distance (" << sqrt(distance.Mag2())
                               << ") smaller than the isolation cut (" << isolation_cut_
                               << "). Staging both tracklets for removal.";
                    unisolatedTracklets.push_back(it0);
                    unisolatedTracklets.push_back(it1);
                }
            }
        }
    }

    // Remove unisolated tracklets
    for(TrackVector::reverse_iterator rit = tracklets.rbegin(); rit != tracklets.rend(); ++rit) {
        if(std::find(unisolatedTracklets.begin(), unisolatedTracklets.end(), --rit.base()) != unisolatedTracklets.end()) {
            // Erase --rit.base(), since (reverse_iterator::base() = iterator + 1)
            LOG(DEBUG) << "Removing unisolated tracklet";
            delete *(--rit.base());
            tracklets.erase(--rit.base());
        }
    }

    // Get timestamp for tracklets
    for(auto& tracklet : tracklets) {
        double tracklet_timestamp = calculate_average_timestamp(tracklet);
        tracklet->setTimestamp(tracklet_timestamp);
    }

    return tracklets;
}

// Filling the histograms for up- & downstream tracklets
void TrackingMultiplet::fill_tracklet_histograms(const streams& stream, TrackVector tracklets) {

    std::string stream_name = stream == upstream ? "upstream" : "downstream";

    trackletMultiplicity[stream]->Fill(static_cast<double>(tracklets.size()));

    if(tracklets.size() > 0) {
        LOG(DEBUG) << "Filling plots for " << stream_name << " tracklets";

        for(auto& tracklet : tracklets) {
            clustersPerTracklet[stream]->Fill(static_cast<double>(tracklet->nClusters()));

            trackletAngleX[stream]->Fill(
                static_cast<double>(Units::convert(tracklet->direction("").X() / tracklet->direction("").Z(), "mrad")));
            trackletAngleY[stream]->Fill(
                static_cast<double>(Units::convert(tracklet->direction("").Y() / tracklet->direction("").Z(), "mrad")));

            trackletPositionAtScattererX[stream]->Fill(tracklet->intercept(scatterer_position_).X());
            trackletPositionAtScattererY[stream]->Fill(tracklet->intercept(scatterer_position_).Y());

            auto trackletClusters = tracklet->clusters();
            for(auto& trackletCluster : trackletClusters) {
                std::string detectorID = trackletCluster->detectorID();
                residualsX[detectorID]->Fill(tracklet->residual(detectorID).X());
                residualsY[detectorID]->Fill(tracklet->residual(detectorID).Y());
            }
        }
    }
}

StatusCode TrackingMultiplet::run(std::shared_ptr<Clipboard> clipboard) {

    LOG(DEBUG) << "Start of event";

    std::map<std::shared_ptr<Detector>, KDTree*> upstream_trees;
    std::map<std::shared_ptr<Detector>, KDTree*> downstream_trees;

    // Store upstream data in KDTrees and define reference detectors
    std::shared_ptr<Detector> reference_up_first = nullptr;
    std::shared_ptr<Detector> reference_up_last = nullptr;
    for(auto& upstream_detector : m_upstream_detectors) {
        auto upstream_detector_ID = upstream_detector->getName();
        LOG(DEBUG) << "Store data for upstream detector " << upstream_detector_ID;

        auto clusters = clipboard->getData<Cluster>(upstream_detector_ID);
        if(clusters == nullptr || clusters->size() == 0) {
            continue;
        }
        LOG(DEBUG) << "Cluster count: " << clusters->size();

        KDTree* clusterTree = new KDTree();
        clusterTree->buildTimeTree(*clusters);
        upstream_trees[upstream_detector] = clusterTree;

        if(reference_up_first == nullptr) {
            reference_up_first = upstream_detector;
        }
        reference_up_last = upstream_detector;
    }

    // Store downstream data in KDTrees and define reference detectors
    std::shared_ptr<Detector> reference_down_first = nullptr;
    std::shared_ptr<Detector> reference_down_last = nullptr;
    for(auto& downstream_detector : m_downstream_detectors) {
        auto downstream_detector_ID = downstream_detector->getName();
        LOG(DEBUG) << "Store data for downstream detector " << downstream_detector_ID;

        auto clusters = clipboard->getData<Cluster>(downstream_detector_ID);
        if(clusters == nullptr || clusters->size() == 0) {
            continue;
        }
        LOG(DEBUG) << "Cluster count: " << clusters->size();

        KDTree* clusterTree = new KDTree();
        clusterTree->buildTimeTree(*clusters);
        downstream_trees[downstream_detector] = clusterTree;

        if(reference_down_first == nullptr) {
            reference_down_first = downstream_detector;
        }
        reference_down_last = downstream_detector;
    }

    // Up- & downstream tracklet finding
    TrackVector upstream_tracklets;
    TrackVector downstream_tracklets;
    if(upstream_trees.size() >= min_hits_upstream_) {
        LOG(DEBUG) << "Reference detectors for upstream tracklet: " << reference_up_first->getName() << " & "
                   << reference_up_last->getName();
        upstream_tracklets = find_multiplet_tracklets(upstream, upstream_trees, reference_up_first, reference_up_last);
    } else {
        LOG(DEBUG) << "Too few hit detectors in upstream arm to find a tracklet";
    }
    if(downstream_trees.size() >= min_hits_downstream_) {
        LOG(DEBUG) << "Reference detectors for downstream tracklet: " << reference_down_first->getName() << " & "
                   << reference_down_last->getName();
        downstream_tracklets =
            find_multiplet_tracklets(downstream, downstream_trees, reference_down_first, reference_down_last);
    } else {
        LOG(DEBUG) << "Too few hit detectors in downstream arm to find a tracklet";
    }

    LOG(DEBUG) << "Found " << upstream_tracklets.size() << " upstream tracklets";
    LOG(DEBUG) << "Found " << downstream_tracklets.size() << " downstream tracklets";

    // Fill histograms for up- and downstream tracklets
    fill_tracklet_histograms(upstream, upstream_tracklets);
    fill_tracklet_histograms(downstream, downstream_tracklets);

    // Multiplet merging
    auto multiplets = std::make_shared<MultipletVector>();
    for(auto& uptracklet : upstream_tracklets) {
        Multiplet* multiplet = nullptr;

        double time_cut_upstream = std::numeric_limits<double>::max();
        for(auto& cluster : uptracklet->clusters()) {
            if(time_cuts_[get_detector(cluster->getDetectorID())] < time_cut_upstream) {
                time_cut_upstream = time_cuts_[get_detector(cluster->getDetectorID())];
            }
        }

        double closestMatchingDistance = scatterer_matching_cut_;
        TrackVector::iterator used_downtracklet;
        for(auto it = downstream_tracklets.begin(); it != downstream_tracklets.end(); ++it) {
            double time_cut_downstream = std::numeric_limits<double>::max();
            for(auto& cluster : (*it)->clusters()) {
                if(time_cuts_[get_detector(cluster->getDetectorID())] < time_cut_downstream) {
                    time_cut_downstream = time_cuts_[get_detector(cluster->getDetectorID())];
                }
            }

            // calculate time cut as the maximum of the minimal time cut of each tracklet.
            double time_cut = std::max(time_cut_upstream, time_cut_downstream);

            if(std::fabs((*it)->timestamp() - uptracklet->timestamp()) > time_cut) {
                LOG(DEBUG) << "Multiplet candidate discarded due to time cut";
                continue;
            }

            auto multipletCandidate = new Multiplet(uptracklet, (*it));
            LOG(DEBUG) << "Got new candidate.";

            multipletCandidate->setScattererPosition(scatterer_position_);
            multipletCandidate->fit();

            double distanceX = multipletCandidate->getOffsetAtScatterer().X();
            double distanceY = multipletCandidate->getOffsetAtScatterer().Y();
            double distance = sqrt(distanceX * distanceX + distanceY * distanceY);

            LOG(DEBUG) << "Multiplet candidate distance (x, y, abs): " << Units::display(distanceX, {"um"}) << "  "
                       << Units::display(distanceY, {"um"}) << "  " << Units::display(distance, {"um"});

            matchingDistanceAtScattererX->Fill(distanceX);
            matchingDistanceAtScattererY->Fill(distanceY);

            if(distance > scatterer_matching_cut_) {
                LOG(DEBUG) << "Multiplet candidate discarded due to high distance at scatterer";
                delete multipletCandidate;
                continue;
            }

            if(distance > closestMatchingDistance) {
                LOG(DEBUG) << "Multiplet candidate discarded - there's a closer match";
                delete multipletCandidate;
                continue;
            }

            LOG(DEBUG) << "Closest multiplet match so far. Proceed as candidate.";
            closestMatchingDistance = distance;
            delete multiplet;
            multiplet = multipletCandidate;
            used_downtracklet = it;
        }

        if(multiplet == nullptr) {
            LOG(DEBUG) << "No matching downstream tracklet found";
            continue;
        }

        LOG(DEBUG) << "Multiplet found";
        multiplet->setTimestamp(
            (multiplet->getUpstreamTracklet()->timestamp() + multiplet->getDownstreamTracklet()->timestamp()) / 2.);

        LOG(DEBUG) << "Deleting downstream tracklet";
        delete *used_downtracklet;
        downstream_tracklets.erase(used_downtracklet);

        multiplets->push_back(multiplet);

        trackChi2->Fill(multiplet->chi2());
        trackChi2ndof->Fill(multiplet->chi2ndof());

        double distanceX = multiplet->getOffsetAtScatterer().X();
        double distanceY = multiplet->getOffsetAtScatterer().Y();

        double kinkX = multiplet->getKinkAtScatterer().X();
        double kinkY = multiplet->getKinkAtScatterer().Y();

        multipletOffsetAtScattererX->Fill(static_cast<double>(Units::convert(distanceX, "um")));
        multipletOffsetAtScattererY->Fill(static_cast<double>(Units::convert(distanceY, "um")));

        multipletKinkAtScattererX->Fill(static_cast<double>(Units::convert(kinkX, "mrad")));
        multipletKinkAtScattererY->Fill(static_cast<double>(Units::convert(kinkY, "mrad")));
    }

    LOG(DEBUG) << "Found " << multiplets->size() << " multiplets";
    multipletMultiplicity->Fill(static_cast<double>(multiplets->size()));

    // Clean up tree and vector objects
    LOG(DEBUG) << "Cleaning up";
    for(auto tree = upstream_trees.cbegin(); tree != upstream_trees.cend();) {
        delete tree->second;
        tree = upstream_trees.erase(tree);
    }
    for(auto tree = downstream_trees.cbegin(); tree != downstream_trees.cend();) {
        delete tree->second;
        tree = downstream_trees.erase(tree);
    }

    for(auto& uptracklet : upstream_tracklets) {
        delete uptracklet;
    }
    for(auto& downtracklet : downstream_tracklets) {
        delete downtracklet;
    }
    upstream_tracklets.clear();
    downstream_tracklets.clear();

    if(multiplets->size() > 0) {
        clipboard->putData(multiplets);
    }

    // Return value telling analysis to keep running
    return StatusCode::Success;
}
