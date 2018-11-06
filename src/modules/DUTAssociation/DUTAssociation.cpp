#include "DUTAssociation.h"

using namespace corryvreckan;
using namespace std;

DUTAssociation::DUTAssociation(Configuration config, std::shared_ptr<Detector> detector)
    : Module(std::move(config), detector), m_detector(detector) {

    timingCut = m_config.get<double>("timingCut", static_cast<double>(Units::convert(200, "ns")));
    spatialCut = m_config.get<XYVector>("spatialCut", 2 * m_detector->pitch());
}

StatusCode DUTAssociation::run(Clipboard* clipboard) {

    // Get the tracks from the clipboard
    Tracks* tracks = reinterpret_cast<Tracks*>(clipboard->get("tracks"));
    if(tracks == nullptr) {
        LOG(DEBUG) << "No tracks on the clipboard";
        return Success;
    }

    // Get the DUT clusters from the clipboard
    Clusters* clusters = reinterpret_cast<Clusters*>(clipboard->get(m_detector->name(), "clusters"));
    if(clusters == nullptr) {
        LOG(DEBUG) << "No DUT clusters on the clipboard";
        return Success;
    }

    // Loop over all tracks
    for(auto& track : (*tracks)) {
        // Loop over all DUT clusters
        for(auto& cluster : (*clusters)) {
            // Check distance between track and cluster
            ROOT::Math::XYZPoint intercept = track->intercept(cluster->globalZ());
            double xdistance = intercept.X() - cluster->globalX();
            double ydistance = intercept.Y() - cluster->globalY();
            if(abs(xdistance) > spatialCut.x() || abs(ydistance) > spatialCut.y()) {
                LOG(DEBUG) << "Discarding DUT cluster with distance (" << abs(xdistance) << "," << abs(ydistance) << ")";
                continue;
            }

            // Check if the cluster is close in time
            if(std::abs(cluster->timestamp() - track->timestamp()) > timingCut) {
                continue;
            }

            LOG(DEBUG) << "Found associated cluster with distance (" << abs(xdistance) << "," << abs(ydistance) << ")";
            track->addAssociatedCluster(cluster);
        }
    }

    // Return value telling analysis to keep running
    return Success;
}
