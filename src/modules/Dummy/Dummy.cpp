#include "Dummy.h"

using namespace corryvreckan;
using namespace std;

Dummy::Dummy(Configuration config, std::vector<Detector*> detectors) : Module(std::move(config), std::move(detectors)) {}

void Dummy::initialise() {

    // Initialise histograms per device
    for(auto& detector : get_detectors()) {
        // Simple histogram per device
        string name = "plotForDevice_" + detector->name();
        plotPerDevice[detector->name()] = new TH2F(name.c_str(), name.c_str(), 256, 0, 256, 256, 0, 256);
    }

    // Initialise single histograms
    string name = "singlePlot";
    singlePlot = new TH1F(name.c_str(), name.c_str(), 1000, 0, 1000);

    // Initialise member variables
    m_eventNumber = 0;
}

StatusCode Dummy::run(Clipboard* clipboard) {

    // Loop over all Timepix3 and make plots
    for(auto& detector : get_detectors()) {

        // Get the pixels
        Pixels* pixels = reinterpret_cast<Pixels*>(clipboard->get(detector->name(), "pixels"));
        if(pixels == nullptr) {
            LOG(DEBUG) << "Detector " << detector->name() << " does not have any pixels on the clipboard";
            continue;
        }

        // Get the clusters
        Clusters* clusters = reinterpret_cast<Clusters*>(clipboard->get(detector->name(), "clusters"));
        if(clusters == nullptr) {
            LOG(DEBUG) << "Detector " << detector->name() << " does not have any clusters on the clipboard";
            continue;
        }

        // Loop over all pixels and make hitmaps
        for(auto& pixel : (*pixels)) {
            // Fill the plots for this device
            plotPerDevice[detector->name()]->Fill(pixel->column(), pixel->row());
        }
    }

    // Fill single histogram
    singlePlot->Fill(m_eventNumber);

    // Increment event counter
    m_eventNumber++;

    // Return value telling analysis to keep running
    return Success;
}

void Dummy::finalise() {

    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}
