#include "EtaCorrection.h"

using namespace corryvreckan;
using namespace std;

EtaCorrection::EtaCorrection(Configuration config, std::vector<Detector*> detectors)
    : Algorithm(std::move(config), std::move(detectors)) {
    m_DUT = m_config.get<std::string>("DUT");
    m_chi2ndofCut = m_config.get<double>("chi2ndofCut", 100.);
}

void EtaCorrection::initialise() {

    // Initialise single histograms
    m_detector = get_detector(m_DUT);
    double pitchX = m_detector->pitchX();
    double pitchY = m_detector->pitchY();
    string name = "etaDistributionX";
    m_etaDistributionX = new TH2F(name.c_str(), name.c_str(), 100., 0., pitchX, 100., 0., pitchY);
    name = "etaDistributionY";
    m_etaDistributionY = new TH2F(name.c_str(), name.c_str(), 100., 0., pitchX, 100., 0., pitchY);
    name = "etaDistributionXprofile";
    m_etaDistributionXprofile = new TProfile(name.c_str(), name.c_str(), 100., 0., pitchX, 0., pitchY);
    name = "etaDistributionYprofile";
    m_etaDistributionYprofile = new TProfile(name.c_str(), name.c_str(), 100., 0., pitchX, 0., pitchY);
    name = "etaDistributionXcorrected";
    m_etaDistributionXcorrected = new TH2F(name.c_str(), name.c_str(), 100., 0., pitchX, 100., 0., pitchY);
    name = "etaDistributionYcorrected";
    m_etaDistributionYcorrected = new TH2F(name.c_str(), name.c_str(), 100., 0., pitchX, 100., 0., pitchY);

    // Initialise member variables
    m_eventNumber = 0;
}

StatusCode EtaCorrection::run(Clipboard* clipboard) {

    // Get the tracks from the clipboard
    Tracks* tracks = (Tracks*)clipboard->get("tracks");
    if(tracks == NULL) {
        LOG(DEBUG) << "No tracks on the clipboard";
        return Success;
    }

    // Loop over all tracks and look at the associated clusters to plot the eta distribution
    for(auto& track : (*tracks)) {

        // Cut on the chi2/ndof
        if(track->chi2ndof() > m_chi2ndofCut) {
            continue;
        }

        // Get the in-pixel track intercept
        PositionVector3D<Cartesian3D<double>> trackIntercept = m_detector->getIntercept(track);
        PositionVector3D<Cartesian3D<double>> trackInterceptLocal = *(m_detector->globalToLocal()) * trackIntercept;
        double pixelInterceptX = m_detector->inPixelX(trackInterceptLocal) / 1000.;
        double pixelInterceptY = m_detector->inPixelY(trackInterceptLocal) / 1000.;
        (pixelInterceptX > m_detector->pitchX() / 2. ? pixelInterceptX -= m_detector->pitchX() / 2.
                                                     : pixelInterceptX += m_detector->pitchX() / 2.);
        (pixelInterceptY > m_detector->pitchY() / 2. ? pixelInterceptY -= m_detector->pitchY() / 2.
                                                     : pixelInterceptY += m_detector->pitchY() / 2.);

        // Look at the associated clusters and plot the eta function
        for(auto& dutCluster : track->associatedClusters()) {

            // Ignore single pixel clusters
            if(dutCluster->size() == 1)
                continue;

            // Get the fraction along the pixel
            double inPixelX = m_detector->pitchX() * (dutCluster->column() - floor(dutCluster->column()));
            double inPixelY = m_detector->pitchY() * (dutCluster->row() - floor(dutCluster->row()));
            if(dutCluster->columnWidth() == 2) {
                m_etaDistributionX->Fill(inPixelX, pixelInterceptX);
                m_etaDistributionXprofile->Fill(inPixelX, pixelInterceptX);
            }
            if(dutCluster->rowWidth() == 2) {
                m_etaDistributionY->Fill(inPixelY, pixelInterceptY);
                m_etaDistributionYprofile->Fill(inPixelY, pixelInterceptY);
            }
        }
    }

    // Increment event counter
    m_eventNumber++;

    // Return value telling analysis to keep running
    return Success;
}

void EtaCorrection::finalise() {

    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}
