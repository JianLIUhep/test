#ifndef CLUSTERING4D_H
#define CLUSTERING4D_H 1

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include "core/module/Module.hpp"
#include "objects/Cluster.h"
#include "objects/Pixel.h"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class Clustering4D : public Module {

    public:
        // Constructors and destructors
        Clustering4D(Configuration config, std::shared_ptr<Detector> detector);
        ~Clustering4D() {}

        // Functions
        void initialise();
        StatusCode run(std::shared_ptr<Clipboard> clipboard);

    private:
        std::shared_ptr<Detector> m_detector;
        static bool sortByTime(Pixel* pixel1, Pixel* pixel2);
        void calculateClusterCentre(Cluster*);
        bool touching(Pixel*, Cluster*);
        bool closeInTime(Pixel*, Cluster*);

        // Cluster histograms
        TH1F* clusterSize;
        TH1F* clusterWidthRow;
        TH1F* clusterWidthColumn;
        TH1F* clusterTot;
        TH2F* clusterPositionGlobal;

        double timingCut;
        int neighbour_radius_row;
        int neighbour_radius_col;
    };
} // namespace corryvreckan
#endif // CLUSTERING4D_H
