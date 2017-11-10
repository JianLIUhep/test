#ifndef DataDump_H
#define DataDump_H 1

#include <dirent.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "core/algorithm/Algorithm.h"
#include "objects/Cluster.h"
#include "objects/Pixel.h"
#include "objects/Track.h"

namespace corryvreckan {
    class DataDump : public Algorithm {

    public:
        // Constructors and destructors
        DataDump(Configuration config, std::vector<Detector*> detectors);
        ~DataDump() {}

        // Functions
        void initialise();
        StatusCode run(Clipboard* clipboard);
        void finalise();

        // Member variables
        int m_eventNumber;
        std::string m_detector;
    };
}
#endif // DataDump_H