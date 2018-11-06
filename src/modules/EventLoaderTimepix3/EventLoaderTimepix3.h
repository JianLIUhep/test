#ifndef TIMEPIX3EVENTLOADER_H
#define TIMEPIX3EVENTLOADER_H 1

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <stdio.h>
#include "core/module/Module.hpp"
#include "objects/Pixel.h"
#include "objects/SpidrSignal.h"

namespace corryvreckan {
    /** @ingroup Modules
     */
    class EventLoaderTimepix3 : public Module {

    public:
        // Constructors and destructors
        EventLoaderTimepix3(Configuration config, std::shared_ptr<Detector> detector);
        ~EventLoaderTimepix3() {}

        // Standard algorithm functions
        void initialise();
        StatusCode run(Clipboard* clipboard);

    private:
        std::shared_ptr<Detector> m_detector;

        // ROOT graphs
        TH1F* pixelToT_beforecalibration;
        TH1F* pixelToT_aftercalibration;
        TH2F* pixelTOTParameterA;
        TH2F* pixelTOTParameterB;
        TH2F* pixelTOTParameterC;
        TH2F* pixelTOTParameterT;
        TH2F* pixelTOAParameterC;
        TH2F* pixelTOAParameterD;
        TH2F* pixelTOAParameterT;
        TH1F* timeshiftPlot;

        bool loadData(Clipboard* clipboard, Pixels*, SpidrSignals*);
        void loadCalibration(std::string path, char delim, std::vector<std::vector<float>>& dat);
        void maskPixels(std::string);

        // configuration paramaters:
        double m_triggerLatency;
        std::string m_inputDirectory;
        int m_minNumberOfPlanes;

        bool temporalSplit;
        size_t m_numberPixelHits;

        bool applyCalibration;
        std::string calibrationPath;
        std::string threshold;

        std::vector<std::vector<float>> vtot;
        std::vector<std::vector<float>> vtoa;

        // Member variables
        std::vector<std::unique_ptr<std::ifstream>> m_files;
        std::vector<std::unique_ptr<std::ifstream>>::iterator m_file_iterator;

        unsigned long long int m_syncTime;
        bool m_clearedHeader;
        long long int m_syncTimeTDC;
        int m_TDCoverflowCounter;

        long long int m_currentEvent;

        unsigned long long int m_prevTime;
        bool m_shutterOpen;
    };
} // namespace corryvreckan
#endif // TIMEPIX3EVENTLOADER_H
