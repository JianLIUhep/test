#include "CLICpixEventLoader.h"

using namespace corryvreckan;
using namespace std;

CLICpixEventLoader::CLICpixEventLoader(Configuration config, std::vector<Detector*> detectors)
    : Module(std::move(config), std::move(detectors)) {
    m_filename = "";
}

void CLICpixEventLoader::initialise() {
    // File structure is RunX/CLICpix/RunX.dat

    // Take input directory from global parameters
    string inputDirectory = m_config.get<std::string>("inputDirectory") + "/CLICpix";

    // Open the root directory
    DIR* directory = opendir(inputDirectory.c_str());
    if(directory == nullptr) {
        LOG(ERROR) << "Directory " << inputDirectory << " does not exist";
        return;
    }
    dirent* entry;

    // Read the entries in the folder
    while((entry = readdir(directory))) {
        // Check for the data file
        string filename = inputDirectory + "/" + entry->d_name;
        if(filename.find(".dat") != string::npos) {
            m_filename = filename;
        }
    }

    // If no data was loaded, give a warning
    if(m_filename.length() == 0)
        LOG(WARNING) << "No data file was found for CLICpix in " << inputDirectory;

    // Open the data file for later
    m_file.open(m_filename.c_str());

    // Make histograms for debugging
    hHitMap = new TH2F("hitMap", "hitMap", 64, 0, 64, 64, 0, 64);
    hPixelToT = new TH1F("pixelToT", "pixelToT", 20, 0, 20);
    hShutterLength = new TH1F("shutterLength", "shutterLength", 3000, 0, 0.3);
    hPixelsPerFrame = new TH1F("pixelsPerFrame", "pixelsPerFrame", 4100, 0, 4100);
}

StatusCode CLICpixEventLoader::run(Clipboard* clipboard) {

    // Assume that the CLICpix is the DUT (if running this algorithm
    auto detector = get_dut();

    // If have reached the end of file, close it and exit program running
    if(m_file.eof()) {
        m_file.close();
        return Failure;
    }

    // Otherwise load a new frame

    // Pixel container, shutter information
    Pixels* pixels = new Pixels();
    double shutterStartTime = 0, shutterStopTime = 0;
    string data;

    int npixels = 0;
    // Read file and load data
    while(getline(m_file, data)) {

        //    LOG(TRACE) <<"Data: "<<data;

        // If line is empty then we have finished this event, stop looping
        if(data.length() < 5)
            break;

        // Check if this is a header/shutter/power info
        if(data.find("PWR_RISE") != string::npos || data.find("PWR_FALL") != string::npos)
            continue;
        if(data.find("SHT_RISE") != string::npos) {
            // Read the shutter start time
            long int timeInt;
            string name;
            istringstream header(data);
            header >> name >> timeInt;
            shutterStartTime = static_cast<double>(timeInt) / (0.04);
            LOG(TRACE) << "Shutter rise time: " << Units::display(shutterStartTime, {"ns", "us", "s"});
            continue;
        }
        if(data.find("SHT_FALL") != string::npos) {
            // Read the shutter stop time
            long int timeInt;
            string name;
            istringstream header(data);
            header >> name >> timeInt;
            shutterStopTime = static_cast<double>(timeInt) / (0.04);
            LOG(TRACE) << "Shutter fall time: " << Units::display(shutterStopTime, {"ns", "us", "s"});
            continue;
        }

        // Otherwise load data
        int row, col, counter, tot(0);
        LOG(TRACE) << "Pixel data: " << data;
        istringstream pixelData(data);
        pixelData >> col >> row >> counter >> tot;
        tot++;
        row = 63 - row;
        LOG(TRACE) << "New pixel: " << col << "," << row << " with tot " << tot;

        // If this pixel is masked, do not save it
        if(detector->masked(col, row))
            continue;
        Pixel* pixel = new Pixel(detector->name(), row, col, tot, 0);
        pixels->push_back(pixel);
        npixels++;
        hHitMap->Fill(col, row);
        hPixelToT->Fill(tot);
    }

    // Now set the event time so that the Timepix3 data is loaded correctly
    clipboard->put_persistent("eventStart", shutterStartTime);
    clipboard->put_persistent("eventEnd", shutterStopTime);
    clipboard->put_persistent("eventLength", (shutterStopTime - shutterStartTime));

    LOG(TRACE) << "Loaded " << npixels << " pixels";
    // Put the data on the clipboard
    if(pixels->size() > 0)
        clipboard->put(detector->name(), "pixels", reinterpret_cast<Objects*>(pixels));

    // Fill histograms
    hPixelsPerFrame->Fill(npixels);
    hShutterLength->Fill(shutterStopTime - shutterStartTime);

    // Return value telling analysis to keep running
    return Success;
}

void CLICpixEventLoader::finalise() {

    LOG(DEBUG) << "Analysed " << m_eventNumber << " events";
}
