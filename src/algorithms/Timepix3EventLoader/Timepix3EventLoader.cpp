#include "Timepix3EventLoader.h"

#include <bitset>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace corryvreckan;
using namespace std;

Timepix3EventLoader::Timepix3EventLoader(Configuration config, std::vector<Detector*> detectors)
    : Algorithm(std::move(config), std::move(detectors)), temporalSplit(false), m_currentTime(0), m_currentEvent(0),
      m_prevTime(0), m_shutterOpen(false) {

    // Take input directory from global parameters
    m_inputDirectory = m_config.get<std::string>("inputDirectory");

    applyTimingCut = m_config.get<bool>("applyTimingCut", false);
    m_timingCut = m_config.get<double>("timingCut", 0.0);
    m_triggerLatency = m_config.get<double>("triggerLatency", 0.0);
    m_minNumberOfPlanes = m_config.get<int>("minNumerOfPlanes", 1);

    // Check whether event length or pixel count should be used to separate events:
    m_eventLength = m_config.get<double>("eventLength", 0.0);
    m_numberPixelHits = m_config.get<int>("number_of_pixelhits", 2000);
}

void Timepix3EventLoader::initialise() {

    if(m_config.has("eventLength")) {
        LOG(INFO) << "Event length set, splitting events by time.";
        temporalSplit = true;
    } else {
        LOG(INFO) << "Splitting events by number of pixel hits on detector plane.";
        temporalSplit = false;
    }

    // File structure is RunX/ChipID/files.dat

    // Open the root directory
    DIR* directory = opendir(m_inputDirectory.c_str());
    if(directory == NULL) {
        throw AlgorithmError("Directory " + m_inputDirectory + " does not exist");
    } else {
        LOG(TRACE) << "Found directory " << m_inputDirectory;
    }
    dirent* entry;
    dirent* file;

    // Buffer for file names:
    std::map<std::string, std::vector<std::string>> detector_files;

    // Read the entries in the folder
    while(entry = readdir(directory)) {

        // Ignore UNIX functional directories:
        if(std::string(entry->d_name).at(0) == '.') {
            continue;
        }

        // If these are folders then the name is the chip ID
        // For some file systems, dirent only returns DT_UNKNOWN - in this case, check the dir. entry starts with "W"
        if(entry->d_type == DT_DIR || (entry->d_type == DT_UNKNOWN && std::string(entry->d_name).at(0) == 'W')) {

            LOG(DEBUG) << "Found directory for detector " << entry->d_name;

            // Open the folder for this device
            string detectorID = entry->d_name;
            string dataDirName = m_inputDirectory + "/" + entry->d_name;
            DIR* dataDir = opendir(dataDirName.c_str());

            // Check if this device has conditions loaded and is a Timepix3
            Detector* detector;
            try {
                detector = get_detector(detectorID);
            } catch(AlgorithmError& e) {
                LOG(WARNING) << e.what();
                continue;
            }

            if(detector->type() != "Timepix3") {
                LOG(WARNING) << "Device with detector ID " << entry->d_name << " is not of type Timepix3.";
                continue;
            }

            // Get all of the files for this chip
            while(file = readdir(dataDir)) {
                string filename = dataDirName + "/" + file->d_name;

                // Check if file has extension .dat
                if(string(file->d_name).find(".dat") != string::npos) {
                    LOG(DEBUG) << "Enqueuing data file for " << detectorID << ": " << filename;
                    detector_files[detectorID].push_back(filename);
                }

                // If not a data file, it might be a trimdac file, with the list of masked pixels etc.
                if(string(file->d_name).find("trimdac") != string::npos) {
                    // Now that we have all of the data files and mask files for this detector, pass the mask file to
                    // parameters
                    LOG(INFO) << "Set mask file " << filename;
                    detector->setMaskFile(filename);

                    // Apply the pixel masking
                    maskPixels(detector, filename);
                }
            }
        }
    }

    // Check that we have files for every detector in the configuration file and sort them correctly:
    for(auto& detector : get_detectors()) {
        std::string detectorID = detector->name();

        if(detector_files.count(detector->name()) == 0) {
            LOG(ERROR) << "No data file found for detector " << detector->name();
            continue;
        }

        // Initialise null values for later
        m_syncTime[detectorID] = 0;
        m_clearedHeader[detectorID] = false;

        // Sort all files by extracting the "serial number" from the file name while ignoring the timestamp:
        std::sort(detector_files[detector->name()].begin(),
                  detector_files[detector->name()].end(),
                  [](std::string a, std::string b) {

                      auto get_serial = [](std::string name) {
                          const auto pos1 = name.find_last_of('-');
                          const auto pos2 = name.find_last_of('.');
                          return name.substr(pos1 + 1, pos2 - pos1 - 1);
                      };
                      return std::stoi(get_serial(a)) < std::stoi(get_serial(b));
                  });

        // Open them:
        for(auto& filename : detector_files[detectorID]) {

            auto new_file = std::make_unique<std::ifstream>(filename);
            if(new_file->is_open()) {
                LOG(DEBUG) << "Opened data file for " << detectorID << ": " << filename;

                // The header is repeated in every new data file, thus skip it for all.
                uint32_t headerID;
                if(!new_file->read(reinterpret_cast<char*>(&headerID), sizeof headerID)) {
                    throw AlgorithmError("Cannot read header ID for " + detectorID + " in file " + filename);
                }
                if(headerID != 1380208723) {
                    throw AlgorithmError("Incorrect header ID for " + detectorID + " in file " + filename + ": " +
                                         std::to_string(headerID));
                }
                LOG(TRACE) << "Header ID: \"" << headerID << "\"";

                // Skip the rest of the file header
                uint32_t headerSize;
                if(!new_file->read(reinterpret_cast<char*>(&headerSize), sizeof headerSize)) {
                    throw AlgorithmError("Cannot read header size for " + detectorID + " in file " + filename);
                }

                // Skip the full header:
                new_file->seekg(headerSize);
                LOG(TRACE) << "Skipped header (" << headerSize << "b)";

                // Store the file in the data vector:
                m_files[detectorID].push_back(std::move(new_file));
            } else {
                throw AlgorithmError("Could not open data file " + filename);
            }
        }
    }

    // Set the file iterator to the first file for every detector:
    for(auto& detector : m_files) {
        m_file_iterator[detector.first] = detector.second.begin();
    }
}

StatusCode Timepix3EventLoader::run(Clipboard* clipboard) {

    // This will loop through each timepix3 registered, and load data from each of them. This can
    // be done in one of two ways: by taking all data in the time interval (t,t+delta), or by
    // loading a fixed number of pixels (ie. 2000 at a time)

    LOG(TRACE) << "== New event";
    int loadedData = 0;

    // End the loop as soon as all detector files are finished:
    StatusCode returnvalue = Failure;

    // Loop through all registered detectors
    for(auto& detector : get_detectors()) {

        // Check if they are a Timepix3
        if(detector->type() != "Timepix3")
            continue;

        string detectorID = detector->name();
        // If all files for this detector have been read, ignore it:
        if(m_file_iterator[detectorID] == m_files[detectorID].end()) {
            continue;
        }

        // This detector has data left:
        returnvalue = Success;

        // Make a new container for the data
        Pixels* deviceData = new Pixels();
        SpidrSignals* spidrData = new SpidrSignals();

        // Load the next chunk of data
        bool data = loadData(clipboard, detector, deviceData, spidrData);

        // If data was loaded then put it on the clipboard
        if(data) {
            loadedData++;
            LOG(DEBUG) << "Loaded " << deviceData->size() << " pixels for device " << detectorID;
            clipboard->put(detectorID, "pixels", (TestBeamObjects*)deviceData);
        }
        clipboard->put(detectorID, "SpidrSignals", (TestBeamObjects*)spidrData);
    }

    // Increment the event time
    clipboard->put_persistent("currentTime", clipboard->get_persistent("currentTime") + m_eventLength);

    // Otherwise tell event loop to keep running
    IFLOG(INFO) {
        if(temporalSplit) {
            LOG_PROGRESS(INFO, "tpx3_loader")
                << "Current time: " << std::setprecision(4) << std::fixed << clipboard->get_persistent("currentTime");
        } else {
            LOG_PROGRESS(INFO, "tpx3_loader") << "Current event: " << m_currentEvent;
        }
    }

    // If no/not enough data in this event then tell the event loop to directly skip to the next event
    if(returnvalue == Success && loadedData < m_minNumberOfPlanes) {
        return NoData;
    }

    return returnvalue;
}

// Function to load the pixel mask file
void Timepix3EventLoader::maskPixels(Detector* detector, string trimdacfile) {

    // Open the mask file
    ifstream trimdacs;
    trimdacs.open(trimdacfile.c_str());

    // Ignore the file header
    string line;
    getline(trimdacs, line);
    int t_col, t_row, t_trim, t_mask, t_tpen;

    // Loop through the pixels in the file and apply the mask
    for(int col = 0; col < 256; col++) {
        for(int row = 0; row < 256; row++) {
            trimdacs >> t_col >> t_row >> t_trim >> t_mask >> t_tpen;
            if(t_mask)
                detector->maskChannel(t_col, t_row);
        }
    }

    // Close the files when finished
    trimdacs.close();
}

// Function to load data for a given device, into the relevant container
bool Timepix3EventLoader::loadData(Clipboard* clipboard, Detector* detector, Pixels* devicedata, SpidrSignals* spidrData) {

    string detectorID = detector->name();

    bool extra = false; // temp

    LOG(DEBUG) << "Loading data for device " << detectorID;
    while(1) {
        // Check if current file is at its end and move to the next:
        if((*m_file_iterator[detectorID])->eof()) {
            m_file_iterator[detectorID]++;
            LOG(INFO) << "Starting to read next file for " << detectorID << ": " << (*m_file_iterator[detectorID]).get();
        }

        // Check if the last file is finished:
        if(m_file_iterator[detectorID] == m_files[detectorID].end()) {
            LOG(INFO) << "EOF for all files of " << detectorID;
            break;
        }

        // Now read the data packets.
        ULong64_t pixdata = 0;
        UShort_t thr = 0;

        // If we can't read data anymore, jump to begin of loop:
        if(!(*m_file_iterator[detectorID])->read(reinterpret_cast<char*>(&pixdata), sizeof pixdata)) {
            LOG(INFO) << "No more data in current file for " << detectorID << ": " << (*m_file_iterator[detectorID]).get();
            continue;
        }

        LOG(TRACE) << "0x" << hex << pixdata << dec << " - " << pixdata;

        // Get the header (first 4 bits) and do things depending on what it is
        // 0x4 is the "heartbeat" signal, 0xA and 0xB are pixel data
        const UChar_t header = ((pixdata & 0xF000000000000000) >> 60) & 0xF;

        // Use header 0x4 to get the long timestamps (called syncTime here)
        if(header == 0x4) {
            LOG(TRACE) << "Found syncTime data";

            // The 0x4 header tells us that it is part of the timestamp
            // There is a second 4-bit header that says if it is the most
            // or least significant part of the timestamp
            const UChar_t header2 = ((pixdata & 0x0F00000000000000) >> 56) & 0xF;

            // This is a bug fix. There appear to be errant packets with garbage data
            // - source to be tracked down.
            // Between the data and the header the intervening bits should all be 0,
            // check if this is the case
            const UChar_t intermediateBits = ((pixdata & 0x00FF000000000000) >> 48) & 0xFF;
            if(intermediateBits != 0x00) {
                LOG(DEBUG) << "Detector " << detectorID << ": intermediateBits error";
                continue;
            }

            // 0x4 is the least significant part of the timestamp
            if(header2 == 0x4) {
                // The data is shifted 16 bits to the right, then 12 to the left in
                // order to match the timestamp format (net 4 right)
                m_syncTime[detectorID] =
                    (m_syncTime[detectorID] & 0xFFFFF00000000000) + ((pixdata & 0x0000FFFFFFFF0000) >> 4);
            }
            // 0x5 is the most significant part of the timestamp
            if(header2 == 0x5) {
                // The data is shifted 16 bits to the right, then 44 to the left in
                // order to match the timestamp format (net 28 left)
                m_syncTime[detectorID] =
                    (m_syncTime[detectorID] & 0x00000FFFFFFFFFFF) + ((pixdata & 0x00000000FFFF0000) << 28);
                if(!m_clearedHeader[detectorID] && (double)m_syncTime[detectorID] / (4096. * 40000000.) < 6.) {
                    m_clearedHeader[detectorID] = true;
                    LOG(TRACE) << "Cleared header";
                }
            }
        }

        // In data taking during 2015 there was sometimes still data left in the buffers at the start of
        // a run. For that reason we keep skipping data until this "header" data has been cleared, when
        // the heart beat signal starts from a low number (~few seconds max)
        if(!m_clearedHeader[detectorID]) {
            LOG(TRACE) << "Header not cleared, skipping data block.";
            continue;
        }

        // Header 0x06 and 0x07 are the start and stop signals for power pulsing
        if(header == 0x0) {

            // Only want to read these packets from the DUT
            if(detectorID != m_config.get<std::string>("DUT")) {
                continue;
            }

            // Get the second part of the header
            const UChar_t header2 = ((pixdata & 0x0F00000000000000) >> 56) & 0xF;

            // New implementation of power pulsing signals from Adrian
            if(header2 == 0x6) {
                LOG(TRACE) << "Found power pulsing - start";

                const uint64_t time((pixdata & 0x0000000FFFFFFFFF) << 12);

                const uint64_t controlbits = ((pixdata & 0x00F0000000000000) >> 52) & 0xF;

                const uint64_t powerOn = ((controlbits & 0x2) >> 1);
                const uint64_t shutterClosed = ((controlbits & 0x1));

                // Stop looking at data if the signal is after the current event window
                // (and rewind the file
                // reader so that we start with this signal next event)
                if(temporalSplit &&
                   ((double)time / (4096. * 40000000.)) > (clipboard->get_persistent("currentTime") + m_eventLength)) {
                    (*m_file_iterator[detectorID])->seekg(-1 * sizeof(pixdata), std::ios_base::cur);
                    LOG(TRACE) << "Signal has a time beyond the current event: " << (double)time / (4096. * 40000000.);
                    break;
                }

                SpidrSignal* powerSignal = (powerOn ? new SpidrSignal("powerOn", time) : new SpidrSignal("powerOff", time));
                spidrData->push_back(powerSignal);
                LOG(DEBUG) << "Power is " << (powerOn ? "on" : "off") << " power! Time: " << std::setprecision(10)
                           << (double)time / (4096. * 40000000.);

                LOG(TRACE) << "Shutter closed: " << hex << shutterClosed << dec;

                SpidrSignal* shutterSignal =
                    (shutterClosed ? new SpidrSignal("shutterClosed", time) : new SpidrSignal("shutterOpen", time));
                if(!shutterClosed) {
                    spidrData->push_back(shutterSignal);
                    m_shutterOpen = true;
                    LOG(TRACE) << "Have opened shutter with signal " << shutterSignal->type() << " at time "
                               << (double)time / (4096. * 40000000.);
                }

                if(shutterClosed && m_shutterOpen) {
                    spidrData->push_back(shutterSignal);
                    m_shutterOpen = false;
                    LOG(TRACE) << "Have closed shutter with signal " << shutterSignal->type() << " at time "
                               << (double)time / (4096. * 40000000.);
                }

                LOG(DEBUG) << "Shutter is " << (shutterClosed ? "closed" : "open") << ". Time: " << std::setprecision(10)
                           << (double)time / (4096. * 40000000.);
            }

            /*
            // 0x6 is power on
            if(header2 == 0x6){
              const uint64_t time( (pixdata & 0x0000000FFFFFFFFF) << 12 );
              SpidrSignal* signal = new SpidrSignal("powerOn",time);
              spidrData->push_back(signal);
              LOG(DEBUG) <<"Turned on power! Time: "<<(double)time/(4096. *
            40000000.);
            }
            // 0x7 is power off
            if(header2 == 0x7){
              const uint64_t time( (pixdata & 0x0000000FFFFFFFFF) << 12 );
              SpidrSignal* signal = new SpidrSignal("powerOff",time);
              spidrData->push_back(signal);
              LOG(DEBUG) <<"Turned off power! Time: "<<(double)time/(4096. *
            40000000.);
            }
             */
        }

        // Header 0x6 indicate trigger data
        if(header == 0x6) {
            const UChar_t header2 = ((pixdata & 0x0F00000000000000) >> 56) & 0xF;
            if(header2 == 0xF) {
                unsigned long long int stamp = (pixdata & 0x1E0) >> 5;
                long long int timestamp_raw = (pixdata & 0xFFFFFFFFE00) >> 9;
                long long int timestamp = 0;
                int triggerNumber = ((pixdata & 0xFFF00000000000) >> 44);
                int intermediate = (pixdata & 0x1F);
                if(intermediate != 0)
                    continue;

                if((m_syncTimeTDC[detectorID] - timestamp_raw) >
                   0x1312d000) // if jump back in time is larger than 1 sec, overflow detected...
                {
                    m_TDCoverflowCounter[detectorID]++;
                }
                m_syncTimeTDC[detectorID] = timestamp_raw;
                timestamp = timestamp_raw + ((unsigned long long int)(m_TDCoverflowCounter[detectorID]) << 35);

                double triggerTime = (timestamp * 25e-9 + stamp * 25e-9 / 12.) / 8.; // 320 MHz clock
                triggerTime -= m_triggerLatency * 1e-9;
                SpidrSignal* triggerSignal = new SpidrSignal("trigger", triggerTime * (4096. * 40000000.));
                spidrData->push_back(triggerSignal);
            }
        }

        // Header 0xA and 0xB indicate pixel data
        if(header == 0xA || header == 0xB) {
            LOG(TRACE) << "Found pixel data";

            // Decode the pixel information from the relevant bits
            const UShort_t dcol = ((pixdata & 0x0FE0000000000000) >> 52);
            const UShort_t spix = ((pixdata & 0x001F800000000000) >> 45);
            const UShort_t pix = ((pixdata & 0x0000700000000000) >> 44);
            const UShort_t col = (dcol + pix / 4);
            const UShort_t row = (spix + (pix & 0x3));

            // Check if this pixel is masked
            if(detector->masked(col, row)) {
                LOG(DEBUG) << "Detector " << detectorID << ": pixel " << col << "," << row << " masked";
                continue;
            }

            // Get the rest of the data from the pixel
            const UShort_t pixno = col * 256 + row;
            const UInt_t data = ((pixdata & 0x00000FFFFFFF0000) >> 16);
            const unsigned int tot = (data & 0x00003FF0) >> 4;
            const uint64_t spidrTime(pixdata & 0x000000000000FFFF);
            const uint64_t ftoa(data & 0x0000000F);
            const uint64_t toa((data & 0x0FFFC000) >> 14);

            // Calculate the timestamp.
            long long int time =
                (((spidrTime << 18) + (toa << 4) + (15 - ftoa)) << 8) + (m_syncTime[detectorID] & 0xFFFFFC0000000000);
            // LOG(DEBUG) <<"Pixel time "<<(double)time/(4096. * 40000000.);
            // LOG(DEBUG) <<"Sync time "<<(double)m_syncTime[detectorID]/(4096. *
            // 40000000.);

            // Add the timing offset from the coniditions file (if any)
            time += (long long int)(detector->timingOffset() * 4096. * 40000000.);

            // The time from the pixels has a maximum value of ~26 seconds. We compare
            // the pixel time
            // to the "heartbeat" signal (which has an overflow of ~4 years) and check
            // if the pixel
            // time has wrapped back around to 0

            // If the counter overflow happens before reading the new heartbeat
            //      while( abs(m_syncTime[detectorID]-time) > 0x0000020000000000 ){
            if(!extra) {
                while(m_syncTime[detectorID] - time > 0x0000020000000000) {
                    time += 0x0000040000000000;
                }
            } else {
                while(m_prevTime - time > 0x0000020000000000) {
                    time += 0x0000040000000000;
                }
            }

            // If events are loaded based on time intervals, take all hits where the
            // time is within this window

            // Ignore pixels if they arrive before the current event window
            //      if(temporalSplit && ((double)time/(4096. *
            //      40000000.)) < (clipboard->get_persistent("currentTime")) ){
            //        continue;
            //      }

            // Stop looking at data if the pixel is after the current event window
            // (and rewind the file
            // reader so that we start with this pixel next event)
            if(temporalSplit &&
               ((double)time / (4096. * 40000000.)) > (clipboard->get_persistent("currentTime") + m_eventLength)) {
                (*m_file_iterator[detectorID])->seekg(-1 * sizeof(pixdata), std::ios_base::cur);
                break;
            }

            // Otherwise create a new pixel object
            Pixel* pixel = new Pixel(detectorID, row, col, (int)tot, time);
            devicedata->push_back(pixel);
            //      bufferedData[detectorID]->push_back(pixel);
            m_prevTime = time;
        }

        // Stop when we reach some large number of pixels (if events not based on
        // time)
        if(!temporalSplit && devicedata->size() >= m_numberPixelHits) {
            break;
        }
    }

    // Now we have data buffered into the temporary storage. We will sort this by time, and then load
    // the data from one event onto it.

    // If no data was loaded, return false
    if(devicedata->empty()) {
        return false;
    }

    // Count events:
    m_currentEvent++;
    return true;
}

void Timepix3EventLoader::finalise() {}
