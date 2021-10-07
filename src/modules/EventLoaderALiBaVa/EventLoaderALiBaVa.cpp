/**
 * @file
 * @brief Implementation of module EventLoaderALiBaVa
 *
 * @copyright Copyright (c) 2020 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <dirent.h>

#include "objects/Pixel.hpp"
#include "ALiBaVa/auxfunctions.h"
#include "EventLoaderALiBaVa.h"
#include "ALiBaVa/HDFRoot.h"

using namespace corryvreckan;

EventLoaderALiBaVa::EventLoaderALiBaVa(Configuration& config, std::shared_ptr<Detector> detector)
    : Module(config, detector) {
      m_detector = detector;
    }

void EventLoaderALiBaVa::initialize() {

  // Take input directory, run number,
  // lower timecut, upper timecut,        --- in nanoseconds
  // first X ALiBaVa events to ignore, channel range,
  // from global parameters in the config.
  // Default values are set first

  config_.setDefault<int>("run", 0);
  config_.setDefault<double>("LowerTimecut", 0);
  config_.setDefault<double>("UpperTimecut", std::numeric_limits<int>::max());
  config_.setDefault<int>("IgnoreEvents", 1);
  config_.setDefault<int>("LowerChannel", 0);
  config_.setDefault<int>("UpperChannel", 255);
  config_.setDefault<double>("Chargecut", std::numeric_limits<int>::max());
  config_.setDefault<bool>("CrosstalkCorrection", true);
  config_.setDefault<double>("CalibrationConstant", 1.0);

  m_inputDirectory = config_.getPath("input_directory");
  m_run = config_.get<int>("run");
  m_timecut_lower = config_.get<int>("LowerTimecut");
  m_timecut_upper = config_.get<int>("UpperTimecut");
  m_ignore_events = config_.get<int>("IgnoreEvents");
  m_lower_channel = config_.get<int>("LowerChannel");
  m_upper_channel = config_.get<int>("UpperChannel");
  m_chargecut = config_.get<int>("Chargecut");
  m_correct_crosstalk = config_.get<bool>("CrosstalkCorrection");
  m_calibration_constant = config_.get<double>("CalibrationConstant");

  // Open the input directory
  DIR* directory = opendir(m_inputDirectory.c_str());
  if(directory == nullptr) {
      LOG(ERROR) << "Directory \'" << m_inputDirectory << "\' does not exist or was not supplied.";
      return;
  }

  // Read the run-files (data, pedestal and calibration) in the folder
  dirent* entry;
  while(entry = readdir(directory)) {
      if(entry->d_type == DT_REG) {
          std::string entryName = entry->d_name;
          if(entryName.find(std::to_string(m_run)+".dat") != std::string::npos || entryName.find("dat_run"+std::to_string(m_run)+".hdf") != std::string::npos){
            m_datafilename = m_inputDirectory + "/" + entryName;
          }
          if(entryName.find(std::to_string(m_run)+".ped") != std::string::npos ||entryName.find("ped_run"+std::to_string(m_run)+".hdf") != std::string::npos){
            m_pedestalfilename = m_inputDirectory + "/" + entryName;
          }
          if(entryName.find(std::to_string(m_run)+".cal") != std::string::npos || entryName.find("cal_run"+std::to_string(m_run)+".hdf") != std::string::npos){
            m_calibrationfilename = m_inputDirectory + "/" + entryName;
          }
      }
  }

  // Log errors in case the files aren't found in the folder.
  // The datafile can also be supplied directly in the config.
  if(m_datafilename.length() == 0) {
      LOG(ERROR) << "No data file was found for ALiBaVa in " << m_inputDirectory;
        return;
  }
  if(m_pedestalfilename.length() == 0) {
      LOG(WARNING) << "No pedestal file was found." << "\n" << "Datafile will be used for pedestal";
  }

  if(m_calibrationfilename.length() == 0) {
      LOG(WARNING) << "No calibration file was found." << "\n" << "Results will be uncalibrated: ADC = charge.";
  }

  // Create a pointer with the data file.
  ALiBaVaPointer = DataFileRoot::OpenFile(m_datafilename.c_str());
  // Call the loader function, which sets everything up.
  // For some reason in ALiBaVa's code, it needs the pointer AND the data file again
  ALiBaVa_loader(ALiBaVaPointer, m_datafilename.c_str(), m_pedestalfilename.c_str(), m_calibrationfilename.c_str(),
      pedestalValues, noiseValues, correctedPedestalValues, correctedNoiseValues, m_lower_channel, m_upper_channel);
  // Set the timecuts
  ALiBaVaPointer->set_timecut(m_timecut_lower, m_timecut_upper);
  // Check how many events the data file contains.
  //nEvents = ALiBaVaPointer->nevents();

// Ignore the first X events to ensure synchronisation, default is X = 1 which ignores the first event.
for(int ievt = 0; ievt < m_ignore_events; ievt++){
    ALiBaVaPointer->read_event();
  }

  // Create histograms
  chargeHist = new TH1F("charge","charge", 50, 0, m_chargecut);
}

StatusCode EventLoaderALiBaVa::run(const std::shared_ptr<Clipboard>& clipboard) {
  // During running, every Corryvreckan event will get one ALiBaVa event
  // Increment the event counter
  iEvent++;
  // Create the pixelvector
  PixelVector pixels;
  // Get the event that already exists on the Corryvreckan clipboard
  auto event = clipboard->getEvent();

  // Read a data event from the ALiBaVa data file
  // Give feedback according to return code
  int return_code = ALiBaVaPointer->read_event();
  if(return_code == -1){
    // LOG(WARNING) << "End of data file reached.";
    return StatusCode::EndRun;
    // Not sure if this is end of run or something else. Need to see difference between HDF5 and binary
  }
  else if(return_code == 0){
    LOG(ERROR) << "There\'s something wrong (0) with the datafile at event number " << iEvent-1;
    LOG(ERROR) << "Terminating run";
    return StatusCode::EndRun;
  }
  else if(return_code == 1){
    // This means the event was read properly.
  }
  else if(return_code == 2 ){
    LOG(ERROR) << "There\'s something wrong (2) with the datafile at event number " << iEvent-1;
    LOG(ERROR) << "Terminating run";
    return StatusCode::EndRun;
  }
  else if(return_code == 3){
    // LOG(ERROR) << "There\'s something wrong (3) with the datafile at event number " << iEvent-1;
    // LOG(ERROR) << "Terminating run";
    // return StatusCode::EndRun;
    // Still need to figure this out... I think it's just end of run. Need to see difference between HDF5 and binary
    // LOG(WARNING) << "End of data file reached.";
    return StatusCode::EndRun;
  }
  else if(return_code == 4){
    LOG(ERROR) << "There\'s something wrong (4) with the HDF5 datafile at event number " << iEvent-1;
    LOG(ERROR) << "Terminating run";
    return StatusCode::EndRun;
  }
  else{
    LOG(ERROR) << "This shouldn\'t happen. Current event number: " << iEvent-1;
    LOG(ERROR) << "Terminating run";
    return StatusCode::EndRun;
  }

  // Process the opened data event, i.e. pedestal correction, common mode noise correction
  ALiBaVaPointer->process_event();
  // This gets the TDC time from the event, allowing timecuts around the event peak
  // The timecut is set in the ALiBaVa_loader() function.
  double TDCTime = ALiBaVaPointer->time();
  if(!ALiBaVaPointer->valid_time(TDCTime)){
    clipboard->putData(pixels, m_detector->getName());
    return StatusCode::NoData;
  }

//// THIS STUFF NEEDS TO MOVE TO SEPARATE MODULES (?)



  if(false){//crosstalk AFTER clustering and calibration
    double channels_CalSig[m_upper_channel-m_lower_channel+1];
    double channels_CalSig_corrected[m_upper_channel-m_lower_channel+1];
    double b_one;
    double b_two;
    if(m_correct_crosstalk){
      for(int chan = m_lower_channel; chan <= m_upper_channel; chan++){
        channels_CalSig[chan-m_lower_channel] = ALiBaVaPointer->get_gain(chan);
      }
      channels_CalSig_corrected[0] = (1+b_one+b_two)*channels_CalSig[0];
      channels_CalSig_corrected[1] = (1+b_one+b_two)*channels_CalSig[1]-b_one*channels_CalSig[0];
      for(int chan = m_lower_channel+2; chan <= m_upper_channel; chan++){
        channels_CalSig_corrected[chan-m_lower_channel] = (1+b_one+b_two)*channels_CalSig[chan-m_lower_channel]-b_one*channels_CalSig[chan-m_lower_channel-1]-b_two*channels_CalSig[chan-m_lower_channel-2];
      }
    }
  }

////

  // This loops over the channels in the current ALiBaVa event
  for(int chan = m_lower_channel; chan <= m_upper_channel; chan++){
    // In order, these are the calibration factor, the calibrated signal, and the ADC signal
    // If a pedestal file is supplied, pedestal and common mode error correction is applied to all
    // double calibration = ALiBaVaPointer->get_gain(chan);
    // double CalSignal = ALiBaVaPointer->signal(chan);
    // double ADCSignal = ALiBaVaPointer->signal(chan)*calibration;
    // double SNRatio = ALiBaVaPointer->sn(chan);
    double CalSignal = ALiBaVaPointer->ADC_signal(chan)*m_calibration_constant;
    double SNRatio = ALiBaVaPointer->sn(chan);

    // The chargecut is applied here
    // std::cout << CalSignal << " smaller than " << m_chargecut << "?\n";
    if(CalSignal < m_chargecut && CalSignal > 0){
    // if(true){
      // Create a pixel for every channel in this event with all the information and put it in the vector.
      // The value in the pixel reserved for the ADC value is used for the S/N ratio multiplied by 100000.
      auto pixel = std::make_shared<Pixel>(m_detector->getName(), chan, 0, SNRatio*100000, CalSignal, 0);
      pixels.push_back(pixel);

      // Fill the histograms
      chargeHist->Fill(CalSignal);
    }
  }

  // Put the created vector of pixels on the clipboard.
  clipboard->putData(pixels, m_detector->getName());

  // If the pixels vector is empty, report this to Corryvreckan
  if(pixels.empty()){
    return StatusCode::NoData;
  }

  // Report the end of this event to Corryvreckan
  return StatusCode::Success;
}


void EventLoaderALiBaVa::finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) {
    ALiBaVaPointer->close();
    delete ALiBaVaPointer;
    LOG(DEBUG) << "Analysed " << iEvent; //<< " of total " << nEvents << " ALiBaVa events";
}
