/**
 * @file
 * @brief Implementation of ROOT data file writer module
 * @copyright Copyright (c) 2019 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * @remarks The implementation of this module is based on the ROOTObjectWriter module of the Allpix Squared project
 */

#include "FileWriter.h"

#include <fstream>
#include <string>
#include <utility>

#include <TBranchElement.h>
#include <TClass.h>

#include "core/utils/file.h"
#include "core/utils/log.h"
#include "core/utils/type.h"

#include "objects/Object.hpp"

using namespace corryvreckan;

FileWriter::FileWriter(Configuration config, std::vector<std::shared_ptr<Detector>> detectors)
    : Module(std::move(config), std::move(detectors)) {}
/**
 * @note Objects cannot be stored in smart pointers due to internal ROOT logic
 */
FileWriter::~FileWriter() {
    // Delete all object pointers
    for(auto& index_data : write_list_) {
        delete index_data.second;
    }
}

void FileWriter::initialise() {

    // Output format
    auto fileFormats = {
        "root",
        "txt",
        "json",
    }; // known output formats - place here or in header?
    output_format_ = m_config.get<std::string>("format", "root");
    if((std::find(std::begin(fileFormats), std::end(fileFormats), output_format_) == std::end(fileFormats))) {
        throw InvalidValueError(m_config, "format", "unknown file format specified");
    }

    // Create output file
    output_file_name_ = createOutputFile(
        corryvreckan::add_file_extension(m_config.get<std::string>("file_name", "data"), output_format_), true);
    if(output_format_ == "root") {
        output_TFile_ = std::make_unique<TFile>(output_file_name_.c_str(), "RECREATE");
        if(output_TFile_->IsZombie()) {
            ModuleError("Creating output file failed, cannot continue");
        }
        output_TFile_->cd();

        // Create event tree:
        event_tree_ = std::make_unique<TTree>("Event", (std::string("Tree of Events").c_str()));
        event_tree_->Bronch("global", "corryvreckan::Event", &event_);
    } else if(output_format_ == "txt") {
        output_file_ = std::make_unique<std::ofstream>(output_file_name_);
        *output_file_ << "# Corryvreckan ASCII data" << std::endl << std::endl;
    } else if(output_format_ == "json") {
        output_file_ = std::make_unique<std::ofstream>(output_file_name_);
        // print first line to create JSON-Array of objects
        *output_file_ << "[" << std::endl;
    } else {
        throw InvalidValueError(m_config, "format", "specified format is not implemented");
    }

    // Read include and exclude list
    if(m_config.has("include") && m_config.has("exclude")) {
        throw InvalidValueError(m_config, "exclude", "include and exclude parameter are mutually exclusive");
    } else if(m_config.has("include")) {
        auto inc_arr = m_config.getArray<std::string>("include");
        include_.insert(inc_arr.begin(), inc_arr.end());
    } else if(m_config.has("exclude")) {
        auto exc_arr = m_config.getArray<std::string>("exclude");
        exclude_.insert(exc_arr.begin(), exc_arr.end());
    }
    m_eventNumber = 0;
}

StatusCode FileWriter::run(std::shared_ptr<Clipboard> clipboard) {

    if(!clipboard->isEventDefined()) {
        ModuleError("No Clipboard event defined, cannot continue");
    }

    if(output_format_ == "root") {
        // Read event from clipboard and write to tree:
        event_ = clipboard->getEvent().get();
        event_tree_->Fill();
        write_cnt_++;
    } else if(output_format_ == "txt") {
        // Print the current event:
        *output_file_ << "=== " << m_eventNumber << " ===" << std::endl;
        *output_file_ << *clipboard->getEvent() << std::endl;
    }

    auto data = clipboard->getAll();
    LOG(DEBUG) << "Clipboard has " << data.size() << " different object types.";

    for(auto& block : data) {
        try {
            auto type_idx = block.first;
            auto class_name = corryvreckan::demangle(type_idx.name());
            auto class_name_full = corryvreckan::demangle(type_idx.name(), true);
            LOG(TRACE) << "Received objects of type \"" << class_name << "\" in " << block.second.size() << " blocks";

            // Check if these objects should be stored
            if((!include_.empty() && include_.find(class_name) == include_.end()) ||
               (!exclude_.empty() && exclude_.find(class_name) != exclude_.end())) {
                LOG(TRACE) << "Ignoring object " << corryvreckan::demangle(type_idx.name())
                           << " because it has been excluded or not explicitly included";
                continue;
            }

            for(auto& detector_block : block.second) {
                // Get the detector name
                std::string detector_name;
                if(!detector_block.first.empty()) {
                    detector_name = detector_block.first;
                }

                if(output_format_ == "root") {
                    auto objects = std::static_pointer_cast<ObjectVector>(detector_block.second);
                    LOG(TRACE) << " - " << detector_name << ": " << objects->size();

                    // Create a new branch of the correct type if this object has not been received before
                    auto index_tuple = std::make_tuple(type_idx, detector_name);
                    if(write_list_.find(index_tuple) == write_list_.end()) {

                        // Add vector of objects to write to the write list
                        write_list_[index_tuple] = new std::vector<Object*>();
                        auto addr = &write_list_[index_tuple];

                        auto new_tree = (trees_.find(class_name) == trees_.end());
                        if(new_tree) {
                            // Create new tree
                            output_TFile_->cd();
                            trees_.emplace(
                                class_name,
                                std::make_unique<TTree>(class_name.c_str(), (std::string("Tree of ") + class_name).c_str()));
                        }

                        std::string branch_name = detector_name.empty() ? "global" : detector_name;

                        trees_[class_name]->Bronch(
                            branch_name.c_str(), (std::string("std::vector<") + class_name_full + "*>").c_str(), addr);

                        if(new_tree) {
                            LOG(DEBUG) << "Pre-filling new tree of " << class_name << " with " << m_eventNumber
                                       << " empty events";
                            for(unsigned int i = 0; i < m_eventNumber; ++i) {
                                trees_[class_name]->Fill();
                            }
                        } else {
                            LOG(DEBUG) << "Pre-filling new branch " << branch_name << " of " << class_name << " with "
                                       << m_eventNumber << " empty events";
                            auto* branch = trees_[class_name]->GetBranch(branch_name.c_str());
                            for(unsigned int i = 0; i < m_eventNumber; ++i) {
                                branch->Fill();
                            }
                        }
                    }

                    // Fill the branch vector
                    for(auto& object : *objects) {
                        ++write_cnt_;
                        write_list_[index_tuple]->push_back(object);
                    }
                } else if(output_format_ == "txt") {
                    *output_file_ << "--- " << (detector_name.empty() ? "<global>" : detector_name) << " ---" << std::endl;
                    auto objects = std::static_pointer_cast<ObjectVector>(detector_block.second);
                    for(auto& object : *objects) {
                        *output_file_ << *object << std::endl;
                    }
                } else if(output_format_ == "json") {
                    auto objects = std::static_pointer_cast<ObjectVector>(detector_block.second);
                    for(auto& object : *objects) {
                        *output_file_ << TBufferJSON::ToJSON(object) << "," << std::endl;
                    }
                }
            }
        } catch(...) {
            LOG(WARNING) << "Cannot process object of type" << corryvreckan::demangle(block.first.name());
            return StatusCode::NoData;
        }
    }

    m_eventNumber++;

    if(output_format_ == "root") {
        LOG(TRACE) << "Writing new objects to tree";
        output_TFile_->cd();

        // Fill the tree with the current received messages
        for(auto& tree : trees_) {
            tree.second->Fill();
        }

        // Clear the current message list
        for(auto& index_data : write_list_) {
            index_data.second->clear();
        }
    }

    return StatusCode::Success;
}

void FileWriter::finalise() {
    if(output_format_ == "root") {
        LOG(TRACE) << "Writing objects to file";
        output_TFile_->cd();

        int branch_count = 0;
        for(auto& tree : trees_) {
            // Update statistics
            branch_count += tree.second->GetListOfBranches()->GetEntries();
        }
        branch_count += event_tree_->GetListOfBranches()->GetEntries();

        // Finish writing to output file
        output_TFile_->Write();

        // Print statistics
        LOG(STATUS) << "Wrote " << write_cnt_ << " objects to " << branch_count << " branches in file:" << std::endl
                    << output_file_name_;

    } else if(output_format_ == "txt") {
        LOG(STATUS) << "Wrote " << m_eventNumber - 1 << " events in file " << output_file_name_;
    } else if(output_format_ == "json") {
        // finalize the JSON Array, add one empty Object to satisfy JSON Rules
        *output_file_ << "{" << std::endl << "}" << std::endl << "]";
        LOG(STATUS) << "Wrote " << m_eventNumber - 1 << " events in file " << output_file_name_;
    }
}
