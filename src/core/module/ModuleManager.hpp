/** @file
 *  @brief Interface to the core framework
 *  @copyright Copyright (c) 2017 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef CORRYVRECKAN_ANALYSIS_H
#define CORRYVRECKAN_ANALYSIS_H

#include <fstream>
#include <map>
#include <vector>

#include <TBrowser.h>
#include <TDirectory.h>
#include <TFile.h>

#include "Module.hpp"
#include "core/clipboard/Clipboard.hpp"
#include "core/config/ConfigManager.hpp"
#include "core/detector/Detector.hpp"

namespace corryvreckan {

    /**
     * @brief Core class of the Corryvreckan analysis framework
     *
     * The analysis class is the core class which allows the event processing to run. It basically contains a vector of
     * modules, each of which is initialised, run on each event and finalised. It does not define what an event is, merely
     * runs each module sequentially and passes the clipboard between them (erasing it at the end of each run sequence). When
     * an module returns a Failure code, the event processing will stop.
     */
    class ModuleManager {
        using ModuleList = std::list<std::shared_ptr<Module>>;

    public:
        // Constructors and destructors
        explicit ModuleManager(std::string config_file_name, std::vector<std::string> options = std::vector<std::string>());
        virtual ~ModuleManager(){};

        // Member functions
        void load();

        void run();
        void timing();
        void initialiseAll();
        void finaliseAll();

        void terminate();
        void reset(){};

        TBrowser* browser;

    protected:
        // Member variables
        std::shared_ptr<Clipboard> m_clipboard;
        Configuration global_config;
        std::vector<std::shared_ptr<Detector>> m_detectors;

    private:
        void load_detectors();
        void load_modules();
        void add_units();

        /**
         * @brief Get a specific detector, identified by its name
         * @param  name Name of the detector to retrieve
         * @return Pointer to the requested detector, nullptr if detector with given name is not found
         */
        std::shared_ptr<Detector> get_detector(std::string name);

        std::shared_ptr<Detector> m_reference;

        // Create string vector for detector types:
        std::vector<std::string> get_type_vector(char* tokens);
        // Log file if specified
        std::ofstream log_file_;

        std::unique_ptr<TFile> m_histogramFile;
        int m_events;
        int m_tracks;

        /**
         * @brief Create unique modules
         * @param library Void pointer to the loaded library
         * @param config Configuration of the module
         * @return An unique module together with its identifier
         */
        std::pair<ModuleIdentifier, Module*> create_unique_module(void* library, Configuration config);

        /**
         * @brief Create detector modules
         * @param library Void pointer to the loaded library
         * @param config Configuration of the module
         * @param dut_only Bollean signalling whether should be instantiated only for DUT detectors
         * @param types List of detector type restrictions imposed by the module itself
         * @return A list of all created detector modules and their identifiers
         */
        std::vector<std::pair<ModuleIdentifier, Module*>>
        create_detector_modules(void* library, Configuration config, bool dut_only, std::vector<std::string> types);

        using IdentifierToModuleMap = std::map<ModuleIdentifier, ModuleList::iterator>;

        ModuleList m_modules;
        IdentifierToModuleMap id_to_module_;

        std::map<std::string, void*> loaded_libraries_;

        std::atomic<bool> m_terminate;
        std::unique_ptr<corryvreckan::ConfigManager> conf_mgr_;

        std::tuple<LogLevel, LogFormat> set_module_before(const std::string&, const Configuration& config);
        void set_module_after(std::tuple<LogLevel, LogFormat> prev);

        std::map<Module*, long double> module_execution_time_;
    };
} // namespace corryvreckan

#endif // CORRYVRECKAN_ANALYSIS_H