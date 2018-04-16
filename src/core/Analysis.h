#ifndef ANALYSIS_H
#define ANALYSIS_H 1

#include <fstream>
#include <map>
#include <vector>

#include "TBrowser.h"
#include "TDirectory.h"
#include "TFile.h"

#include "Detector.h"
#include "clipboard/Clipboard.hpp"
#include "config/ConfigManager.hpp"
#include "config/Configuration.hpp"
#include "module/Module.hpp"

//-------------------------------------------------------------------------------
// The analysis class is the core class which allows the event processing to
// run. It basically contains a vector of modules, each of which is
// initialised,
// run on each event and finalised. It does not define what an event is, merely
// runs each module sequentially and passes the clipboard between them
// (erasing
// it at the end of each run sequence). When an module returns a 0, the event
// processing will stop.
//-------------------------------------------------------------------------------

namespace corryvreckan {
    class Analysis {

    public:
        // Constructors and destructors
        explicit Analysis(std::string config_file_name, std::vector<std::string> options = std::vector<std::string>());
        virtual ~Analysis(){};

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
        Clipboard* m_clipboard;
        Configuration global_config;
        std::vector<Detector*> detectors;

    private:
        void load_detectors();
        void load_modules();
        void add_units();

        // Log file if specified
        std::ofstream log_file_;

        TFile* m_histogramFile;
        TDirectory* m_directory;
        int m_events;
        int m_tracks;

        std::vector<Module*> m_modules;
        std::map<std::string, void*> loaded_libraries_;

        std::atomic<bool> m_terminate;
        std::unique_ptr<corryvreckan::ConfigManager> conf_mgr_;

        Module* create_module(void* library, corryvreckan::Configuration config);
        std::tuple<LogLevel, LogFormat> set_module_before(const std::string&, const Configuration& config);
        void set_module_after(std::tuple<LogLevel, LogFormat> prev);
    };
} // namespace corryvreckan
#endif // ANALYSIS_H
