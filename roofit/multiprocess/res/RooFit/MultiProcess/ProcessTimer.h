#ifndef ROOT_ROOFIT_MultiProcess_ProcessTimer
#define ROOT_ROOFIT_MultiProcess_ProcessTimer

#include <chrono>
#include <string>
#include <sys/types.h> // pid_t
#include <map>
#include <list>
#include <tuple>

using namespace std;

class ProcessTimer{

    // Map of a list of timepoints, even timepoints are start times, uneven timepoints are end times
    static map<string, list<chrono::time_point<chrono::steady_clock>>> durations;
    static chrono::time_point<chrono::steady_clock> begin;
    static pid_t process;

public:

    static void setup(pid_t proc) { ProcessTimer::process = proc; ProcessTimer::begin = chrono::steady_clock::now();};

    static pid_t get_process() {return ProcessTimer::process; };

    static void start_timer(string section_name);

    static void end_timer(string section_name, bool write_now = false);

    static void print_durations(string to_print = "all");

    static void print_timestamps();

    static void write_file();

};

#endif // ROOT_ROOFIT_MultiProcess_ProcessTimer
