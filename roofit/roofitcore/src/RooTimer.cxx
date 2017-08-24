#include "RooTimer.h"
#include <iostream>

double RooTimer::timing_s() {
  return _timing_s;
}

void RooTimer::set_timing_s(double timing_s) {
  _timing_s = timing_s;
}

void RooTimer::store_object_timing(const std::string &name, const std::string &category) {
  RooTimer::objectTiming[category][name] = _timing_s;  // subscript operator overwrites existing values, insert does not
}

std::vector<RooJsonListFile> RooTimer::timing_outfiles;
std::vector<double> RooTimer::timings;


RooWallTimer::RooWallTimer() {
  start();
}

void RooWallTimer::start() {
  _timing_begin = std::chrono::high_resolution_clock::now();
}

void RooWallTimer::stop() {
  _timing_end = std::chrono::high_resolution_clock::now();
  set_timing_s(std::chrono::duration_cast<std::chrono::nanoseconds>(_timing_end - _timing_begin).count() / 1.e9);
}


RooCPUTimer::RooCPUTimer() {
  start();
}

void RooCPUTimer::start() {
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &_timing_begin);
}

namespace RooFit {
  long long nanoseconds_from_timespec(timespec time_point) {
    long long timing_nsec = time_point.tv_nsec + 1000000000 * time_point.tv_sec; // 1'000'000'000 in c++14
    return timing_nsec;
  }

  double diff_seconds_timespecs(timespec begin, timespec end) {
    long long timing_end_nsec = nanoseconds_from_timespec(end);
    long long timing_begin_nsec = nanoseconds_from_timespec(begin);
    return (timing_end_nsec - timing_begin_nsec) / 1.e9;
  }
}

void RooCPUTimer::stop() {
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &_timing_end);
  set_timing_s(RooFit::diff_seconds_timespecs(_timing_begin, _timing_end));
}


////////////////////////////////////////////////////////////////////////////////
/// static map used to store timings of named objects
std::map< std::string, std::map< std::string, double> > RooTimer::objectTiming;

////////////////////////////////////////////////////////////////////////////////
/// Flag to switch on timers used for performance benchmarks. Use at your own peril!
///
/// 1: not used in code, used in scripts to time the entire minimization
/// 2: timing_RATS_evaluate_full
/// 3: timing_RATS_evaluate_mpmaster_perCPU
/// 4: timing_RRMPFE_evaluate_full
/// 5: timing_wall_RRMPFE_evaluate_client
/// 6: timing_cpu_RRMPFE_evaluate_client
/// 7: timing_RRMPFE_calculate_initialize
/// 8: timing_RRMPFE_serverloop_p
/// 9: timing_RRMPFE_serverloop_while_p
/// 10: time communication overhead
int RooTimer::_timing_flag = 0;

Bool_t RooTimer::time_numInts() {
  return _time_numInts;
}

int RooTimer::timing_flag() {
  return _timing_flag;
}

void RooTimer::set_timing_flag(int flag) {
  _timing_flag = flag;
}

void RooTimer::set_time_numInts(Bool_t flag) {
  _time_numInts = flag;
}

Bool_t RooTimer::time_evaluate_partition() {
  return _time_evaluate_partition;
}

void RooTimer::set_time_evaluate_partition(Bool_t flag) {
  std::cout << "WARNING: RooTimer::set_time_evaluate_partition() is best set before using RooRealMPFE to fork to multiple processes. When resetting the flag value after forking, the new setting is not synchronized to other processes." << std::endl;
  _time_evaluate_partition = flag;
}

Bool_t RooTimer::time_MPFE_forks() {
  return _time_MPFE_forks;
}

void RooTimer::set_time_MPFE_forks(Bool_t flag) {
  std::cout << "WARNING: RooTimer::set_time_MPFE_forks() is best set before using RooRealMPFE to fork to multiple processes. When resetting the flag value after forking, the new setting is not synchronized to other processes." << std::endl;
  _time_MPFE_forks = flag;
}

Bool_t RooTimer::_time_numInts = kFALSE;
Bool_t RooTimer::_time_evaluate_partition = kFALSE;
Bool_t RooTimer::_time_MPFE_forks = kFALSE;

// Output streamer; stops time, streams the timing_s value and restarts the timer afterwards
std::ostream& operator<<(std::ostream& out, RooTimer &timer){
  timer.stop();
  return out << timer.timing_s();
  timer.start();
}
