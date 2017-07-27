#ifndef ROO_TIMER
#define ROO_TIMER

#include <chrono>
#include <ctime>
#include <string>
#include <vector>
#include <map>
#include "RooJsonListFile.h"
#include "Rtypes.h"

namespace RooFit {
  long long nanoseconds_from_timespec(timespec time_point);
  double diff_seconds_timespecs(timespec begin, timespec end);
}

class RooTimer {
public:
  virtual void start() = 0;
  virtual void stop() = 0;
  double timing_s();
  void set_timing_s(double timing_s);
  void store_object_timing(const std::string &name, const std::string &category);

//  static std::map<std::string,double> objectTiming;
  static std::map< std::string, std::map< std::string, double> > objectTiming;

  static std::vector<RooJsonListFile> timing_outfiles;
  static std::vector<double> timings;

  static int timing_flag();
  static void set_timing_flag(int flag);

  static Bool_t time_numInts();
  static void set_time_numInts(Bool_t flag);

  static Bool_t time_evaluate_partition();
  static void set_time_evaluate_partition(Bool_t flag);

  static Bool_t time_MPFE_forks();
  static void set_time_MPFE_forks(Bool_t flag);

private:
  double _timing_s;

  static int _timing_flag;
  static Bool_t _time_numInts;
  static Bool_t _time_evaluate_partition;
  static Bool_t _time_MPFE_forks;
};

class RooWallTimer: public RooTimer {
public:
  RooWallTimer();
  virtual void start();
  virtual void stop();

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> _timing_begin, _timing_end;
};

/// @class RooCPUTimer
/// Measures the CPU time on the local process. Note that for multi-process runs,
/// e.g. when using RooRealMPFE, the child process CPU times are not included!
/// Use a separate timer in child processes to measure their CPU timing.
class RooCPUTimer: public RooTimer {
public:
  RooCPUTimer();
  virtual void start();
  virtual void stop();

private:
  struct timespec _timing_begin, _timing_end;
};

#endif
