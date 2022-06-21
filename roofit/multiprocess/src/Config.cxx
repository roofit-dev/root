/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2021, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#include "RooFit/MultiProcess/Config.h"
#include "RooFit/MultiProcess/JobManager.h"
#include "RooFit/MultiProcess/ProcessTimer.h"

#include <thread> // std::thread::hardware_concurrency()
#include <cstdio>

namespace RooFit {
namespace MultiProcess {

/** \class Config
 *
 * \brief Configuration for MultiProcess infrastructure
 *
 * This class offers user-accessible configuration of the MultiProcess infrastructure.
 * Since the rest of the MultiProcess classes are only accessible at compile time, a
 * separate class is needed to set configuration. Currently, the only configurable
 * part is the number of workers to be deployed.
 *
 * The default number of workers is set using 'std::thread::hardware_concurrency()'.
 * To change it, use 'Config::setDefaultNWorkers()' to set it to a different value
 * before creation of a new JobManager instance. Note that it cannot be set to zero
 * and also cannot be changed after JobManager has been instantiated.
 *
 * Use Config::getDefaultNWorkers() to access the current value.
 */

void Config::setDefaultNWorkers(unsigned int N_workers)
{
   if (JobManager::is_instantiated()) {
      printf("Warning: Config::setDefaultNWorkers cannot set number of workers after JobManager has been instantiated!\n");
   } else if (N_workers == 0) {
      printf("Warning: Config::setDefaultNWorkers cannot set number of workers to zero.\n");
   } else {
      defaultNWorkers_ = N_workers;
   }
}

unsigned int Config::getDefaultNWorkers()
{
   return defaultNWorkers_;
}

/// Set a write interval for the timer, 0 (default) means times are not 
/// written intermittently any other (positive) value means times are written
/// every `write_interval` secondss
void Config::setTimerWriteInterval(int write_interval)
{
  ProcessTimer::set_write_interval(write_interval);
}

// initialize static member
unsigned int Config::defaultNWorkers_ = std::thread::hardware_concurrency();

} // namespace MultiProcess
} // namespace RooFit
