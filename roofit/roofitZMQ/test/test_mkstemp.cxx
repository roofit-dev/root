/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2016-2021, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */

#include "RooFit_ZMQ/ZeroMQSvc.h"

#include "gtest/gtest.h"

#include <cstdlib>  // mkstemp
#include <unistd.h> // for "legacy" systems

TEST(BindToTmpFile, mkstemp)
{
   char filename[] = "/tmp/roofitMP_XXXXXX";
   while (mkstemp(filename) >= 0) {}
   auto socket = zmqSvc().socket(zmq::socket_type::push);
   char address [50];
   sprintf(address, "ipc://%s", filename);
   try {
      socket.bind(address);
   } catch (const zmq::error_t &) {
      printf("caught an exception\n");
   }
}