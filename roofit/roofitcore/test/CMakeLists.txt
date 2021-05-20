# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# @author Patrick Bos, NL eScience Center, 2018
# @author Danilo Piparo CERN, 2018

# Only Google Test tests are added here, the other tests are added elsewhere
ROOT_ADD_GTEST(simple simple.cxx LIBRARIES RooFitCore)
ROOT_ADD_GTEST(testRooGradMinimizer RooGradMinimizer.cxx LIBRARIES RooFitCore)
ROOT_ADD_GTEST(testBidirMMapPipe testBidirMMapPipe.cxx LIBRARIES RooFitCore m)
ROOT_ADD_GTEST(testWorkspace testWorkspace.cxx LIBRARIES RooFitCore RooFit RooStats)
if(NOT MSVC OR win_broken_tests)
  ROOT_ADD_GTEST(testRooDataHist testRooDataHist.cxx LIBRARIES RooFitCore)
endif()
ROOT_ADD_GTEST(testRooWrapperPdf testRooWrapperPdf.cxx LIBRARIES Gpad RooFitCore)
ROOT_ADD_GTEST(testGenericPdf testGenericPdf.cxx LIBRARIES RooFitCore)
ROOT_ADD_GTEST(testRooProxy testRooProxy.cxx LIBRARIES RooFitCore WILLFAIL)
ROOT_ADD_GTEST(testRooAbsPdf testRooAbsPdf.cxx LIBRARIES RooFitCore)
ROOT_ADD_GTEST(testRooDataSet testRooDataSet.cxx LIBRARIES Tree RooFitCore)
ROOT_ADD_GTEST(testRooFormula testRooFormula.cxx LIBRARIES RooFitCore)

ROOT_ADD_GTEST(testMPFEnll MultiProcess/MPFEnll.cpp LIBRARIES RooFitCore)

find_package(ZeroMQ REQUIRED)

ROOT_ADD_GTEST(testLikelihoodGradientJob TestStatistics/LikelihoodGradientJob.cpp LIBRARIES RooFitMultiProcess RooFitCore RooFit RooStats m ${ZeroMQ_LIBRARY})
ROOT_ADD_GTEST(testRooRealL TestStatistics/RooRealL.cpp LIBRARIES RooFitCore RooFit)
