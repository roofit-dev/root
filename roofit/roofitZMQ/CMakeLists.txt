############################################################################
# CMakeLists.txt file for building ROOT roofitcore/ZMQ package
# @author Patrick Bos, Netherlands eScience Center
############################################################################

ROOT_STANDARD_LIBRARY_PACKAGE(RooFitZMQ
    HEADERS
        RooFit_ZMQ/Utility.h
        RooFit_ZMQ/ZeroMQSvc.h
        RooFit_ZMQ/functions.h
        RooFit_ZMQ/ZeroMQPoller.h
        RooFit_ZMQ/ppoll.h
    SOURCES
        src/ZeroMQSvc.cpp
        src/ZeroMQPoller.cpp
        src/functions.cpp
        src/ppoll.cpp
    DICTIONARY_OPTIONS
        -writeEmptyRootPCM
    DEPENDENCIES
        Core # TVersionCheck, RegisterModule, these are always called from G__RooFitZMQ.cxx.o
    )

target_link_libraries(RooFitZMQ PUBLIC ZeroMQ cppzmq)
ROOT_ADD_TEST_SUBDIRECTORY(test)
