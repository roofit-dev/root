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
    LIBRARIES
        ZeroMQ
        cppzmq
    )

target_link_libraries(RooFitZMQ PRIVATE ZeroMQ)
target_include_directories(RooFitZMQ PUBLIC $<BUILD_INTERFACE:${ZeroMQ_INCLUDE_DIR}>)
target_link_libraries(RooFitZMQ PRIVATE cppzmq)
target_include_directories(RooFitZMQ PRIVATE $<BUILD_INTERFACE:${cppzmq_INCLUDE_DIR}>)
if (NOT MSVC AND NOT APPLE)
    target_link_options(RooFitZMQ PRIVATE "LINKER:-Bsymbolic")
endif()
ROOT_ADD_TEST_SUBDIRECTORY(test)
