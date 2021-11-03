############################################################################
# CMakeLists.txt file for building ROOT roofitcore/ZMQ package
# @author Patrick Bos, Netherlands eScience Center
############################################################################

ROOT_LINKER_LIBRARY(RooFitZMQ
        src/ZeroMQSvc.cpp
        src/ZeroMQPoller.cpp
        src/functions.cpp
        src/ppoll.cpp
    )

target_link_libraries(RooFitZMQ PUBLIC libzmq cppzmq)
ROOT_ADD_TEST_SUBDIRECTORY(test)

ROOT_INSTALL_HEADERS()
