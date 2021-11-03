############################################################################
# CMakeLists.txt file for building ROOT RooFitMultiProcess package
# @author Patrick Bos, Netherlands eScience Center
############################################################################

ROOT_LINKER_LIBRARY(RooFitMultiProcess
        src/worker.cxx
        src/Messenger.cxx
        src/ProcessManager.cxx
        src/util.cxx
        src/Queue.cxx
        src/JobManager.cxx
        src/Job.cxx
    DEPENDENCIES
        RooFitZMQ
)

target_link_libraries(RooFitMultiProcess PUBLIC RooFitZMQ)

ROOT_ADD_TEST_SUBDIRECTORY(test)

ROOT_INSTALL_HEADERS()
