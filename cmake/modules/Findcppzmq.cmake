if (NOT ZeroMQ_FOUND)
    message(FATAL_ERROR "Search for libzmq first!")
endif ()

find_path(CPPZMQ_INCLUDE_DIRS "zmq.hpp"
        HINTS "${ZeroMQ_INCLUDE_DIR}"
        PATHS "${CPPZMQ_DIR}" "$ENV{CPPZMQ_DIR}"
        PATH_SUFFIXES "include" "cppzmq"
        )
mark_as_advanced (CPPZMQ_INCLUDE_DIRS)

if (CPPZMQ_INCLUDE_DIRS)
    add_library (cppzmq INTERFACE IMPORTED)
    set_target_properties(cppzmq PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${CPPZMQ_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES libzmq
            )
    set (CPPZMQ_LIBRARIES cppzmq)
endif ()

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (CPPZMQ DEFAULT_MSG CPPZMQ_INCLUDE_DIRS)