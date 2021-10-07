include(ExternalProject)

set(ZeroMQ_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/BUILTIN_ZeroMQ-prefix)
set(ZeroMQ_LIBNAME ${CMAKE_STATIC_LIBRARY_PREFIX}zmq${CMAKE_STATIC_LIBRARY_SUFFIX})

set(ZeroMQ_LIBRARY ${ZeroMQ_PREFIX}/lib/${ZeroMQ_LIBNAME} CACHE INTERNAL "" FORCE)
set(ZeroMQ_LIBRARIES ${ZeroMQ_LIBRARY} CACHE INTERNAL "" FORCE)

if(NOT WIN32)
    # 2021-10-07: libzmq's minimum CMake version does not yet support CXX_VISIBILITY_PRESET, so we need to pass the flag manually
    set(cxx_visibility_flag "-fvisibility=hidden")
else()
    set(cxx_visibility_flag "")
endif()

ExternalProject_Add(BUILTIN_ZeroMQ
    #    URL https://github.com/zeromq/libzmq/archive/5d8d857540323e2d85c64a7edde1ad5280cad04b.tar.gz
    #    URL_HASH SHA256=750bc1c67a3740a5438da8762931983cb08dc649e1a658cfc34de91c2cf9222c
    GIT_REPOSITORY https://github.com/egpbos/libzmq.git
    GIT_SHALLOW TRUE
    GIT_TAG origin/disable_export
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DWITH_PERF_TOOL=OFF
        -DZMQ_BUILD_TESTS=OFF
        -DENABLE_DRAFTS=ON
        -DENABLE_NO_EXPORT=ON
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        # TODO: remove -w flag when updating to v4.3.5 or higher, current warnings should be gone by then
        -DCMAKE_CXX_FLAGS=${cxx_visibility_flag}\ ${CMAKE_CXX_FLAGS}\ -w
    BUILD_BYPRODUCTS ${ZeroMQ_LIBRARIES}
    )

set(ZeroMQ_FOUND TRUE CACHE BOOL "" FORCE)

set(ZeroMQ_INCLUDE_DIR  ${ZeroMQ_PREFIX}/include CACHE INTERNAL "" FORCE)
set(ZeroMQ_INCLUDE_DIRS ${ZeroMQ_PREFIX}/include CACHE INTERNAL "" FORCE)

add_library(ZeroMQ INTERFACE)
target_include_directories(ZeroMQ INTERFACE $<BUILD_INTERFACE:${ZeroMQ_INCLUDE_DIR}>)
target_link_libraries(ZeroMQ INTERFACE $<BUILD_INTERFACE:${ZeroMQ_LIBRARIES}>)
add_dependencies(ZeroMQ BUILTIN_ZeroMQ)

# necessary for cppzmq:
add_library(libzmq ALIAS ZeroMQ)

set(ZeroMQ_DIR ${ZeroMQ_PREFIX}/share/cmake/ZeroMQ CACHE INTERNAL "" FORCE)
