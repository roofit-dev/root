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
    URL https://github.com/zeromq/libzmq/archive/7c2df78b49a3aa63e654b3f3526adf71ed091534.tar.gz
    URL_HASH SHA256=fcc1b0648afa5d92e0ff0e6e93beb28cbbe008a5f98c228ff97144ba6e4a6c3e
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DWITH_PERF_TOOL=OFF
        -DZMQ_BUILD_TESTS=OFF
        -DENABLE_DRAFTS=ON
        -DENABLE_NO_EXPORT=ON
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
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
