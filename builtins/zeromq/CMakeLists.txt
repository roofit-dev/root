include(ExternalProject)

set(ZeroMQ_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/BUILTIN_ZeroMQ-prefix)
set(ZeroMQ_LIBNAME ${CMAKE_STATIC_LIBRARY_PREFIX}zmq${CMAKE_STATIC_LIBRARY_SUFFIX})

set(ZeroMQ_LIBRARY ${ZeroMQ_PREFIX}/lib/${ZeroMQ_LIBNAME} CACHE INTERNAL "" FORCE)
set(ZeroMQ_LIBRARIES ${ZeroMQ_LIBRARY} CACHE INTERNAL "" FORCE)

ExternalProject_Add(BUILTIN_ZeroMQ
    URL https://github.com/zeromq/libzmq/archive/5d8d857540323e2d85c64a7edde1ad5280cad04b.tar.gz
    URL_HASH SHA256=750bc1c67a3740a5438da8762931983cb08dc649e1a658cfc34de91c2cf9222c
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DWITH_PERF_TOOL=OFF
        -DZMQ_BUILD_TESTS=OFF
        -DENABLE_DRAFTS=ON
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        # TODO: remove -w flag when updating to v4.3.5 or higher, current warnings should be gone by then
        -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}\ -w
    BUILD_BYPRODUCTS ${ZeroMQ_LIBRARIES}
    )

ExternalProject_Get_Property(BUILTIN_ZeroMQ install_dir)

set(ZeroMQ_FOUND TRUE CACHE BOOL "" FORCE)

set(ZeroMQ_INCLUDE_DIR  ${ZeroMQ_PREFIX}/include CACHE INTERNAL "" FORCE)
set(ZeroMQ_INCLUDE_DIRS ${ZeroMQ_PREFIX}/include CACHE INTERNAL "" FORCE)

add_library(ZeroMQ INTERFACE)
target_include_directories(ZeroMQ INTERFACE $<BUILD_INTERFACE:${ZeroMQ_INCLUDE_DIR}>)
target_link_libraries(ZeroMQ INTERFACE $<BUILD_INTERFACE:${ZeroMQ_LIBRARIES}>)
add_dependencies(ZeroMQ BUILTIN_ZeroMQ)

add_library(ZeroMQ::ZeroMQ ALIAS ZeroMQ)
add_library(libzmq ALIAS ZeroMQ)

set(ZeroMQ_DIR ${ZeroMQ_PREFIX}/share/cmake/ZeroMQ CACHE INTERNAL "" FORCE)

ExternalProject_Add(BUILTIN_cppzmq
    URL https://github.com/zeromq/cppzmq/archive/refs/tags/v4.8.1.tar.gz
    URL_HASH SHA256=7a23639a45f3a0049e11a188e29aaedd10b2f4845f0000cf3e22d6774ebde0af
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${install_dir}
        -DCPPZMQ_BUILD_TESTS=OFF
    DEPENDS BUILTIN_ZeroMQ
    )
