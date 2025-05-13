add_executable(ctrecon-math.app)

if (MSVC)
    target_compile_options(ctrecon-math.app
        PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:/utf-8>
            $<$<COMPILE_LANGUAGE:CXX>:/Zc:preprocessor>
            $<$<COMPILE_LANGUAGE:CXX>:/std:c++23preview>
    )
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(ctrecon-math.app
        PRIVATE
            $<$<COMPILE_LANGUAGE:CXX>:-Wall>
            $<$<COMPILE_LANGUAGE:CXX>:-Wextra>
            $<$<COMPILE_LANGUAGE:CXX>:-Wpedantic>
            $<$<COMPILE_LANGUAGE:CXX>:-std=c++2b>
            $<$<COMPILE_LANGUAGE:CXX>:-finput-charset=UTF-8>
            $<$<COMPILE_LANGUAGE:CXX>:-fexec-charset=UTF-8>
    )
endif()


target_include_directories(ctrecon-math.app
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
)

target_sources(ctrecon-math.app
    PRIVATE
        app.cpp
)




target_link_libraries(ctrecon-math.app
    PRIVATE
        ctrecon-math
        ctrecon-math.static
)

find_package(fmt REQUIRED)
target_link_libraries(ctrecon-math.app
    PRIVATE
        fmt::fmt
)