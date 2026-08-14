cmake_minimum_required(VERSION 3.16)

project(crtbp_cpp_port LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(crtbp_cpp
    src/ode.cpp
    src/crtbp.cpp
    src/hill4bp.cpp
)

target_include_directories(crtbp_cpp
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_options(crtbp_cpp
    PRIVATE
        $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall>
        $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wextra>
        $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wpedantic>
)

add_executable(hill4bp_poincare
    apps/hill4bp_poincare.cpp
)

target_link_libraries(hill4bp_poincare
    PRIVATE
        crtbp_cpp
)

target_compile_options(hill4bp_poincare
    PRIVATE
        $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall>
        $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wextra>
        $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wpedantic>
)

include(CTest)
if(BUILD_TESTING)
    add_executable(capd_decimal_literal_tests
        tests/capd_decimal_literal_tests.cpp
    )
    target_compile_options(capd_decimal_literal_tests
        PRIVATE
            $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall>
            $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wextra>
            $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wpedantic>
    )
    add_test(NAME capd_decimal_literal COMMAND capd_decimal_literal_tests)

    add_executable(hill4bp_regression_tests
        tests/hill4bp_regression_tests.cpp
    )
    target_link_libraries(hill4bp_regression_tests
        PRIVATE
            crtbp_cpp
    )
    target_compile_options(hill4bp_regression_tests
        PRIVATE
            $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall>
            $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wextra>
            $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wpedantic>
    )
    add_test(NAME hill4bp_regression COMMAND hill4bp_regression_tests)
    set_tests_properties(hill4bp_regression PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    add_test(
        NAME neck_sweep_shell_syntax
        COMMAND bash -n "${CMAKE_CURRENT_SOURCE_DIR}/apps/run_neck_sweep_capd.sh"
    )
    add_test(
        NAME neck_sweep_contract
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_neck_sweep.sh"
    )
    add_test(
        NAME capd_build_script_shell_syntax
        COMMAND bash -n "${CMAKE_CURRENT_SOURCE_DIR}/apps/build_capd_apps.sh"
    )
    add_test(
        NAME professional_neck_study_shell_syntax
        COMMAND bash -n
            "${CMAKE_CURRENT_SOURCE_DIR}/apps/run_professional_neck_study.sh"
    )
    add_test(
        NAME professional_neck_study_contract
        COMMAND bash
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_professional_neck_study.sh"
    )
endif()
