cmake_minimum_required(VERSION 3.15)
project(cec2026clsp VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CORE_SRC
    src/instance_loader.cpp
    src/model.cpp
    src/evo_ops.cpp
    src/nsgaii.cpp
    src/nsgaiii.cpp
    src/moead.cpp
    src/algo_clsp.cpp
    src/experiment_runner.cpp
)

add_executable(clsp_demo src/main.cpp ${CORE_SRC})
target_include_directories(clsp_demo PRIVATE ${CMAKE_SOURCE_DIR}/src)
if(NOT DEFINED INSTANCE_ROOT)
    set(INSTANCE_ROOT "${CMAKE_SOURCE_DIR}")
endif()
target_compile_definitions(clsp_demo PRIVATE INSTANCE_ROOT="${INSTANCE_ROOT}")
