project(EnvelopeExtraction)

cmake_minimum_required(VERSION 3.2.0)

# C++14 support
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++14")

# Sources
set (EnvelopeExtraction_SRC
        "src/main.cpp"
        "src/HilbertTransform.cpp"
        )

add_executable(EnvelopeExtraction ${EnvelopeExtraction_SRC})

# EnvelopeExtraction includes
set(FFTW_LIB_DIR "/usr/local/lib")
set(FFTW_SRC_DIR "/usr/local/include")
#set(FFTW_SRC_DIR "" CACHE PATH "")
#set(FFTW_LIB_DIR "" CACHE PATH "")
message("EnvelopeExtraction source directory: ${FFTW_SRC_DIR}")
message("EnvelopeExtraction lib directory: ${FFTW_LIB_DIR}")
include_directories(${FFTW_SRC_DIR})

target_link_libraries(EnvelopeExtraction
    ${FFTW_LIB_DIR}/libfftw3.a
    )

# OpenCV
#find_package(OpenCV REQUIRED)
#if (OpenCV_FOUND)
#  target_link_libraries(burgercpp ${OpenCV_LIBS})
#endif()

# Include other stuff, like Units
include_directories("${CMAKE_SOURCE_DIR}/include")
