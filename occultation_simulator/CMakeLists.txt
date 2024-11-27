add_library(${LIBRARY_NAME} STATIC numpy.c diffraction.c slcio.c)


find_package(GSL REQUIRED)
if (GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    target_link_libraries(${LIBRARY_NAME} PUBLIC GSL::gsl GSL::gslcblas)
endif()

# Find FFTW3
find_package(FFTW3 REQUIRED COMPONENTS DOUBLE_LIB)
if (FFTW3_FOUND)
    include_directories(${FFTW3_INCLUDE_DIRS})
    target_link_libraries(${LIBRARY_NAME} PUBLIC fftw3 fftw3_threads)
endif()

# Find HDF5
set(HDF5_NO_FIND_PACKAGE_CONFIG_FILE TRUE)
set(HDF5_USE_STATIC_LIBRARIES OFF)
set(HDF5_PREFER_PARALLEL FALSE)
find_package(HDF5 REQUIRED COMPONENTS C HL)

if (HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
    target_link_libraries(${LIBRARY_NAME} PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
endif()


target_include_directories(${LIBRARY_NAME} PUBLIC "./")
