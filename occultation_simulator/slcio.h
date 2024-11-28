//
// Created by cest on 3/7/24.
// simulation light curve input output

#ifndef OCCULTATION_LIGHT_CURVES_HDF5_H
#define OCCULTATION_LIGHT_CURVES_HDF5_H
#include <gsl/gsl_matrix.h>
#include <hdf5.h>

typedef struct HDF5Wrapper HDF5Wrapper;

typedef void (*HDF5WriteDatasetFunc)(HDF5Wrapper* self, const char* datasetName, gsl_matrix *gslMatrix);
typedef gsl_matrix* (*HDF5ReadDatasetFunc)(HDF5Wrapper* self, const char* datasetName, int rows, int cols);
typedef void (*HDF5FreeFunc)(HDF5Wrapper* self);
typedef double (*HDF5ReadScalar)(HDF5Wrapper* self, const char* datasetName);
typedef void (*HDF5WriteScalar)(HDF5Wrapper* self, const char* datasetName, double value);

struct HDF5Wrapper {
    hid_t file_id;
    HDF5WriteDatasetFunc writeDataset;
    HDF5ReadDatasetFunc readDataset;
    HDF5ReadScalar  readScalar;
    HDF5WriteScalar writeScalar;
    HDF5FreeFunc free;
};

HDF5Wrapper* NewHDF5Wrapper(const char* filename);


#endif //OCCULTATION_LIGHT_CURVES_HDF5_H
