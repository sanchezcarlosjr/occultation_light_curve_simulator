//
// Created by cest on 3/7/24.
//

#include "slcio.h"

void HDF5Wrapper_writeDataset(HDF5Wrapper* self, const char* datasetName, gsl_matrix *gslMatrix) {
    hid_t dataset_id = -1, dataspace_id = -1;
    hsize_t dims[2] = {gslMatrix->size1, gslMatrix->size2};
    dataspace_id = H5Screate_simple(2, dims, NULL);
    if (dataspace_id < 0) {
        fprintf(stderr, "Error creating dataspace\n");
        goto cleanup;
    }
    dataset_id = H5Dcreate2(self->file_id, datasetName, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (dataset_id < 0) {
        fprintf(stderr, "Error creating dataset\n");
        goto cleanup;
    }

    // Write the matrix data to the dataset
    if (H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, gslMatrix->data) < 0) {
        fprintf(stderr, "Error writing to dataset\n");
    }

    cleanup:
    // Attempt to close the dataspace if it was successfully created
    if (dataspace_id >= 0) {
        H5Sclose(dataspace_id);
    }
    // Attempt to close the dataset if it was successfully created
    if (dataset_id >= 0) {
        H5Dclose(dataset_id);
    }
}

gsl_matrix* HDF5Wrapper_readDataset(HDF5Wrapper* self, const char* datasetName, int rows, int cols) {
    hid_t dataset_id;
    gsl_matrix* matrix = gsl_matrix_alloc(rows, cols);

    if ((dataset_id = H5Dopen2(self->file_id, datasetName, H5P_DEFAULT)) < 0) {
        fprintf(stderr, "unable to open the dataset in file");
        H5Fclose(self->file_id);
        gsl_matrix_free(matrix);
        return NULL;
    }

    if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix->data) < 0) {
        fprintf(stderr, "unable to read the dataset from file. ");
        H5Dclose(dataset_id);
        H5Fclose(self->file_id);
        gsl_matrix_free(matrix);
        return NULL;
    }

    H5Dclose(dataset_id);

    return matrix;
}

void HDF5Wrapper_free(HDF5Wrapper* self) {
    H5Fclose(self->file_id);
    free(self);
}

// Constructor-like function implementation
HDF5Wrapper* NewHDF5Wrapper(const char* filename) {
    HDF5Wrapper* wrapper = malloc(sizeof(HDF5Wrapper));
    if (wrapper != NULL) {
        // Open the HDF5 file; create if it does not exist
        wrapper->file_id = H5Fopen(filename, H5F_ACC_RDWR | H5F_ACC_CREAT, H5P_DEFAULT);
        if (wrapper->file_id < 0) {
            free(wrapper);
            return NULL;
        }
        wrapper->writeDataset = HDF5Wrapper_writeDataset;
        wrapper->readDataset = HDF5Wrapper_readDataset;
        wrapper->free = HDF5Wrapper_free;
    }
    return wrapper;
}