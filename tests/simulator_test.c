#include <unity.h>
#include "diffraction.h"
#include "numpy.h"
#include "slcio.h"
#include <stdio.h>
#include <gsl/gsl_matrix.h>

void print_gsl_matrix(const gsl_matrix* m) {
    size_t rows = m->size1;
    size_t cols = m->size2;

    printf("[\n");
    for (size_t i = 0; i < rows; i++) {
        printf(" [");
        for (size_t j = 0; j < cols; j++) {
            printf("%8.4f", gsl_matrix_get(m, i, j));
            if (j < cols - 1)
                printf(", ");
        }
        printf("]%s\n", (i < rows - 1) ? "," : "");
    }
    printf("]\n");
}


double M=2048;
double lambda=600e-9;
double nLamb=10;
double d=3000;
double ua=45;


void setUp(void) {
    // Initialize your test setup here
}

void tearDown(void) {
    // Clean up your test setup here
}

void test_it_calculates_the_optimal_plane_size(void)
{
    HDF5Wrapper* hdf5 = NewHDF5Wrapper("testing_simulation.h5");
    double expected = hdf5->readScalar(hdf5, "total_plane_size");
    double D = calcPlano(d, lambda, ua);
    TEST_ASSERT_FLOAT_WITHIN(0.000000000001, D, expected);
}

void test_it_generates_a_circular_obstruction(void) {
    HDF5Wrapper* hdf5 = NewHDF5Wrapper("testing_simulation.h5");
    gsl_matrix* expected_matrix = hdf5->readDataset(hdf5, "circular_pupil", M, M);
    TEST_ASSERT_EQUAL(expected_matrix->size1, M);
    TEST_ASSERT_EQUAL(expected_matrix->size2, M);
    TEST_ASSERT_EQUAL(gsl_matrix_get(expected_matrix, 0, 0), 1);
    double D = calcPlano(d, lambda, ua);
    gsl_matrix* actual_matrix = pupilCO(M, D, d);
    int result = gsl_matrix_equal_with_tolerance(expected_matrix, actual_matrix, 0.000000000001);
    TEST_ASSERT_EQUAL(1, result);
    gsl_matrix_free(expected_matrix);
    gsl_matrix_free(actual_matrix);
}

void test_it_generates_a_contact_binary(void) {
    HDF5Wrapper* hdf5 = NewHDF5Wrapper("testing_simulation.h5");
    gsl_matrix* expected_matrix = hdf5->readDataset(hdf5, "O2", M, M);
    TEST_ASSERT_EQUAL(expected_matrix->size1, M);
    TEST_ASSERT_EQUAL(expected_matrix->size2, M);
    TEST_ASSERT_EQUAL(gsl_matrix_get(expected_matrix, 0, 0), 1);
    gsl_matrix* actual_matrix = pupilDoble(M, calcPlano(d, lambda, ua), d);
    TEST_ASSERT_EQUAL(gsl_matrix_get(actual_matrix, 0, 0), 1);
    TEST_ASSERT_EQUAL(gsl_matrix_get(actual_matrix, 1024, 1024), 0);
    int result = gsl_matrix_equal_with_tolerance(expected_matrix, actual_matrix, 0);
    TEST_ASSERT_EQUAL(1, result);
    gsl_matrix_free(expected_matrix);
    gsl_matrix_free(actual_matrix);
}


void test_it_calculates_distance_of_object_in_meters(void) {
    HDF5Wrapper* hdf5 = NewHDF5Wrapper("testing_simulation.h5");
    double expected = hdf5->readScalar(hdf5, "object_distance");
    double actual = astronomicalUnitsToMeters(ua);
    TEST_ASSERT_FLOAT_WITHIN(0.000000000001, actual, expected);
}

void test_calculates_monochromatic_diffraction_pattern(void) {
    HDF5Wrapper* hdf5 = NewHDF5Wrapper("testing_simulation.h5");
    gsl_matrix* expected_matrix = hdf5->readDataset(hdf5, "monochromatic_diffraction_pattern", M, M);
    TEST_ASSERT_EQUAL(expected_matrix->size1, M);
    TEST_ASSERT_EQUAL(expected_matrix->size2, M);

    TEST_ASSERT_EQUAL(gsl_matrix_get(expected_matrix, 0, 0), 1);


    double total_plane_size = hdf5->readScalar(hdf5, "total_plane_size");
    double object_distance = hdf5->readScalar(hdf5, "object_distance");
    gsl_matrix* circular_pupil = hdf5->readDataset(hdf5, "circular_pupil", M, M);


    gsl_matrix* actual_matrix = fresnel(circular_pupil, M, total_plane_size, object_distance, lambda);

    int result = gsl_matrix_equal_with_tolerance(expected_matrix, actual_matrix, 0.7);
    TEST_ASSERT_EQUAL(1, result);


    gsl_matrix_free(expected_matrix);
    hdf5->free(hdf5);
}

void test_calculates_chromatic_diffraction_pattern(void) {
    HDF5Wrapper* hdf5 = NewHDF5Wrapper("testing_simulation.h5");
    gsl_matrix* expected_matrix = hdf5->readDataset(hdf5, "chromatic_diffraction_pattern", M, M);
    TEST_ASSERT_EQUAL(expected_matrix->size1, M);
    TEST_ASSERT_EQUAL(expected_matrix->size2, M);

    TEST_ASSERT_EQUAL(gsl_matrix_get(expected_matrix, 0, 0), 1);


    double total_plane_size = hdf5->readScalar(hdf5, "total_plane_size");
    double object_distance = hdf5->readScalar(hdf5, "object_distance");
    double nEst = 1;
    gsl_matrix* circular_pupil = hdf5->readDataset(hdf5, "circular_pupil", M, M);

    gsl_matrix* chromatic_diffraction_pattern = spectra(circular_pupil, M, total_plane_size, object_distance, nEst, nLamb);

    int result = gsl_matrix_equal_with_tolerance(expected_matrix, chromatic_diffraction_pattern, 0.7);
    TEST_ASSERT_EQUAL(1, result);


    gsl_matrix_free(expected_matrix);
    hdf5->free(hdf5);
}

void test_calculates_pattern_for_extended_source(void) {
    HDF5Wrapper* hdf5 = NewHDF5Wrapper("testing_simulation.h5");
    hdf5->readText(hdf5, "star_type");

}





int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_calculates_pattern_for_extended_source);
    RUN_TEST(test_it_calculates_the_optimal_plane_size);
    RUN_TEST(test_it_generates_a_circular_obstruction);
    RUN_TEST(test_it_calculates_distance_of_object_in_meters);
    RUN_TEST(test_calculates_monochromatic_diffraction_pattern);
    RUN_TEST(test_calculates_chromatic_diffraction_pattern);

    return UNITY_END();
}
