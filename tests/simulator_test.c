#include <unity.h>
#include "diffraction.h"
#include "numpy.h"
#include "slcio.h"

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
    TEST_ASSERT_EQUAL(result, 1);
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
    TEST_ASSERT_EQUAL(result, 1);
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
    
    
}



int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_it_calculates_the_optimal_plane_size);
    RUN_TEST(test_it_generates_a_circular_obstruction);
    RUN_TEST(test_it_calculates_distance_of_object_in_meters);

    return UNITY_END();
}
