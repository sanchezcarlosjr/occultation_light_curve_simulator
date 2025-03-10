//
// Created by cest on 1/31/24.
//

#include "numpy.h"
#include <gsl/gsl_matrix.h>

int gsl_matrix_equal_with_tolerance(const gsl_matrix* matrixA, const gsl_matrix* matrixB, double tolerance) {
    if (matrixA->size1 != matrixB->size1 || matrixA->size2 != matrixB->size2) {
        // Matrices have different dimensions
        return 0;
    }

    for (size_t i = 0; i < matrixA->size1; ++i) {
        for (size_t j = 0; j < matrixA->size2; ++j) {
            double diff = fabs(gsl_matrix_get(matrixA, i, j) - gsl_matrix_get(matrixB, i, j));
            if (diff > tolerance) {
                printf("Diference in element (row=%zu, column=%zu, matrixA=%f, matrixB=%f, difference=%f)\n", i, j, gsl_matrix_get(matrixA, i, j), gsl_matrix_get(matrixB, i, j), diff);
                return 0;
            }
        }
    }

    // All elements are equal within the specified tolerance
    return 1;
}

gsl_matrix* gsl_zeros(int rows, int cols) {
    // Allocate a matrix of size rows x cols
    gsl_matrix* m = gsl_matrix_alloc(rows, cols);

    // Check if the matrix was successfully allocated
    if (m != NULL) {
        // Set all elements of the matrix to zero
        gsl_matrix_set_zero(m);
    }

    return m;  // Return the pointer to the newly created matrix
}

void gsl_print_matrix(const gsl_matrix* m) {
    for (size_t i = 0; i < m->size1; i++) {  // size1 is the number of rows
        for (size_t j = 0; j < m->size2; j++) {  // size2 is the number of columns
            printf("%g ", gsl_matrix_get(m, i, j));
        }
        printf("\n");
    }
}

void gsl_linspace(double start, double end, size_t num, gsl_vector* vector) {
    double step = (end - start) / (num - 1);
    for (size_t i = 0; i < num; i++) {
        gsl_vector_set(vector, i, start + step * i);
    }
}

void linspace(double *arr, double start, double end, int n) {
    double step = (end - start) / (n - 1);
    for(int i = 0; i < n; i++) {
        arr[i] = start + i * step;
    }
}

void meshgrid(double *x, double *y, int n, double **X, double **Y) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            X[i][j] = x[j];
            Y[i][j] = y[i];
        }
    }
}

Matrix zeros(int rows, int cols) {
    double** matrix = (double**) malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*) calloc(cols, sizeof(double));  // Use calloc to initialize to 0
    }
    return matrix;
}

Matrix full(int rows, int cols, int value) {
    double** matrix = (double**) malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = value;
        }
    }
    return matrix;
}

Matrix ones(int rows, int cols) {
    return full(rows, cols, 1);
}

void freeMatrix(Matrix matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

Array arange(int start, int stop, int step) {
    Array array;
    array.length = (int)((stop - start) / step) + 1;
    array.array = (double *)malloc(array.length * sizeof(double));
    for (int i = 0; i < array.length; i++) {
        array.array[i] = step * (i + 1);
    }
    return array;
}


// Function to initialize a 2D array of complex numbers with zeros
void zerosComplex(ComplexMatrix array, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array[i][j] = 0.0 + 0.0 * I; // Initialize to complex zero
        }
    }
}

ComplexMatrix allocateComplex2DArray(int rows, int cols) {
    ComplexMatrix array = (ComplexMatrix) malloc(rows * sizeof(double complex *));
    if (!array) return NULL;

    for (int i = 0; i < rows; i++) {
        array[i] = (double complex *) malloc(cols * sizeof(double complex));
        if (!array[i]) {
            // Free previously allocated memory in case of failure
            for (int j = 0; j < i; j++) free(array[j]);
            free(array);
            return NULL;
        }
    }
    return array;
}

// Function to free a dynamically allocated 2D array of complex numbers
void freeComplex2DArray(ComplexMatrix array, int rows) {
    for (int i = 0; i < rows; i++) free(array[i]);
    free(array);
}

typedef int (*MathFunc)(float, int);

int do_math(float arg1, int arg2) {
    return arg2;
}

int call_a_func(MathFunc call_this) {
    int output = call_this(5.5, 7);

    return output;
}
