//
// Created by cest on 1/31/24.
//

#include "diffraction.h"
#include "numpy.h"
#include <string.h>

#include <libgen.h> // For dirname function
typedef double complex** DiffractionPattern;


void getLibDir(char *libdir, const char *filepath) {
    strcpy(libdir, filepath);
    char *dir = dirname(libdir); // Extract directory name
    strcpy(libdir, dir); // Copy the directory path back into libdir
}

/**
 * Convert Cartesian coordinates to polar coordinates.
 *
 * @param x The x-coordinate.
 * @param y The y-coordinate.
 * @param phi Pointer to the angle variable to be set.
 * @param rho Pointer to the radius variable to be set.
 */
void cart2pol(double x, double y, double *phi, double *rho) {
    *phi = atan2(y, x);          // Calculate the angle in radians
    *rho = sqrt(x * x + y * y);  // Calculate the radius
}

void pol2cart(double rho, double phi, double *x, double *y) {
    *x = rho * cos(phi);  // Calculate the x coordinate
    *y = rho * sin(phi);  // Calculate the y coordinate
}


/* Transform astronomical units to meters
 *  @param double ua Astronomical units
 *
 */
double astronomicalUnitsToMeters(double ua) {
    return ua * 1.496e11;
}


/**
 * Generate a circular obstruction.
 *
 * @param M The size of the matrix in pixels.
 * @param D The size of the matrix in meters.
 * @param d The central dimming in meters.
 */
gsl_matrix *pupilCO(int M, double D, double d) {
    gsl_matrix *P = gsl_matrix_alloc(M, M);
    gsl_vector *m = gsl_vector_alloc(M);

    // Populate the vector m with linearly spaced values
    for (int i = 0; i < M; i++) {
        gsl_vector_set(m, i, -D / 2 + i * D / (M - 1));
    }

    double phi, rho;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            double x = gsl_vector_get(m, i);
            double y = gsl_vector_get(m, j);
            cart2pol(x, y, &phi, &rho);
            gsl_matrix_set(P, i, j, (rho >= d / 2) ? 1.0 : 0.0);
        }
    }

    gsl_vector_free(m);
    return P;
}


/**
 * Translates a diffraction pattern represented by a GSL complex matrix.
 * The translation is applied by shifting the matrix in both X and Y directions
 * with wrap-around (circular shift).
 *
 * @param P The input GSL complex matrix representing the diffraction pattern.
 * @param dx The translation distance in the X direction (columns).
 * @param dy The translation distance in the Y direction (rows).
 * @return A new GSL complex matrix representing the translated diffraction pattern.
 */
gsl_matrix* translate(const gsl_matrix* P, int dx, int dy) {
    // Allocate a new matrix of the same size as P, initialized to zeros
    gsl_matrix* result = gsl_matrix_calloc(P->size1, P->size2);

    // Get the dimensions of the matrix
    size_t rows = P->size1;
    size_t cols = P->size2;

    // Loop over each element in the matrix P
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            // Calculate the new position after translation
            int new_i = i + dy;
            int new_j = j + dx;

            // Check if the new position is within the bounds of the matrix
            if (new_i >= 0 && new_i < rows && new_j >= 0 && new_j < cols) {
                // Get the value from the original matrix
                double val = gsl_matrix_get(P, i, j);
                // Set the value in the new position of the result matrix
                gsl_matrix_set(result, new_i, new_j, val);
            }
        }
    }

    // Return the translated matrix
    return result;
}

/**
 * Create a dynamic circular obstruction within a square matrix.
 *
 * @param M Size of the square matrix in pixels.
 * @param D Size of the matrix in meters.
 * @param d Central dimming in meters, treated as if it were circular.
 * @param sepX Horizontal separation between the centers of the two obstructions.
 * @param sepY Vertical separation between the centers of the two obstructions.
 * @param r1 Radius of the larger circular obstruction.
 * @param r2 Radius of the smaller circular obstruction.
 * @return A GSL matrix representing the combined circular obstructions.
 */
gsl_matrix* pupilDoble(int M, double D, double d) {
    double r1 = (d / 2) * 0.65;
    double r2 = sqrt(pow(d / 2, 2) - pow(r1, 2));
    double d1 = r1 * 2;
    double d2 = r2 * 2;
    double Dx = 0.45 * d1 + 0.45 * d2;  // Orientation in X
    double Dy = 0;
    int sepX = (int)((Dx / 2) / D * M);
    int sepY = (int)((Dy / 2) / D * M);

    gsl_matrix* P1 = gsl_matrix_alloc(M, M);
    gsl_matrix* P2 = gsl_matrix_alloc(M, M);

    // Generate obstructions
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            double x = -D / 2 + D * i / M;
            double y = -D / 2 + D * j / M;
            double r = sqrt(x * x + y * y);
            gsl_matrix_set(P1, i, j, r >= r1 ? 1.0 : 0.0);  // Large obstruction
            gsl_matrix_set(P2, i, j, r >= r2 ? 1.0 : 0.0);  // Small obstruction
        }
    }

    // Translate and combine obstructions using the 'translate' function
    gsl_matrix* translatedP1 = translate(P1, -sepX, sepY);
    gsl_matrix* translatedP2 = translate(P2, sepX, sepY);
    gsl_matrix* P = gsl_matrix_alloc(M, M);
    gsl_matrix_set_zero(P);

    // Combine and binarize
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            double val1 = gsl_matrix_get(translatedP1, i, j);
            double val2 = gsl_matrix_get(translatedP2, i, j);

            gsl_matrix_set(P, i, j, val1+val2 == 2 ? 1.0 : 0.0);
        }
    }

    gsl_matrix_free(P1);
    gsl_matrix_free(P2);
    gsl_matrix_free(translatedP1);
    gsl_matrix_free(translatedP2);

    return P;
}

/**
 * Wrapper function to create a double circular obstruction with a custom separator.
 *
 * @param M Size of the square matrix in pixels.
 * @param D Size of the matrix in meters.
 * @param d Central dimming in meters, treated as if it were circular.
 * @param sepX Horizontal separation between the centers of the two obstructions.
 * @param sepY Vertical separation between the centers of the two obstructions.
 * @return A GSL matrix representing the double circular obstructions with a custom separator.
 */
//gsl_matrix_complex* pupilDobleWithCustomSeparator(int M, double D, double d, int sepX, int sepY) {
//    double r1 = (d / 2) * 0.65;
//    double r2 = sqrt(pow(d / 2, 2) - pow(r1, 2));
//    double d1 = r1 * 2;
//    double d2 = r2 * 2;
//    double Dx = 0.45 * d1 + 0.45 * d2;
//    double Dy = 0;
//    return dynamicPupilDoble(M, D, d, (int) ((Dx / 2) / D * M) + sepX, (int) ((Dy / 2) / D * M) + sepY, r1, r2);
//}
//
///**
// * Wrapper function to create a double circular obstruction with no separator.
// *
// * @param M Size of the square matrix in pixels.
// * @param D Size of the matrix in meters.
// * @param d Central dimming in meters, treated as if it were circular.
// * @return A GSL matrix representing the double circular obstructions with no separator.
// */
//gsl_matrix_complex* pupilDoble(int M, double D, double d) {
//    return pupilDobleWithCustomSeparator(M, D, d, 0, 0);
//}


// Function to generate a circular aperture
// TODO: Ellipse
ComplexMatrix pupilCA(int M, double D, double d) {
    ComplexMatrix P = allocateComplex2DArray(M, M);
    zerosComplex(P, M, M);
    double step = D / M;
    double phi, rho;

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            double x = -D / 2 + i * step;  // Calculate x coordinate
            double y = -D / 2 + j * step;  // Calculate y coordinate
            cart2pol(x, y, &phi, &rho);    // Convert to polar coordinates
            P[i][j] = (rho <= d / 2) ? 1.0 : 0.0;  // Mark points within the aperture
        }
    }

    return P;
}

// Function to generate a square obstruction
double **pupilSO(int M, double D, double d) {
    double **P = ones(M, M);
    int t = (int) (M * d / D);  // Size of the obstruction in pixels
    int c = M / 2;  // Center of the matrix

    // Calculate the bounds of the obstruction
    int start = c - t / 2;
    int end = c + t / 2;

    // Create the square obstruction
    for (int i = start; i < end; i++) {
        for (int j = start; j < end; j++) {
            if (i >= 0 && i < M && j >= 0 && j < M) {
                P[i][j] = 0.0;  // Mark the obstruction area with zeros
            }
        }
    }

    return P;
}


void pupilSA(ComplexMatrix P, int M, double D, double d) {
    int t = M * d / D; // Calculate the size of the central obscuration in pixels
    int c = M / 2;     // Center of the array

    // Initialize the entire array to 0
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            P[i][j] = 0.0;
        }
    }

    // Set the central square area to 1
    for (int i = c - t / 2; i < c + t / 2; i++) {
        for (int j = c - t / 2; j < c + t / 2; j++) {
            if (i >= 0 && i < M && j >= 0 && j < M) { // Check boundaries
                P[i][j] = 1.0;
            }
        }
    }
}

void fresnel(double complex *U0, int nx, int ny, int M, double plano, double z, double lmda, double complex *output_intensity) {
    double k = 2 * M_PI / lmda;
    double x = (plano / M) * nx; // Normally nx = M, so x = plano in meters
    double y = (plano / M) * ny;
    double fx = 1 / x; // Spatial frequency in m**-1
    double fy = 1 / y;

    // Allocate memory for u, v, O, H, and U arrays
    double complex *u = (double complex *) fftw_malloc(sizeof(double complex) * nx);
    double complex *v = (double complex *) fftw_malloc(sizeof(double complex) * ny);
    double complex *O = (double complex *) fftw_malloc(sizeof(double complex) * nx * ny);
    double complex *H = (double complex *) fftw_malloc(sizeof(double complex) * nx * ny);
    double complex *U = (double complex *) fftw_malloc(sizeof(double complex) * nx * ny);

    // Compute u and v arrays
    for (int i = 0; i < nx; i++) {
        u[i] = (i - nx / 2) * fx;
    }
    for (int j = 0; j < ny; j++) {
        v[j] = (j - ny / 2) * fy;
    }

    // Apply FFT to U0 using FFTW
    fftw_plan plan_forward = fftw_plan_dft_2d(nx, ny, U0, O, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    // Compute H array
    // np.exp(1j*k*z)*np.exp(-1j*np.pi*(lmda*z)*(u**2+v**2))
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double u_squared = creal(u[i]) * creal(u[i]);
            double v_squared = creal(v[j]) * creal(v[j]);
            H[i * ny + j] = cexp(I * k * z) * cexp(-I * M_PI * lmda * z * (u_squared + v_squared));
        }
    }

    // Multiply O and H element-wise
    for (int i = 0; i < nx * ny; i++) {
        U[i] = O[i] * H[i];
    }

    // Apply inverse FFT to U using FFTW
    fftw_plan plan_backward = fftw_plan_dft_2d(nx, ny, U, output_intensity, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    // Compute the intensity pattern
    for (int i = 0; i < nx * ny; i++) {
        output_intensity[i] = cpow(cabs(output_intensity[i]), 2);
    }

    // Free memory and destroy FFTW plans
    fftw_free(u);
    fftw_free(v);
    fftw_free(O);
    fftw_free(H);
    fftw_free(U);
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
}

/* listadat.txt--> A0=1;A1=2;A2=3;A3=4;A4=5;A5=6;A7=7;F0=8;F2=9;F3=10;F5=11;F6=12;F7=13;F8=14;
 * G0=15;G1=16;G2=17;G5=18;G8=19;K0=20;K1=21;K2=22;K3=23;K4=24;K5=25;K7=26;
 * M0=27;M1=28;M2=29;M3=30;M4=31;M5=32;M6=33;M7=34;M8=35
*/
void spectra(double complex *U0, int nx, int ny, int M, double plano, double z, int nEst, int nLmdas, double complex *acc) {
    char libdir[1024];
    getLibDir(libdir, __FILE__);
    // Define file paths - you'll need to define how libdir is set
    char refFilePath[1024];
    sprintf(refFilePath, "%s/listadat.txt", libdir);

    FILE *refFile = fopen(refFilePath, "r");
    if (!refFile) {
        perror("Failed to open reference file");
        return;
    }

    // Read the nEst-th line from the reference file to get the data file name
    char line[256];
    char dataFileName[256];
    int currentLine = 0;
    while (fgets(line, sizeof(line), refFile)) {
        currentLine++;
        if (currentLine == nEst) {
            strncpy(dataFileName, line, strlen(line) - 1);  // Remove newline character
            break;
        }
    }
    fclose(refFile);

    // Construct the path to the data file
    char dataFilePath[1024];
    sprintf(dataFilePath, "%s/spectra/%s", libdir, dataFileName);

    // Open the data file
    FILE *dataFile = fopen(dataFilePath, "r");
    if (!dataFile) {
        perror("Failed to open data file");
        return;
    }

    int U0Size = nx * ny;

    // Read data from the file and perform calculations
    double lamda, peso;
    double complex *output_intensity = (double complex *) malloc(U0Size * sizeof(double complex));
    int count = 0;

    // Assuming the data file has two columns: lamda and peso
    while (fscanf(dataFile, "%lf,%lf", &lamda, &peso) == 2) {
        fresnel(U0, nx, ny, M, plano, z, lamda * 1e-10, output_intensity);

        // Weight the result by 'peso' and accumulate
        for (int i = 0; i < U0Size; i++) {
            acc[i] += output_intensity[i] * peso;
        }

        count++;
        if (count >= nLmdas) break;  // Stop after nLmdas iterations
    }

    fclose(dataFile);
    free(output_intensity);

    // Normalize and finalize the result
    for (int i = 0; i < U0Size; i++) {
        acc[i] /= acc[0];  // Normalize by the first element
    }


}

typedef struct {
    char tipo[10];  // Spectral type
    double T;       // Temperature
    double M;       // Absolute magnitude
    double L;       // Luminosity relative to the Sun
} Star;

/*
 * Function to calculate the apparent radii of stars
% mV --> Apparent magnitude
% nEst --> star number
% ua --> Distance to the object in astronomical units
% Absolute magnitudes in order from A0 type stars to M8 type
* Star[] -> The caller must read the data, for example, from "estrellas.csv", and convert it into a Star struct.
% M0 = [1.5 1.7 1.8 2.0 2.1 2.2 2.4 3.0 3.3 3.5 3.7 4.0 4.3 4.4 4.7 4.9 5.0...
%     5.2 2.6 6.0 6.2 6.4 6.7 7.1 7.4 8.1 8.7 9.4 10.1 10.7 11.2 12.3 13.4...
%     13.9 14.4];
OUT --> type, R_star: spectral type chosen and calculated star radius, respectively
 */
void calcRstar(double mV, int nEst, double ua, Star stars[], int numStars, char *tipo, double *R_star) {
    // Constants
    double ua_meters = 1.496e11 * ua;  // Distance in meters
    double Tsol = 5780;                // Solar temperature in Kelvin
    double Rsol = 6.96e8;              // Solar radius in meters

    if (nEst < 1 || nEst > numStars) {
        printf("Invalid star number.\n");
        return;
    }

    // Select the star
    Star star = stars[nEst - 1];

    // Calculations
    double d1 = pow(10, (mV - star.M + 5) / 5);
    double d = 3.085e16 * d1;  // Convert from parsecs to meters
    double Rst = sqrt(star.L) / pow(star.T / Tsol, 2);  // Star radius in Rsol
    double alfa = (Rsol * Rst) / d;  // Angular size of the star in radians
    *R_star = alfa * ua_meters;  // Size of the star in meters

    // Copy spectral type to output
    snprintf(tipo, 10, "%s", star.tipo);
}

/**
 * Calculates the average point diffraction by simulating the diffraction pattern translation
 * and intensity accumulation across multiple steps.
 *
 * @param diffractionPattern A GSL complex matrix representing the initial diffraction pattern.
 * @param R_star The apparent radius of the star in meters.
 * @param plano The size of the screen (plane) in meters.
 * @param M The size of the matrix in pixels (assuming a square matrix).
 * @param d The diameter of the object causing the diffraction in meters.
 * @param intensityOut A GSL complex matrix to store the output intensity pattern.
 */
void promedioPD(gsl_matrix_complex *diffractionPattern, double R_star, double plano, int M, double d, gsl_matrix_complex *intensityOut) {

}


void extraerPerfil(double **I0, int M, double D, double T, double b, double *x, double *y) {

}

/**
 * Function to calculate the optimal plane size (for both the object and diffraction) for small objects (<10km),
 * avoiding the FFT scaling issue.
 *
 * @param d Diameter of the object in meters.
 * @param lmda Wavelength in meters.
 * @param ua Distance of the object in Astronomical Units (AU).
 *
 * @return Plane size in meters (one dimension).
 */
double calcPlano(double d, double lmda, double ua) {
    const double AU_TO_METERS = 1.496e11; // Conversion factor from AU to meters
    double z = ua * AU_TO_METERS; // Distance in meters
    double fscale = sqrt(lmda * z / 2); // Fresnel scale
    double Rho = d / (2 * fscale);
    double plano = (50 * d) / Rho; // Size of the plane in meters
    return plano;
}

double SNR_TAOS2(double mV) {
    // Polynomial coefficients
    double p1 = 1.5792;
    double p2 = -57.045;
    double p3 = 515.04;

    // Polynomial calculation for SNR
    double SNR = p1 * mV * mV + p2 * mV + p3;

    return SNR;
}


/*
 * Add Poisson noise to an image
diffractionPattern --> image matrix
mV --> apparent magnitude of the star
OUT --> In: matrix with added noise, assuming NOISE=1/SNR calculated from TAOS-II
 */
void addNoise(DiffractionPattern diffractionPattern, int M, double mV) {
    double noise = 1 / SNR_TAOS2(mV);

    // Initialize GSL random number generator
    const gsl_rng_type *T;
    gsl_rng *GSLRandomNumberGenerator;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    GSLRandomNumberGenerator = gsl_rng_alloc(T);

    double mean = 0.0; // For calculating the mean of the noise mask

    // Add Poisson noise to each pixel and calculate mean of the noise mask
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            double poisson_noise = gsl_ran_poisson(GSLRandomNumberGenerator, diffractionPattern[i][j]);
            diffractionPattern[i][j] += poisson_noise; // Add Poisson noise to the original image
            mean += poisson_noise;
        }
    }
    mean /= (M * M); // Calculate mean of the noise mask

    // Normalize the noise mask and add weighted noise according to TAOS-II
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            double normalized_noise = (diffractionPattern[i][j] / mean) * noise - noise;
            diffractionPattern[i][j] += normalized_noise; // Add the normalized and weighted noise to the image
        }
    }

    gsl_rng_free(GSLRandomNumberGenerator); // Free the GSL random number generator
}

/*
 * Function to sample the diffraction profile by obtaining the average point
lc --> diffraction profile or light curve
D --> size of the plane in meters
vr --> object velocity ~5000 m/s (positive if it goes against Earth's velocity)
fps --> frames per second of the camera, 20 for TAOS-2
toff --> Time offset within the sampling period
vE --> translational velocity of the Earth == 29800 m/s
opangle --> angle from the object's opposition: O, S, E
ua --> Distance of the object in Astronomical Units
OUT --> s_lin, lc_lin, s_pun, lc_pun: time vectors for lines, sample in lines, time at points, and sample at points, RESPECTIVELY
 */
void sampling(double *lc, int tam, double D, double vr, double fps, double toff, double vE, double opangle, double ua,
              double **s_lin, double **lc_lin, double **s_pun, double **lc_pun, int *size_s_lin, int *size_lc_lin, int *size_s_pun, int *size_lc_pun) {
    double T = 1 / fps;  // Exposure time
    double OA = opangle * M_PI / 180;  // Opposition angle in radians

    // Tangential velocity of the object relative to Earth
    double Vt = vE * (cos(OA) - sqrt((1 / ua) * (1 - (1 / (ua * ua)) * sin(OA) * sin(OA)))) + vr;

    double t = D / Vt;  // Visibility of the plane in seconds
    int Nm = (int)(t / T);  // Total number of samples in the observation plane
    int dpix = tam / Nm;
    int pixoffset = (int)(toff * fps);  // Convert time offset to pixel offset

    // Allocate memory for output arrays
    *size_s_lin = tam;  // Adjust the size according to your sampling
    *size_lc_lin = tam;
    *size_s_pun = Nm;  // Adjust the size according to your sampling
    *size_lc_pun = Nm;
    *s_lin = (double *)malloc(*size_s_lin * sizeof(double));
    *lc_lin = (double *)malloc(*size_lc_lin * sizeof(double));
    *s_pun = (double *)malloc(*size_s_pun * sizeof(double));
    *lc_pun = (double *)malloc(*size_lc_pun * sizeof(double));

    // Sample the diffraction profile
    // Assuming the sampling logic is implemented similar to the Python version
    // You need to translate the logic of sampling `lc` into `lc_lin` and `lc_pun` here

    // Fill the time vectors `s_lin` and `s_pun`
    for (int i = 0; i < *size_s_lin; i++) {
        (*s_lin)[i] = -t / 2 + i * (t / (*size_s_lin - 1));
    }
    for (int i = 0; i < *size_s_pun; i++) {
        (*s_pun)[i] = -t / 2 + i * (t / (*size_s_pun - 1));
    }
}


// Function to find peaks in the light curve using the derivative method
// x, y: Arrays containing the occultation data (distance and amplitude)
// n: Number of elements in the x and y arrays
// D: Diameter of the object in meters
// fil: Threshold value for identifying peaks, default is 0.005
// peaks: Array to store the indices of the peaks
// peakValues: Array to store the values of the peaks
// Returns the number of peaks found
int searchPeaks(double *x, double *y, int n, double D, double fil, int **peaks, double **peakValues) {
    double *yp = (double *)malloc((n - 1) * sizeof(double));  // Array for the derivative of y
    int *tempPeaks = (int *)malloc(n * sizeof(int));          // Temporary array for peak indices
    double *tempPeakValues = (double *)malloc(n * sizeof(double));  // Temporary array for peak values
    int count = 0;  // Counter for the number of peaks

    // Calculate the derivative of y
    for (int i = 0; i < n - 1; i++) {
        yp[i] = y[i + 1] - y[i];
    }

    // Find indices where the derivative is close to 0 and within the region of interest
    for (int i = 0; i < n - 1; i++) {
        if (fabs(yp[i]) < fil && fabs(x[i]) < (D / 2)) {
            tempPeaks[count] = i;
            tempPeakValues[count] = y[i];
            count++;
        }
    }

    // Allocate memory for the output arrays with the exact number of peaks
    *peaks = (int *)malloc(count * sizeof(int));
    *peakValues = (double *)malloc(count * sizeof(double));

    // Copy data to the output arrays, filtering out duplicates
    int lastIdx = -1;  // Last index added to the peaks array
    for (int i = 0; i < count; i++) {
        if (i == 0 || fabs(tempPeakValues[i] - tempPeakValues[lastIdx]) > fil) {
            (*peaks)[lastIdx + 1] = tempPeaks[i];
            (*peakValues)[lastIdx + 1] = tempPeakValues[i];
            lastIdx++;
        }
    }

    // Free temporary arrays
    free(yp);
    free(tempPeaks);
    free(tempPeakValues);

    return lastIdx + 1;  // Return the number of peaks found
}
