#include "numpy.h"
#include "cmdline.h"
#include "diffraction.h"
#include "slcio.h"
#include <stdlib.h>

int main(int argc, char *argv[]) {
    struct gengetopt_args_info args_info;

    if (cmdline_parser(argc, argv, &args_info) != 0) {
        exit(1);
    }

    HDF5Wrapper* hdf5 = NewHDF5Wrapper(args_info.output_arg);
    if (hdf5 == NULL) {
        fprintf(stderr, "Failed to initialize HDF5Wrapper.\n");
        return 1;
    }

    hdf5->writeScalar(hdf5, "mesh_size_pixels",      (double) args_info.M_arg);
    hdf5->writeScalar(hdf5, "wavelength_meters",     args_info.lamb_arg);

    hdf5->writeScalar(hdf5, "earth_translation_speed_mps", (double) args_info.vE_arg);
    hdf5->writeScalar(hdf5, "body_speed_mps",        (double) args_info.vr_arg);
    hdf5->writeScalar(hdf5, "angle_degrees",         (double) args_info.ang_arg);
    hdf5->writeScalar(hdf5, "frames_per_second",     (double) args_info.fps_arg);
    hdf5->writeScalar(hdf5, "apparent_magnitude",    (double) args_info.mV_arg);
    hdf5->writeScalar(hdf5, "spectral_type_index",   (double) args_info.nEst_arg);
    hdf5->writeScalar(hdf5, "num_wavelengths",       (double) args_info.nLamb_arg);

    hdf5->writeScalar(hdf5, "object_diameter_meters",(double) args_info.d_arg);
    hdf5->writeScalar(hdf5, "distance_au",           args_info.ua_arg);
    hdf5->writeScalar(hdf5, "time_offset_pixels",    (double) args_info.toffset_arg);
    hdf5->writeScalar(hdf5, "reading_direction_deg", (double) args_info.T_arg);
    hdf5->writeScalar(hdf5, "impact_parameter_m",    (double) args_info.b_arg);

    hdf5->writeScalar(hdf5, "eccentricity",          (double) args_info.eccentricity_arg);

    double D = calcPlano(args_info.d_arg, args_info.lamb_arg, args_info.ua_arg);
    hdf5->writeScalar(hdf5, "total_plane_size", D);
    gsl_matrix* O1 = pupilCO(args_info.M_arg, D, args_info.d_arg);
    hdf5->writeDataset(hdf5, "pupil_obstruction", O1);


    double object_distance = astronomicalUnitsToMeters(args_info.ua_arg);
    double total_plane_size = D;
    gsl_matrix* monochromatic_diffraction_pattern = fresnel(O1, args_info.M_arg, total_plane_size, object_distance, args_info.lamb_arg);
    hdf5->writeDataset(hdf5, "monochromatic_diffraction_pattern", monochromatic_diffraction_pattern);

    gsl_matrix* chromatic_diffraction_pattern = spectra(O1, args_info.M_arg, total_plane_size, object_distance, args_info.nEst_arg, args_info.nLamb_arg);
    hdf5->writeDataset(hdf5, "chromatic_diffraction_pattern", chromatic_diffraction_pattern);

    


    hdf5->free(hdf5);

    gsl_matrix_free(O1);

    cmdline_parser_free(&args_info);
    return 0;
}
