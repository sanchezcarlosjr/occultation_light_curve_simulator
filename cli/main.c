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

    double D = calcPlano(args_info.d_arg, args_info.lamb_arg, args_info.ua_arg);
    hdf5->writeScalar(hdf5, "total_plane_size", D);
    gsl_matrix* O1 = pupilCO(args_info.M_arg, D, args_info.d_arg);
    hdf5->writeDataset(hdf5, "pupil_obstruction", O1);

    hdf5->free(hdf5);

    gsl_matrix_free(O1);

    cmdline_parser_free(&args_info);
    return 0;
}
