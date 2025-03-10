# Contributing
Contributions are welcomed. Open a pull-request or an issue.

To contribute to this repository, you must install CMake, a C99-compatible compiler, GSL, and FFTW. Additionally, modifications to the command line interface (CLI) may necessitate the installation of gengetopt.
As a good practice, we provide you with a test suite through Unity. This is a ready-to-go Git Submodule.

# Installation
1. Clone the repository
```
git clone --recurse-submodules -j8 https://github.com/sanchezcarlosjr/occultation_light_curve_simulator.git
```

2. Install global dependencies.

a. Install FFTW from this repository.
```
  cd external/fftw-3.3.10
  ./configure
  make
  make install
```

b. Install GSL, Unity, and hdf5 from the official docs.

2. Cmake
```
cmake -B build && cmake --build build && cd build/bin/
```

3. Pull data to test purposes.
```
git lfs pull
```

# Task runner
https://taskfile.dev/


## Code of conduct
This project adheres to the [Open Code of Conduct][code-of-conduct]. By participating, you are expected to honor this code.

[code-of-conduct]: https://github.com/spotify/code-of-conduct/blob/master/code-of-conduct.md

