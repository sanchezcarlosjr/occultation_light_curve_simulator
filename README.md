# Occultation light curve simulator

The Trans-Neptunian Automated Occultation Survey II (TAOS II) event simulator is a valuable tool for modeling and predicting occultation events caused by small objects in the Kuiper Belt and beyond. Integrating this simulator into machine learning (ML) pipelines can indeed provide significant insights. Generating synthetic data or simulating rare events, can help train ML models to better detect and analyze actual astronomical observations, enhancing the accuracy of object detection and characterization in space.

![Generate light curves and save it with HDF5](./assets/generate-light-curves-and-save-it-with-hdf5.png)

# Playground
Learn how to install and use the program through this [playground](https://colab.research.google.com/drive/1GCPLfTBvZLvwUEgk9O1yfWWH1MQAXUHs?usp=sharing).

# Installation
1. Install the shared libraries.

For example, on Arch Linux, you would run
```
sudo pacman -S fftw gsl hdf5
```

2. Download the [latest version](https://github.com/sanchezcarlosjr/occultation_light_curves/releases/latest/download/slc) from GitHub Releases, so you don't need to compile the repository; it just works.

# Usage
Start your simulation with default parameters and save your data into HDF5.

```bash
slc -o polymele.h5
```

## Available commands
```slc``` supports several options, each accessible through the ```slc``` command and through our library. For help on individual commands, add --help following the command name. The commands are available on [here](./cli/cli.ggo).


# HDF5 Viewers
https://myhdf5.hdfgroup.org/


# Arch Linux (AUR)

1. Install the AUR package
```
yay -S occultation_light_curve_simulator
```

# Contribution with Arch Linux

1. Clone the repository
```
git clone https://github.com/sanchezcarlosjr/occultation_light_curve_simulator.git
```

Build the package to ensure it works correctly.
2. Install with PKGBUILD
```
makepkg -csi
```
Explanation:
-c: Cleans up work files after build.
-s: Installs missing dependencies.
-i: Installs the package after building.


# Contribution
To contribute to this repository, you must install CMake, a C99-compatible compiler, GSL, and FFTW. Additionally, modifications to the command line interface (CLI) may necessitate the installation of gengetopt.
As a good practice, we provide you with a test suite through Unity. This is a ready-to-go Git Submodule.

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
