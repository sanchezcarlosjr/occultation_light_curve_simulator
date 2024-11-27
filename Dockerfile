# Use the specified version of Arch Linux
FROM archlinux

# Update the system and install base packages
RUN pacman -Syu --noconfirm \
    base-devel \
    git \
    wget \
    cmake \
    gcc \
    make \
    gsl \
    hdf5

# Install FFTW3
RUN pacman -S --noconfirm fftw

# Clone the repository with submodules
WORKDIR /root
RUN git clone --recurse-submodules -j8 https://github.com/sanchezcarlosjr/occultation_light_curve_simulator.git

# Build the project using CMake
WORKDIR /root/occultation_light_curve_simulator

RUN cmake -B build . && cmake --build build

RUN cp build/bin/slc /usr/local/bin/

WORKDIR /data

ENTRYPOINT ["slc"]
