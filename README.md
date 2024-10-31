# NS-EOF: Navier Stokes - Ernst-Otto-Fischer Teaching Code

## General Information

This code has been tested under Ubuntu 22.04 or higher. Other Linux distributions will probably work.
However, we do not recommend using Windows or MacOS.
If you do not have Linux installed on your computer, please use WSL, Virtual Machine, Docker or similar.

## Dependencies

### Build Essentials

On Ubuntu, start with the installation of all required build essentials:
`apt install build-essential gcc g++ pkg-config`

### MPI (recommended OpenMPI)

* Under Ubuntu you can simply run `apt install libopenmpi-dev` to install OpenMPI.

### PETSc

* Please install [PETSc](https://petsc.org/release/) on your system.
* We recommend to use `apt install petsc-dev` to install.

### CMake

* Run `apt install cmake` to install CMake on your system.

## Docker

A prebuilt Docker image is available in [Dockerhub](https://hub.docker.com/r/tumi5/ns-eof).

To use the prebuilt image, we do recommended that you first clone NS-EOF, then navigate into the NS-EOF checkout
and run the Docker container interactively by mapping the NS-EOF directory into the container (here, we use `work` as our mapping point):

```shell
docker run -it -v ${PWD}:/work --rm --privileged tumi5/ns-eof /bin/bash
```

Navigate into the `work` directory and continue with the steps below.

## Build and Test

As build system configurator we use CMake. To compile the code execute the following commands in this directory:

* Create a build directory: `mkdir build`. You can also choose any other name for your build directory.
* Switch to this directory: `cd build`
* Run CMake: `cmake ..` (for an overview of all available options, use `ccmake ..`)
* For a `Debug` build, run `cmake .. -DCMAKE_BUILD_TYPE=Debug`
* Run Make: `make` (or `make -j` to compile with multiple cores).
* Run Tests: Some basic unit tests have been implemented (`make test`). Feel free to add your own test cases inside the `Tests` folder.

## Running a Simulation

* Run the code in serial via `./NS-EOF-Runner path/to/your/configuration`
  * Example: `./NS-EOF-Runner ../ExampleCases/Cavity2D.xml`
* Run the code in parallel via `mpirun -np nproc ./NS-EOF-Runner path/to/your/configuration`
  * Example: `mpirun -np 4 ./NS-EOF-Runner ../ExampleCases/Cavity2DParallel.xml`

## Adding New Source Files

You can add new source files by just creating them somewhere within the `Source` folder. CMake automatically detects these files and adds them to the build.

## Development Hints & FAQ

### It does not compile and everything seems fine?

Make sure to use `make clean` before you use `make`. Sometimes there are build artifacts from previous build processes that spoil your current compilation process. `make clean` takes care of deleting everything that should not be there and allows the compiler to start from scratch.

Sometimes it is also helpful to delete the `build` folder and create a new one, following the steps from the compilation section above.

### How can I see all the compiler flags the generated Makefile is using?

Instead of using `make`, run `VERBOSE=1 make`. You can also run `make -n` to invoke a dry run where you see what the Makefile would do in case of compilation.

### How can I see the test output?

Instead of using `make test`, run `ctest --verbose`.

### WARNING! There are options you set that were not used! There are 2 unused database options.

Just ignore these warnings. If you run the code in serial, PETSc does not require the two parameters for parallel runs and vice versa.
