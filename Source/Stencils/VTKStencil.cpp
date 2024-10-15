#include "StdAfx.hpp"

#include "VTKStencil.hpp"

Stencils::VTKStencil::VTKStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters),
  written_(false),
  prefix_(parameters.vtk.prefix) {

  if (parameters_.parallel.rank == 0) {
    const std::string outputFolder = "Output/" + prefix_;

#ifdef _MSC_VER
    const int success = _mkdir(outputFolder.data());
    if (success != 0) {
      if (errno != EEXIST) {
        spdlog::error("Cannot create folder {}", prefix_);
        throw std::runtime_error("Error while creating folder for VTK output");
      }
    }
#else
    char str_mkdir[256] = "mkdir -p ";
    strcat(str_mkdir, outputFolder.data());
    const int success = system(str_mkdir);
    if (success != 0) {
      spdlog::error("Cannot create folder {}", prefix_);
      throw std::runtime_error("Error while creating folder for VTK output");
    }
#endif
  }

  std::ofstream gitignore;
  gitignore.open(prefix_ + "/.gitignore");
  gitignore << "*";
  gitignore.close();
}

void Stencils::VTKStencil::writeVTKHeader(std::ostream& file) const {
  file << "# vtk DataFile Version 2.0" << std::endl << "NS-EOF" << std::endl << "ASCII" << std::endl << std::endl;
}

void Stencils::VTKStencil::writePoints(std::ostream& file, RealType simulationTime) const {
  // Number of points in every direction
  int px = parameters_.parallel.localSize[0] + 1;
  int py = parameters_.parallel.localSize[1] + 1;
  int pz = parameters_.geometry.dim == 2 ? 1 : parameters_.parallel.localSize[2] + 1;

  std::string grid;
  char        buffer[256];

  grid.reserve((file.precision() + 6) * px * py * pz * 3);

  sprintf(
    buffer,
    "DATASET STRUCTURED_GRID\nFIELD FieldData 1\nTIME 1 1 double\n%f\nDIMENSIONS %d %d %d\nPOINTS %d float\n",
    simulationTime,
    px,
    py,
    pz,
    px * py * pz
  );
  grid.append(buffer);

  if (parameters_.geometry.dim == 3) {
    for (int k = 2; k < 2 + pz; k++) {
      for (int j = 2; j < 2 + py; j++) {
        for (int i = 2; i < 2 + px; i++) {
          // Determine positions of grid points at lower/left/front corner of the respective grid cell (i,j,k) -> use
          // meshsize-ptr
          sprintf(
            buffer,
            "%f %f %f\n",
            parameters_.meshsize->getPosX(i, j, k),
            parameters_.meshsize->getPosY(i, j, k),
            parameters_.meshsize->getPosZ(i, j, k)
          );
          grid.append(buffer);
        }
      }
    }
  } else {
    for (int j = 2; j < 2 + py; j++) {
      for (int i = 2; i < 2 + px; i++) {
        sprintf(buffer, "%f %f 0.0\n", parameters_.meshsize->getPosX(i, j), parameters_.meshsize->getPosY(i, j));
        grid.append(buffer);
      }
    }
  }
  grid.append("\n");
  file << grid;
}

void Stencils::VTKStencil::apply(FlowField& flowField, int i, int j) {
  ASSERTION(FieldStencil<FlowField>::parameters_.geometry.dim == 2);

  RealType pressure    = 0.0;
  RealType velocity[2] = {0.0, 0.0};

  if ((flowField.getFlags().getValue(i, j) & OBSTACLE_SELF) == 0) {
    flowField.getPressureAndVelocity(pressure, velocity, i, j);

    pressureStream_ << pressure << std::endl;
    velocityStream_ << velocity[0] << " " << velocity[1] << " 0" << std::endl;
  } else {
    pressureStream_ << "0.0" << std::endl;
    velocityStream_ << "0.0 0.0 0.0" << std::endl;
  }
}

void Stencils::VTKStencil::apply(FlowField& flowField, int i, int j, int k) {
  ASSERTION(FieldStencil<FlowField>::parameters_.geometry.dim == 3);

  RealType pressure    = 0.0;
  RealType velocity[3] = {0.0, 0.0, 0.0};

  if ((flowField.getFlags().getValue(i, j, k) & OBSTACLE_SELF) == 0) {
    flowField.getPressureAndVelocity(pressure, velocity, i, j, k);

    pressureStream_ << pressure << std::endl;
    velocityStream_ << velocity[0] << " " << velocity[1] << " " << velocity[2] << std::endl;
  } else {
    pressureStream_ << "0.0" << std::endl;
    velocityStream_ << "0.0 0.0 0.0" << std::endl;
  }
}

void Stencils::VTKStencil::openFile(int timeStep, RealType simulationTime) {
  written_ = false;
  std::stringstream namestream;
  std::string       name;
  namestream.precision(4);
  namestream
    << "Output"
    << "/" << prefix_ << "/" << prefix_ << "." << parameters_.parallel.rank << "." << timeStep << ".vtk";
  name = namestream.str();
  ofile_.open(name.c_str());
  namestream.str("");

  writeVTKHeader(ofile_);
  writePoints(ofile_, simulationTime);
}

void Stencils::VTKStencil::write(FlowField& flowField, int timeStep, RealType simulationTime) {
  openFile(timeStep, simulationTime);

  if (FieldStencil<FlowField>::parameters_.geometry.dim == 2) {
    // Write pressure
    ofile_
      << "CELL_DATA " << flowField.getNx() * flowField.getNy() << std::endl
      << "SCALARS pressure float 1" << std::endl
      << "LOOKUP_TABLE default" << std::endl;
    ofile_ << pressureStream_.str() << std::endl;
    pressureStream_.str("");

    // Write velocity
    ofile_ << "VECTORS velocity float" << std::endl;
    ofile_ << velocityStream_.str() << std::endl;
    velocityStream_.str("");
  }

  if (FieldStencil<FlowField>::parameters_.geometry.dim == 3) {
    // Write pressure
    ofile_
      << "CELL_DATA " << flowField.getNx() * flowField.getNy() * flowField.getNz() << std::endl
      << "SCALARS pressure float 1" << std::endl
      << "LOOKUP_TABLE default" << std::endl;
    ofile_ << pressureStream_.str() << std::endl;
    pressureStream_.str("");

    // Write velocity
    ofile_ << "VECTORS velocity float" << std::endl;
    ofile_ << velocityStream_.str() << std::endl;
    velocityStream_.str("");
  }

  written_ = true;
  closeFile();
}

void Stencils::VTKStencil::closeFile() { ofile_.close(); }
