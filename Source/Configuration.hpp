#pragma once

#include "Definitions.hpp"
#include "Parameters.hpp"

class Configuration {
private:
  int         dim_;
  std::string filename_;

public:
  Configuration();
  Configuration(const std::string& filename);
  ~Configuration() = default;

  void setFileName(const std::string& filename);
  void loadParameters(Parameters& parameters, const MPI_Comm& communicator = PETSC_COMM_WORLD);
};
