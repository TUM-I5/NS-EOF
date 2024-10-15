#include "StdAfx.hpp"

#include "DataStructures.hpp"

ScalarField::ScalarField(int Nx, int Ny):
  Field<RealType>(Nx, Ny, 1, 1) {

  initialize();
}

ScalarField::ScalarField(int Nx, int Ny, int Nz):
  Field<RealType>(Nx, Ny, Nz, 1) {

  initialize();
}

RealType& ScalarField::getScalar(int i, int j, int k) { return data_[index2array(i, j, k)]; }

void ScalarField::show(const std::string title) {
  std::cout << std::endl << "--- " << title << " ---" << std::endl;
  for (int k = 0; k < sizeZ_; k++) {
    for (int j = sizeY_ - 1; j > -1; j--) {
      for (int i = 0; i < sizeX_; i++) {
        std::cout << getScalar(i, j, k) << "\t";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void ScalarField::initialize() {
  for (int i = 0; i < size_; i++) {
    data_[i] = 0.0;
  }
}

VectorField::VectorField(int Nx, int Ny):
  Field<RealType>(Nx, Ny, 1, 2) {

  initialize();
}

VectorField::VectorField(int Nx, int Ny, int Nz):
  Field<RealType>(Nx, Ny, Nz, 3) {

  initialize();
}

RealType* VectorField::getVector(int i, int j, int k) { return &data_[index2array(i, j, k)]; }

void VectorField::show(const std::string title) {
  std::cout << std::endl << "--- " << title << " ---" << std::endl;
  std::cout << "Component 1" << std::endl;
  for (int k = 0; k < sizeZ_; k++) {
    for (int j = sizeY_ - 1; j > -1; j--) {
      for (int i = 0; i < sizeX_; i++) {
        std::cout << getVector(i, j, k)[0] << "\t";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  std::cout << "Component 2" << std::endl;
  for (int k = 0; k < sizeZ_; k++) {
    for (int j = sizeY_ - 1; j > -1; j--) {
      for (int i = 0; i < sizeX_; i++) {
        std::cout << getVector(i, j, k)[1] << "\t";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void VectorField::initialize() {
  for (int i = 0; i < size_; i++) {
    data_[i] = 0.0;
  }
}

IntScalarField::IntScalarField(int Nx, int Ny):
  Field<int>(Nx, Ny, 1, 1) {

  initialize();
}

IntScalarField::IntScalarField(int Nx, int Ny, int Nz):
  Field<int>(Nx, Ny, Nz, 1) {

  initialize();
}

void IntScalarField::initialize() {
  for (int i = 0; i < size_; i++) {
    data_[i] = 0;
  }
}

int& IntScalarField::getValue(int i, int j, int k) { return data_[index2array(i, j, k)]; }

void IntScalarField::show(const std::string title) {
  std::cout << std::endl << "--- " << title << " ---" << std::endl;
  for (int k = 0; k < sizeZ_; k++) {
    for (int j = sizeY_ - 1; j > -1; j--) {
      for (int i = 0; i < sizeX_; i++) {
        std::cout << getValue(i, j, k) << "\t";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}
