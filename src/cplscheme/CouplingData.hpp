#pragma once

#include <Eigen/Core>
#include "mesh/Data.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace cplscheme {

struct CouplingData {  // @todo: should be a class from a design standpoint. See https://github.com/precice/precice/pull/865#discussion_r495825098
  using DataMatrix = Eigen::MatrixXd;

  Eigen::VectorXd &values()
  {
    data->values();
  }

  void extrapolateData(double t)
  {
    data->values() = waveform.sample(t);
  }

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &values() const
  {
    return couplingValues;
  }

  /// Data values of previous iteration (1st col) and previous time windows.
  DataMatrix oldValues;

  mesh::PtrData data;

  mesh::PtrMesh mesh;

  Eigen::VectorXd couplingValues;

  ///  True, if the data values if this CouplingData requires to be initialized by a participant.
  bool requiresInitialization;

  int getDimensions()
  {
    PRECICE_ASSERT(data != nullptr);
    return data->getDimensions();
  }

  int getSize()
  {
    PRECICE_ASSERT(couplingValues.size() == data->values().size());
    return couplingValues.size();
  }

  /**
   * @brief Default constructor, not to be used!
   *
   * Necessary when compiler creates template code for std::map::operator[].
   */
  CouplingData()
  {
    PRECICE_ASSERT(false);
  }

  CouplingData(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      bool          requiresInitialization)
      : data(data),
        mesh(mesh),
        requiresInitialization(requiresInitialization)
  {
    PRECICE_ASSERT(data != nullptr);
    PRECICE_ASSERT(mesh != nullptr);
    PRECICE_ASSERT(mesh.use_count() > 0);
    std::cout << "CouplingData::CouplingData calls copyDataFromMesh:" << std::endl;
    copyDataFromMesh();
  }

  void copyDataFromMesh(){
    PRECICE_ASSERT(data != nullptr);
    std::cout << "copy:\n" << data->values() << "\n---" << std::endl;
    if(data->values().size() > 0)
      couplingValues = data->values();
    std::cout << "stored:\n" << couplingValues << "\n---" << std::endl;
  }

  void copyDataToMesh(){
    PRECICE_ASSERT(data != nullptr);
    std::cout << "copy:\n" << couplingValues << "\n---" << std::endl;
    if(couplingValues.size() > 0)
      data->values() = couplingValues;
  }
};

} // namespace cplscheme
} // namespace precice
