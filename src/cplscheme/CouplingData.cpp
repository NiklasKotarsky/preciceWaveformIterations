#include "cplscheme/CouplingData.hpp"

#include <utility>

#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice::cplscheme {

CouplingData::CouplingData(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          requiresInitialization,
    int           extrapolationOrder)
    : requiresInitialization(requiresInitialization),
      _data(std::move(data)),
      _mesh(std::move(mesh)),
      _extrapolation(extrapolationOrder)
{
  PRECICE_ASSERT(_data != nullptr);
  _previousIteration = Eigen::VectorXd::Zero(getSize());
  PRECICE_ASSERT(_mesh != nullptr);
  PRECICE_ASSERT(_mesh.use_count() > 0);
}

int CouplingData::getDimensions() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->getDimensions();
}

int CouplingData::getSize() const
{
  PRECICE_ASSERT(_data != nullptr);
  return values().size();
}

Eigen::VectorXd &CouplingData::values()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->values();
}

const Eigen::VectorXd &CouplingData::values() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->values();
}

Eigen::MatrixXd &CouplingData::gradientValues()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->gradientValues();
}

const Eigen::MatrixXd &CouplingData::gradientValues() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->gradientValues();
}

bool CouplingData::hasGradient() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->hasGradient();
}

int CouplingData::meshDimensions() const
{
  return _mesh->getDimensions();
}

void CouplingData::storeIteration()
{
  _previousIteration = this->values();
}

const Eigen::VectorXd CouplingData::previousIteration() const
{
  return _previousIteration;
}

int CouplingData::getPreviousIterationSize() const
{
  return previousIteration().size();
}

int CouplingData::getMeshID()
{
  return _mesh->getID();
}

int CouplingData::getDataID()
{
  return _data->getID();
}

std::string CouplingData::getDataName()
{
  return _data->getName();
}

std::vector<int> CouplingData::getVertexOffsets()
{
  return _mesh->getVertexOffsets();
}

void CouplingData::initializeExtrapolation()
{
  _extrapolation.initialize(getSize());
  storeIteration();
}

void CouplingData::moveToNextWindow()
{
  _extrapolation.moveToNextWindow();
  values() = _extrapolation.getInitialGuess();
}

void CouplingData::storeExtrapolationData()
{
  _extrapolation.store(values());
}

Eigen::VectorXd CouplingData::getStoredTimesAscending()
{
  return _timeStepsStorage.getTimes();
}

void CouplingData::clearTimeStepsStorage()
{
  _timeStepsStorage.clear();
}

void CouplingData::storeDataAtTime(Eigen::VectorXd data, double relativeDt, bool override)
{
  // @todo move logic inside time::Storage? Similar to what also happens in time::Waveform.
  double endOfWindow = 1.0;
  PRECICE_ASSERT(relativeDt > 0.0);
  PRECICE_ASSERT(relativeDt <= endOfWindow);
  if (override && math::greaterEquals(_timeStepsStorage.maxStoredNormalizedDt(), relativeDt)) {
    PRECICE_ASSERT(math::equals(endOfWindow, _timeStepsStorage.maxStoredNormalizedDt()), _timeStepsStorage.maxStoredNormalizedDt());
    _timeStepsStorage.clear();
  } else {
    PRECICE_ASSERT(math::greater(relativeDt, _timeStepsStorage.maxStoredNormalizedDt()), relativeDt, _timeStepsStorage.maxStoredNormalizedDt());
  }
  _timeStepsStorage.setValueAtTime(relativeDt, data);
}

void CouplingData::overrideDataAtEndWindowTime(Eigen::VectorXd data)
{
  _timeStepsStorage.overrideDataAtEndWindowTime(data);
}

Eigen::VectorXd CouplingData::getDataAtTime(double relativeDt)
{
  return _timeStepsStorage.getValueAtTime(relativeDt);
}

Eigen::VectorXd CouplingData::getSerialized()
{
  int  nValues        = getSize();
  int  nTimeSteps     = _timeStepsStorage.nTimes();
  auto serializedData = Eigen::VectorXd(nTimeSteps * nValues);
  auto timesAndValues = _timeStepsStorage.getTimesAndValues();
  auto times          = timesAndValues.first;
  auto values         = timesAndValues.second;

  for (int timeId = 0; timeId < nTimeSteps; timeId++) {
    auto time  = times(timeId);
    auto slice = values.col(timeId);
    for (int valueId = 0; valueId < nValues; valueId++) {
      serializedData(valueId * nTimeSteps + timeId) = slice(valueId);
    }
  }
  return serializedData;
}

void CouplingData::storeFromSerialized(Eigen::VectorXd timesAscending, Eigen::VectorXd serializedData)
{
  PRECICE_ASSERT(timesAscending.size() * getSize() == serializedData.size());
  for (int timeId = 0; timeId < timesAscending.size(); timeId++) {
    auto slice = Eigen::VectorXd(getSize());
    for (int valueId = 0; valueId < slice.size(); valueId++) {
      slice(valueId) = serializedData(valueId * timesAscending.size() + timeId);
    }
    auto time = timesAscending(timeId);
    PRECICE_ASSERT(time > 0.0 && time <= 1.0); // time <= 0 or time > 1 is not allowed.
    this->storeDataAtTime(slice, time);
  }
  this->values() = this->getDataAtTime(_timeStepsStorage.maxStoredNormalizedDt()); // store data in values to make this non-breaking.
}

} // namespace precice::cplscheme
