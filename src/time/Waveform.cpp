#include "time/Waveform.hpp"
#include <algorithm>
#include <unsupported/Eigen/Splines>
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "math/differences.hpp"
#include "mesh/Data.hpp"
#include "time/Time.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice::time {

Waveform::Waveform(const int degree)
    : _degree(degree)
{
  PRECICE_ASSERT(Time::MIN_WAVEFORM_DEGREE <= _degree && _degree <= Time::MAX_WAVEFORM_DEGREE);
}

int Waveform::getDegree() const
{
  return _degree;
}

time::Storage &Waveform::timeStepsStorage()
{
  return _timeStepsStorage;
}

const time::Storage &Waveform::timeStepsStorage() const
{
  return _timeStepsStorage;
}

Eigen::VectorXd Waveform::sample(double normalizedDt) const
{

  PRECICE_ASSERT(math::equals(this->_timeStepsStorage.maxStoredNormalizedDt(), time::Storage::WINDOW_END), this->_timeStepsStorage.maxStoredNormalizedDt()); // sampling is only allowed, if a window is complete.

  return this->_timeStepsStorage.sampleAt(double normalizedDt)
}

} // namespace precice::time
