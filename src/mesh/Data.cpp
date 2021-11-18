#include "Data.hpp"
#include <algorithm>
#include <utility>

#include "precice/types.hpp"
#include "time/Waveform.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace mesh {

size_t Data::_dataCount = 0;

Data::Data()
    : _name(""),
      _id(-1),
      _dimensions(0)
{
  PRECICE_ASSERT(false);
}

Data::Data(
    std::string name,
    DataID      id,
    int         dimensions)
    : _values(),
      _name(std::move(name)),
      _id(id),
      _dimensions(dimensions)
{
  PRECICE_ASSERT(dimensions > 0, dimensions);
  _ptrWaveform = time::PtrWaveform(new time::Waveform(Data::EXTRAPOLATION_ORDER, Data::INTERPOLATION_ORDER));
  _dataCount++;
}

Data::~Data()
{
  _dataCount--;
}

Eigen::VectorXd &Data::values()
{
  return _values;
}

const Eigen::VectorXd &Data::values() const
{
  return _values;
}

const std::string &Data::getName() const
{
  return _name;
}

DataID Data::getID() const
{
  return _id;
}

void Data::toZero()
{
  auto begin = _values.data();
  auto end   = begin + _values.size();
  std::fill(begin, end, 0.0);
}

int Data::getDimensions() const
{
  return _dimensions;
}

size_t Data::getDataCount()
{
  return _dataCount;
}

void Data::resetDataCount()
{
  _dataCount = 0;
}

time::PtrWaveform Data::waveform()
{
  return _ptrWaveform;
}

void Data::setExtrapolationOrder(int extrapolationOrder)
{
  PRECICE_CHECK((extrapolationOrder == 0) || (extrapolationOrder == 1) || (extrapolationOrder == 2),
                "Extrapolation order has to be  0, 1, or 2.");
  _ptrWaveform->setExtrapolationOrder(extrapolationOrder);
}

} // namespace mesh
} // namespace precice
