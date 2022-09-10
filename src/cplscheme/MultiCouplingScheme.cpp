#include "MultiCouplingScheme.hpp"
#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
#include <ostream>
#include <type_traits>
#include <utility>
#include "acceleration/Acceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace cplscheme {

MultiCouplingScheme::MultiCouplingScheme(
    double                             maxTime,
    int                                maxTimeWindows,
    double                             timeWindowSize,
    int                                validDigits,
    const std::string &                localParticipant,
    std::map<std::string, m2n::PtrM2N> m2ns,
    constants::TimesteppingMethod      dtMethod,
    CouplingMode                       cplMode,
    const std::string &                controller,
    int                                maxIterations,
    int                                extrapolationOrder)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, localParticipant, maxIterations, cplMode, dtMethod, extrapolationOrder),
      _controller(controller), _isController(controller == localParticipant)
{
  _m2ns = m2ns;
  if (_m2ns.size() > 1) {
    // @todo implement MultiCouplingScheme for explicit coupling
    PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  }
  // Controller participant never does the first step, because it is never the first participant
  setDoesFirstStep(!_isController);
  PRECICE_DEBUG("MultiCoupling scheme is created for {}.", localParticipant);
}

void MultiCouplingScheme::exchangeInitialData()
{
  if (_isController) {
    if (receivesInitializedData()) {
      for (auto &receiveExchange : _receiveDataVector) {
        receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
      }
      checkDataHasBeenReceived();
    }
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.first], sendExchange.second);
      }
    }
  } else {
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.first], sendExchange.second);
      }
    }
    if (receivesInitializedData()) {
      for (auto &receiveExchange : _receiveDataVector) {
        receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
      }
      checkDataHasBeenReceived();
    }
  }
  PRECICE_DEBUG("Initial data is exchanged in MultiCouplingScheme");
}

bool MultiCouplingScheme::exchangeDataAndAccelerate()
{
  PRECICE_DEBUG("Computed full length of iteration");

  bool convergence = true;

  if (_isController) {
    PRECICE_DEBUG("Receiving data...");
    for (auto &receiveExchange : _receiveDataVector) {
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Perform acceleration (only second participant)...");
      convergence = doImplicitStep();
      for (const auto &m2nPair : _m2ns) {
        sendConvergence(m2nPair.second, convergence);
      }
    }
    PRECICE_DEBUG("Sending data...");
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
  } else {
    PRECICE_DEBUG("Sending data...");
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
    PRECICE_DEBUG("Receiving data...");
    if (isImplicitCouplingScheme()) {
      convergence = receiveConvergence(_m2ns[_controller]);
    }
    for (auto &receiveExchange : _receiveDataVector) {
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();
  }
  return convergence;
}

typedef std::map<int, PtrCouplingData> DataMap;

const DataMap MultiCouplingScheme::getAccelerationData()
{
  DataMap accelerationData;
  for (auto &data : allCouplingData()) {
    PRECICE_ASSERT(accelerationData.count(data->getDataID()) == 0);
    accelerationData[data->getDataID()] = data;
  }
  return accelerationData;
}

} // namespace cplscheme
} // namespace precice
