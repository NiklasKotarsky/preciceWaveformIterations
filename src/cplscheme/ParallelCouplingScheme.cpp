#include "ParallelCouplingScheme.hpp"

#include <utility>

#include "cplscheme/BiCouplingScheme.hpp"
#include "logging/LogMacros.hpp"

namespace precice {
namespace cplscheme {

ParallelCouplingScheme::ParallelCouplingScheme(
    double                        maxTime,
    int                           maxTimeWindows,
    double                        timeWindowSize,
    int                           validDigits,
    const std::string &           firstParticipant,
    const std::string &           secondParticipant,
    const std::string &           localParticipant,
    m2n::PtrM2N                   m2n,
    constants::TimesteppingMethod dtMethod,
    CouplingMode                  cplMode,
    int                           maxIterations,
    int                           extrapolationOrder)
    : BiCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, firstParticipant,
                       secondParticipant, localParticipant, std::move(m2n), maxIterations, cplMode, dtMethod, extrapolationOrder) {}

bool ParallelCouplingScheme::exchangeDataAndAccelerate()
{
  bool convergence = true;

  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Sending data...");
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
    PRECICE_DEBUG("Receiving data...");
    if (isImplicitCouplingScheme()) {
      convergence = receiveConvergence(_m2ns[_otherParticipant]);
    }
    for (auto &receiveExchange : _receiveDataVector) {
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();
  } else { // second participant
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
  }
  return convergence;
}

typedef std::map<int, PtrCouplingData> DataMap;

const DataMap ParallelCouplingScheme::getAccelerationData()
{
  PRECICE_ASSERT(!doesFirstStep(), "Only the second participant should do the acceleration.");
  DataMap accelerationData;
  for (auto &data : allCouplingData()) {
    PRECICE_ASSERT(accelerationData.count(data->getDataID()) == 0);
    accelerationData[data->getDataID()] = data;
  }
  return accelerationData;
}

} // namespace cplscheme
} // namespace precice
