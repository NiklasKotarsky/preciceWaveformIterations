#include "SerialCouplingScheme.hpp"
#include <cmath>
#include <memory>
#include <ostream>
#include <utility>

#include <vector>
#include "acceleration/Acceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/BiCouplingScheme.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/M2N.hpp"
#include "math/differences.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace cplscheme {

SerialCouplingScheme::SerialCouplingScheme(
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
    int                           extrapolationOrder,
    bool                          experimentalAPI)
    : BiCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, firstParticipant,
                       secondParticipant, localParticipant, std::move(m2n), maxIterations, cplMode, dtMethod, extrapolationOrder),
      _experimental(experimentalAPI)
{
  if (dtMethod == constants::FIRST_PARTICIPANT_SETS_TIME_WINDOW_SIZE) {
    if (doesFirstStep()) {
      _participantSetsTimeWindowSize = true;
      setTimeWindowSize(UNDEFINED_TIME_WINDOW_SIZE);
    } else {
      _participantReceivesTimeWindowSize = true;
    }
  }
}

void SerialCouplingScheme::receiveAndSetTimeWindowSize()
{
  PRECICE_TRACE();
  if (_participantReceivesTimeWindowSize) {
    double dt = UNDEFINED_TIME_WINDOW_SIZE;
    getM2N()->receive(dt);
    PRECICE_DEBUG("Received time window size of {}.", dt);
    PRECICE_ASSERT(not math::equals(dt, UNDEFINED_TIME_WINDOW_SIZE));
    PRECICE_ASSERT(not doesFirstStep(), "Only second participant can receive time window size.");
    setTimeWindowSize(dt);
  }
}

void SerialCouplingScheme::initializeImplementation()
{
  // determine whether initial data needs to be communicated
  determineInitialSend(getSendData());
  determineInitialReceive(getReceiveData());

  // If the second participant initializes data, the first receive for the
  // second participant is done in initializeData() instead of initialize().
  if (not doesFirstStep() && not sendsInitializedData() && isCouplingOngoing()) {
    PRECICE_DEBUG("Receiving data");
    receiveAndSetTimeWindowSize();
    receiveData(getM2N(), getReceiveData());
    checkDataHasBeenReceived();
  }
}

void SerialCouplingScheme::exchangeInitialData()
{
  if (doesFirstStep()) {
    if (sendsInitializedData()) {
      PRECICE_ASSERT(isImplicitCouplingScheme() && _experimental, "First participant cannot send data during initialization, if experimental=\"false\".");
      sendData(getM2N(), getSendData());
    }
    if (receivesInitializedData()) {
      receiveData(getM2N(), getReceiveData());
      checkInitialDataHasBeenReceived();
    }
  } else { // second participant
    if (receivesInitializedData()) {
      PRECICE_ASSERT(isImplicitCouplingScheme() && _experimental, "Only first participant can receive data during initialization, if experimental=\"false\".");
      receiveData(getM2N(), getReceiveData()); // @todo has to be stored in waveform, if this gets triggered. Otherwise zero.
      checkInitialDataHasBeenReceived();
      moveToNextWindow(); // @todo need to do much more here: not only extrapolation, but also persist current data in some buffer.
    }
    if (sendsInitializedData()) {
      // The second participant sends the initialized data to the first participant
      // here, which receives the data on call of initialize().
      sendData(getM2N(), getSendData());
      receiveAndSetTimeWindowSize();
      // This receive replaces the receive in initialize().
      receiveData(getM2N(), getReceiveData());
      if (not receivesInitializedData()) {
        checkDataHasBeenReceived();
      }
    }
  }
}

bool SerialCouplingScheme::exchangeDataAndAccelerate()
{
  bool convergence = true;

  if (doesFirstStep()) { // first participant
    PRECICE_DEBUG("Sending data...");
    if (_participantSetsTimeWindowSize) {
      PRECICE_DEBUG("sending time window size of {}", getComputedTimeWindowPart());
      getM2N()->send(getComputedTimeWindowPart());
    }
    sendData(getM2N(), getSendData());
    if (isImplicitCouplingScheme()) {
      convergence = receiveConvergence(getM2N());
    }
    PRECICE_DEBUG("Receiving data...");
    receiveData(getM2N(), getReceiveData());
    checkDataHasBeenReceived();
  } else { // second participant
    if (isImplicitCouplingScheme()) {
      PRECICE_DEBUG("Test Convergence and accelerate...");
      convergence = doImplicitStep();
      sendConvergence(getM2N(), convergence);
    }
    PRECICE_DEBUG("Sending data...");
    sendData(getM2N(), getSendData());
    // the second participant does not want new data in the last iteration of the last time window
    if (isCouplingOngoing() || (isImplicitCouplingScheme() && not convergence)) {
      if (_participantReceivesTimeWindowSize) {
        receiveAndSetTimeWindowSize();
      }
      PRECICE_DEBUG("Receiving data...");
      receiveData(getM2N(), getReceiveData());
      checkDataHasBeenReceived();
    }
  }

  return convergence;
}

} // namespace cplscheme
} // namespace precice
