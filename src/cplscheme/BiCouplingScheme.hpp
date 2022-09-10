#pragma once

#include <string>
#include <vector>

#include "BaseCouplingScheme.hpp"
#include "cplscheme/Constants.hpp"
#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace cplscheme {
class CouplingData;

/**
 * @brief Abstract base class for coupling schemes with two participants.
 *
 * ! General description
 * A BiCouplingScheme couples two participants. It is a specialization of
 * BaseCouplingScheme.
 *
 */
class BiCouplingScheme : public BaseCouplingScheme {

public:
  BiCouplingScheme(
      double                        maxTime,
      int                           maxTimeWindows,
      double                        timeWindowSize,
      int                           validDigits,
      std::string                   firstParticipant,
      std::string                   secondParticipant,
      const std::string &           localParticipant,
      m2n::PtrM2N                   m2n,
      int                           maxIterations,
      CouplingMode                  cplMode,
      constants::TimesteppingMethod dtMethod,
      int                           extrapolationOrder);

  /**
   * @returns true, if coupling scheme has sendData with given DataID
   */
  bool hasSendData(DataID dataID)
  {
    return getSendData(dataID) != nullptr;
  }

protected:
  /// name of the other participant
  std::string _otherParticipant = "unknown";

  /**
   * @brief Exchanges data, if it has to be initialized.
   */
  void exchangeInitialData() override final;

private:
  mutable logging::Logger _log{"cplscheme::BiCouplingScheme"};

  /// First participant name.
  std::string _firstParticipant = "unknown";

  /// Second participant name.
  std::string _secondParticipant = "unknown";
};

} // namespace cplscheme
} // namespace precice
