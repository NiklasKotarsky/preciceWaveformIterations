#pragma once

#include "cplscheme/MultiCouplingScheme.hpp"

namespace precice {
namespace testing {
/*
 * @brief A fixture that is used to access private functions of the ParallelCouplingScheme class.
 *
 * The fixture can be used to call private functions for individual testing.
 */
struct ParallelCouplingSchemeFixture {
  static bool isImplicitCouplingScheme(cplscheme::MultiCouplingScheme &cplscheme);

  static cplscheme::CouplingData *getReceiveData(cplscheme::MultiCouplingScheme &cplscheme, int dataID);

  static cplscheme::CouplingData *getSendData(cplscheme::MultiCouplingScheme &cplscheme, int dataID);
};
} // namespace testing
} // namespace precice
