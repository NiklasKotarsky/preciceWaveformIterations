#include "testing/ParallelCouplingSchemeFixture.hpp"

namespace precice {
namespace testing {

bool ParallelCouplingSchemeFixture::isImplicitCouplingScheme(cplscheme::MultiCouplingScheme &cplscheme)
{
  return cplscheme.isImplicitCouplingScheme();
}

cplscheme::CouplingData *ParallelCouplingSchemeFixture::getReceiveData(cplscheme::MultiCouplingScheme &cplscheme, int dataID)
{
  return cplscheme.getReceiveData(dataID);
}

cplscheme::CouplingData *ParallelCouplingSchemeFixture::getSendData(cplscheme::MultiCouplingScheme &cplscheme, int dataID)
{
  return cplscheme.getSendData(dataID);
}
} // namespace testing
} // namespace precice
