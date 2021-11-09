#include <Eigen/Core>
#include "logging/Logger.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "time/config/TimeConfiguration.hpp"

using namespace precice;
using namespace precice::time;

BOOST_AUTO_TEST_SUITE(TimeTests)

BOOST_AUTO_TEST_SUITE(Configuration)

BOOST_AUTO_TEST_CASE(dummy)
{
  PRECICE_TEST(1_rank);
  BOOST_TEST(false);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
