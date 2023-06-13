#include <Eigen/Core>
#include <memory>
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "m2n/M2N.hpp"
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/IntraComm.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(TestingTests) // Use name of the module, e.g. subdirectory below src/, suffixed with Tests

BOOST_AUTO_TEST_SUITE(Examples) // If your file contains multiple tests, put them in a test suite

/// This test runs on a single processor.
BOOST_AUTO_TEST_CASE(SingleProcessor)
{
  PRECICE_TEST(1_rank);
  /* Do not use DEBUG, TRACE, INFO calls inside tests, if you need to log a message use
     BOOST_TEST_MESSAGE("I have done that " << whatIHaveDone);

     From Boost Test Documentation:
     "Messages generated by this tool do not appear in test log output with default value of the
     active log level threshold. For these messages to appear the active log level threshold has to
     be set to a value below or equal to "message"."

     In order to get this output to the terminal, use testprecice with --log-level=message.
  */

  BOOST_TEST(0 == 0); // Always use BOOST_TEST

  Eigen::Vector3d one(1, 2, 3);
  Eigen::Vector3d two(1, 2, 3);
  // Use testing::equals instead of math::equals when comparing in tests.
  // This gives you a report which coordinates fail to compare.
  BOOST_TEST(testing::equals(one, two));
}

/// Test with a modified numerical tolerance
BOOST_AUTO_TEST_CASE(NumericalTolerance,
                     *boost::unit_test::tolerance(1e-4))
{
  PRECICE_TEST(1_rank);
  // Default tolerance is 1e-9, it can be changed for the entire case or even suite
  // using the decorator above
  BOOST_TEST(1.0 == 1.0001);

  // Or on a per test basis
  BOOST_TEST(1.0 == 1.01, boost::test_tools::tolerance(0.1));
}

/// Use testing::Deleted to unconditionally delete the test
BOOST_AUTO_TEST_CASE(Deleted,
                     *testing::Deleted())
{
  PRECICE_TEST(1_rank);
  BOOST_TEST(false);
}

/// Test that requires 4 processors.
/*
 * If less than 4 procs are available the test will fail.
 */
BOOST_AUTO_TEST_CASE(FourProcTests)
{
  PRECICE_TEST(4_ranks);
  // Don't copy over that line, it's for testing the example
  BOOST_TEST(context.size == 4);
  BOOST_TEST(context.hasSize(4));
}

/// Test that runs on 2 processors.
BOOST_AUTO_TEST_CASE(TwoProcTests)
{
  PRECICE_TEST(2_ranks);

  // Put your test code here
  BOOST_TEST(context.hasSize(2));
}

#ifndef PRECICE_NO_MPI
/// Test that requires 4 processors and an intra-participant communication
/*
 * For some primary tests, you might need an intra-participant communication. This example shows how to set one up.
 * Please note: Such tests always need to be excluded for compilation without MPI (PRECICE_NO_MPI).
 */
BOOST_AUTO_TEST_CASE(FourProcTestsWithPrimaryCommmunication)
{
  // The short syntax won't work here. You have to name the context
  PRECICE_TEST(""_on(4_ranks).setupIntraComm())
  // In this test you can use an intra-participant communication, here is an example how:
  BOOST_TEST(context.hasSize(4));
  BOOST_TEST(utils::IntraComm::getCommunication()->isConnected());
}

/// Test that requires 2 participants "A" on 1 rank and "B" on 2 ranks
BOOST_AUTO_TEST_CASE(NamedContexts)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(2_ranks));

  if (context.isNamed("A")) {
    BOOST_TEST(context.hasSize(1));
  } else if (context.isNamed("B")) {
    BOOST_TEST(context.hasSize(2));
  } else {
    BOOST_TEST(false);
  }
}

/// Tests that requires an m2n communication
/*
 * For some primary tests, you might need an m2n communication (e.g. partition or cplscheme).
 * This example shows how to set up one. Call .connectPrimary() on the context and pass the participants to be connected.
 * M2N requires Events, thus you also need to list it as a requirement using Require::Events.
 * Please note: Such tests always need to be excluded for compilation without MPI (PRECICE_NO_MPI).
 */
BOOST_AUTO_TEST_CASE(TwoProcTestsWithM2NCommunication)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  BOOST_TEST(context.hasSize(1));
  BOOST_TEST(context.isRank(0));
  BOOST_TEST(context.isPrimary());

  auto m2n = context.connectPrimaryRanks("A", "B");

  //This is how you can access the m2n communication
  BOOST_TEST(m2n->getPrimaryRankCommunication()->isConnected());

  // Automatically finalizes Events
}

#endif // PRECICE_NO_MPI

BOOST_AUTO_TEST_CASE(TwoProcTestsWithPETSc)
{
  PRECICE_TEST(2_ranks, Require::PETSc); // implies Require::Events
  BOOST_TEST(context.hasSize(2));

  // Automatically finalizes PETSc and Events
}

/// Integration tests with two participants.
/*
 * For integration tests (tests that directly use the preCICE API), often, you need two participants
 * where each participant uses it own communicator, i.e. each participant should not see that he is
 * part of a test.
 * In this case, you can simply create the participants and create a Participant object.
 * The context-object is of type TestContext and provides access to the name of the current context and the rank and size of its communicator.
 */
BOOST_AUTO_TEST_CASE(IntegrationTestsWithTwoParticipants)
{
  PRECICE_TEST("Solid"_on(2_ranks), "Fluid"_on(2_ranks));

  if (context.isNamed("Solid")) {
    // This is the participant Solid
    BOOST_TEST(context.hasSize(2));

    // You can now create a Participant object for your first participant
    // You can use context.name, context.rank, context.size
  } else {
    // This is the participant Fluid
    BOOST_TEST(context.hasSize(2));

    // You can now create a Participant object for your second participant
    // You can use context.name, context.rank, context.size
  }
}

BOOST_AUTO_TEST_SUITE_END() // Examples
BOOST_AUTO_TEST_SUITE_END() // TestingTests
