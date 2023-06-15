#ifndef PRECICE_NO_MPI

#include <precice/SolverInterface.hpp>
#include <vector>
#include "testing/Testing.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief Test to ensure that the acceleration is applied to every timestep.
 *
 *
 */
BOOST_AUTO_TEST_CASE(WaveformSubcyclingWithConstantAcceleration)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank));

  using Eigen::Vector2d;

  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
  Vector2d                 vertex{0.0, 0.0};

  typedef double (*DataFunction)(double);
  DataFunction dataOneFunction = [](double t) -> double {
    return (double) (2 + t);
  };

  if (context.isNamed("A")) {
    const precice::MeshID meshID   = interface.getMeshID("A-Mesh");
    int                   vertexID = interface.setMeshVertex(meshID, vertex.data());
    int                   dataID   = interface.getDataID("Data", meshID);

    if (interface.requiresInitialData()) {
      double writeData = dataOneFunction(0);
      interface.writeScalarData(dataID, vertexID, writeData);
    }

    double maxDt           = interface.initialize();
    double dt              = maxDt / 5; // 5 substeps per iteration
    double time            = 0;
    double checkpoint_time = 0;
    while (interface.isCouplingOngoing()) {

      if (interface.requiresWritingCheckpoint())
        checkpoint_time = time;

      interface.writeScalarData(dataID, vertexID, dataOneFunction(time));
      time += dt;

      interface.advance(dt);
      if (interface.requiresReadingCheckpoint())
        time = checkpoint_time;
    }

    interface.finalize();

  } else {

    BOOST_TEST(context.isNamed("B"));
    const precice::MeshID meshID          = interface.getMeshID("B-Mesh");
    int                   vertexID        = interface.setMeshVertex(meshID, vertex.data());
    int                   dataID          = interface.getDataID("Data", meshID);
    double                maxDt           = interface.initialize();
    double                dt              = maxDt / 10; // subcycling with waveform iteration
    double                time            = 0;
    double                checkpoint_time = 0;

    while (interface.isCouplingOngoing()) {

      if (interface.requiresWritingCheckpoint())
        checkpoint_time = time;

      double value = -1.0;
      interface.readScalarData(dataID, vertexID, dt, value);
      BOOST_TEST(value == 0.1 * dataOneFunction(time)); // due to constant acceleration
      interface.advance(dt);
      time += dt;
      if (interface.requiresWritingCheckpoint())
        checkpoint_time = time;
    }
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // ParallelCoupling

#endif // PRECICE_NO_MPI
