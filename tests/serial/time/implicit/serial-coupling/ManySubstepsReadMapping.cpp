#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief Test to run a simple coupling with third order waveform subcycling.
 */
BOOST_AUTO_TEST_CASE(ManySubstepsReadMapping)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  Participant precice(context.name, context.config(), 0, 1);

  std::string meshName;
  std::string writeDataName;
  std::string readDataName;

  typedef double (*DataFunction)(double);

  // use third degree function to test third order waveforms
  DataFunction dataOneFunction = [](double t) -> double {
    return (double) (2 + t * t * t);
  };
  DataFunction dataTwoFunction = [](double t) -> double {
    return (double) (10 + t * t * t);
  };
  DataFunction writeFunction;
  DataFunction readFunction;

  if (context.isNamed("SolverOne")) {
    meshName      = "MeshOne";
    writeDataName = "DataOne";
    writeFunction = dataOneFunction;
    readDataName  = "DataTwo";
    readFunction  = dataTwoFunction;
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshName      = "MeshTwo";
    writeDataName = "DataTwo";
    writeFunction = dataTwoFunction;
    readDataName  = "DataOne";
    readFunction  = dataOneFunction;
  }

  const int dim       = 100;
  int       nSubsteps = 500; // perform subcycling on solvers. 4 steps happen in each window.

  Eigen::Matrix<double, 2 * dim, 1> v0 = Eigen::VectorXd::LinSpaced(Eigen::Sequential, 2 * dim, 0.0, 2.0);
  VertexID                          vertexID[dim];
  precice.setMeshVertices(meshName, v0, vertexID);

  double writeData[dim] = {1.0};
  double readData[dim]  = {0.0};
  int    nWindows       = 1; // perform 5 windows.
  int    timestep       = 0;
  double time           = 0;

  if (precice.requiresInitialData()) {
    precice.writeData(meshName, writeDataName, vertexID, writeData);
  }

  precice.initialize();
  double maxDt              = precice.getMaxTimeStepSize();
  double windowDt           = maxDt;
  int    timestepCheckpoint = 0;
  double dt                 = windowDt / nSubsteps; // Timestep length desired by solver. E.g. 4 steps  with size 1/4
  double currentDt          = dt;                   // Timestep length used by solver
  double timeCheckpoint     = 0.0;
  int    iterations         = 0;

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      timeCheckpoint     = time;
      timestepCheckpoint = timestep;
      iterations         = 0;
    }

    precice.readData(meshName, readDataName, vertexID, currentDt, readData);

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    timestep++;
    precice.writeData(meshName, writeDataName, vertexID, writeData);
    precice.advance(currentDt);
    double maxDt = precice.getMaxTimeStepSize();
    if (precice.requiresReadingCheckpoint()) {
      time     = timeCheckpoint;
      timestep = timestepCheckpoint;
      iterations++;
    }
    currentDt = dt > maxDt ? maxDt : dt;
  }

  precice.finalize();
  BOOST_TEST(timestep == nWindows * nSubsteps);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI
