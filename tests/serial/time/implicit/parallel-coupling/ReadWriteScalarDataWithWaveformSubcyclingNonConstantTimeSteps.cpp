#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(ParallelCoupling)

/**
 * @brief Test to run a simple coupling with first order waveform subcycling.
 *
 * Provides a dt argument to the read function. A first order waveform is used.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSubcyclingNonConstantTimeSteps)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  SolverInterface precice(context.name, context.config(), 0, 1);

  MeshID meshID;
  DataID writeDataID;
  DataID readDataID;

  typedef double (*DataFunction)(double);

  DataFunction dataOneFunction = [](double t) -> double {
    return (double) (2 + t);
  };
  DataFunction dataTwoFunction = [](double t) -> double {
    return (double) (10 + t);
  };
  DataFunction writeFunction;
  DataFunction readFunction;

  if (context.isNamed("SolverOne")) {
    meshID        = precice.getMeshID("MeshOne");
    writeDataID   = precice.getDataID("DataOne", meshID);
    writeFunction = dataOneFunction;
    readDataID    = precice.getDataID("DataTwo", meshID);
    readFunction  = dataTwoFunction;

  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    meshID        = precice.getMeshID("MeshTwo");
    writeDataID   = precice.getDataID("DataTwo", meshID);
    writeFunction = dataTwoFunction;
    readDataID    = precice.getDataID("DataOne", meshID);
    readFunction  = dataOneFunction;
  }

  double   writeData, readData;
  VertexID vertexID;

  vertexID = precice.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());
  int    timestepCheckpoint;
  double time = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeScalarData(writeDataID, vertexID, writeData);
  }

  double              maxDt = precice.initialize();
  std::vector<double> dt_v; // the different dt used
  if (context.isNamed("SolverOne")) {
    dt_v = {maxDt / 5, 2 * maxDt / 5, 2 * maxDt / 5, 0}; // initialize the different timesteps and one extra timestep
  } else {
    dt_v = {maxDt / 6, 2 * maxDt / 6, maxDt / 6, 2 * maxDt / 6, 0}; // initialize the different timesteps ratios and one extra timestep
  }

  double currentDt = dt_v[0];
  double timeCheckpoint;
  int    iterations;
  int    nbr_timestep = 0;

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      nbr_timestep   = 0;
      currentDt      = dt_v[0];
      timeCheckpoint = time;
      iterations     = 0;
    }
    precice.readScalarData(readDataID, vertexID, currentDt, readData);

    if (iterations == 0) { // in the first iteration of each window, we only have one sample of data. Therefore constant interpolation
      BOOST_TEST(readData == readFunction(timeCheckpoint));
    } else { // in the following iterations we have two samples of data. Therefore linear interpolation
      BOOST_TEST(readData == readFunction(time + currentDt));
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    writeData = writeFunction(time);
    precice.writeScalarData(writeDataID, vertexID, writeData);
    maxDt = precice.advance(currentDt);
    nbr_timestep += 1;
    currentDt = dt_v[nbr_timestep];
    currentDt = currentDt > maxDt ? maxDt : currentDt;

    if (precice.requiresReadingCheckpoint()) {

      nbr_timestep = 0;
      currentDt    = dt_v[0];
      time         = timeCheckpoint;
      iterations++;
    }
  }

  precice.finalize();
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // ParallelCoupling

#endif // PRECICE_NO_MPI
