#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_SUITE(ParallelCoupling)

/**
 * @brief Test to run a simple coupling with a non uniform time grid and waveform subcycling.
 *
 * Provides a dt argument to the read function. A first order waveform is used.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSubcyclingNonConstantTimeSteps)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));

  Participant precice(context.name, context.config(), 0, 1);

  typedef double (*DataFunction)(double);

  DataFunction dataOneFunction = [](double t) -> double {
    return (double) (2 + t);
  };
  DataFunction dataTwoFunction = [](double t) -> double {
    return (double) (10 + t);
  };
  DataFunction writeFunction;
  DataFunction readFunction;

  std::string meshName, writeDataName, readDataName;
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

  double   writeData, readData;
  VertexID vertexID;
  double   v0[] = {0, 0, 0};
  vertexID      = precice.setMeshVertex(meshName, v0);
  int    timestepCheckpoint;
  double time = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
  }

  precice.initialize();
  double              maxDt = precice.getMaxTimeStepSize();
  std::vector<double> dt_v; // the different dt used
  if (context.isNamed("SolverOne")) {
    dt_v = {maxDt / 5, 2 * maxDt / 5, 2 * maxDt / 5, 0}; // initialize the different timesteps and one extra timestep
  } else {
    dt_v = {maxDt / 6, 2 * maxDt / 6, maxDt / 6, 2 * maxDt / 6, 0}; // initialize the different timesteps and one extra timestep
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
    precice.readData(meshName, readDataName, {&vertexID, 1}, currentDt, {&readData, 1});

    if (iterations == 0) { // in the first iteration of each window, we only have one sample of data. Therefore constant interpolation
      BOOST_TEST(readData == readFunction(timeCheckpoint));
    } else { // in the following iterations we have two samples of data. Therefore linear interpolation
      BOOST_TEST(readData == readFunction(time + currentDt));
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
    precice.advance(currentDt);
    maxDt = precice.getMaxTimeStepSize();
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
