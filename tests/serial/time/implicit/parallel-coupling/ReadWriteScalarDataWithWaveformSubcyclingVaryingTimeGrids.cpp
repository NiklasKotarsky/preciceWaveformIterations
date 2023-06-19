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
 * @brief Test to run a simple coupling with time grids that change with the waveform iteration and subcycling.
 *
 * Provides a dt argument to the read function. A first order waveform is used.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSubcyclingVaryingTimeGrids)
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

  double v0[] = {0, 0, 0};
  vertexID    = precice.setMeshVertex(meshName, v0);
  int    timestepCheckpoint;
  double time = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
  }

  precice.initialize();
  double time_window_size = precice.getMaxTimeStepSize();
  double maxDt            = time_window_size;
  double dt               = maxDt;
  double timeCheckpoint;
  int    iterations;
  int    nbr_timestep = 0;

  while (precice.isCouplingOngoing()) {
    if (precice.requiresWritingCheckpoint()) {
      dt             = maxDt;
      timeCheckpoint = time;
      iterations     = 0;
    }
    precice.readData(meshName, readDataName, {&vertexID, 1}, dt, {&readData, 1});

    if (iterations == 0) { // in the first iteration of each window, we only have one sample of data. Therefore constant interpolation
      BOOST_TEST(readData == readFunction(timeCheckpoint));
    } else { // in the following iterations we have two samples of data. Therefore linear interpolation
      BOOST_TEST(readData == readFunction(time + dt));
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += dt;
    writeData = writeFunction(time);
    precice.writeData(meshName, writeDataName, {&vertexID, 1}, {&writeData, 1});
    precice.advance(dt);
    maxDt = precice.getMaxTimeStepSize();
    dt    = dt > maxDt ? maxDt : dt;

    if (precice.requiresReadingCheckpoint()) {
      nbr_timestep = 0;
      dt *= 0.5;
      time = timeCheckpoint;
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
