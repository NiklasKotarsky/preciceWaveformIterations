#ifndef PRECICE_NO_MPI

#include <precice/SolverInterface.hpp>
#include <vector>
#include "testing/Testing.hpp"

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
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSubcyclingVaryingTimeGrids)
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

  double time_window_size = precice.initialize();
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
    precice.readScalarData(readDataID, vertexID, dt, readData);

    if (iterations == 0) { // in the first iteration of each window, we only have one sample of data. Therefore constant interpolation
      BOOST_TEST(readData == readFunction(timeCheckpoint));
    } else { // in the following iterations we have two samples of data. Therefore linear interpolation
      BOOST_TEST(readData == readFunction(time + dt));
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += dt;
    writeData = writeFunction(time);
    precice.writeScalarData(writeDataID, vertexID, writeData);
    maxDt = precice.advance(dt);
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
