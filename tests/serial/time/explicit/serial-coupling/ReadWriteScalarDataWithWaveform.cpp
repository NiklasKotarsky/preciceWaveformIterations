#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Time)
BOOST_AUTO_TEST_SUITE(Explicit)
BOOST_AUTO_TEST_SUITE(SerialCoupling)

/**
 * @brief Test to run a simple coupling with subcycling.
 *
 * Ensures that each time step provides its own data, but preCICE only exchanges data at the end of the window.
 */
BOOST_AUTO_TEST_CASE(ReadWriteScalarDataWithWaveformSecond)
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
  VertexID vertexID = precice.setMeshVertex(meshID, Eigen::Vector3d(0.0, 0.0, 0.0).data());

  int    nSubsteps  = 4; // perform subcycling on solvers. 4 steps happen in each window.
  int    nWindows   = 5; // perform 5 windows.
  int    timestep   = 0;
  int    timewindow = 0;
  double time       = 0;

  if (precice.requiresInitialData()) {
    writeData = writeFunction(time);
    precice.writeScalarData(writeDataID, vertexID, writeData);
  }

  double maxDt = precice.initialize();
  BOOST_TEST(maxDt == 2.0); // use window size != 1.0 to be able to detect more possible bugs
  double windowDt  = maxDt;
  double dt        = windowDt / (nSubsteps - 0.5); // Solver always tries to do a timestep of fixed size.
  double currentDt = dt > maxDt ? maxDt : dt;      // determine actual timestep length; must fit into remaining time in window
  double timeCheckpoint;

  while (precice.isCouplingOngoing()) {

    precice.readScalarData(readDataID, vertexID, currentDt, readData);

    if (context.isNamed("SolverOne")) { // in the first iteration of each window, we only have one sample of data. Therefore constant interpolation
      BOOST_TEST(readData == readFunction(timeCheckpoint));
    } else { // in the following iterations we have two samples of data. Therefore linear interpolation
      BOOST_TEST(readData == readFunction(time + currentDt));
    }

    precice.readScalarData(readDataID, vertexID, currentDt / 2, readData);

    if (context.isNamed("SolverOne")) { // in the first iteration of each window, we only have one sample of data. Therefore constant interpolation
      BOOST_TEST(readData == readFunction(timeCheckpoint));
    } else { // in the following iterations we have two samples of data. Therefore linear interpolation
      BOOST_TEST(readData == readFunction(time + currentDt / 2));
    }

    // solve usually goes here. Dummy solve: Just sampling the writeFunction.
    time += currentDt;
    // Need to keep track of the timewindows to check correct order on the first participant!
    timestep++;
    if (timestep % nSubsteps == 0) {
      timeCheckpoint += windowDt;
    }

    writeData = writeFunction(time);
    precice.writeScalarData(writeDataID, vertexID, writeData);
    maxDt = precice.advance(currentDt);

    currentDt = dt > maxDt ? maxDt : dt;
  }

  BOOST_TEST(timestep == nWindows * nSubsteps);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Time
BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // SerialCoupling

#endif // PRECICE_NO_MPI
