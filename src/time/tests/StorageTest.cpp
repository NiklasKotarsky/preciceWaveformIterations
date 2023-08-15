#include <Eigen/Core>
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "time/Storage.hpp"

using namespace precice;
using namespace precice::time;

BOOST_AUTO_TEST_SUITE(TimeTests)
BOOST_AUTO_TEST_SUITE(StorageTests)

// create storage and test for correct initial values.
BOOST_AUTO_TEST_CASE(testInitialize)
{
  PRECICE_TEST(1_rank);
  auto storage = Storage(0);
  int  nValues = 3;
  BOOST_TEST(storage.nTimes() == 0);
  storage.initialize(time::Sample{1, Eigen::VectorXd::Ones(nValues)});
  BOOST_TEST(storage.nDofs() == nValues);
  BOOST_TEST(storage.nTimes() == 2);
  for (int i = 0; i < nValues; i++) {
    BOOST_TEST(storage.getValuesAtOrAfter(0)(i) == 1);
    BOOST_TEST(storage.getValuesAtOrAfter(0.5)(i) == 1);
    BOOST_TEST(storage.getValuesAtOrAfter(1)(i) == 1);
  }
}

BOOST_AUTO_TEST_CASE(testInitializationVector)
{
  PRECICE_TEST(1_rank);

  const int       interpolationDegree = 0;
  const int       valuesSize          = 3;
  Eigen::VectorXd value(valuesSize);
  auto            storage = Storage(0);
  value << 0, 0, 0;
  storage.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(storage.nDofs() == valuesSize);

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(storage.sample(0.0)(i), 0.0));
    BOOST_TEST(testing::equals(storage.sample(1.0)(i), 0.0));
  }
}

// create storage and trim it.
BOOST_AUTO_TEST_CASE(testClear)
{
  PRECICE_TEST(1_rank);
  auto storage = Storage(0);
  int  nValues = 3;
  BOOST_TEST(storage.nTimes() == 0);
  storage.initialize(time::Sample{1, Eigen::VectorXd::Ones(nValues)});
  BOOST_TEST(storage.nDofs() == nValues);
  BOOST_TEST(storage.nTimes() == 2);
  BOOST_TEST(storage.maxStoredNormalizedDt() == 1.0);
  storage.trim();
  BOOST_TEST(storage.nDofs() == nValues);
  BOOST_TEST(storage.nTimes() == 1);
  BOOST_TEST(storage.maxStoredNormalizedDt() == 0.0);
}

// create storage, add some values and then move to next window.
BOOST_AUTO_TEST_CASE(testMove)
{
  PRECICE_TEST(1_rank);
  auto storage = Storage(0);
  int  nValues = 3;
  BOOST_TEST(storage.nTimes() == 0);
  storage.initialize(time::Sample{1, Eigen::VectorXd::Ones(nValues)});
  BOOST_TEST(storage.nDofs() == nValues);
  BOOST_TEST(storage.nTimes() == 2);
  BOOST_TEST(storage.maxStoredNormalizedDt() == 1.0);
  storage.trim();
  BOOST_TEST(storage.nTimes() == 1);
  storage.setSampleAtTime(0.5, time::Sample{1, Eigen::VectorXd::Ones(nValues)});
  BOOST_TEST(storage.nTimes() == 2);
  BOOST_TEST(storage.maxStoredNormalizedDt() == 0.5);
  storage.setSampleAtTime(1.0, time::Sample{1, Eigen::VectorXd::Zero(nValues)});
  BOOST_TEST(storage.nTimes() == 3);
  BOOST_TEST(storage.maxStoredNormalizedDt() == 1.0);
  for (int i = 0; i < nValues; i++) {
    BOOST_TEST(storage.getValuesAtOrAfter(0)(i) == 1);
    BOOST_TEST(storage.getValuesAtOrAfter(0.5)(i) == 1);
    BOOST_TEST(storage.getValuesAtOrAfter(1)(i) == 0);
  }
  storage.move();
  BOOST_TEST(storage.nDofs() == nValues);
  BOOST_TEST(storage.nTimes() == 2);
  BOOST_TEST(storage.maxStoredNormalizedDt() == 1.0);
  for (int i = 0; i < nValues; i++) {
    BOOST_TEST(storage.getValuesAtOrAfter(0)(i) == 0);
    BOOST_TEST(storage.getValuesAtOrAfter(1)(i) == 0);
  }
}

// get times and values
BOOST_AUTO_TEST_CASE(testGetTimesAndValues)
{
  PRECICE_TEST(1_rank);
  auto storage = Storage(0);
  int  nValues = 3;
  storage.initialize(time::Sample{1, Eigen::VectorXd::Ones(nValues)});
  storage.trim();
  storage.setSampleAtTime(0.5, time::Sample{1, Eigen::VectorXd::Ones(nValues)});
  storage.setSampleAtTime(1.0, time::Sample{1, Eigen::VectorXd::Zero(nValues)});
  auto times = storage.getTimes();
  BOOST_TEST(times[0] == 0.0);
  BOOST_TEST(times[1] == 0.5);
  BOOST_TEST(times[2] == 1.00);
  auto timesAndValues = storage.getTimesAndValues();
  BOOST_TEST(timesAndValues.first[0] == 0.0);
  BOOST_TEST(timesAndValues.first[1] == 0.5);
  BOOST_TEST(timesAndValues.first[2] == 1.00);
  for (int i = 0; i < nValues; i++) {
    BOOST_TEST(timesAndValues.second.col(0)(i) == 1);
    BOOST_TEST(timesAndValues.second.col(1)(i) == 1);
    BOOST_TEST(timesAndValues.second.col(2)(i) == 0);
  }
}

BOOST_AUTO_TEST_SUITE(InterpolationTests)

BOOST_AUTO_TEST_CASE(testInterpolateDataZerothDegree)
{
  PRECICE_TEST(1_rank);

  // Test zeroth degree interpolation
  const int       interpolationDegree = 0;
  const int       valuesSize          = 1;
  Eigen::VectorXd value(valuesSize);
  Storage         waveform(interpolationDegree);
  value(0) = 0.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value(0) = 1.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(waveform.nDofs() == valuesSize);
  BOOST_TEST(waveform.nTimes() == 2);

  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 1.0));

  value(0) = 2.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  waveform.move();
  BOOST_TEST(waveform.nTimes() == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  value(0) = 3.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 3.0));
}

BOOST_AUTO_TEST_CASE(testInterpolateDataFirstDegree)
{
  PRECICE_TEST(1_rank);

  // Test first degree interpolation
  const int       interpolationDegree = 1;
  const int       valuesSize          = 1;
  Eigen::VectorXd value(valuesSize);
  Storage         waveform(interpolationDegree);
  value(0) = 0.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value(0) = 1.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(waveform.nDofs() == valuesSize);
  BOOST_TEST(waveform.nTimes() == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 0.5));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 1.0));

  value(0) = 2.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  waveform.move();
  BOOST_TEST(waveform.nTimes() == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  value(0) = 3.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.5));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 3.0));
}

// Remove or modify this feature? Creating a second degree interpolant by using data from previous windows is difficult, because this would require several pieces of data during initialization. What would be useful: Generating a second degree interpolant from multiple samples in a single window (if available). This would go into the least-squares direction
BOOST_AUTO_TEST_CASE(testInterpolateDataSecondDegree)
{
  PRECICE_TEST(1_rank);

  // Test second degree interpolation, but there are not enough samples. Therefore, always only first degree.
  const int       interpolationDegree = 2;
  const int       valuesSize          = 1;
  Eigen::VectorXd value(valuesSize);
  Storage         waveform(interpolationDegree);
  value(0) = 0.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value(0) = 1.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(waveform.nDofs() == valuesSize);
  BOOST_TEST(waveform.nTimes() == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 0.5));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 1.0));

  waveform.trim();

  value(0) = 2.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  waveform.move();
  BOOST_TEST(waveform.nTimes() == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 2.0));

  value(0) = 8.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 5.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 8.0));

  waveform.trim();

  value(0) = 4.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 4.0));

  waveform.move();
  BOOST_TEST(waveform.nTimes() == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 4.0));

  value(0) = 8.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(testing::equals(waveform.sample(0.0)(0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(0.5)(0), 6.0));
  BOOST_TEST(testing::equals(waveform.sample(1.0)(0), 8.0));
}

BOOST_AUTO_TEST_CASE(testInterpolateDataFirstDegreeVector)
{
  PRECICE_TEST(1_rank);

  // Test first degree interpolation
  const int       interpolationDegree = 1;
  const int       valuesSize          = 3;
  Eigen::VectorXd value(valuesSize);
  Storage         waveform(interpolationDegree);
  value << 0, 0, 0;
  waveform.setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value << 1, 2, 3;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(waveform.nDofs() == valuesSize);
  BOOST_TEST(waveform.nTimes() == 2);

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), 0 * value[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i), 0.5 * value[i]));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i), value[i]));
  }

  waveform.trim();

  value << 2, 4, 2;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), 0 * value[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i), 0.5 * value[i]));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i), value[i]));
  }

  waveform.move();
  BOOST_TEST(waveform.nTimes() == 2);

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), value[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i), value[i]));
    BOOST_TEST(testing::equals(waveform.sample(1.0)(i), value[i]));
  }

  Eigen::VectorXd value0 = value;
  value << 1, 2, 3;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  for (int i = 0; i < valuesSize; i++) {
    BOOST_TEST(testing::equals(waveform.sample(0.0)(i), value0[i]));
    BOOST_TEST(testing::equals(waveform.sample(0.5)(i), 0.5 * value0[i] + 0.5 * value[i]));
  }
}

BOOST_AUTO_TEST_CASE(testPiecewiseInterpolateDataZerothDegree)
{
  PRECICE_TEST(1_rank);

  // Test zeroth degree interpolation
  const int       interpolationDegree = 0;
  const int       valuesSize          = 1;
  Eigen::VectorXd value(valuesSize);
  Storage         waveform(interpolationDegree);
  value(0) = 0.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value(0) = 0.5;
  waveform.setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});
  value(0) = 1.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(waveform.nDofs() == valuesSize);
  BOOST_TEST(waveform.nTimes() == 3);
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.5));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.5));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 1.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 1.0));

  waveform.trim();

  value(0) = 1.5;
  waveform.setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});

  value(0) = 2.0;
  waveform.setSampleAtTime(1.0 * time::Storage::WINDOW_END, time::Sample{1, value});
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 1.5));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 1.5));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.0));

  waveform.move();
  BOOST_TEST(waveform.nTimes() == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.0));

  value(0) = 3.0;
  waveform.setSampleAtTime(1.0 * time::Storage::WINDOW_END, time::Sample{1, value});
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 3.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 3.0));

  waveform.trim();

  value(0) = 1.5;
  waveform.setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});

  value(0) = 4.0;
  waveform.setSampleAtTime(1.0 * time::Storage::WINDOW_END, time::Sample{1, value});
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.0));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 1.5));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 1.5));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 4.0));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 4.0));
}

BOOST_AUTO_TEST_CASE(testPiecewiseInterpolateDataFirstDegree)
{
  PRECICE_TEST(1_rank);

  // Test zeroth degree interpolation
  const int       interpolationDegree = 1;
  const int       valuesSize          = 1;
  Eigen::VectorXd value(valuesSize);
  Storage         waveform(interpolationDegree);
  value(0) = 0.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value(0) = 0.5;
  waveform.setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});
  value(0) = 1.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(waveform.nDofs() == valuesSize);
  BOOST_TEST(waveform.nTimes() == 3);
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.25));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.50));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 0.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 1.00));

  value(0) = 1.5;
  waveform.setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});

  value(0) = 2.0;
  waveform.setSampleAtTime(1.0 * time::Storage::WINDOW_END, time::Sample{1, value});
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 0.75));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 1.50));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 1.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.00));

  waveform.move();
  BOOST_TEST(waveform.nTimes() == 2);

  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.00));

  value(0) = 3.0;
  waveform.setSampleAtTime(1.0 * time::Storage::WINDOW_END, time::Sample{1, value});
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 2.25));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 2.50));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 3.00));

  waveform.trim();

  value(0) = 1.5;
  waveform.setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});

  value(0) = 4.0;
  waveform.setSampleAtTime(1.0 * time::Storage::WINDOW_END, time::Sample{1, value});
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 2.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), 1.75));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 1.50));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 2.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 4.00));
}

BOOST_AUTO_TEST_CASE(testPiecewiseInterpolateDataSecondDegree)
{
  PRECICE_TEST(1_rank);

  // Test zeroth degree interpolation
  const int       interpolationDegree = 2;
  const int       valuesSize          = 1;
  Storage         waveform(interpolationDegree);
  Eigen::VectorXd value(valuesSize);
  value(0) = 0.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_START, time::Sample{1, value});
  value(0) = 0.0;
  waveform.setSampleAtTime(0.5 * time::Storage::WINDOW_END, time::Sample{1, value});
  value(0) = 2.0;
  waveform.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{1, value});

  BOOST_TEST(waveform.nDofs() == valuesSize);
  BOOST_TEST(waveform.nTimes() == 3);
  BOOST_TEST(testing::equals(waveform.sample(0.00)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.25)(0), -0.25));
  BOOST_TEST(testing::equals(waveform.sample(0.50)(0), 0.00));
  BOOST_TEST(testing::equals(waveform.sample(0.75)(0), 0.75));
  BOOST_TEST(testing::equals(waveform.sample(1.00)(0), 2.00));
}

BOOST_AUTO_TEST_CASE(testPiecewiseInterpolateDataThirdDegree)
{
  PRECICE_TEST(1_rank);

  // Test zeroth degree interpolation
  const int interpolationDegree = 3;
  const int valuesSize          = 1;
  Storage   waveform(interpolationDegree);

  // linearly increasing values
  Eigen::VectorXd value(valuesSize);
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    value(0) = t;
    waveform.setSampleAtTime(t * time::Storage::WINDOW_END, time::Sample{1, value});
  }

  BOOST_TEST(waveform.nDofs() == valuesSize);
  BOOST_TEST(waveform.nTimes() == 5);

  for (double t : std::vector<double>{0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t));
  }

  waveform.trim();

  // quadratically increasing values
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    value(0) = t * t;
    waveform.setSampleAtTime(t * time::Storage::WINDOW_END, time::Sample{1, value});
  }

  // interpolates given values
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t));
  }

  // introduces no approximation error w.r.t function
  for (double t : std::vector<double>{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t));
  }

  waveform.trim();

  // cubically increasing values
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    value(0) = t * t * t;
    waveform.setSampleAtTime(t * time::Storage::WINDOW_END, time::Sample{1, value});
  }

  // interpolates given values
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t));
  }

  // introduces no approximation error w.r.t function
  for (double t : std::vector<double>{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t));
  }

  waveform.trim();

  // cubically increasing values, but with non-uniform spacing
  for (double t : std::vector<double>{0, 0.01, 0.1, 0.2, 1}) {
    value(0) = t * t * t;
    waveform.setSampleAtTime(t * time::Storage::WINDOW_END, time::Sample{1, value});
  }

  // interpolates given values
  for (double t : std::vector<double>{0, 0.01, 0.1, 0.2, 1}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t));
  }

  // introduces no approximation error w.r.t function
  for (double t : std::vector<double>{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t));
  }

  waveform.trim();

  // quadratically increasing values, but with non-uniform spacing
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    value(0) = t * t * t * t;
    waveform.setSampleAtTime(t * time::Storage::WINDOW_END, time::Sample{1, value});
  }

  // interpolates given values
  for (double t : std::vector<double>{0, 0.25, 0.5, 0.75, 1}) {
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t * t));
  }

  // introduces approximation error w.r.t function
  for (double t : std::vector<double>{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
    double tol = 0.015625; // error < h**3 = 0.015625
    BOOST_TEST(testing::equals(waveform.sample(t)(0), t * t * t * t, tol));
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
