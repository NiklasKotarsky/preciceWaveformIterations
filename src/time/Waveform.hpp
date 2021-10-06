#pragma once

#include <Eigen/Core>
#include "logging/Logger.hpp"

namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
struct WaveformFixture;
} // namespace testing

namespace time {

class Waveform {
  friend struct testing::WaveformFixture; // Make the fixture friend of this class
public:
  /// To be used, when the interpolation order is not defined for this Waveform.
  static const int UNDEFINED_INTERPOLATION_ORDER;

  /**
   * @brief Waveform object which stores data of current and past time windows for performing extrapolation.
   * @param initializedNumberOfData defines how many pieces of data one sample in time consists of
   * @param extrapolatioOrder defines the maximum extrapolation order supported by this Waveform and reserves storage correspondingly
   * @param interpolationOrder defines the maximum interpolation order supported by this Waveform and reserves storage correspondingly
   */
  Waveform(const int initializedNumberOfData,
           const int extrapolationOrder,
           const int interpolationOrder);

  /**
   * @brief resizes _timeWindows to store more data. Used for already created waveforms.
   * @param numberOfData defines how many pieces of data one sample in time consists of
   */
  void resizeData(int numberOfData);

  /**
   * @brief Updates entry in _timeWindows corresponding to this window with given data
   * @param data new sample for this time window
   */
  void store(const Eigen::VectorXd &data);

  /**
   * @brief Updates entry in _timeWindows corresponding to a given column ID with given data
   * @param data new sample for this time window
   * @param columnID ID of column to be updated
   */
  void storeAt(const Eigen::VectorXd data, int columnID);

  /**
   * @brief Called, when moving to the next time window. All entries in _timeWindows are shifted. The new entry is initialized as the value from the last window (= constant extrapolation)
   */
  void moveToNextWindow();

  /**
   * @brief sample Waveform. Uses interpolation with Waveform's interpolation order, if necessary
   * @param normalizedDt time where the sampling inside the window happens. 0 refers to the beginning of the window and 1 to the end.
   */
  Eigen::VectorXd sample(const double normalizedDt);

  /**
   * @brief getter for Eigen::MatrixXd containing data of current and past time windows. Each column represents a sample in time, with col(0)
   * being the current time window.
   */
  const Eigen::MatrixXd &lastTimeWindows();

  /// @todo try to make this private!
  /**
   * @brief returns number of data per sample in time stored by this waveform
   */
  int numberOfData(); // @todo bad naming, consider renaming. See https://github.com/precice/precice/pull/1094#pullrequestreview-771715472

  /// @todo try to make this private!
  /**
   * @brief returns number of samples in time stored by this waveform
   */
  int numberOfSamples();

private:
  /// Data values of time windows.
  Eigen::MatrixXd _timeWindows;

  /// extrapolation order for this waveform
  const int _extrapolationOrder;

  /// interpolation order for this waveform
  const int _interpolationOrder;

  /// number of valid samples in _timeWindows
  int _numberOfValidSamples;

  /**
   * @brief returns number of valid samples in time stored by this waveform
   */
  int numberOfValidSamples();

  mutable logging::Logger _log{"time::Waveform"};

  /**
   * @brief Extrapolates data _timeWindows using an extrapolation scheme of given order. 
   * 
   * If the order condition cannot be satisfied, since there are not enough samples available, the order is automatically reduced.
   * If order two is required, but only two samples are available, the extrapolation order is automatically reduced to one.
   */
  Eigen::VectorXd extrapolateData();

  /**
   * @brief Interpolates data inside current time window using an interpolation scheme of the order of this Waveform.
   *
   * @param normalizedDt time where the sampling inside the window happens. 0 refers to the beginning of the window and 1 to the end.
   */
  Eigen::VectorXd interpolateData(const double normalizedDt);
};

} // namespace time
} // namespace precice
