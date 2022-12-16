#include "MultiCouplingScheme.hpp"
#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
#include <ostream>
#include <type_traits>
#include <utility>
#include "acceleration/Acceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"

namespace precice::cplscheme {

MultiCouplingScheme::MultiCouplingScheme(
    double                             maxTime,
    int                                maxTimeWindows,
    double                             timeWindowSize,
    int                                validDigits,
    const std::string &                localParticipant,
    std::map<std::string, m2n::PtrM2N> m2ns,
    constants::TimesteppingMethod      dtMethod,
    const std::string &                controller,
    int                                maxIterations,
    int                                extrapolationOrder)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, localParticipant, maxIterations, Implicit, dtMethod, extrapolationOrder),
      _m2ns(std::move(m2ns)), _controller(controller), _isController(controller == localParticipant)
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  // Controller participant never does the first step, because it is never the first participant
  setDoesFirstStep(!_isController);
  PRECICE_DEBUG("MultiCoupling scheme is created for {}.", localParticipant);
}

void MultiCouplingScheme::determineInitialDataExchange()
{
  for (auto &sendExchange : _sendDataVector) {
    determineInitialSend(sendExchange.second);
  }
  for (auto &receiveExchange : _receiveDataVector) {
    determineInitialReceive(receiveExchange.second);
  }
}

std::vector<std::string> MultiCouplingScheme::getCouplingPartners() const
{
  std::vector<std::string> partnerNames;

  for (const auto &m2nPair : _m2ns) {
    partnerNames.emplace_back(m2nPair.first);
  }

  return partnerNames;
}

void MultiCouplingScheme::clearTimeStepSendStorage()
{
  for (auto &sendExchange : _sendDataVector) {
    for (const DataMap::value_type &pair : sendExchange.second) {
      pair.second->clearTimeStepsStorage();
    }
  }
}

void MultiCouplingScheme::exchangeInitialData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");

  if (_isController) {
    if (receivesInitializedData()) {
      for (auto &receiveExchange : _receiveDataVector) {
        for (const DataMap::value_type &pair : receiveExchange.second) {
          pair.second->clearTimeStepsStorage();
        }
        receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
      }
      checkDataHasBeenReceived();
    } else {
      for (auto &receiveExchange : _receiveDataVector) {
        initializeZeroReceiveData(receiveExchange.second);
      }
    }
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.first], sendExchange.second);
      }
    }
  } else {
    if (sendsInitializedData()) {
      for (auto &sendExchange : _sendDataVector) {
        sendData(_m2ns[sendExchange.first], sendExchange.second);
      }
    }
    if (receivesInitializedData()) {
      for (auto &receiveExchange : _receiveDataVector) {
        for (const DataMap::value_type &pair : receiveExchange.second) {
          pair.second->clearTimeStepsStorage();
        }
        receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
      }
      checkDataHasBeenReceived();
    } else {
      for (auto &receiveExchange : _receiveDataVector) {
        initializeZeroReceiveData(receiveExchange.second);
      }
    }
  }
  PRECICE_DEBUG("Initial data is exchanged in MultiCouplingScheme");
}

void MultiCouplingScheme::moveToNextWindow()
{
  for (auto &receiveExchange : _receiveDataVector) {
    for (auto &receiveData : receiveExchange.second) {
      receiveData.second->moveToNextWindow();
    }
  }
}

void MultiCouplingScheme::storeTimeStepSendData(double relativeDt)
{
  PRECICE_ASSERT(relativeDt > time::Storage::WINDOW_START);
  PRECICE_ASSERT(relativeDt <= time::Storage::WINDOW_END);
  for (auto &sendExchange : _sendDataVector) {
    for (auto &sendData : sendExchange.second) {
      auto values = sendData.second->values();
      sendData.second->storeValuesAtTime(relativeDt, values);
    }
  }
}

void MultiCouplingScheme::storeTimeStepReceiveData(double relativeDt)
{
  PRECICE_ASSERT(math::greaterEquals(relativeDt, time::Storage::WINDOW_START), relativeDt);
  PRECICE_ASSERT(math::greaterEquals(time::Storage::WINDOW_END, relativeDt), relativeDt);
  for (auto &receiveExchange : _receiveDataVector) {
    for (auto &receiveData : receiveExchange.second) {
      bool mustOverride = true;
      receiveData.second->storeValuesAtTime(relativeDt, receiveData.second->values(), mustOverride);
    }
  }
}

void MultiCouplingScheme::retreiveTimeStepReceiveData(double relativeDt)
{
  // @todo breaks, if different receiveData live on different time-meshes. This is a realistic use-case for multi coupling! Should use a different signature here to individually retreiveTimeStepData from receiveData. Would also be helpful for mapping.
  PRECICE_ASSERT(math::greaterEquals(relativeDt, time::Storage::WINDOW_START), relativeDt);
  PRECICE_ASSERT(math::greaterEquals(time::Storage::WINDOW_END, relativeDt), relativeDt);
  for (auto &receiveExchange : _receiveDataVector) {
    for (auto &receiveData : receiveExchange.second) {
      receiveData.second->values() = receiveData.second->getValuesAtTime(relativeDt);
      //retreiveTimeStepForData(relativeDt, receiveData.second->getDataID());  // @todo: Causes problems. Why? Because DataID is not unique for MultiCouplingScheme's allData.
    }
  }
}

std::vector<double> MultiCouplingScheme::getReceiveTimes()
{
  //@todo stub implementation. Should walk over all receive data, get times and ensure that all times vectors actually hold the same times (since otherwise we would have to get times individually per data)
  //@todo subcycling is not supported for MultiCouplingScheme, because this needs a complicated interplay of picking the right data in time and mapping this data. This is hard to realize with the current implementation.
  auto times = std::vector<double>({time::Storage::WINDOW_END});
  return times;
}

typedef std::map<int, PtrCouplingData> DataMap;

const DataMap MultiCouplingScheme::getAllData()
{
  DataMap allData;
  // @todo use C++17 std::map::merge
  for (auto &sendData : _sendDataVector) {
    allData.insert(sendData.second.begin(), sendData.second.end());
  }
  for (auto &receiveData : _receiveDataVector) {
    allData.insert(receiveData.second.begin(), receiveData.second.end());
  }
  return allData;
}

void MultiCouplingScheme::performReceiveOfFirstAdvance()
{
  return; // no action needed.
}

void MultiCouplingScheme::exchangeFirstData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  // @todo implement MultiCouplingScheme for explicit coupling

  PRECICE_DEBUG("Computed full length of iteration");

  if (_isController) {
    for (auto &receiveExchange : _receiveDataVector) {
      for (const DataMap::value_type &pair : receiveExchange.second) {
        pair.second->clearTimeStepsStorage();
      }
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();

  } else {
    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
  }
}

void MultiCouplingScheme::exchangeSecondData()
{
  PRECICE_ASSERT(isImplicitCouplingScheme(), "MultiCouplingScheme is always Implicit.");
  // @todo implement MultiCouplingScheme for explicit coupling

  if (_isController) {
    doImplicitStep();
    for (const auto &m2nPair : _m2ns) {
      sendConvergence(m2nPair.second);
    }

    for (auto &sendExchange : _sendDataVector) {
      sendData(_m2ns[sendExchange.first], sendExchange.second);
    }
  } else {
    receiveConvergence(_m2ns[_controller]);

    for (auto &receiveExchange : _receiveDataVector) {
      for (const DataMap::value_type &pair : receiveExchange.second) {
        pair.second->clearTimeStepsStorage();
      }
      receiveData(_m2ns[receiveExchange.first], receiveExchange.second);
    }
    checkDataHasBeenReceived();
  }
  if (hasConverged()) {
    for (const DataMap::value_type &data : getAllData()) {
      data.second->moveTimeStepsStorage();
    }
  }
}

void MultiCouplingScheme::addDataToSend(
    const mesh::PtrData &data,
    mesh::PtrMesh        mesh,
    bool                 initialize,
    const std::string &  to)
{
  int id = data->getID();
  PRECICE_DEBUG("Configuring send data to {}", to);
  PtrCouplingData ptrCplData(new CouplingData(data, std::move(mesh), initialize, getExtrapolationOrder()));
  _sendDataVector[to].emplace(id, ptrCplData);
}

void MultiCouplingScheme::addDataToReceive(
    const mesh::PtrData &data,
    mesh::PtrMesh        mesh,
    bool                 initialize,
    const std::string &  from)
{
  int id = data->getID();
  PRECICE_DEBUG("Configuring receive data from {}", from);
  PtrCouplingData ptrCplData(new CouplingData(data, std::move(mesh), initialize, getExtrapolationOrder()));
  _receiveDataVector[from].emplace(id, ptrCplData);
}

} // namespace precice::cplscheme
