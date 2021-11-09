#pragma once

#include <map>
#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace time {

class TimeConfiguration : public xml::XMLTag::Listener {
public:
  TimeConfiguration();

  /// Callback method required when using xml::XMLTag.
  virtual void xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag);

  /// Callback method required when using xml::XMLTag.
  virtual void xmlEndTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag);

private:
  logging::Logger _log{"time::TimeConfiguration"};
};
} // namespace time
} // namespace precice
