#include "ReadDataContext.hpp"

#include "time/Waveform.hpp"

namespace precice::impl {

logging::Logger ReadDataContext::_log{"impl::ReadDataContext"};

ReadDataContext::ReadDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    int           interpolationOrder)
    : DataContext(data, mesh)
{
  _waveform = std::make_shared<time::Waveform>(interpolationOrder);
}

void ReadDataContext::appendMappingConfiguration(const MappingContext &mappingContext, const MeshContext &meshContext)
{
  PRECICE_ASSERT(!hasReadMapping(), "The read data context must be unique. Otherwise we would have an ambiguous read data operation on the user side.")
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData data = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(data != _providedData, "Data the read mapping is mapping from needs to be different from _providedData");
  appendMapping(mappingContext, data, _providedData);
  PRECICE_ASSERT(hasReadMapping());
}

int ReadDataContext::getInterpolationOrder() const
{
  return _waveform->getInterpolationOrder();
}

void ReadDataContext::storeDataInWaveform(double relativeDt)
{
  _waveform->store(_providedData->values(), relativeDt); // store mapped or received _providedData in the _waveform
}

DataID ReadDataContext::getProvidedDataID()
{
  PRECICE_ASSERT(_providedData);
  return _providedData->getID();
}

std::vector<DataID> ReadDataContext::getFromDataIDs()
{
  PRECICE_ASSERT(_mappingContexts.size() > 0, "No mapping context. Use getProvidedDataID.");
  std::vector<DataID> dataIds;
  for (unsigned int i = 0; i < _mappingContexts.size(); ++i) {
    const DataID fromDataID = getFromDataID(i);
    dataIds.emplace_back(fromDataID);
  }
  return dataIds;
}

void ReadDataContext::mapFromData(DataID id)
{
  // @todo needs refactoring. See mapData()
  // Execute the mapping
  for (unsigned int i = 0; i < _mappingContexts.size(); ++i) {
    const DataID fromDataID = getFromDataID(i);
    if (fromDataID == id) {
      const DataID toDataID = getToDataID(i);
      // Reset the toData before executing the mapping
      _toData[i]->toZero();
      _mappingContexts[i].mapping->map(fromDataID, toDataID);
    }
  }
}

Eigen::VectorXd ReadDataContext::sampleWaveformAt(double normalizedDt)
{
  return _waveform->sample(normalizedDt);
}

void ReadDataContext::initializeWaveform()
{
  PRECICE_ASSERT(not hasWriteMapping(), "Write mapping does not need waveforms.");
  _waveform->initialize(_providedData->values());
}

void ReadDataContext::moveToNextWindow()
{
  _waveform->moveToNextWindow();
}

} // namespace precice::impl
