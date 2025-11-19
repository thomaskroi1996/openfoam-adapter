#include "PressureTemporalDerivative.H"
#include "coupledPressureFvPatchField.H"

using namespace Foam;

preciceAdapter::FF::PressureTemporalDerivative::PressureTemporalDerivative(
    const Foam::fvMesh& mesh,
    const std::string namePTD)
:
mesh_(mesh),
p_(const_cast<volScalarField*>(&mesh.lookupObject<volScalarField>(namePTD))),
pOld_(const_cast<volScalarField*>(&mesh.lookupObject<volScalarField>(namePTD))),
dpdt_(const_cast<volScalarField*>(&mesh.lookupObject<volScalarField>(namePTD)))
{
    dataType_ = scalar;
}

std::size_t preciceAdapter::FF::PressureTemporalDerivative::write(double* buffer,
                                                                  bool meshConnectivity,
                                                                  const unsigned int dim)
{
    
    int bufferIndex = 0;
    
    if (this->locationType_ == LocationType::volumeCenters)
    {
        if (cellSetNames_.empty())
        {
            for (const auto& cell : dpdt_->internalField())
            {
                // when implementing more schemes, move rhs to function, so that its more readable
                buffer[bufferIndex++] = (p_->internalField()[cell] - pOld_->internalField()[cell]) / mesh_.time().deltaTValue();
                pOld_->internalField()[cell] = p_->internalField()[cell];
                // buffer[bufferIndex++] = cell;
            }
        }
        else
        {
            for (const auto& cellSetName : cellSetNames_)
            {
                cellSet overlapRegion(dpdt_->mesh(), cellSetName);
                const labelList& cells = overlapRegion.toc();
                
                for (const auto& currentCell : cells)
                {
                    buffer[bufferIndex++] = (p_->internalField()[currentCell] - pOld_->internalField()[currentCell]) / mesh_.time().deltaTValue();
                    pOld_->internalField()[currentCell] = p_->internalField()[currentCell];
                }
            }
        }
    }
    
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);
        
        forAll(dpdt_->boundaryFieldRef()[patchID], i)
        {
            buffer[bufferIndex++] = (p_->boundaryFieldRef()[patchID][i] - pOld_->boundaryFieldRef()[patchID][i]) / mesh_.time().deltaTValue();
            pOld_->boundaryFieldRef()[patchID][i] = p_->boundaryFieldRef()[patchID][i];
        }
    }

    return bufferIndex;
}

void preciceAdapter::FF::PressureTemporalDerivative::read(double* buffer, const unsigned int dim)
{
    int bufferIndex = 0;

    if (this->locationType_ == LocationType::volumeCenters)
    {
        if (cellSetNames_.empty())
        {
            for (auto& cell : dpdt_->ref())
            {
                cell = buffer[bufferIndex++];
            }
        }
        else
        {
            for (const auto& cellSetName : cellSetNames_)
            {
                cellSet overlapRegion(dpdt_->mesh(), cellSetName);
                const labelList& cells = overlapRegion.toc();

                for (const auto& currentCell : cells)
                {
                    dpdt_->ref()[currentCell] = buffer[bufferIndex++];
                }
            }
        }
    }

    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        scalarField* valuePatchPtr = &dpdt_->boundaryFieldRef()[patchID];
        if (isA<coupledPressureFvPatchField>(dpdt_->boundaryFieldRef()[patchID]))
        {
            valuePatchPtr = &refCast<coupledPressureFvPatchField>(
                                 dpdt_->boundaryFieldRef()[patchID])
                                 .refValue();
        }
        scalarField& valuePatch = *valuePatchPtr;

        forAll(dpdt_->boundaryFieldRef()[patchID], i)
        {
            valuePatch[i] =
                buffer[bufferIndex++];
        }
    }
}

bool preciceAdapter::FF::PressureTemporalDerivative::isLocationTypeSupported(const bool meshConnectivity) const
{
    if (meshConnectivity)
    {
        return (this->locationType_ == LocationType::faceCenters);
    }
    else
    {
        return (this->locationType_ == LocationType::faceCenters || this->locationType_ == LocationType::volumeCenters);
    }
}

std::string preciceAdapter::FF::PressureTemporalDerivative::getDataName() const
{
    return "PressureTemporalDerivative";
}