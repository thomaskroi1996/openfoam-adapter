#include "PressureTemporalDerivative.H"
#include "coupledPressureFvPatchField.H"

using namespace Foam;

preciceAdapter::FF::PressureTemporalDerivative::PressureTemporalDerivative(
    const Foam::fvMesh& mesh,
    const std::string namePTD)
:
mesh_(mesh),
p_(const_cast<volScalarField*>(&mesh.lookupObject<volScalarField>(namePTD))),
firstStep_(true)
{
    dataType_ = scalar;
}

// this function gets called in Interface.C line 574-ish
// double* buffer points to a memory location that has a size of: dim_ * numDataLocations_ (vector) or numDataLocations_ (scalar)
// we can just create a new buffer of the same size for pOld_
// cleaner would be to somehow use the openFoam API for using writeable objects for pOld_
// because pOld_->internalField()[cell] = p_->internalField()[cell] doesnt work when we instantiate with lookupObject or lookupObjectRef
std::size_t preciceAdapter::FF::PressureTemporalDerivative::write(double* buffer,
                                                                  bool meshConnectivity,
                                                                  const unsigned int dim)
{

    std::cout << "We are now in PressureTemporalDerivative::write()" << std::endl;
    int bufferIndex = 0;


    // can't be private member because we dont have access to buffer
    // so only create if firstStep_ = true
    if (firstStep_)
    {
        std::cout << "firstStep_: " << firstStep_ << std::endl;
        size_t nCells = p_->internalField().size();
        std::cout << "nCells: " << nCells << std::endl;
        pOld_.resize(80200, 0.0); // JUST FOR TESTING OBVIOUSLY, it works that way! nCells = 80000
        firstStep_ = false;
    }
    
    if (this->locationType_ == LocationType::volumeCenters)
    {
        if (cellSetNames_.empty())
        {
            for (const auto& cell : p_->internalField())
            {
                // when implementing more schemes, move rhs to function, so that its more readable
                // std::cout << mesh_.time().deltaTValue() << std::endl; deltaTValue is 5e-5
                // if ((p_->internalField()[cell] - pOld_[bufferIndex]) > 1e10)
                // {
                //     std::cout << "p_: " << p_->internalField()[cell] << std::endl;
                //     std::cout << "pOld_: " << pOld_[bufferIndex] << std::endl;
                //     std::cout << "Subtraction part of derivative: " << (p_->internalField()[cell] - pOld_[bufferIndex]) << std::endl;
                // }
                buffer[bufferIndex] = (cell - pOld_[bufferIndex]) / mesh_.time().deltaTValue(); // 
                pOld_[bufferIndex] = cell;
                bufferIndex++;
            }
        }
        else
        {
            for (const auto& cellSetName : cellSetNames_)
            {
                cellSet overlapRegion(p_->mesh(), cellSetName);
                const labelList& cells = overlapRegion.toc();
                
                for (const auto& currentCell : cells)
                {
                    buffer[bufferIndex] = (p_->internalField()[currentCell] - pOld_[bufferIndex]) / mesh_.time().deltaTValue();
                    pOld_[bufferIndex] = p_->internalField()[currentCell];
                    bufferIndex++;                
                }
            }
        }
    }
    
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);
        
        forAll(p_->boundaryFieldRef()[patchID], i)
        {
            buffer[bufferIndex] = (p_->boundaryFieldRef()[patchID][i] - pOld_[bufferIndex]) / mesh_.time().deltaTValue();
            pOld_[bufferIndex] = p_->boundaryFieldRef()[patchID][i];
            bufferIndex++;
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
            for (auto& cell : p_->ref())
            {
                cell = buffer[bufferIndex++];
            }
        }
        else
        {
            for (const auto& cellSetName : cellSetNames_)
            {
                cellSet overlapRegion(p_->mesh(), cellSetName);
                const labelList& cells = overlapRegion.toc();

                for (const auto& currentCell : cells)
                {
                    p_->ref()[currentCell] = buffer[bufferIndex++];
                }
            }
        }
    }

    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        scalarField* valuePatchPtr = &p_->boundaryFieldRef()[patchID];
        if (isA<coupledPressureFvPatchField>(p_->boundaryFieldRef()[patchID]))
        {
            valuePatchPtr = &refCast<coupledPressureFvPatchField>(
                                 p_->boundaryFieldRef()[patchID])
                                 .refValue();
        }
        scalarField& valuePatch = *valuePatchPtr;

        forAll(p_->boundaryFieldRef()[patchID], i)
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