#include "PressureTemporalDerivative.H"
#include "coupledPressureFvPatchField.H"

using namespace Foam;

preciceAdapter::FF::PressureTemporalDerivative::PressureTemporalDerivative(
    const Foam::fvMesh& mesh,
    const std::string name) 
:   
    mesh_(mesh),
    p_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(name))),
    pOld_(IOobject
        (
            "p", // Field name
            runTime.timeName(), // Current time
            mesh, // The mesh object
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh),
    dpdt_(IOobject
        (
            "dpdt", // Field name
            runTime.timeName(), // Current time
            mesh, // The mesh object
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh)
    
{
    dataType_ = scalar;
}

void preciceAdapter::FF::PressureTemporalDerivative::compute()
{
    std::cout << "Computing dpdt..." << std::endl;
    if (firstStep_)
    {
        pOld_ = p_; 
        firstStep_ = false;
    }

    std::cout << dpdt_.internalField()

    dpdt_.internalField() = (p_.internalField() - pOld_.internalField()) / mesh_.time().deltaTValue();

    pOld_ = p_; // update previous pressure
}

void preciceAdapter::FF::PressureTemporalDerivative::compute(double* buffer, double* pOldFromBuffer_){
    //probably not useful, but just in case
}

const Foam::volScalarField& preciceAdapter::FF::PressureTemporalDerivative::field() const
{
    return dpdt_;
}

std::size_t preciceAdapter::FF::PressureTemporalDerivative::write(double* buffer,
                                              bool meshConnectivity,
                                              const unsigned int dim)
{
    //we don't want to compute dpdt from the buffer, we just use member vars
    compute();

    std::cout << "dpdt computed." << std::endl;

    int bufferIndex = 0;

    if (this->locationType_ == LocationType::volumeCenters)
    {
        if (cellSetNames_.empty())
        {
            for (const auto& cell : dpdt_->internalField())
            {
                buffer[bufferIndex++] = cell;
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
                    // Copy dpdt into the buffer
                    buffer[bufferIndex++] = dpdt_->internalField()[currentCell];

                }
            }
        }
    }

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(dpdt_->boundaryFieldRef()[patchID], i)
        {
            // Copy dpdt into the buffer
            buffer[bufferIndex++] =
                dpdt_->boundaryFieldRef()[patchID][i];
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
                    // Copy the pressure into the buffer
                    dpdt_->ref()[currentCell] = buffer[bufferIndex++];
                }
            }
        }
    }

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the pressure value boundary patch
        scalarField* valuePatchPtr = &dpdt_->boundaryFieldRef()[patchID];
        if (isA<coupledPressureFvPatchField>(dpdt_->boundaryFieldRef()[patchID]))
        {
            valuePatchPtr = &refCast<coupledPressureFvPatchField>(
                                 dpdt_->boundaryFieldRef()[patchID])
                                 .refValue();
        }
        scalarField& valuePatch = *valuePatchPtr;

        // For every cell of the patch
        forAll(dpdt_->boundaryFieldRef()[patchID], i)
        {
            // Set the pressure as the buffer value
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