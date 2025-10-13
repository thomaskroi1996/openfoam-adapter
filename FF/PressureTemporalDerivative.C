#include "PressureTemporalDerivative.H"
#include "coupledPressureFvPatchField.H"

using namespace Foam;

preciceAdapter::FF::PressureTemporalDerivative::PressureTemporalDerivative(
    const Foam::fvMesh& mesh,
    const std::string name_dpdt) 
:   
    mesh_(mesh),
    p_(
    const_cast<volScalarField*>(
        &mesh.lookupObject<volScalarField>(name_dpdt))), //issue here, we should probably have name_dpdt be "p", since there is no dp/dt in the data yet
    pOld_(p_.internalField().size()),
    dpdt_(p_.internalField().size())
{}

void preciceAdapter::FF::PressureTemporalDerivative::compute()
{
    if (firstStep_)
    {
        pOld_ = p_;
        firstStep_ = false;
    }

    dpdt_.internalField() = (p_.internalField() - pOld_.internalField()) / mesh_.time().deltaTValue();

    pOld_ = p_; // update previous pressure
}

const Foam::volScalarField& preciceAdapter::FF::PressureTemporalDerivative::field() const
{
    return dpdt_;
}

std::size_t preciceAdapter::FF::PressureTemporalDerivative::write(double* buffer,
                                              bool meshConnectivity,
                                              const unsigned int dim)
{
    //JUST COMPUTE DP/DT HERE??? compute() hmmmmmmmm

    int bufferIndex = 0;

    if (meshConnectivity)
        return 0;

    for (const auto& cell : p_->internalField())
        {
            buffer[bufferIndex++] = cell;
        }
    return bufferIndex;
}

std::string preciceAdapter::FF::PressureTemporalDerivative::getDataName() const
{
    return "PressureTemporalDerivative";
}