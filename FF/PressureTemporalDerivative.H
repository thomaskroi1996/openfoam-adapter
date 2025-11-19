#ifndef FF_PRESSURE_TEMPORAL_DERIVATIVE_H
#define FF_PRESSURE_TEMPORAL_DERIVATIVE_H

#include "CouplingDataUser.H"
#include "fvCFD.H"

namespace preciceAdapter
{
namespace FF
{

class PressureTemporalDerivative : public CouplingDataUser
{
private:
    const  Foam::fvMesh& mesh_; // needed to access mesh outside of constructor (for deltaTValue() method, dt)
    Foam::volScalarField* p_;
    Foam::volScalarField* pOld_;
    Foam::volScalarField* dpdt_;
    std::vector<double> pOldFromBuffer;
    bool firstStep_ = true;

public:
    //Constructor
    PressureTemporalDerivative(
        const Foam::fvMesh& mesh,
        const std::string namePTD);

    //Write Temporal Derivative of Pressure into the buffer
    std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) final;

    //Read the Temporal Derivative of Pressure from the buffer
    void read(double* buffer, const unsigned int dim) final;

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;
};

}

}

#endif

