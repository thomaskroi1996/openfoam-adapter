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
        Foam::fvMesh mesh_;
        Foam::volScalarField p_;
        Foam::volScalarField pOld_;
        Foam::volScalarField dpdt_;
        bool firstStep_ = true;

    public:
        //Constructor
        PressureTemporalDerivative(
        const Foam::fvMesh& mesh,
        const std::string name_dpdt);

        //Compute dpdt field
        void compute();
            
        //Return dpdt field
        const Foam::volScalarField& field() const;

        //Write Temporal Derivative of Pressure into the buffer
        std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) final;

        //Read the Temporal Derivative of Pressure from the buffer
        void read(double* buffer, bool meshConnectivity) final;

        bool isLocationTypeSupported(const bool meshConnectivity) const final;

        //- Get the name of the current data field
        std::string getDataName() const final;

};

}

}

#endif