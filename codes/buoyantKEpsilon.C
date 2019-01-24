template<class BasicTurbulenceModel>
tmp<volScalarField>
buoyantKEpsilon<BasicTurbulenceModel>::Gcoef() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");
    return
        (Cg_*this->Cmu_)*this->alpha_*this->k_*(g & fvc::grad(this->rho_))
       /(this->epsilon_ + this->epsilonMin_);
}

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
buoyantKEpsilon<BasicTurbulenceModel>::kSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");
    if (mag(g.value()) > SMALL)
    {
        return -fvm::SuSp(Gcoef(), this->k_);
    }
    else
    {
        return kEpsilon<BasicTurbulenceModel>::kSource();
    }
}

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
buoyantKEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");
    if (mag(g.value()) > SMALL)
    {
        vector gHat(g.value()/mag(g.value()));
        volScalarField v(gHat & this->U_);
        volScalarField u
        (
            mag(this->U_ - gHat*v)
          + dimensionedScalar("SMALL", dimVelocity, SMALL)
        );
        return -fvm::SuSp(this->C1_*tanh(mag(v)/u)*Gcoef(), this->epsilon_);
    }
    else
    {
        return kEpsilon<BasicTurbulenceModel>::epsilonSource();
    }
}

