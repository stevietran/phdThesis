{
    // Reading user defined Sc number
    scalar Sct = runTime.controlDict().lookupOrDefault<scalar>("Sc", 1.0);
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);
    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];
            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff()/Sct, Yi)
             ==
                reaction->R(Yi)
              + fvOptions(rho, Yi)
            );
            YiEqn.relax();
            fvOptions.constrain(YiEqn);
            YiEqn.solve(mesh.solver("Yi"));
            fvOptions.correct(Yi);
            Yi.max(0.0);
            Yt += Yi;
        }
    }
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
