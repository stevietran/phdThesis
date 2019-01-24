// Local references
const alphaField& alpha = this->alpha_;
const rhoField& rho = this->rho_;
const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
const volVectorField& U = this->U_;
volScalarField& nut = this->nut_;

volScalarField::Internal divU
(
    fvc::div(fvc::absolute(this->phi(), U))().v()
);
tmp<volTensorField> tgradU = fvc::grad(U);
volScalarField::Internal G
(
    this->GName(),
    nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
);
tgradU.clear();
// Update epsilon and G at the wall
epsilon_.boundaryFieldRef().updateCoeffs();
// Dissipation equation
tmp<fvScalarMatrix> epsEqn
(
    fvm::ddt(alpha, rho, epsilon_)
  + fvm::div(alphaRhoPhi, epsilon_)
  - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
 ==
    C1_*alpha()*rho()*G*epsilon_()/k_()
  - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilon_)
  - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
  + epsilonSource()
  + fvOptions(alpha, rho, epsilon_)
);
epsEqn.ref().relax();
epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
solve(epsEqn);
bound(epsilon_, this->epsilonMin_);
// Turbulent kinetic energy equation
tmp<fvScalarMatrix> kEqn
(
    fvm::ddt(alpha, rho, k_)
  + fvm::div(alphaRhoPhi, k_)
  - fvm::laplacian(alpha*rho*DkEff(), k_)
 ==
    alpha()*rho()*G
  - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
  - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
  + kSource()
  + fvOptions(alpha, rho, k_)
);
kEqn.ref().relax();
solve(kEqn);
bound(k_, this->kMin_);
correctNut();

