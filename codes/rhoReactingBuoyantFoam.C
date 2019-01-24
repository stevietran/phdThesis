#include "postProcess.H"
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createControl.H"
#include "createTimeControls.H"
#include "initContinuityErrs.H"
#include "createFields.H"
#include "createFieldRefs.H"
#include "createFvOptions.H"
turbulence->validate();
while (runTime.run())
{
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setDeltaT.H"
    runTime++;
    Info<< "Time = " << runTime.timeName() << nl << endl;
    #include "rhoEqn.H"
    // Pressure-velocity PIMPLE corrector loop
    while (pimple.loop())
    {
        #include "UEqn.H"
        #include "YEqn.H"
        #include "EEqn.H"
        // Pressure corrector loop
        while (pimple.correct())
        {
            #include "pEqn.H"
        }
        if (pimple.turbCorr())
        {
            turbulence->correct();
        }
    }
    rho = thermo.rho();
    runTime.write();
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
}
Info<< "End\n" << endl;
return 0;

