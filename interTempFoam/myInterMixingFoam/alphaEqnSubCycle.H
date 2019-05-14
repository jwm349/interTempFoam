if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", rhoPhi.dimensions(), 0)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

{
    volScalarField rhoNew(alpha1*rho1 + alpha2*rho2 + alpha3*rho3);

    //solve(fvm::ddt(rho) + fvc::div(rhoPhi));
    //Info<< "density error = "
    //    << max((mag(rho - rhoNew)/mag(rhoNew))().primitiveField()) << endl;

    rho == rhoNew;
//ADDED---------------------->
	volScalarField rhoCpNew(alpha1*rho1*cp1 + alpha2*rho2*cp2 + alpha3*rho3*cp3);
	rhoCp == rhoCpNew;
//<--------------------ADDED
}
