{
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

    /*volScalarField hOperator(UEqn.H1());

    forAll(phiHbyA, facei)
    {
        if(ibmGeometryData::geoDataTable()[geoDataName].vofGrad()[facei] =! 0)
        {
            Info << ibmGeometryData::geoDataTable()[geoDataName].vofGrad()[facei] << ", " << hOperator[mesh.faceOwner()[facei]] << ", " << hOperator[mesh.faceNeighbour()[facei]] << endl;
        }
    }*/

    // phiHbyA überall 0, wo Solid oder Interface liegt.
    phiHbyA *= geometryMask;

    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p);

    // phiHbyA überall 0, wo Solid oder Interface liegt.
    phiHbyA *= geometryMask;

    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        rAtU = 1.0/(1.0/rAU - UEqn.H1());
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();

        // phiHbyA überall 0, wo Solid oder Interface liegt.
        phiHbyA *= geometryMask;

        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    tUEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU(), MRF);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();

            // phiHbyA überall 0, wo Solid oder Interface liegt.
            phi *= geometryMask;
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U = HbyA - rAtU()*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
