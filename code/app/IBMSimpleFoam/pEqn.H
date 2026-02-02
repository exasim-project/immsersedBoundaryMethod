{
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    
    // Wenn phiHbyA = 0 sein soll in und um IBM:
    //phiHbyA *= geometryMask_face;

    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p); //Behandelt boundary patches und macht phiHbyA kompatibel zu den BCs.

    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        rAtU = 1.0/(1.0/rAU - UEqn.H1());
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf(); //Sollte bei zeroGradient eine 0 addieren.

        // Wenn phiHbyA = 0 sein soll in und um IBM:
        //phiHbyA *= geometryMask_face;

        HbyA -= (rAU - rAtU())*fvc::grad(p); //Sollte bei zeroGradient eine 0 addieren.
    }

    tUEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU(), MRF); //Bringt p RBs in Matrix ein.

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

            // Wenn phiHbyA = 0 sein soll in und um IBM:
            //phi *= geometryMask_face;
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
