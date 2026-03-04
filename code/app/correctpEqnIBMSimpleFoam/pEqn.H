{
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    
    // Wenn phiHbyA = 0 sein soll in und um IBM:
    phiHbyA *= geometryMask_face;

    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p); //Behandelt boundary patches und macht phiHbyA kompatibel zu den BCs.

    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        rAtU = 1.0/(1.0/rAU - UEqn.H1());
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf(); //fvc::snGrad(p) = 0 bei zeroGradient.

        // Wenn phiHbyA = 0 sein soll in und um IBM:
        phiHbyA *= geometryMask_face;

        HbyA -= (rAU - rAtU())*fvc::grad(p); //fvc::grad(.) teilt durch Zellvolumen. Sieht also die ganze Zelle anstatt Fluidanteil.
    }

    tUEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU(), MRF); //Bringt p RBs in Matrix ein.

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA) //fvc::div(.) teilt durch Zellvolumen. Sieht also die ganze Zelle anstatt Fluidanteil.
        );

        scalarField& diag   = pEqn.diag();
        scalarField& source = pEqn.source();
        scalarField& upper  = pEqn.upper();
        scalarField& lower  = pEqn.lower();

        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();

        // 1) Remove coupling ONLY across interface faces
        forAll(nei, facei) // internal faces only: neighbour list length
        {
            const label co = own[facei];
            const label cn = nei[facei];

            if (fluid[co] != fluid[cn])
            {
                const scalar u = upper[facei];
                const scalar l = lower[facei];

                // remove their contribution from diagonal (ldu-consistent)
                diag[co] += u;
                diag[cn] += l;

                // kill coupling
                upper[facei] = 0.0;
                lower[facei] = 0.0;
            }
        }

        // 2) Enforce p = pRefValue in all zone cells (turn rows into identity)
        //    Optional but usually what you want if zone should be "constant pressure volume"
        const cellList& cells = mesh.cells();

        for (const label c : solidCellId)
        {
            // Remove any remaining couplings incident to c (also inside-zone couplings),
            // otherwise the row won't be a clean identity.
            for (const label facei : cells[c])
            {
                if (!mesh.isInternalFace(facei)) continue;

                const label co = own[facei];
                const label cn = nei[facei];

                if (co == c)
                {
                    const scalar u = upper[facei];
                    diag[co] += u;      // undo diag contribution
                    upper[facei] = 0.0;

                    // also remove reverse coupling in neighbour eqn
                    const scalar l = lower[facei];
                    diag[cn] += l;
                    lower[facei] = 0.0;
                }
                else if (cn == c)
                {
                    const scalar l = lower[facei];
                    diag[cn] += l;
                    lower[facei] = 0.0;

                    const scalar u = upper[facei];
                    diag[co] += u;
                    upper[facei] = 0.0;
                }
            }

            // Identity row: p_c = pRefValue
            diag[c]   = 1.0;
            source[c] = pRefValue;
        }

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux(); //pEqn.flux() = 0 bei zeroGradient.

            // Wenn phiHbyA = 0 sein soll in und um IBM:
            phi *= geometryMask_face;
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U = HbyA - rAtU()*fvc::grad(p); //fvc::grad(.) teilt durch Zellvolumen. Sieht also die ganze Zelle anstatt Fluidanteil.
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
