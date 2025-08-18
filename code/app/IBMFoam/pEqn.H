volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

// phiHbyA auf Surface zu 0 setzen, wenn Interface zwischen Zellen liegt!
forAll(phiHbyA, facei){
    if(ibmGeometryData::geoDataTable()[geoDataName].vofGrad()[facei] != 0){
        phiHbyA[facei] = 0;
    }
}

if (pimple.ddtCorr())
{
    phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf));
}
else
{
    phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU));
}

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

//surfaceScalarField phirAtU erstellen welches auf den Surface 0 sind bei dem Interface zwischen Zellen liegt!
surfaceScalarField rAtUFace("rAtUFace", fvc::interpolate(rAtU()));
forAll(rAtUFace, facei){
    if(ibmGeometryData::geoDataTable()[geoDataName].vofGrad()[facei] != 0){
        if(ibmGeometryData::geoDataTable()[geoDataName].ownerCorrection()[facei]){
            rAtUFace[facei] = rAtU()[mesh.faceOwner()[facei]]; //face wert wird auf 1/A[korrigierte Zelle] gesetzt (darf nicht zu klein sein)
        }else{
            rAtUFace[facei] = rAtU()[mesh.faceNeighbour()[facei]];
        }
    }
}

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    /*fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );*/

    // pEqn mit korrigiertem rAtU
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtUFace, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);
