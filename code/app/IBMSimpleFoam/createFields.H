Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

const dynamicFvMesh& dynMesh = static_cast<const dynamicFvMesh&>(mesh);
const word geoDataName = mesh.time().controlDict().get<word>("geometry");
if(!ibmGeometryData::geoDataTable().found(geoDataName))
{
    ibmGeometryData::geoDataTable().insert(geoDataName, ibmGeometryData(dynMesh));
}

//  #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   //
// Hier phiMask-Feld erzeugen!
surfaceScalarField geometryMask
(
    IOobject
    (
        "geometryMask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("one", dimless, 1.0)
);

forAll(phi, facei)
{
    if( ibmGeometryData::geoDataTable()[geoDataName].vofField()[mesh.faceOwner()[facei]] == 0 ||
        ibmGeometryData::geoDataTable()[geoDataName].vofField()[mesh.faceNeighbour()[facei]] == 0)
    {   // Alle faces die entweder vollständig im Solid liegen oder ein Interface dazwischen haben.
        geometryMask[facei] = 0.0;
    }
}
//  #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   //


autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
