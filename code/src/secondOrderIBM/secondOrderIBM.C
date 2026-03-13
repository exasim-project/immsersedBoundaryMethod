/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "secondOrderIBM.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "gravityMeshObject.H"

#include "fvm.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(secondOrderIBM, 0);
    addToRunTimeSelectionTable(option, secondOrderIBM, dictionary);
}
}

Foam::tmp<Foam::fv::laplacianScheme<Foam::scalar, Foam::scalar>> Foam::fv::secondOrderIBM::ibmSchemeScalar()
{
    // Read interpolation scheme from fvSchemes dictionary
    const dictionary& surfInterpDict =
        mesh().schemesDict().subDict("interpolationSchemes");

    const word key("IBM");

    // Look up entry for IBM, fall back to "default"
    ITstream& interpData =
        surfInterpDict.found(key)
        ? surfInterpDict.lookup(key)
        : surfInterpDict.lookup("default");

    // Construct interpolation scheme for scalars
    tmp<surfaceInterpolationScheme<scalar>> tinterp =
        surfaceInterpolationScheme<scalar>::New(mesh(), interpData);
    

    // Read snGrad scheme from fvSchemes dictionary
    const dictionary& snGradDict = mesh().schemesDict().subDict("snGradSchemes");

    // Look up scheme name for IBM, fall back to "default"
    ITstream& schemeData =
        snGradDict.found(key)
        ? snGradDict.lookup(key)
        : snGradDict.lookup("default");

    // Construct snGrad scheme from dictionary (runtime selection)
    tmp<snGradScheme<Foam::scalar>> tsngrad =
        snGradScheme<Foam::scalar>::New(mesh(), schemeData);

    
    // Construct Laplacian scheme with both
    tmp<fv::laplacianScheme<scalar, scalar>> scheme
    (
        new fv::ibmGaussLaplacianScheme<scalar, scalar>
        (
            mesh(),
            tinterp,
            tsngrad
        )
    );

    static_cast<fv::ibmGaussLaplacianScheme<scalar, scalar>&>
    (
        scheme.ref()).setGeometry(dict_.get<word>("geometry")
    );
    return scheme;
}

Foam::tmp<Foam::fv::laplacianScheme<Foam::vector, Foam::scalar>> Foam::fv::secondOrderIBM::ibmSchemeVector()
{
    // Read interpolation scheme from fvSchemes dictionary
    const dictionary& surfInterpDict =
        mesh().schemesDict().subDict("interpolationSchemes");

    const word key("IBM");

    // Look up entry for IBM, fall back to "default"
    ITstream& interpData =
        surfInterpDict.found(key)
        ? surfInterpDict.lookup(key)
        : surfInterpDict.lookup("default");

    // Construct interpolation scheme for scalars
    tmp<surfaceInterpolationScheme<scalar>> tinterp =
        surfaceInterpolationScheme<scalar>::New(mesh(), interpData);
    

    // Read snGrad scheme from fvSchemes dictionary
    const dictionary& snGradDict = mesh().schemesDict().subDict("snGradSchemes");

    // Look up scheme name for IBM, fall back to "default"
    ITstream& schemeData =
        snGradDict.found(key)
        ? snGradDict.lookup(key)
        : snGradDict.lookup("default");

    // Construct snGrad scheme from dictionary (runtime selection)
    tmp<snGradScheme<Foam::vector>> tsngrad =
        snGradScheme<Foam::vector>::New(mesh(), schemeData);

    
    // Construct Laplacian scheme with both
    tmp<fv::laplacianScheme<vector, scalar>> scheme
    (
        new fv::ibmGaussLaplacianScheme<vector, scalar>
        (
            mesh(),
            tinterp,
            tsngrad
        )
    );

    static_cast<fv::ibmGaussLaplacianScheme<vector, scalar>&>
    (
        scheme.ref()).setGeometry(dict_.get<word>("geometry")
    );
    return scheme;
}

Foam::tmp<Foam::volScalarField> Foam::fv::secondOrderIBM::nu()
{
    return tmp<Foam::volScalarField>
    (
        new volScalarField 
        (
            IOobject
            (
                "nu",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            nu_ 
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::secondOrderIBM::secondOrderIBM
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(sourceName, modelType, dict, mesh),
    nu_(dict.get<dimensionedScalar>("nu")),
    penalty_(dict.getOrDefault<scalar>("penalty", 1e12)),
    geometryName_(dict.get<word>("geometry"))
{
    coeffs_.readEntry("fields", fieldNames_);

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::secondOrderIBM::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    eqn += ibmSchemeScalar()->fvmLaplacian(nu(), eqn.psi());
}

void Foam::fv::secondOrderIBM::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    eqn += ibmSchemeVector()->fvmLaplacian(nu(), eqn.psi());
}

void Foam::fv::secondOrderIBM::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    eqn += rho*ibmSchemeScalar()->fvmLaplacian(nu(), eqn.psi());
}

void Foam::fv::secondOrderIBM::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    eqn += rho*ibmSchemeVector()->fvmLaplacian(nu(), eqn.psi());
}

void Foam::fv::secondOrderIBM::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if(!ibmGeometryData::geoDataTable().found(geometryName_))
    { 
        Info << geometryName_ << " not found in " << ibmGeometryData::geoDataTable().toc() << endl;
        return; 
    }

    volScalarField& fluid = ibmGeometryData::geoDataTable()[geometryName_].vofField();
    eqn.diag() += max(eqn.diag())*(1-fluid)*penalty_;
}

void Foam::fv::secondOrderIBM::correct
(
    volVectorField& U
)
{
 
}

bool Foam::fv::secondOrderIBM::read(const dictionary& dict)
{
    NotImplemented;

    return false;
}


// ************************************************************************* //
