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
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    eqn += fvm::laplacian(nu_, eqn.psi(), "laplacian(IBM," + eqn.psi().name() + ")");
}


void Foam::fv::secondOrderIBM::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    eqn += rho*fvm::laplacian(nu_, eqn.psi(), "laplacian(IBM," + eqn.psi().name() + ")");
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

    Info << "Correcting U and phi." << endl;

    volScalarField& fluid = ibmGeometryData::geoDataTable()[geometryName_].vofField();
    eqn.diag() += max(eqn.diag())*(1-fluid)*100;   // Faktor muss erhöht werden (Hier: 100)
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
