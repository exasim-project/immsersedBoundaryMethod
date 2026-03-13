/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::IBMlinearP

Group
    grpFvSurfaceInterpolationSchemes

Description
    Central-differencing interpolation scheme class

SourceFiles
    IBMlinearP.C

\*---------------------------------------------------------------------------*/

#ifndef IBMlinearP_H
#define IBMlinearP_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

#include "ibmGeometryData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IBMlinearP Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IBMlinearP
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- No copy assignment
        void operator=(const IBMlinearP&) = delete;


public:

    //- Runtime type information
    TypeName("IBMlinearP");


    // Constructors

        //- Construct from mesh
        IBMlinearP(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        IBMlinearP(const fvMesh& mesh, Istream&)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        IBMlinearP
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().surfaceInterpolation::weights();
        }

        //- Return the face-interpolate of the given cell field
        //  with the given weighting factors dotted with given field Sf
        template<class SFType>
        static tmp< GeometricField<typename innerProduct<typename SFType::value_type, Type>::type,fvsPatchField,surfaceMesh>>
        dotInterpolate
        (
            const SFType& Sf,
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const tmp<surfaceScalarField>& tlambdas
        );
        //- Return the face-interpolate of the given cell field
        //  with the given weighting factors
        static tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const tmp<surfaceScalarField>&
        );

        //- Return the face-interpolate of the given cell field
        //  with explicit correction dotted with given field Sf
        virtual
        tmp<GeometricField<typename innerProduct<vector, Type>::type,fvsPatchField,surfaceMesh>>
        dotInterpolate
        (
            const surfaceVectorField& Sf,
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        virtual
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<Type,fvPatchField, volMesh>&
        ) const;

        //- Return the face-interpolate of the given tmp cell field
        //  with explicit correction dotted with given field Sf
        tmp<GeometricField<typename innerProduct<vector, Type>::type,fvsPatchField,surfaceMesh>>
        dotInterpolate
        (
            const surfaceVectorField& Sf,
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&
        ) const;
        //- Return the face-interpolate of the given tmp cell field
        //  with explicit correction
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&
        ) const;
};

template<>
tmp<GeometricField<typename innerProduct<vector, scalar>::type,fvsPatchField,surfaceMesh>>
IBMlinearP<scalar>::dotInterpolate
(
    const surfaceVectorField& Sf,
    const GeometricField<scalar, fvPatchField, volMesh>&
) const
{
    NotImplemented;

    return
        tmp
        <
            GeometricField
            <
                typename innerProduct<vector, scalar>::type,
                fvsPatchField,
                surfaceMesh
            >
        >
        (
            GeometricField
            <
                typename innerProduct<vector, scalar>::type,
                fvsPatchField,
                surfaceMesh
            >::null()
        );
}

/*
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
linearInterpolate(const GeometricField<Type, fvPatchField, volMesh>& vf)
{
    return surfaceInterpolationScheme<Type>::interpolate
    (
        vf,
        vf.mesh().surfaceInterpolation::weights()
    );
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
linearInterpolate(const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tinterp =
        linearInterpolate(tvf());
    tvf.clear();
    return tinterp;
}*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
