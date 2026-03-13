/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Central-differencing interpolation scheme class

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "IBMlinearP.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceInterpolationScheme(IBMlinearP)
}

template<class Type>
template<class SFType>
Foam::tmp<Foam::GeometricField<typename Foam::innerProduct<typename SFType::value_type, Type>::type,Foam::fvsPatchField,Foam::surfaceMesh>>
Foam::IBMlinearP<Type>::dotInterpolate
(
    const SFType& Sf,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<surfaceScalarField>& tlambdas
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "Interpolating "
            << vf.type() << " "
            << vf.name()
            << " from cells to faces without explicit correction"
            << endl;
    }

    typedef typename Foam::innerProduct<typename SFType::value_type, Type>::type
        RetType;

    const surfaceScalarField& lambdas = tlambdas();

    const Field<Type>& vfi = vf;
    const scalarField& lambda = lambdas;

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<RetType, fvsPatchField, surfaceMesh>> tsf
    (
        new GeometricField<RetType, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            Sf.dimensions()*vf.dimensions()
        )
    );
    GeometricField<RetType, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    Field<RetType>& sfi = sf.primitiveFieldRef();

    const typename SFType::Internal& Sfi = Sf();

    //IBM:
    Field<Type> newVfi(vfi);
    const dynamicFvMesh& dynMesh = static_cast<const dynamicFvMesh&>(mesh);
    const word geoDataName = mesh.time().controlDict().get<word>("geometry");
    if(!ibmGeometryData::geoDataTable().found(geoDataName))
    {
        ibmGeometryData::geoDataTable().insert(geoDataName, ibmGeometryData(dynMesh));
    }
    Foam::ibmGeometryData& geoData = ibmGeometryData::geoDataTable()[geoDataName];
    //surfaceScalarField& deltaQuotient = geoData.deltaQuotient();
    volScalarField& vofField = geoData.vofField();
    surfaceScalarField& vofGrad = geoData.vofGrad();
    
    for (label fi=0; fi<P.size(); fi++)
    {
        Type vP = newVfi[P[fi]];
        Type vN = newVfi[N[fi]];

        // Flussinterpolation, mit zeroGradient RB auf IBM face:
        if(vofGrad[fi] != 0)
        {
            if(vofField[P[fi]] == 1) //P ist Fluid, N ist Solid
            {
                vN = vP;                    
            }else{
                vP = vN;
            }
        }
        //  #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   //

        sfi[fi] = Sfi[fi] & (lambda[fi]*(vP - vN) + vN);

        // Fluss in IBM zu 0 setzen:
        if( vofField[P[fi]] == 0 && vofField[N[fi]] == 0){
            sfi[fi] *= 0;
        }
        //  #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   //
    }

    // Interpolate across coupled patches using given lambdas

    typename GeometricField<RetType, fvsPatchField, surfaceMesh>::
        Boundary& sfbf = sf.boundaryFieldRef();

    forAll(lambdas.boundaryField(), pi)
    {
        const fvsPatchScalarField& pLambda = lambdas.boundaryField()[pi];
        const typename SFType::Patch& pSf = Sf.boundaryField()[pi];
        fvsPatchField<RetType>& psf = sfbf[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            psf =
                pSf
              & (
                    pLambda*vf.boundaryField()[pi].patchInternalField()
                  + (1.0 - pLambda)*vf.boundaryField()[pi].patchNeighbourField()
                );
        }
        else
        {
            psf = pSf & vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();

//    tsf.ref().oriented() = Sf.oriented();

    return tsf;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::IBMlinearP<Type>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<surfaceScalarField>& tlambdas
)
{
    return dotInterpolate(geometricOneField(), vf, tlambdas);
}

template<class Type>
Foam::tmp<Foam::GeometricField<typename Foam::innerProduct<Foam::vector, Type>::type,Foam::fvsPatchField,Foam::surfaceMesh>>
Foam::IBMlinearP<Type>::dotInterpolate
(
    const surfaceVectorField& Sf,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "Interpolating "
            << vf.type() << " "
            << vf.name()
            << " from cells to faces"
            << endl;
    }
    
    tmp
    <
        GeometricField
        <
            typename Foam::innerProduct<Foam::vector, Type>::type,
            fvsPatchField,
            surfaceMesh
        >
    > tsf = dotInterpolate(Sf, vf, weights(vf));
    
    tsf.ref().oriented() = Sf.oriented();

    /*if (corrected())
    {
        Info << "Bin in corrected" << endl;
        tsf.ref() += Sf & correction(vf);
    }*/

    return tsf;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::IBMlinearP<Type>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "Interpolating "
            << vf.type() << " "
            << vf.name()
            << " from cells to faces"
            << endl;
    }

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
        = interpolate(vf, weights(vf));

    /*if (corrected())
    {
        tsf.ref() += correction(vf);
    }*/

    return tsf;
}

template<class Type>
Foam::tmp<Foam::GeometricField<typename Foam::innerProduct<Foam::vector, Type>::type,Foam::fvsPatchField,Foam::surfaceMesh>>
Foam::IBMlinearP<Type>::dotInterpolate
(
const surfaceVectorField& Sf,
const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
) const
{
tmp
<
    GeometricField
    <
        typename Foam::innerProduct<Foam::vector, Type>::type,
        fvsPatchField,
        surfaceMesh
    >
> tSfDotinterpVf = dotInterpolate(Sf, tvf());

tvf.clear();
return tSfDotinterpVf;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::IBMlinearP<Type>::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
) const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tinterpVf
        = interpolate(tvf());
    tvf.clear();
    return tinterpVf;
}

// ************************************************************************* //
