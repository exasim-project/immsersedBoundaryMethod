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
 
\*---------------------------------------------------------------------------*/
 
#include "ibmGaussLaplacianScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvMatrices.H"
#include "unitConversion.H"
#include "searchableSurfaces.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace fv
{

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
template<class Type, class GType>
void ibmGaussLaplacianScheme<Type, GType>::createGeoData(const fvMesh& mesh)
{
    const dynamicFvMesh& dynMesh = static_cast<const dynamicFvMesh&>(mesh);

    name_ = mesh.time().controlDict().get<word>("geometry");

    if(!ibmGeometryData::geoDataTable().found(name_))
    {
        ibmGeometryData::geoDataTable().insert(name_, ibmGeometryData(dynMesh));
    }
}


template<class Type, class GType>
ibmGeometryData& ibmGaussLaplacianScheme<Type, GType>::geoData()
{
    return ibmGeometryData::geoDataTable()[name_];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type>>
ibmGaussLaplacianScheme<Type, GType>::fvmLaplacianUncorrected
(
    const surfaceScalarField& gammaMagSf,
    const surfaceScalarField& deltaCoeffs,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    surfaceScalarField surfCentreCorrection(gammaMagSf*(1.0-geoData().deltaQuotient())*mag(geoData().vofGrad())); 
    surfaceScalarField faceCorrection(gammaMagSf*mag(geoData().vofGrad()));

    forAll(faceCorrection, facei){    // Wenn Unterschied zwischen upper und lower gemacht werden soll.
        if(geoData().ownerCorrection()[facei]){    //Wenn Owner korrigiert wird, muss faceCorrection in lower.
            fvm.lower()[facei] = -faceCorrection.primitiveField()[facei];
        }else{
            fvm.upper()[facei] = -faceCorrection.primitiveField()[facei];
        }
    }
    //fvm.upper() = -faceCorrection.primitiveField(); // Wenn kein Unterschied zwischen upper und lower gemacht werden soll.

    volScalarField volCentreCorrection
    (
        IOobject
        (
            "volCentreCorrection",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("volCentreCorrection", dimless, 0.0)
    );

    forAll(mesh.C(), celli){
        forAll(mesh.cells()[celli], i){
            const label& facei = mesh.cells()[celli][i];
            if (mesh.isInternalFace(facei)) //Damit keine Werte auf den Rändern berücksichtigt werden.
            {
                volCentreCorrection[celli] += surfCentreCorrection[facei] * geoData().vofField()[celli];
            }
        }
    }

    if(mesh.time().outputTime())
    {
        geoData().vofField().write();
        volCentreCorrection.write();
    }

    fvm.diag() = volCentreCorrection.primitiveField();


    forAll(vf.boundaryField(), patchi)  // Wenn auskommentiert: RB werden am Gitterrand nicht mehr angewendet
    {
        const fvPatchField<Type>& pvf = vf.boundaryField()[patchi];
        const fvsPatchScalarField& pGamma = gammaMagSf.boundaryField()[patchi];
        const fvsPatchScalarField& pDeltaCoeffs =
            deltaCoeffs.boundaryField()[patchi];
 
        if (pvf.coupled())
        {
            fvm.internalCoeffs()[patchi] =
                pGamma*pvf.gradientInternalCoeffs(pDeltaCoeffs);
            fvm.boundaryCoeffs()[patchi] =
               -pGamma*pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
        }
        else
        {
            fvm.internalCoeffs()[patchi] = pGamma*pvf.gradientInternalCoeffs();
            fvm.boundaryCoeffs()[patchi] = -pGamma*pvf.gradientBoundaryCoeffs();
        }
    }
 
    return tfvm;
}
 

template<class Type, class GType>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
ibmGaussLaplacianScheme<Type, GType>::gammaSnGradCorr
(
    const surfaceVectorField& SfGammaCorr,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();
 
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tgammaSnGradCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "gammaSnGradCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            SfGammaCorr.dimensions()
           *vf.dimensions()*mesh.deltaCoeffs().dimensions()
        )
    );
    tgammaSnGradCorr.ref().oriented() = SfGammaCorr.oriented();
 
    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        tgammaSnGradCorr.ref().replace
        (
            cmpt,
            fvc::dotInterpolate(SfGammaCorr, fvc::grad(vf.component(cmpt)))
        );
    }
 
    return tgammaSnGradCorr;
}
 
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ibmGaussLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{   
    const fvMesh& mesh = this->mesh();
 
    tmp<GeometricField<Type, fvPatchField, volMesh>> tLaplacian
    (
        fvc::div(this->tsnGradScheme_().snGrad(vf)*mesh.magSf())
    );
 
    tLaplacian.ref().rename("laplacian(" + vf.name() + ')');
 
    return tLaplacian;
}
 
 
template<class Type, class GType>
tmp<fvMatrix<Type>>
ibmGaussLaplacianScheme<Type, GType>::fvmLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();
 
    const surfaceVectorField Sn(mesh.Sf()/mesh.magSf());
 
    const surfaceVectorField SfGamma(mesh.Sf() & gamma);
    const GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn
    (
        SfGamma & Sn
    );
    const surfaceVectorField SfGammaCorr(SfGamma - SfGammaSn*Sn);
 
    tmp<fvMatrix<Type>> tfvm = fvmLaplacianUncorrected
    (
        SfGammaSn,
        this->tsnGradScheme_().deltaCoeffs(vf),
        vf
    );
    fvMatrix<Type>& fvm = tfvm.ref();
    
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tfaceFluxCorrection
        = gammaSnGradCorr(SfGammaCorr, vf);
    
    if (this->tsnGradScheme_().corrected())
    {
        tfaceFluxCorrection.ref() +=
            SfGammaSn*this->tsnGradScheme_().correction(vf);
    }
 
    fvm.source() -= mesh.V()*fvc::div(tfaceFluxCorrection())().primitiveField();
 
    if (mesh.fluxRequired(vf.name()))
    {
        fvm.faceFluxCorrectionPtr() = tfaceFluxCorrection.ptr();
    }
 
    return tfvm;
}
 

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh>>
ibmGaussLaplacianScheme<Type, GType>::fvcLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = this->mesh();
 
    const surfaceVectorField Sn(mesh.Sf()/mesh.magSf());
    const surfaceVectorField SfGamma(mesh.Sf() & gamma);
    const GeometricField<scalar, fvsPatchField, surfaceMesh> SfGammaSn
    (
        SfGamma & Sn
    );
    const surfaceVectorField SfGammaCorr(SfGamma - SfGammaSn*Sn);
 
    tmp<GeometricField<Type, fvPatchField, volMesh>> tLaplacian
    (
        fvc::div
        (
            SfGammaSn*this->tsnGradScheme_().snGrad(vf)
          + gammaSnGradCorr(SfGammaCorr, vf)
        )
    );
 
    tLaplacian.ref().rename
    (
        "ibmLaplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
 
    return tLaplacian;
}
 
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
} // End namespace fv
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
} // End namespace Foam
 
// ************************************************************************* //