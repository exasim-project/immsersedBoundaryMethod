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
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace fv
{
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type>>
ibmGaussLaplacianScheme<Type, GType>::fvmLaplacianUncorrected
(
    const surfaceScalarField& gammaMagSf,
    const surfaceScalarField& deltaCoeffs,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    //Info << "Bin im ibmGaussLapacianScheme" << endl;
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions() //gammaMagSf = gamma*(surface area)
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.upper() = deltaCoeffs.primitiveField()*gammaMagSf.primitiveField();
    fvm.negSumDiag();

    const fvMesh& mesh = this->mesh();
 
    const surfaceVectorField Sn(mesh.Sf()/mesh.magSf()); //warum keine Referenz?
    
    volScalarField vofField_    //später mit Indicator definieren
    (
        IOobject
        (
            "vof",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("vof", dimless, 1.0)
    );

    forAll(vofField_,i){
        if(i<250*5 || i>=(vofField_.size()-250*5)){vofField_[i] = 0.0;}
    }

    const scalar& delta_quotient_ = 0.01;
    
    surfaceScalarField vofGrad = fvc::snGrad(vofField_);

    //Korrekturterm:

    //Annahme: vofGrad ist immer entweder 0 oder 1/delx
    // Falsche Berechnung der Korrekturen
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    /*surfaceScalarField surfCentreCorrection(gammaMagSf*(1-delta_quotient_)*mag(vofGrad));
    surfaceScalarField faceCorrection(gammaMagSf*mag(vofGrad)); //deltaCoeffs*gammaMagSf*(vofGrad/deltaCoeffs)

    volScalarField volCentreCorrection(vofField_);
    forAll(mesh.C(), celli){
        volCentreCorrection[celli] = 0.0;
        forAll(mesh.cells()[celli], i){
            const label& facei = mesh.cells()[celli][i];
            volCentreCorrection[celli] += surfCentreCorrection[facei];
            //Info << mesh.cells()[celli][i] << endl;
        }
    }
    fvm.diag() += volCentreCorrection.primitiveField()*vofField_.primitiveField();*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //Ansonsten:
    // Korrekte Berechnung

    surfaceScalarField surfCentreCorrection(deltaCoeffs*gammaMagSf*(1-delta_quotient_));
    surfaceScalarField faceCorrection(deltaCoeffs*gammaMagSf);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    forAll(mesh.Sf(), facei){
        if(mag(vofGrad[facei])<1e-8){
            surfCentreCorrection[facei] = 0.0;
            faceCorrection[facei] = 0.0;
        }else{
            if(vofField_[mesh.faceOwner()[facei]]==1){ //Korrektur Fluid oder Solid zuweisen?
                fvm.diag()[mesh.faceOwner()[facei]] += surfCentreCorrection[facei];
            }else{
                fvm.diag()[mesh.faceNeighbour()[facei]] += surfCentreCorrection[facei];
            }
        }
    }
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //Oder ohne zweite if-Abfrage:
    // Simulation unterbricht mit: [stack trace] Floating point exception
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    /*forAll(mesh.Sf(), facei){
        if(mag(vofGrad[facei])<1e-8){
            surfCentreCorrection[facei] = 0.0;
            faceCorrection[facei] = 0.0;
        }
    }

    volScalarField volCentreCorrection(vofField_);

    forAll(mesh.C(), celli){
        volCentreCorrection[celli] = 0.0;
        forAll(mesh.cells()[celli], facei){
            const label& facei = mesh.cells()[celli][i];
            volCentreCorrection[celli] += surfCentreCorrection[facei];
        }
    }
    
    fvm.diag() += volCentreCorrection.primitiveField()*vofField_;*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fvm.upper() -= faceCorrection.primitiveField();
    
    Info << "IBM Matrix berechnet!" << endl;

    /*
    Owner ist bei einem face immer der mit dem niedrigeren Index, neighbor der mit dem höheren.
    Gehe faces durch und schaue, ob zwischen owner und neighbor flüssige zu feste Phase wechselt.
    Wenn ja:
        fvm.upper[face] -= deltaCoeffs[face].value()*gammaMagSf[face].value();

        fvm.diag() += deltaCoeffs[face].value()*gammaMagSf[face].value()* (1-delta_quotient);
    */

    //  Berechnung hinter negSumDiag(): 
    /*for (register label face=0; face<l.size(); face++)
    {
        Diag[l[face]] -= Lower[face];
        Diag[u[face]] -= Upper[face];
    }*/
    
    forAll(vf.boundaryField(), patchi)
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