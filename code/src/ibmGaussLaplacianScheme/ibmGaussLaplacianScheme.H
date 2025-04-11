/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::fv::ibmGaussLaplacianScheme
 
Group
    grpFvLaplacianSchemes
 
Description
    Basic second-order laplacian using face-gradients and Gauss' theorem.
 
SourceFiles
    ibmGaussLaplacianScheme.C
 
\*---------------------------------------------------------------------------*/
 
#ifndef ibmGaussLaplacianScheme_H
#define ibmGaussLaplacianScheme_H
 
#include "laplacianScheme.H"
#include "indicator.H"
#include "advectionSchemes.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace fv
{
 
/*---------------------------------------------------------------------------*\
                       Class ibmGaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/
 
template<class Type, class GType>
class ibmGaussLaplacianScheme
:
    public fv::laplacianScheme<Type, GType>
{

private:

    //-
        volScalarField& alpha1_;

        volVectorField& U_;

        surfaceScalarField& phi_;

        autoPtr<advectionSchemes> advector_;

        const volVectorField& interfaceNormal_;

        const volVectorField& interfaceCentre_;

        const boolList& isInterfaceCell_;

        volScalarField vofField_;

        surfaceScalarField vofGrad_;

        surfaceScalarField delta_quotient_;

protected:

    // Private Member Functions
 
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
 
        //- No copy construct
        ibmGaussLaplacianScheme(const ibmGaussLaplacianScheme&) = delete;
 
        //- No copy assignment
        void operator=(const ibmGaussLaplacianScheme&) = delete;
 
 
public:
 
    //- Runtime type information
    TypeName("IBMGauss");
 
 
    // Constructors
 
        //- Construct null
        //ibmGaussLaplacianScheme(const fvMesh& mesh);
        ibmGaussLaplacianScheme(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh),
            alpha1_(mesh.lookupObjectRef<volScalarField>("alpha.water")),
            U_(mesh.lookupObjectRef<volVectorField>("U")),
            phi_(mesh.lookupObjectRef<surfaceScalarField>("phi")),
            advector_(advectionSchemes::New(alpha1_,phi_,U_)),
            interfaceNormal_(advector_->surf().normal()),
            interfaceCentre_(advector_->surf().centre()),
            isInterfaceCell_(advector_->surf().interfaceCell()),
            vofField_
            (
                IOobject
                (
                    "vof",
                    alpha1_.instance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("vof", dimless, -1.0)
            ),
            vofGrad_(fvc::snGrad(vofField_)),
            delta_quotient_
            (
                IOobject
                (
                    "delta_quotient",
                    phi_.instance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("delta_quotient", dimless, 1)
            )
        {
            //advector_.reset(advectionSchemes::New(alpha1_,phi_,U_));
            advector_->surf().reconstruct();

            forAll(mesh.C(), celli){        // binäres Feld
                if(alpha1_[celli] < 0.5){
                    vofField_[celli] = 0;
                }else{
                    vofField_[celli] = 1;
                }
            }
            
            vofGrad_ = fvc::snGrad(vofField_);
        }
 
        
        //- Construct from Istream
        //ibmGaussLaplacianScheme(const fvMesh& mesh, Istream& is);
        ibmGaussLaplacianScheme(const fvMesh& mesh, Istream& is)
        :
            laplacianScheme<Type, GType>(mesh, is),
            alpha1_(mesh.lookupObjectRef<volScalarField>("alpha.water")),
            U_(mesh.lookupObjectRef<volVectorField>("U")),
            phi_(mesh.lookupObjectRef<surfaceScalarField>("phi")),
            advector_(advectionSchemes::New(alpha1_,phi_,U_)),
            interfaceNormal_(advector_->surf().normal()),
            interfaceCentre_(advector_->surf().centre()),
            isInterfaceCell_(advector_->surf().interfaceCell()),
            vofField_
            (
                IOobject
                (
                    "vof",
                    alpha1_.instance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("vof", dimless, -1.0)
            ),
            vofGrad_(fvc::snGrad(vofField_)),
            delta_quotient_
            (
                IOobject
                (
                    "delta_quotient",
                    phi_.instance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("delta_quotient", dimless, 1)
            )
        {
            //advector_.reset(advectionSchemes::New(alpha1_,phi_,U_));
            advector_->surf().reconstruct();

            forAll(mesh.C(), celli){        // binäres Feld
                if(alpha1_[celli] < 0.5){
                    vofField_[celli] = 0;
                }else{
                    vofField_[celli] = 1;
                }
            }

            vofGrad_ = fvc::snGrad(vofField_);
        }
 
        //- Construct from mesh, interpolation and snGradScheme schemes
        ibmGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs),
            alpha1_(mesh.lookupObjectRef<volScalarField>("alpha.water")),
            U_(mesh.lookupObjectRef<volVectorField>("U")),
            phi_(mesh.lookupObjectRef<surfaceScalarField>("phi")),
            advector_(advectionSchemes::New(alpha1_,phi_,U_)),
            interfaceNormal_(advector_->surf().normal()),
            interfaceCentre_(advector_->surf().centre()),
            isInterfaceCell_(advector_->surf().interfaceCell()),
            vofField_
            (
                IOobject
                (
                    "vof",
                    alpha1_.instance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("vof", dimless, -1.0)
            ),
            vofGrad_(fvc::snGrad(vofField_)),
            delta_quotient_
            (
                IOobject
                (
                    "delta_quotient",
                    phi_.instance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("delta_quotient", dimless, 1)
            )
        {
            //advector_.reset(advectionSchemes::New(alpha1_,phi_,U_));
            advector_->surf().reconstruct();

            forAll(mesh.C(), celli){        // binäres Feld
                if(alpha1_[celli] < 0.5){
                    vofField_[celli] = 0;
                }else{
                    vofField_[celli] = 1;
                }
            }

            vofGrad_ = fvc::snGrad(vofField_);
        }
 
 
    //- Destructor
    virtual ~ibmGaussLaplacianScheme() = default;
 
 
    // Member Functions
        
        vector bisektion
        (
            const vector a,
            const vector b,
            const int max_iter,
            const vector centre,
            const scalar radius
        );

        scalar intersection
        (
            const vector& c_1, 
            const vector& c_2, 
            const vector& c_E, 
            const vector& normal
        );

        tmp<fvMatrix<Type>> fvmLaplacianUncorrected
        (
            const surfaceScalarField& gammaMagSf,
            const surfaceScalarField& phi_,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
 
        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );
 
        tmp<fvMatrix<Type>> fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
 
        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
};
 
 
// Use macros to emulate partial-specialisation of the Laplacian functions
// for scalar diffusivity gamma
 
#define defineFvmLaplacianScalarGamma(Type)                                    \
                                                                               \
template<>                                                                     \
tmp<fvMatrix<Type>> ibmGaussLaplacianScheme<Type, scalar>::fvmLaplacian           \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,                 \
    const GeometricField<Type, fvPatchField, volMesh>&                         \
);                                                                             \
                                                                               \
template<>                                                                     \
tmp<GeometricField<Type, fvPatchField, volMesh>>                               \
ibmGaussLaplacianScheme<Type, scalar>::fvcLaplacian                               \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,                 \
    const GeometricField<Type, fvPatchField, volMesh>&                         \
);
 
 
defineFvmLaplacianScalarGamma(scalar);
defineFvmLaplacianScalarGamma(vector);
defineFvmLaplacianScalarGamma(sphericalTensor);
defineFvmLaplacianScalarGamma(symmTensor);
defineFvmLaplacianScalarGamma(tensor);
 
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
} // End namespace fv
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
} // End namespace Foam
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#ifdef NoRepository
    #include "ibmGaussLaplacianScheme.C"
#endif
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#endif
 
// ************************************************************************* //