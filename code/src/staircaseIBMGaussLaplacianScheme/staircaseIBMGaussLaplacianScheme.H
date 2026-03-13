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
    Foam::fv::staircaseIBMGaussLaplacianScheme
 
Group
    grpFvLaplacianSchemes
 
Description
    Basic second-order laplacian using face-gradients and Gauss' theorem.
 
SourceFiles
    staircaseIBMGaussLaplacianScheme.C
 
\*---------------------------------------------------------------------------*/
 
#ifndef staircaseIBMGaussLaplacianScheme_H
#define staircaseIBMGaussLaplacianScheme_H
 
#include "laplacianScheme.H"
#include "advectionSchemes.H"

#include "ibmGeometryData.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace fv
{
 
/*---------------------------------------------------------------------------*\
                       Class staircaseIBMGaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/
 
template<class Type, class GType>
class staircaseIBMGaussLaplacianScheme
:
    public fv::laplacianScheme<Type, GType>
{

private:

protected:

    // Private Member Functions
 
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> gammaSnGradCorr
        (
            const surfaceVectorField& SfGammaCorr,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
 
        //- No copy construct
        staircaseIBMGaussLaplacianScheme(const staircaseIBMGaussLaplacianScheme&) = delete;
 
        //- No copy assignment
        void operator=(const staircaseIBMGaussLaplacianScheme&) = delete;

        word name_;
        word* geometryName_;

        static HashTable<ibmGeometryData> geoData_;

        void createGeoData(const fvMesh& mesh);
        
        ibmGeometryData& geoData();
 
 
public:
 
    //- Runtime type information
    TypeName("staircaseIBMGauss");
 
 
    // Constructors
 
        //- Construct null
        //staircaseIBMGaussLaplacianScheme(const fvMesh& mesh);
        staircaseIBMGaussLaplacianScheme(const fvMesh& mesh)
        :
            laplacianScheme<Type, GType>(mesh)
        {}
 
        
        //- Construct from Istream
        //staircaseIBMGaussLaplacianScheme(const fvMesh& mesh, Istream& is);
        staircaseIBMGaussLaplacianScheme(const fvMesh& mesh, Istream& is)
        :
            laplacianScheme<Type, GType>(mesh, is)
        {}
 
        //- Construct from mesh, interpolation and snGradScheme schemes
        staircaseIBMGaussLaplacianScheme
        (
            const fvMesh& mesh,
            const tmp<surfaceInterpolationScheme<GType>>& igs,
            const tmp<snGradScheme<Type>>& sngs
        )
        :
            laplacianScheme<Type, GType>(mesh, igs, sngs)
        {}
 
 
    //- Destructor
    virtual ~staircaseIBMGaussLaplacianScheme() = default;
 
 
    // Member Functions
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

        void setGeometry(const word& name){ geometryName_ = new word(name); }
};
 
 
// Use macros to emulate partial-specialisation of the Laplacian functions
// for scalar diffusivity gamma
 
#define defineFvmLaplacianScalarGamma(Type)                                    \
                                                                               \
template<>                                                                     \
tmp<fvMatrix<Type>> staircaseIBMGaussLaplacianScheme<Type, scalar>::fvmLaplacian           \
(                                                                              \
    const GeometricField<scalar, fvsPatchField, surfaceMesh>&,                 \
    const GeometricField<Type, fvPatchField, volMesh>&                         \
);                                                                             \
                                                                               \
template<>                                                                     \
tmp<GeometricField<Type, fvPatchField, volMesh>>                               \
staircaseIBMGaussLaplacianScheme<Type, scalar>::fvcLaplacian                               \
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
    #include "staircaseIBMGaussLaplacianScheme.C"
#endif
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#endif
 
// ************************************************************************* //