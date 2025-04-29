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
    Foam::fv::ibmGeometryData
 
Group
 
Description
 
SourceFiles
    ibmGeometryData.C
 
\*---------------------------------------------------------------------------*/
 
#ifndef ibmGeometryData_H
#define ibmGeometryData_H
 
#include "laplacianScheme.H"
#include "indicator.H"
#include "advectionSchemes.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                       Class ibmGeometryData Declaration
\*---------------------------------------------------------------------------*/
 
class ibmGeometryData
{

private:

protected:
    const dynamicFvMesh& mesh_;
    label timeIndex_;

    static HashTable<ibmGeometryData> geoData_;

    bool updated();
    void update();
    void correct();

    tmp<scalarField> computeDistance();
    tmp<volScalarField> computeAlpha(tmp<scalarField> distance);

    void computeDeltaQuotient
    (
        const volVectorField& interfaceNormal,
        const volVectorField& interfaceCentre,
        const boolList& isInterfaceCell
    );

    autoPtr<surfaceScalarField> deltaQuotientPtr_;
    autoPtr<surfaceScalarField> vofGradPtr_;
    autoPtr<volScalarField> vofFieldPtr_;

    scalar intersection
    (
        const vector& c_1, 
        const vector& c_2, 
        const vector& c_E, 
        const vector& normal
    );
 
public:
 
 
    // Constructors
        ibmGeometryData(const dynamicFvMesh& mesh);
 
    //- Destructor
        virtual ~ibmGeometryData() = default;
    
    surfaceScalarField& deltaQuotient() { update(); return deltaQuotientPtr_(); }
    volScalarField& vofField()          { update(); return vofFieldPtr_(); }
    surfaceScalarField& vofGrad()       { update(); return vofGradPtr_(); }
    
    static HashTable<ibmGeometryData>& geoDataTable(){ return geoData_; }
};
 
} // End namespace Foam
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 
#endif
 
// ************************************************************************* //