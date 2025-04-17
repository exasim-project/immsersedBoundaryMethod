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
#include "cutCellIso.H"
#include "cutFaceIso.H"
#include "searchableSurfaces.H"
#include "foamTool.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace fv
{
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    // Bisektionsverfahren:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*template<class Type, class GType>
vector ibmGaussLaplacianScheme<Type, GType>::bisektion
(
    vector a,
    vector b,
    int max_iter,
    vector centre,
    scalar radius
)
{    
    vector c;
    
    for (int i = 0; i < max_iter; ++i) {
        c = (a + b) / 2;
        if(
            (mag(c-centre) < radius &&
            mag(a-centre) < radius) ||
            (mag(c-centre) > radius &&
            mag(a-centre) > radius)
        ){
            a = c;
        }else{
            b = c;
        }
    }

    return c;
}*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Line-Plane Intersection:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class Type, class GType>
scalar ibmGaussLaplacianScheme<Type, GType>::intersection
(
    const vector& c_1, 
    const vector& c_2, 
    const vector& c_E, 
    const vector& normal
)
{
    const scalar& dividend = normal & (c_E - c_1);
    const scalar& divisor = normal & (c_2 - c_1);
    
    return dividend/(divisor + SMALL);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Funktion aus setAlphaField:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
template<class Type, class GType>
void ibmGaussLaplacianScheme<Type, GType>::setAlpha
(
    volScalarField& alpha1,
    DynamicList<List<point>>& facePts,
    scalarField& f,
    const bool writeOBJ // immer false
)
{
    const fvMesh& mesh = alpha1.mesh();
    cutCellIso cutCell(mesh, f);
    cutFaceIso cutFace(mesh, f);

    forAll(alpha1, cellI)
    {
        cutCell.calcSubCell(cellI, 0.0);

        alpha1[cellI] = max(min(cutCell.VolumeOfFluid(), 1), 0); // hier könnte problem liegen

        if (writeOBJ && (mag(cutCell.faceArea()) >= 1e-14))
        {
            facePts.append(cutCell.facePoints());
        }
    }

    // Setting boundary alpha1 values
    forAll(mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].size() > 0)
        {
            const label start = mesh.boundary()[patchi].patch().start();
            scalarField& alphap = alpha1.boundaryFieldRef()[patchi];
            const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];

            forAll(alphap, patchFacei)
            {
                const label facei = patchFacei + start;
                cutFace.calcSubFace(facei, 0.0);
                alphap[patchFacei] =
                    mag(cutFace.subFaceArea())/magSfp[patchFacei];
            }
        }
    }
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
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.upper() = deltaCoeffs.primitiveField()*gammaMagSf.primitiveField();
    fvm.negSumDiag();

    const fvMesh& mesh = this->mesh();

    // cylinder channel:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    /*const vector& centre = vector(0, 0, 0);
    const scalar& radius = 0.5;
    forAll(vofField_, i){
        if(
            mag(mesh.C()[i]-centre) < radius
        ){
            vofField_[i] = 0.0;
        }
    }
    if(mesh.time().outputTime()){ vofField_.write(); }*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // sin channel in phase:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    /*forAll(vofField_, i){
        if(
            mesh.C()[i].component(1) < 0.6 + 0.5*sin(mesh.C()[i].component(0)*4/5*M_PI) ||
            mesh.C()[i].component(1) > 1.4 + 0.5*sin(mesh.C()[i].component(0)*4/5*M_PI)
        ){
            vofField_[i] = 0.0;
        }
    }*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // sin channel not in phase:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    /*forAll(vofField_, i){
        if(
            mesh.C()[i].component(1) < 0.35 + 0.25*cos(mesh.C()[i].component(0)*4/5*M_PI) ||
            mesh.C()[i].component(1) > 1.65 - 0.25*cos(mesh.C()[i].component(0)*4/5*M_PI)
        ){
            vofField_[i] = 0.0;
        }
    }*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // channel with angle:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    /*const scalar& alpha = 35.0; // Winkel in [deg]
    forAll(vofField_, i){
        if(
            mesh.C()[i].component(1) < tan(Foam::degToRad(alpha))*mesh.C()[i].component(0) + 0.1 ||
            mesh.C()[i].component(1) > tan(Foam::degToRad(alpha))*mesh.C()[i].component(0) + 0.9
        ){
            vofField_[i] = 0.0;
        }
    }*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // gridStudy:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    /*const scalar& k = 0;
    const scalar cellLenght = pow(mesh.V()[0]*10, 1/2); // delta_z = 0.1
    forAll(vofField_, i){
        if(
            mesh.C()[i].component(1) < 0.1+k*cellLenght ||
            mesh.C()[i].component(1) > 0.9+k*cellLenght // fünftel der feinsten Zelle, k:1...5
        ){
            vofField_[i] = 0.0;
        }
    }*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Erweiterung auf stl-Datei als input:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    const dynamicFvMesh& dynMesh = static_cast<const dynamicFvMesh&>(mesh);
    dictionary dict = foamTools::readDict
    (
        dynMesh,
        foamTools::dictFile(dynMesh, "distances/stlDict")
    );
    autoPtr<searchableSurfaces> geomPtr_;
    geomPtr_.reset(nullptr);
    geomPtr_.reset
    (
        new searchableSurfaces
        (
            IOobject
            (
                "abc",                             // dummy name
                mesh.time().constant(),           // directory
                "triSurface",                      // instance
                mesh.time(),                      // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict.subDict("geometry"),
            true                // allow single-region shortcut
        )
    );

    const pointField& pc = mesh.points();
    scalarField distance(mesh.nPoints(), Zero);    // distance ist f Field

    labelList surfaces;
    List<pointIndexHit> nearestInfo;
    geomPtr_().findNearest
    (
        pc,
        scalarField(pc.size(), GREAT),
        surfaces,
        nearestInfo
    );

    forAll(nearestInfo, i)
    {
        distance[i] = mag(nearestInfo[i].hitPoint()-pc[i]);
    }

    /*volScalarField distanceCell
    (
        IOobject
        (
            "distanceCell",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("distanceCell", dimLength, 0)
    );
    
    scalar n;
    forAll(mesh.C(), celli){
        n = 0;
        forAll(mesh.cellPoints()[celli], i){
            n += 1;
            const label& pointi = mesh.cellPoints()[celli][i];
            distanceCell[celli] += distance[pointi];
        }
        distanceCell[celli] = distanceCell[celli]/n;
    }
    distanceCell.write();*/

        // Hierdurch wird Fehler verursacht:
    List<volumeType> volType;
    forAll(geomPtr_(), sID)
    {
        geomPtr_()[sID].getVolumeType
        (
            pc,
            volType
        );
        forAll(volType, pointi)
        {
            if(volType[pointi] == volumeType::INSIDE)
            {
                distance[pointi] = -1*distance[pointi];
            }
        }
    }

    /*volScalarField distanceCell
    (
        IOobject
        (
            "distanceCell",
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("distanceCell", dimLength, 0)
    );
    
    scalar n;
    forAll(mesh.C(), celli){
        n = 0;
        forAll(mesh.cellPoints()[celli], i){
            n += 1;
            const label& pointi = mesh.cellPoints()[celli][i];
            distanceCell[celli] += distance[pointi];
        }
        distanceCell[celli] = distanceCell[celli]/n;
    }
    distanceCell.write();*/


    // Dinge die in Constructor verschoben wurden:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    autoPtr<advectionSchemes> advector_;
    //volScalarField& alpha1 = mesh.lookupObjectRef<volScalarField>("alpha.water"); // Zeile, wenn alpha.water verwendet wird
    volScalarField alpha1     // Zeile, wenn geomPtr verwendet wird
    (
        IOobject
        (
            "alpha.water",
            vf.instance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha.water", dimless, 0)
    );
    volVectorField& U = mesh.lookupObjectRef<volVectorField>("U");
    surfaceScalarField& phi = mesh.lookupObjectRef<surfaceScalarField>("phi");
    
    DynamicList<List<point>> facePts;
    setAlpha(alpha1, facePts, distance, false);

    advector_.reset(advectionSchemes::New(alpha1,phi,U));
    advector_->surf().reconstruct();

    //volScalarField& vofField_ = alpha1;

    volScalarField vofField_    // Später mit Indicator definieren
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
        dimensionedScalar("vof", dimless, -1.0)
    );
    
    forAll(mesh.C(), celli){
        if(alpha1[celli] < 0.5){
            vofField_[celli] = 0;
        }else{
            vofField_[celli] = 1;
        }
    }
    
    //const scalar& delta_quotient_ = 2; // Werte zwischen 1 und unendlich
    surfaceScalarField delta_quotient_    // in Constructor verlegt
    (
        IOobject
        (
            "delta_quotient",
            deltaCoeffs.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("delta_quotient", dimless, 1)
    );

    const volVectorField& interfaceNormal_ = advector_->surf().normal();
    const volVectorField& interfaceCentre_ = advector_->surf().centre();
    const boolList& isInterfaceCell_ = advector_->surf().interfaceCell();
    
    surfaceScalarField vofGrad_ = fvc::snGrad(vofField_);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Distance calculation:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    volScalarField interfaceDistanceX_
    (
        IOobject
        (
            "interfaceDistance_X",
            vofField_.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("interfaceDistance_X", dimLength, -1)
    );
    volScalarField interfaceDistanceY_
    (
        IOobject
        (
            "interfaceDistance_Y",
            vofField_.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("interfaceDistance_Y", dimLength, -1)
    );    

    const vector& x = vector(1, 0, 0);
    const vector& y = vector(0, 1, 0);
    
    scalar dist_ownerInterface;
    scalar dist_neighbourInterface;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    forAll(delta_quotient_, i){

        if(
            vofGrad_[i] != 0
        ){
            const label& owneri = mesh.faceOwner()[i];
            const vector& cOwneri = mesh.C()[owneri];
            const vector& interfaceNormalOwner = interfaceNormal_[owneri];
            const vector& interfaceCentreOwner = interfaceCentre_[owneri];

            const label& neighbouri = mesh.faceNeighbour()[i];
            const vector& cNeighbouri = mesh.C()[neighbouri];
            const vector& interfaceNormalNeighbour = interfaceNormal_[neighbouri];
            const vector& interfaceCentreNeighbour = interfaceCentre_[neighbouri];

            dist_ownerInterface = 1e10;
            dist_neighbourInterface = 1e10;

            const scalar& isYdir = (cOwneri-cNeighbouri) & x;
            const scalar& isXdir = (cOwneri-cNeighbouri) & y;

            if(
                vofField_[owneri]   ==  1
            ){  // Dann muss Owner korrigiert werden.

                if(isInterfaceCell_[owneri]){
                    dist_ownerInterface = intersection(cOwneri, cNeighbouri, interfaceCentreOwner, interfaceNormalOwner);
                }
                if(isInterfaceCell_[neighbouri]){
                    dist_neighbourInterface = intersection(cOwneri, cNeighbouri, interfaceCentreNeighbour, interfaceNormalNeighbour);
                }

                if(
                    dist_ownerInterface <= dist_neighbourInterface &&
                    dist_ownerInterface >= 0
                ){
                    delta_quotient_[i] = 1/(dist_ownerInterface + SMALL);

                    if(mag(isYdir) < 1e-7){
                        interfaceDistanceY_[owneri] = dist_ownerInterface / deltaCoeffs[i];
                    }
                    if(mag(isXdir) < 1e-7){
                        interfaceDistanceX_[owneri] = dist_ownerInterface / deltaCoeffs[i];
                    }
                }else{
                    delta_quotient_[i] = 1/(dist_neighbourInterface + SMALL);

                    if(mag(isYdir) < 1e-7){
                        interfaceDistanceY_[owneri] = dist_neighbourInterface / deltaCoeffs[i];
                    }
                    if(mag(isXdir) < 1e-7){
                        interfaceDistanceX_[owneri] = dist_neighbourInterface / deltaCoeffs[i];
                    }
                }
            }else{  // Dann muss Neighbour korrigiert werden.

                if(isInterfaceCell_[owneri]){
                    dist_ownerInterface = intersection(cNeighbouri, cOwneri, interfaceCentreOwner, interfaceNormalOwner);
                }
                if(isInterfaceCell_[neighbouri]){
                    dist_neighbourInterface = intersection(cNeighbouri, cOwneri, interfaceCentreNeighbour, interfaceNormalNeighbour);
                }
                
                if(
                    dist_ownerInterface <= dist_neighbourInterface &&
                    dist_ownerInterface >= 0
                ){
                    delta_quotient_[i] = 1/(dist_ownerInterface + SMALL);

                    if(mag(isYdir) < 1e-7){
                        interfaceDistanceY_[neighbouri] = dist_ownerInterface / deltaCoeffs[i];
                    }
                    if(mag(isXdir) < 1e-7){
                        interfaceDistanceX_[neighbouri] = dist_ownerInterface / deltaCoeffs[i];
                    }
                }else{
                    delta_quotient_[i] = 1/(dist_neighbourInterface + SMALL);
                    
                    if(mag(isYdir) < 1e-7){
                        interfaceDistanceY_[neighbouri] = dist_neighbourInterface / deltaCoeffs[i];
                    }
                    if(mag(isXdir) < 1e-7){
                        interfaceDistanceX_[neighbouri] = dist_neighbourInterface / deltaCoeffs[i];
                    }
                }
            }
        }
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    /*scalar del_x_1 = 0.0;   // x-Abstand zu Cellcentre des unteren Randes
    scalar del_x_2 = 0.0;   // x-Abstand zu Cellcentre des oberen Randes
    scalar del_y_1 = 0.0;   // y-Abstand zu Cellcentre des unteren Randes
    scalar del_y_2 = 0.0;   // y-Abstand zu Cellcentre des oberen Randes*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // cylinder channel:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    /*forAll(delta_quotient_, i){

        if(vofGrad_[i] != 0.0){                          // nur faces, die auch Abstand benötigen

            const label& faceOwneri = mesh.faceOwner()[i];
            const vector& cOwneri = mesh.C()[faceOwneri];
            const label& faceNeighbouri = mesh.faceNeighbour()[i];
            const vector& cNeighbouri = mesh.C()[faceNeighbouri];
            const scalar& delX = 1/deltaCoeffs[i];      // Abstand benachbarter Cellcentre. 1/deltaCoeffs, da deltaCoeffs als 1/(x_i-x_i-1) auf orthogonalen Gitter definiert ist.
            
            vector cylinderIntersection = bisektion(cOwneri, cNeighbouri, 10, centre, radius);

            if(vofField_[faceOwneri] == 1.0){           // muss Abstand von owner genommen werden?

                delta_quotient_[i] = delX/mag(cOwneri - cylinderIntersection);

            }else{                                      // Abstand von neighbour!
                
                delta_quotient_[i] = delX/mag(cNeighbouri - cylinderIntersection);
            }
        }
    }*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // channel with angle:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    /*forAll(delta_quotient_, i){

        if(vofGrad_[i] != 0.0){                          // nur faces, die auch Abstand benötigen

            const label& faceOwneri = mesh.faceOwner()[i];
            const vector& cOwneri = mesh.C()[faceOwneri];
            const label& faceNeighbouri = mesh.faceNeighbour()[i];
            const vector& cNeighbouri = mesh.C()[faceNeighbouri];
            const scalar& delX = 1/deltaCoeffs[i];      // Abstand benachbarter Cellcentre. 1/deltaCoeffs, da deltaCoeffs als 1/(x_i-x_i-1) auf orthogonalen Gitter definiert ist.

            if(vofField_[faceOwneri] == 1.0){           // muss Abstand von owner genommen werden?

                if(mesh.Sf()[i].component(1) == 0.0){   // Abstand in x-Richtung? (y-Richtung ändert sich nicht)

                    del_x_1 = mag(cOwneri.component(0) - (cOwneri.component(1) - 0.1)/tan(Foam::degToRad(alpha)));
                    del_x_2 = mag(cOwneri.component(0) - (cOwneri.component(1) - 0.9)/tan(Foam::degToRad(alpha)));

                    if(del_x_1*del_x_2 == 0){
                        delta_quotient_[i] = 1e24;
                    }else if(del_x_1 < del_x_2){              // Welche Wand ist die relevante?
                        delta_quotient_[i] = delX/del_x_1;
                    }else{
                        delta_quotient_[i] = delX/del_x_2;
                    }

                }else{                                  // Abstand in y-Richtung!

                    del_y_1 = mag(cOwneri.component(1) - (tan(Foam::degToRad(alpha))*cOwneri.component(0) + 0.1));
                    del_y_2 = mag(cOwneri.component(1) - (tan(Foam::degToRad(alpha))*cOwneri.component(0) + 0.9));

                    if(del_y_1*del_y_2 == 0){
                        delta_quotient_[i] = 1e24;
                    }else if(del_y_1 < del_y_2){              // Welche Wand ist die Relevante?
                        delta_quotient_[i] = delX/del_y_1;
                    }else{
                        delta_quotient_[i] = delX/del_y_2;
                    }
                }

            }else{                                      // Abstand von neighbour!
                if(mesh.Sf()[i].component(1) == 0.0){   // Abstand in x-Richtung?

                    del_x_1 = mag(cNeighbouri.component(0) - (cNeighbouri.component(1) - 0.1)/tan(Foam::degToRad(alpha)));
                    del_x_2 = mag(cNeighbouri.component(0) - (cNeighbouri.component(1) - 0.9)/tan(Foam::degToRad(alpha)));
                    
                    if(del_x_1*del_x_2 == 0){
                        delta_quotient_[i] = 1e24;
                    }else if(del_x_1 < del_x_2){              // Welche Wand ist die relevante?
                        delta_quotient_[i] = delX/del_x_1;
                    }else{
                        delta_quotient_[i] = delX/del_x_2;
                    }

                }else{                                  // Abstand in y-Richtung!

                    del_y_1 = mag(cNeighbouri.component(1) - (tan(Foam::degToRad(alpha))*cNeighbouri.component(0) + 0.1));
                    del_y_2 = mag(cNeighbouri.component(1) - (tan(Foam::degToRad(alpha))*cNeighbouri.component(0) + 0.9));

                    if(del_y_1*del_y_2 == 0){
                        delta_quotient_[i] = 1e24;
                    }else if(del_y_1 < del_y_2){              // Welche Wand ist die relevante?
                        delta_quotient_[i] = delX/del_y_1;
                    }else{
                        delta_quotient_[i] = delX/del_y_2;
                    }
                }
            }
        }
    }*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // sin Channel in phase:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    /*forAll(delta_quotient_, i){

        if(vofGrad_[i] != 0.0){                          // nur faces, die auch Abstand benötigen

            const label& faceOwneri = mesh.faceOwner()[i];
            const vector& cOwneri = mesh.C()[faceOwneri];
            const label& faceNeighbouri = mesh.faceNeighbour()[i];
            const vector& cNeighbouri = mesh.C()[faceNeighbouri];
            const scalar& delX = 1/deltaCoeffs[i];      // Abstand benachbarter Cellcentre. 1/deltaCoeffs, da deltaCoeffs als 1/(x_i-x_i-1) auf orthogonalen Gitter definiert ist.

            if(vofField_[faceOwneri] == 1.0){           // muss Abstand von owner genommen werden?

                if(mesh.Sf()[i].component(1) == 0.0){   // Abstand in x-Richtung? (y-Richtung ändert sich nicht)

                    if(
                        cOwneri.component(1) >= 0.1 && 
                        cOwneri.component(1) < 0.9
                    ){
                        del_x_1 = mag(cOwneri.component(0) - (5/(4*M_PI)*asin(2*(cOwneri.component(1)-0.6))));
                        if(del_x_1 == 0){
                            delta_quotient_[i] = 1e24;
                        }else{
                            delta_quotient_[i] = delX/del_x_1;
                        }
                    }else if(
                        cOwneri.component(1) >= 0.9 && 
                        cOwneri.component(1) <= 1.1
                    ){
                        del_x_1 = mag(cOwneri.component(0) - (5/(4*M_PI)*asin(2*(cOwneri.component(1)-0.6))));
                        del_x_2 = mag(cOwneri.component(0) - (5/(4*M_PI)*asin(2*(cOwneri.component(1)-1.4))));

                        if(del_x_1*del_x_2 == 0){
                            delta_quotient_[i] = 1e24;
                        }else if(del_x_1 < del_x_2){              // Welche Wand ist die relevante?
                            delta_quotient_[i] = delX/del_x_1;
                        }else{
                            delta_quotient_[i] = delX/del_x_2;
                        }
                    }else if(
                        cOwneri.component(1) > 1.1 && 
                        cOwneri.component(1) <= 1.9
                    ){
                        del_x_2 = mag(cOwneri.component(0) - (5/(4*M_PI)*asin(2*(cOwneri.component(1)-1.4))));
                        if(del_x_2 == 0){
                            delta_quotient_[i] = 1e24;
                        }else{
                            delta_quotient_[i] = delX/del_x_2;
                        }
                    }

                }else{                                  // Abstand in y-Richtung!

                    del_y_1 = mag(cOwneri.component(1) - (0.6 + 0.5*sin(cOwneri.component(0)*4/5*M_PI)));
                    del_y_2 = mag(cOwneri.component(1) - (1.4 + 0.5*sin(cOwneri.component(0)*4/5*M_PI)));

                    if(del_y_1*del_y_2 == 0){
                        delta_quotient_[i] = 1e24;
                    }else if(del_y_1 < del_y_2){              // Welche Wand ist die Relevante?
                        delta_quotient_[i] = delX/del_y_1;
                    }else{
                        delta_quotient_[i] = delX/del_y_2;
                    }
                }

            }else{                                      // Abstand von neighbour!
                if(mesh.Sf()[i].component(1) == 0.0){   // Abstand in x-Richtung?

                    if(
                        cNeighbouri.component(1) >= 0.1 && 
                        cNeighbouri.component(1) < 0.9
                    ){
                        del_x_1 = mag(cNeighbouri.component(0) - (5/(4*M_PI)*asin(2*(cNeighbouri.component(1)-0.6))));
                        if(del_x_1 == 0){
                            delta_quotient_[i] = 1e24;
                        }else{
                            delta_quotient_[i] = delX/del_x_1;
                        }
                    }else if(
                        cNeighbouri.component(1) >= 0.9 && 
                        cNeighbouri.component(1) <= 1.1
                    ){
                        del_x_1 = mag(cNeighbouri.component(0) - (5/(4*M_PI)*asin(2*(cNeighbouri.component(1)-0.6))));
                        del_x_2 = mag(cNeighbouri.component(0) - (5/(4*M_PI)*asin(2*(cNeighbouri.component(1)-1.4))));

                        if(del_x_1*del_x_2 == 0){
                            delta_quotient_[i] = 1e24;
                        }else if(del_x_1 < del_x_2){              // Welche Wand ist die relevante?
                            delta_quotient_[i] = delX/del_x_1;
                        }else{
                            delta_quotient_[i] = delX/del_x_2;
                        }
                    }else if(
                        cNeighbouri.component(1) > 1.1 && 
                        cNeighbouri.component(1) <= 1.9
                    ){
                        del_x_2 = mag(cNeighbouri.component(0) - (5/(4*M_PI)*asin(2*(cNeighbouri.component(1)-1.4))));
                        if(del_x_2 == 0){
                            delta_quotient_[i] = 1e24;
                        }else{
                            delta_quotient_[i] = delX/del_x_2;
                        }
                    }

                }else{                                  // Abstand in y-Richtung!

                    del_y_1 = mag(cNeighbouri.component(1) - (0.6 + 0.5*sin(cNeighbouri.component(0)*4/5*M_PI)));
                    del_y_2 = mag(cNeighbouri.component(1) - (1.4 + 0.5*sin(cNeighbouri.component(0)*4/5*M_PI)));

                    if(del_y_1*del_y_2 == 0){
                        delta_quotient_[i] = 1e24;
                    }else if(del_y_1 < del_y_2){              // Welche Wand ist die relevante?
                        delta_quotient_[i] = delX/del_y_1;
                    }else{
                        delta_quotient_[i] = delX/del_y_2;
                    }
                }
            }
        }
    }*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // sin Channel not in phase:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    /*forAll(delta_quotient_, i){

        if(vofGrad_[i] != 0.0){                          // nur faces, die auch Abstand benötigen

            const label& faceOwneri = mesh.faceOwner()[i];
            const vector& cOwneri = mesh.C()[faceOwneri];
            const label& faceNeighbouri = mesh.faceNeighbour()[i];
            const vector& cNeighbouri = mesh.C()[faceNeighbouri];
            const scalar& delX = 1/deltaCoeffs[i];      // Abstand benachbarter Cellcentre. 1/deltaCoeffs, da deltaCoeffs als 1/(x_i-x_i-1) auf orthogonalen Gitter definiert ist.

            if(vofField_[faceOwneri] == 1.0){           // muss Abstand von owner genommen werden?

                if(mesh.Sf()[i].component(1) == 0.0){   // Abstand in x-Richtung? (y-Richtung ändert sich nicht)

                    if(
                        cOwneri.component(1) >= 0.1 && 
                        cOwneri.component(1) <= 0.6
                    ){
                        del_x_1 = mag(cOwneri.component(0) - (5/(4*M_PI)*acos(4*(cOwneri.component(1))-0.35)));
                        if(del_x_1 == 0){
                            delta_quotient_[i] = 1e24;
                        }else{
                            delta_quotient_[i] = delX/del_x_1;
                        }
                    }else if(
                        cOwneri.component(1) >= 1.4 && 
                        cOwneri.component(1) <= 1.9
                    ){
                        del_x_2 = mag(cOwneri.component(0) - (5/(4*M_PI)*acos(4*(1.65-cOwneri.component(1)))));
                        if(del_x_2 == 0){
                            delta_quotient_[i] = 1e24;
                        }else{
                            delta_quotient_[i] = delX/del_x_2;
                        }
                    }

                }else{                                  // Abstand in y-Richtung!

                    del_y_1 = mag(cOwneri.component(1) - (0.35 + 0.25*cos(cOwneri.component(0)*4/5*M_PI)));
                    del_y_2 = mag(cOwneri.component(1) - (1.65 - 0.25*cos(cOwneri.component(0)*4/5*M_PI)));

                    if(del_y_1*del_y_2 == 0){
                        delta_quotient_[i] = 1e24;
                    }else if(del_y_1 < del_y_2){              // Welche Wand ist die Relevante?
                        delta_quotient_[i] = delX/del_y_1;
                    }else{
                        delta_quotient_[i] = delX/del_y_2;
                    }
                }

            }else{                                      // Abstand von neighbour!
                if(mesh.Sf()[i].component(1) == 0.0){   // Abstand in x-Richtung?

                    if(
                        cNeighbouri.component(1) >= 0.1 && 
                        cNeighbouri.component(1) <= 0.6
                    ){
                        del_x_1 = mag(cNeighbouri.component(0) - (5/(4*M_PI)*acos(4*(cNeighbouri.component(1))-0.35)));
                        if(del_x_1 == 0){
                            delta_quotient_[i] = 1e24;
                        }else{
                            delta_quotient_[i] = delX/del_x_1;
                        }
                    }else if(
                        cNeighbouri.component(1) >= 1.4 && 
                        cNeighbouri.component(1) <= 1.9
                    ){
                        del_x_2 = mag(cNeighbouri.component(0) - (5/(4*M_PI)*acos(4*(1.65-cNeighbouri.component(1)))));
                        if(del_x_2 == 0){
                            delta_quotient_[i] = 1e24;
                        }else{
                            delta_quotient_[i] = delX/del_x_2;
                        }
                    }

                }else{                                  // Abstand in y-Richtung!

                    del_y_1 = mag(cNeighbouri.component(1) - (0.35 + 0.25*cos(cNeighbouri.component(0)*4/5*M_PI)));
                    del_y_2 = mag(cNeighbouri.component(1) - (1.65 - 0.25*cos(cNeighbouri.component(0)*4/5*M_PI)));

                    if(del_y_1*del_y_2 == 0){
                        delta_quotient_[i] = 1e24;
                    }else if(del_y_1 < del_y_2){              // Welche Wand ist die relevante?
                        delta_quotient_[i] = delX/del_y_1;
                    }else{
                        delta_quotient_[i] = delX/del_y_2;
                    }
                }
            }
        }
    }*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // gridStudy:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    /*forAll(delta_quotient_, i){
        if(
            mesh.Cf()[i].component(1) < 0.5
        ){
            delta_quotient_[i] = 2;
        }else{
            delta_quotient_[i] = 2;
        }
    }*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //Korrekturterm:
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    surfaceScalarField surfCentreCorrection(gammaMagSf*(1.0-delta_quotient_)*mag(vofGrad_)); // benötigt binäres Feld
    surfaceScalarField faceCorrection(gammaMagSf*mag(vofGrad_)); // benötigt binäres Feld

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
                volCentreCorrection[celli] += surfCentreCorrection[facei];
            }
        }
    }

    if(mesh.time().outputTime())
    {
        vofField_.write();
        alpha1.write();
        interfaceDistanceX_.write();
        interfaceDistanceY_.write(); 
        volCentreCorrection.write();
    }

    fvm.diag() += volCentreCorrection.primitiveField() + (1-vofField_.primitiveField())*1e24;
    // Geschwindigkeit wird im Solid zu 0 gesetzt (benötigt binäres Feld)
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fvm.upper() -= faceCorrection.primitiveField();


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