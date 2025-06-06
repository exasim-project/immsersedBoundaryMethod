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
 
#include "ibmGeometryData.H"

#include "cutCellIso.H"
#include "cutFaceIso.H"
#include "searchableSurfaces.H"
 

Foam::HashTable<Foam::ibmGeometryData> Foam::ibmGeometryData::geoData_;

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
bool Foam::ibmGeometryData::updated()
{
    if(mesh_.time().timeIndex() == timeIndex_) return true;

    timeIndex_ = mesh_.time().timeIndex();
    return false;
}

void Foam::ibmGeometryData::update()
{
    if(!updated() && mesh_.changing()) 
    {
        deltaQuotientPtr_.clear();
        vofGradPtr_.clear();
        vofFieldPtr_.clear();
        correct();
    }
}

void Foam::ibmGeometryData::correct()
{
    volScalarField alpha(computeAlpha(computeDistance()));

    if(!vofFieldPtr_.valid())
    {
        vofFieldPtr_.reset
        (
            new volScalarField( "vofField", pos(alpha - 0.5 - 1e-6) ) // Achtung: Extrem wichtig, dass hier ein "-" und kein "+" vor der 1e-6 steht!!!
        );
    }

    if(!deltaQuotientPtr_.valid())
    {
        deltaQuotientPtr_.reset
        (
            new surfaceScalarField 
            (
                IOobject
                (
                    "delta_quotient",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("delta_quotient", dimless, 1)
            )
        );
    }

    if(!vofGradPtr_.valid())
    {
        vofGradPtr_.reset(new surfaceScalarField(fvc::snGrad( vofFieldPtr_() )));
    }

    volVectorField& U = mesh_.lookupObjectRef<volVectorField>("U");
    surfaceScalarField& phi = mesh_.lookupObjectRef<surfaceScalarField>("phi");

    advector_ = advectionSchemes::New(alpha,phi,U);
    advector_->surf().reconstruct();
    
    computeDeltaQuotient
    (
        advector_->surf().normal(),
        advector_->surf().centre(),
        advector_->surf().interfaceCell()
    );
}

Foam::tmp<Foam::scalarField> Foam::ibmGeometryData::computeDistance()
{
    IOdictionary dict
    (
        IOobject
        (
            "stlDict",
            mesh_.time().system(),
            mesh_.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<searchableSurfaces> geomPtr_(nullptr);
    geomPtr_.reset
    (
        new searchableSurfaces
        (
            IOobject
            (
                "abc",                            
                mesh_.time().constant(),         
                "triSurface",                 
                mesh_.time(),                   
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict.subDict("geometry"),
            true             
        )
    );

    const pointField& pc = mesh_.points();
    tmp<scalarField> distancePtr( new scalarField(mesh_.nPoints(), Zero)); 
    scalarField& distance = distancePtr.ref();

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

    return distancePtr;
}


Foam::tmp<Foam::volScalarField> Foam::ibmGeometryData::computeAlpha(tmp<scalarField> distance)
{
    tmp<volScalarField> alphaPtr;
    alphaPtr.reset
    ( 
        new volScalarField 
        (
            IOobject
            (
                "alpha.ibm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("alpha.ibm", dimless, 0)
        )
    );
    volScalarField& alpha = alphaPtr.ref();

    cutCellIso cutCell(mesh_, distance.ref());
    cutFaceIso cutFace(mesh_, distance.ref());

    forAll(alpha, cellI)
    {
        cutCell.calcSubCell(cellI, 0.0);
        alpha[cellI] = max(min(cutCell.VolumeOfFluid(), 1), 0); 
    }

    // Setting boundary alpha values
    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].size() > 0)
        {
            const label start = mesh_.boundary()[patchi].patch().start();
            scalarField& alphap = alpha.boundaryFieldRef()[patchi];
            const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];

            forAll(alphap, patchFacei)
            {
                const label facei = patchFacei + start;
                cutFace.calcSubFace(facei, 0.0);
                alphap[patchFacei] =
                    mag(cutFace.subFaceArea())/magSfp[patchFacei];
            }
        }
    }
    return alphaPtr;
}


void Foam::ibmGeometryData::computeDeltaQuotient
(
    const volVectorField& interfaceNormal,
    const volVectorField& interfaceCentre,
    const boolList& isInterfaceCell
)
{
    surfaceScalarField& delta_quotient_ = deltaQuotientPtr_();
    surfaceScalarField& vofGrad_ = vofGradPtr_();
    volScalarField& vofField_ = vofFieldPtr_();

    //################################################
    volScalarField isInterfacecell
    (
        IOobject
        (
            "isInterfacecell",
            vofField_.instance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("isInterfacecell", dimless, -1)
    );
    forAll(mesh_.C(), celli){
        if(isInterfaceCell[celli]){
            isInterfacecell[celli] = 1;
        }else{
            isInterfacecell[celli] = 0;
        }
    }
    isInterfacecell.write();

    volScalarField delQuoKlEins
    (
        IOobject
        (
            "delQuoKlEins",
            vofField_.instance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("delQuoKlEins", dimless, 0)
    );
    //################################################
    
    forAll(delta_quotient_, i){

        if( vofGrad_[i] != 0 )
        {
            const label& owneri = mesh_.faceOwner()[i];
            const vector& cOwneri = mesh_.C()[owneri];

            const label& neighbouri = mesh_.faceNeighbour()[i];
            const vector& cNeighbouri = mesh_.C()[neighbouri];

            scalar dist_ownerInterface = GREAT;
            scalar dist_neighbourInterface = GREAT;

            if( vofField_[owneri]   ==  1 )
            {  // Dann muss Owner korrigiert werden.

                if(isInterfaceCell[owneri])
                {
                    dist_ownerInterface = intersection(cOwneri, cNeighbouri, interfaceCentre[owneri], interfaceNormal[owneri]);
                }
                if(isInterfaceCell[neighbouri]) //Es ist korrekt, dass hier keine else if() ist, da der Fall, dass beide Zellen ein Interface haben, abgedeckt sein muss.
                {
                    dist_neighbourInterface = intersection(cOwneri, cNeighbouri, interfaceCentre[neighbouri], interfaceNormal[neighbouri]);
                }

                if(
                    !isInterfaceCell[owneri] && !isInterfaceCell[neighbouri]    //Das ist der Fall, wenn das Interface exakt auf einer Gitterfläche liegt.
                ){
                    delta_quotient_[i] = 1/intersection(cOwneri, cNeighbouri, mesh_.Cf()[i], mesh_.Sf()[i]/mesh_.magSf()[i]);
                }
                else if(
                    dist_ownerInterface <= dist_neighbourInterface &&
                    dist_ownerInterface >= 0
                )
                {
                    delta_quotient_[i] = 1/(dist_ownerInterface + SMALL); //Evtl. SMALL in einen Wert ändern, der im Bereich des num. Fehlers liegt, um negative dist_* zu umgehen.
                }
                else
                {
                    delta_quotient_[i] = 1/(dist_neighbourInterface + SMALL);
                }
            }
            else{  // Dann muss Neighbour korrigiert werden.

                if(isInterfaceCell[owneri])
                {
                    dist_ownerInterface = intersection(cNeighbouri, cOwneri, interfaceCentre[owneri], interfaceNormal[owneri]);
                }
                if(isInterfaceCell[neighbouri])
                {
                    dist_neighbourInterface = intersection(cNeighbouri, cOwneri, interfaceCentre[neighbouri], interfaceNormal[neighbouri]);
                }
                
                if(
                    !isInterfaceCell[owneri] && !isInterfaceCell[neighbouri]    //Das ist der Fall, wenn das Interface exakt auf einem Gitter face liegt.
                ){
                    delta_quotient_[i] = 1/intersection(cNeighbouri, cOwneri, mesh_.Cf()[i], mesh_.Sf()[i]/mesh_.magSf()[i]);
                }
                else if(
                    dist_ownerInterface <= dist_neighbourInterface &&
                    dist_ownerInterface >= 0
                )
                {
                    delta_quotient_[i] = 1/(dist_ownerInterface + SMALL);
                }
                else
                {
                    delta_quotient_[i] = 1/(dist_neighbourInterface + SMALL);
                }
            }

            /*if(delta_quotient_[i] < 1){
                if(vofField_[owneri]   ==  1){
                    delQuoKlEins[owneri] = delta_quotient_[i];
                }else{
                    delQuoKlEins[neighbouri] = delta_quotient_[i];
                }
            }*/

        }
    }
    //delQuoKlEins.write();
}

Foam::scalar Foam::ibmGeometryData::intersection
(
    const vector& c_1, 
    const vector& c_2, 
    const vector& c_E, 
    const vector& normal
)
{
    const scalar& dividend = normal & (c_E - c_1);
    const scalar& divisor = normal & (c_2 - c_1);
    
//clip funktion verwenden void Foam::isoAdvection::applyBruteForceBounding()

    if(dividend/(divisor + SMALL) >= 1){
        return (1 - SMALL);
    }else if(dividend/(divisor + SMALL) < 1e-3){
        return 1e-3;
    }else{
        return dividend/(divisor + SMALL);
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Foam::ibmGeometryData::ibmGeometryData(const dynamicFvMesh& mesh)
:
mesh_(mesh),
timeIndex_(mesh.time().timeIndex())
{
    correct();
}

 
// ************************************************************************* //