/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "IOPosition.H"

#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "processorPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wallPolyPatch.H"
#include "wedgePolyPatch.H"
#include "meshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TrackData>
void Foam::barycentricParticle::prepareForParallelTransfer
(
    const label patchi,
    TrackData& td
)
{
    // Convert the face index to be local to the processor patch
    facei_ = patchFace(patchi, facei_);
}


template<class TrackData>
void Foam::barycentricParticle::correctAfterParallelTransfer
(
    const label patchi,
    TrackData& td
)
{
    const coupledPolyPatch& ppp =
        refCast<const coupledPolyPatch>(mesh_.boundaryMesh()[patchi]);

    celli_ = ppp.faceCells()[facei_];

    // Have patch transform the position
    //ppp.transformPosition(position_, facei_); // <-- !!!

    // Transform the properties
    if (!ppp.parallel())
    {
        const tensor& T =
        (
            ppp.forwardT().size() == 1
          ? ppp.forwardT()[0]
          : ppp.forwardT()[facei_]
        );
        transformProperties(T);
    }
    else if (ppp.separated())
    {
        const vector& s =
        (
            (ppp.separation().size() == 1)
          ? ppp.separation()[0]
          : ppp.separation()[facei_]
        );
        transformProperties(-s);
    }

    tetFacei_ = facei_ + ppp.start();

    // Faces either side of a coupled patch have matched base indices,
    // tetPtI is specified relative to the base point, already and
    // opposite circulation directions by design, so if the vertices
    // are:
    // source:
    // face    (a b c d e f)
    // fPtI     0 1 2 3 4 5
    //            +
    // destination:
    // face    (a f e d c b)
    // fPtI     0 1 2 3 4 5
    //                  +
    // where a is the base point of the face are matching , and we
    // have fPtI = 1 on the source processor face, i.e. vertex b, then
    // this because of the face circulation direction change, vertex c
    // is the characterising point on the destination processor face,
    // giving the destination fPtI as:
    //     fPtI_d = f.size() - 1 - fPtI_s = 6 - 1 - 1 = 4
    // This relationship can be verified for other points and sizes of
    // face.

    tetPtI_ = mesh_.faces()[tetFacei_].size() - 1 - tetPtI_;

    // Reset the face index for the next tracking operation
    if (stepFraction_ > (1.0 - SMALL))
    {
        stepFraction_ = 1.0;
        facei_ = -1;
    }
    else
    {
        facei_ += ppp.start();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::barycentricParticle::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    IOobject procIO(c.fieldIOobject("origProcId", IOobject::MUST_READ));

    if (procIO.headerOk())
    {
        IOField<label> origProcId(procIO);
        c.checkFieldIOobject(c, origProcId);
        IOField<label> origId(c.fieldIOobject("origId", IOobject::MUST_READ));
        c.checkFieldIOobject(c, origId);

        label i = 0;
        forAllIter(typename CloudType, c, iter)
        {
            barycentricParticle& p = iter();

            p.origProc_ = origProcId[i];
            p.origId_ = origId[i];
            i++;
        }
    }
}


template<class CloudType>
void Foam::barycentricParticle::writeFields(const CloudType& c)
{
    // Write the cloud position file
    IOPosition<CloudType> ioP(c);
    ioP.write();

    label np =  c.size();

    IOField<label> origProc
    (
        c.fieldIOobject("origProcId", IOobject::NO_READ),
        np
    );
    IOField<label> origId(c.fieldIOobject("origId", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename CloudType, c, iter)
    {
        origProc[i] = iter().origProc_;
        origId[i] = iter().origId_;
        i++;
    }

    origProc.write();
    origId.write();
}


template<class TrackData>
Foam::label Foam::barycentricParticle::track
(
    const vector& endPosition,
    TrackData& td
)
{
    facei_ = -1;

    // Tracks to endPosition or stop on boundary
    while (!onBoundary() && stepFraction_ < 1.0 - SMALL)
    {
        stepFraction_ += trackToFace(endPosition, td)*(1.0 - stepFraction_);
    }

    return facei_;
}


template<class TrackData>
Foam::scalar Foam::barycentricParticle::trackToFace
(
    const vector& endPosition,
    TrackData& td
)
{
    scalar trackFraction = 0.0;

    td.tetTriI = -1;

    facei_ = -1;

    // Tracking loop
    do
    {
        trackFraction += trackToTri(endPosition, td)*(1 - trackFraction);

        if (td.tetTriI == -1)
        {
            // The track is complete
            return trackFraction;
        }
        else if (td.tetTriI == 0)
        {
            // Set the face so that the loop will exit for face/patch processing
            facei_ = tetFacei_;
        }
        else
        {
            // Move to the next tet and continue the track
            changeTet(td);
        }

    } while(facei_ == -1);

    // Face/patch processing
    if (internalFace(facei_))
    {
        changeCell(td);
    }
    else
    {
        NotImplemented;
    }

    return trackFraction;
}


template<class TrackData>
Foam::scalar Foam::barycentricParticle::trackToTri
(
    const vector& endPosition,
    TrackData& td
)
{
    const vector x0 = position();
    // !!! <-- We really need the input as a displacement, not an end position
    const vector x1 = endPosition - x0;
    const Barycentric<scalar> y0 = barycentric_;

    if (debug)
    {
        Info<< "Tracking from " << x0 << " to " << x0 + x1 << endl;
    }

    // Get the tet geometry
    vector centre;
    scalar detA;
    tensor T;
    tetReverseTransform(centre, detA, T);

    if (debug)
    {
        vector o, b, v1, v2;
        tetGeometry(o, b, v1, v2);
        Info<< "Tet points o=" << o << ", b=" << b
            << ", v1=" << v1 << ", v2=" << v2 << endl
            << "Tet determinant = " << detA << endl
            << "Start local coordinates = " << y0 << endl;
    }

    // Calculate the local tracking displacement
    const Barycentric<scalar> Tx1 = toBarycentric(T, x1);

    if (debug)
    {
        Info<< "Local displacement = " << Tx1 << "/" << detA << endl;
    }

    // Calculate the hit fraction
    label iH = -1;
    scalar muH = detA == 0 ? VGREAT : mag(1/detA);
    for (label i = 0; i < 4; ++ i)
    {
        if (Tx1[i] < 0)
        {
            scalar mu = - y0[i]/Tx1[i];

            if (debug)
            {
                Info<< "Hit on tet face " << i << " at local coordinate "
                    << y0 + mu*Tx1 << ", " << mu*detA*100 << "\% of the "
                    << "way along the track" << endl;
            }

            if (0 <= mu && mu < muH)
            {
                iH = i;
                muH = mu;
            }
        }
    }

    // Set the new coordinates
    Barycentric<scalar> yH = y0 + muH*Tx1;

    // Remove tolerance issues in the event of a hit
    if (iH != -1)
    {
        yH.replace(iH, 0);
    }

    // Set the new position and hit index
    barycentric_ = yH;
    td.tetTriI = iH;

    if (debug)
    {
        if (iH != -1)
        {
            Info<< "Track hit tet face " << iH << " first" << endl;
        }
        else
        {
            Info<< "Track hit no tet faces" << endl;
        }
        Info<< "End local coordinates = " << yH << endl
            << "End global coordinates = " << position() << endl
            << "Tracking displacement = " << position() - x0 << endl
            << muH*detA*100 << "\% of the track completed" << endl;
    }

    // Return the proportion of the track that has been completed
    return muH*detA;
}


template<class TrackData>
void Foam::barycentricParticle::changeTet(const TrackData& td)
{
    const bool isOwner = mesh_.faceOwner()[tetFacei_] == celli_;

    const label firstTetPtI = 1;
    const label lastTetPtI = mesh_.faces()[tetFacei_].size() - 2;

    if (td.tetTriI == 1)
    {
        changeFace(td);
    }
    else if (td.tetTriI == 2)
    {
        if (isOwner)
        {
            if (tetPtI_ == lastTetPtI)
            {
                changeFace(td);
            }
            else
            {
                reflect();
                tetPtI_ += 1;
            }
        }
        else
        {
            if (tetPtI_ == firstTetPtI)
            {
                changeFace(td);
            }
            else
            {
                reflect();
                tetPtI_ -= 1;
            }
        }
    }
    else if (td.tetTriI == 3)
    {
        if (isOwner)
        {
            if (tetPtI_ == firstTetPtI)
            {
                changeFace(td);
            }
            else
            {
                reflect();
                tetPtI_ -= 1;
            }
        }
        else
        {
            if (tetPtI_ == lastTetPtI)
            {
                changeFace(td);
            }
            else
            {
                reflect();
                tetPtI_ += 1;
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Changing tet without changing cell should only happen when the "
            << "track is on triangle 1, 2 or 3."
            << exit(FatalError);
    }
}


template<class TrackData>
void Foam::barycentricParticle::changeFace(const TrackData& td)
{
    // Get the tet topology
    label basei, vertex1i, vertex2i;
    tetMeshIndices(basei, vertex1i, vertex2i);

    // Get the shared edge and the pre-rotation
    edge sharedEdge;
    if (td.tetTriI == 1)
    {
        sharedEdge = edge(vertex1i, vertex2i);
    }
    else if (td.tetTriI == 2)
    {
        sharedEdge = edge(vertex2i, basei);
    }
    else if (td.tetTriI == 3)
    {
        sharedEdge = edge(basei, vertex1i);
    }
    else
    {
        FatalErrorInFunction
            << "Changing face without changing cell should only happen when the "
            << "track is on triangle 1, 2 or 3."
            << exit(FatalError);
    }

    // Find the face in the same cell that shares the edge, and the
    // corresponding tetrahedra point
    tetPtI_ = -1;
    const label nFaces = mesh_.cells()[celli_].size();
    for (label cellFaceI = 0; tetPtI_ == -1 && cellFaceI < nFaces; ++ cellFaceI)
    {
        const label newFaceI = mesh_.cells()[celli_][cellFaceI];

        // Exclude the current face
        if (tetFacei_ == newFaceI)
        {
            continue;
        }

        const class face& newFace = mesh_.faces()[newFaceI];

        // Loop over the edges, looking for the shared one
        label edgeComp = 0;
        label edgeI = 0;
        for (; edgeComp == 0 && edgeI < newFace.size(); ++ edgeI)
        {
            edgeComp = edge::compare(sharedEdge, newFace.faceEdge(edgeI));
        }

        // If the face does not contain the edge, then move on to the next face
        if (edgeComp == 0)
        {
            continue;
        }

        // Correct the edge index based on whether the face is owned or
        // neighbours the current cell, and whether the comparison was in order
        // or not
        const bool isOwner = mesh_.faceOwner()[newFaceI] == celli_;
        if (isOwner && edgeComp == 1)
        {
            edgeI = newFace.prevLabel(edgeI);
        }
        else if (!isOwner && edgeComp == -1)
        {
            edgeI = newFace.nextLabel(edgeI);
        }

        // Set the new face and tet point. The loop will now exit as the tet
        // point has a value.
        const label newBaseI = max(0, mesh_.tetBasePtIs()[newFaceI]);
        tetFacei_ = newFaceI;
        tetPtI_ = (edgeI + newFace.size() - newBaseI) % newFace.size();
    }

    if (tetPtI_ == -1)
    {
        FatalErrorInFunction
            << "The search for an edge-connected face and tet-point failed."
            << exit(FatalError);
    }

    // Pre-rotation puts the shared edge opposite the base of the tetrahedron
    if (sharedEdge.otherVertex(vertex1i) == -1)
    {
        rotate(false);
    }
    else if (sharedEdge.otherVertex(vertex2i) == -1)
    {
        rotate(true);
    }

    // Update the new tet topology
    tetMeshIndices(basei, vertex1i, vertex2i);

    // Reflect to account for the change of triangle orientation on the new face
    reflect();

    // Post rotation puts the shared edge back in the correct location
    if (sharedEdge.otherVertex(vertex1i) == -1)
    {
        rotate(true);
    }
    else if (sharedEdge.otherVertex(vertex2i) == -1)
    {
        rotate(false);
    }
}


template<class TrackData>
void Foam::barycentricParticle::changeCell(const TrackData& td)
{
    // Set the cell to be the one on the other side of the face
    const label ownerCellI = mesh_.faceOwner()[tetFacei_];
    const bool isOwner = celli_ == ownerCellI;
    celli_ = isOwner ? mesh_.faceNeighbour()[tetFacei_] : ownerCellI;

    // Reflect to account for the change of triangle orientation in the new cell
    reflect();
}


template<class TrackData>
void Foam::barycentricParticle::hitFace(TrackData&)
{}


template<class TrackData>
bool Foam::barycentricParticle::hitPatch
(
    const polyPatch&,
    TrackData&,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


template<class TrackData>
void Foam::barycentricParticle::hitWedgePatch
(
    const wedgePolyPatch& wpp,
    TrackData&
)
{
    FatalErrorInFunction
        << "Hitting a wedge patch should not be possible."
        << abort(FatalError);

    vector nf = normal();
    nf /= mag(nf);

    transformProperties(I - 2.0*nf*nf);
}


template<class TrackData>
void Foam::barycentricParticle::hitSymmetryPlanePatch
(
    const symmetryPlanePolyPatch& spp,
    TrackData&
)
{
    vector nf = normal();
    nf /= mag(nf);

    transformProperties(I - 2.0*nf*nf);
}


template<class TrackData>
void Foam::barycentricParticle::hitSymmetryPatch
(
    const symmetryPolyPatch& spp,
    TrackData&
)
{
    vector nf = normal();
    nf /= mag(nf);

    transformProperties(I - 2.0*nf*nf);
}


template<class TrackData>
void Foam::barycentricParticle::hitCyclicPatch
(
    const cyclicPolyPatch& cpp,
    TrackData& td
)
{
    facei_ = cpp.transformGlobalFace(facei_);

    celli_ = mesh_.faceOwner()[facei_];

    tetFacei_ = facei_;

    // See note in correctAfterParallelTransfer for tetPtI_ addressing.
    tetPtI_ = mesh_.faces()[tetFacei_].size() - 1 - tetPtI_;

    const cyclicPolyPatch& receiveCpp = cpp.neighbPatch();
    label patchFacei = receiveCpp.whichFace(facei_);

    // Now the particle is on the receiving side

    // Have patch transform the position
    //receiveCpp.transformPosition(position_, patchFacei); // <-- !!!

    // Transform the properties
    if (!receiveCpp.parallel())
    {
        const tensor& T =
        (
            receiveCpp.forwardT().size() == 1
          ? receiveCpp.forwardT()[0]
          : receiveCpp.forwardT()[patchFacei]
        );
        transformProperties(T);
    }
    else if (receiveCpp.separated())
    {
        const vector& s =
        (
            (receiveCpp.separation().size() == 1)
          ? receiveCpp.separation()[0]
          : receiveCpp.separation()[patchFacei]
        );
        transformProperties(-s);
    }
}


template<class TrackData>
void Foam::barycentricParticle::hitCyclicAMIPatch
(
    const cyclicAMIPolyPatch& cpp,
    TrackData& td,
    const vector& direction
)
{
    const cyclicAMIPolyPatch& receiveCpp = cpp.neighbPatch();

    // Patch face index on sending side
    label patchFacei = facei_ - cpp.start();

    // Patch face index on receiving side - also updates position
    //patchFacei = cpp.pointFace(patchFacei, direction, position_); // <-- !!!

    if (patchFacei < 0)
    {
        FatalErrorInFunction
            << "Particle lost across " << cyclicAMIPolyPatch::typeName
            << " patches " << cpp.name() << " and " << receiveCpp.name()
            << " at position " << position() << abort(FatalError);
    }

    // Convert face index into global numbering
    facei_ = patchFacei + receiveCpp.start();

    celli_ = mesh_.faceOwner()[facei_];

    tetFacei_ = facei_;

    // See note in correctAfterParallelTransfer for tetPtI_ addressing.
    tetPtI_ = mesh_.faces()[tetFacei_].size() - 1 - tetPtI_;

    // Now the particle is on the receiving side

    // Transform the properties
    if (!receiveCpp.parallel())
    {
        const tensor& T =
        (
            receiveCpp.forwardT().size() == 1
          ? receiveCpp.forwardT()[0]
          : receiveCpp.forwardT()[patchFacei]
        );
        transformProperties(T);
    }
    else if (receiveCpp.separated())
    {
        const vector& s =
        (
            (receiveCpp.separation().size() == 1)
          ? receiveCpp.separation()[0]
          : receiveCpp.separation()[patchFacei]
        );
        transformProperties(-s);
    }
}


template<class TrackData>
void Foam::barycentricParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    TrackData&
)
{}


template<class TrackData>
void Foam::barycentricParticle::hitWallPatch
(
    const wallPolyPatch&,
    TrackData&,
    const tetIndices&
)
{}


template<class TrackData>
void Foam::barycentricParticle::hitPatch(const polyPatch&, TrackData&)
{}


// ************************************************************************* //
