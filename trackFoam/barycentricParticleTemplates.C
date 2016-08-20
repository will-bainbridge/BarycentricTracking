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
    NotImplemented;

    return 0.0;
}


template<class TrackData>
Foam::scalar Foam::barycentricParticle::trackToTri
(
    const vector& endPosition,
    TrackData& td
)
{
    NotImplemented;

    return 0.0;
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
