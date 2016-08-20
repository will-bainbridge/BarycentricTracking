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

#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::barycentricParticle::getNewParticleID() const
{
    label id = particleCount_++;

    if (id == labelMax)
    {
        WarningInFunction
            << "Particle counter has overflowed. This might cause problems"
            << " when reconstructing particle tracks." << endl;
    }

    return id;
}


inline const Foam::polyMesh& Foam::barycentricParticle::mesh() const
{
    return mesh_;
}

inline const Foam::Barycentric<Foam::scalar>&
Foam::barycentricParticle::barycentric() const
{
    return barycentric_;
}


inline Foam::Barycentric<Foam::scalar>&
Foam::barycentricParticle::barycentric()
{
    return barycentric_;
}


inline Foam::label Foam::barycentricParticle::cell() const
{
    return celli_;
}


inline Foam::label& Foam::barycentricParticle::cell()
{
    return celli_;
}


inline Foam::label Foam::barycentricParticle::tetFace() const
{
    return tetFacei_;
}


inline Foam::label& Foam::barycentricParticle::tetFace()
{
    return tetFacei_;
}


inline Foam::label Foam::barycentricParticle::tetPt() const
{
    return tetPtI_;
}


inline Foam::label& Foam::barycentricParticle::tetPt()
{
    return tetPtI_;
}


inline Foam::tetIndices Foam::barycentricParticle::currentTetIndices() const
{
    return tetIndices(celli_, tetFacei_, tetPtI_, mesh_);
}


inline Foam::tetPointRef Foam::barycentricParticle::currentTet() const
{
    return currentTetIndices().tet(mesh_);
}


inline Foam::vector Foam::barycentricParticle::normal() const
{
    return currentTetIndices().faceTri(mesh_).normal();
}


inline Foam::vector Foam::barycentricParticle::oldNormal() const
{
    return currentTetIndices().oldFaceTri(mesh_).normal();
}


inline Foam::label Foam::barycentricParticle::face() const
{
    return facei_;
}


inline Foam::label& Foam::barycentricParticle::face()
{
    return facei_;
}


inline void Foam::barycentricParticle::initCellFacePt()
{
    if (celli_ != -1)
    {
        mesh_.findTetFacePt(celli_, position(), tetFacei_, tetPtI_);
    }

    if (celli_ == -1 || tetFacei_ == -1 || tetPtI_ == -1)
    {
        mesh_.findCellFacePt(position(), celli_, tetFacei_, tetPtI_);

        if (celli_ == -1)
        {
            FatalErrorInFunction
                << "cell, tetFace and tetPt search failure at position "
                << position() << abort(FatalError);
        }
    }
}


inline bool Foam::barycentricParticle::hasCellFacePt() const
{
    // This is a bit of a mess. The particle class is very fundamentally built
    // around the concept of staged construction. After construction, the
    // tetrahedra may be uninitialised. Therefore we have to consider the stored
    // coordinates as cartesian up until the point at which the tet is fully
    // specified. This function provides that logic.
    return celli_ != -1 && tetFacei_ != -1 && tetPtI_ != -1;
}


inline bool Foam::barycentricParticle::onBoundary() const
{
    return facei_ != -1 && facei_ >= mesh_.nInternalFaces();
}


inline Foam::scalar& Foam::barycentricParticle::stepFraction()
{
    return stepFraction_;
}


inline Foam::scalar Foam::barycentricParticle::stepFraction() const
{
    return stepFraction_;
}


inline Foam::label Foam::barycentricParticle::origProc() const
{
    return origProc_;
}


inline Foam::label& Foam::barycentricParticle::origProc()
{
    return origProc_;
}


inline Foam::label Foam::barycentricParticle::origId() const
{
    return origId_;
}


inline Foam::label& Foam::barycentricParticle::origId()
{
    return origId_;
}


inline bool Foam::barycentricParticle::softImpact() const
{
    return false;
}


inline Foam::scalar Foam::barycentricParticle::currentTime() const
{
    return
        mesh_.time().value()
      + stepFraction_*mesh_.time().deltaTValue();
}


inline bool Foam::barycentricParticle::internalFace(const label facei) const
{
    return mesh_.isInternalFace(facei);
}


bool Foam::barycentricParticle::boundaryFace(const label facei) const
{
    return !internalFace(facei);
}


inline Foam::label Foam::barycentricParticle::patch(const label facei) const
{
    return mesh_.boundaryMesh().whichPatch(facei);
}


inline Foam::label Foam::barycentricParticle::patchFace
(
    const label patchi,
    const label facei
) const
{
    return mesh_.boundaryMesh()[patchi].whichFace(facei);
}


inline Foam::label Foam::barycentricParticle::faceInterpolation() const
{
    return facei_;
}


// ************************************************************************* //
