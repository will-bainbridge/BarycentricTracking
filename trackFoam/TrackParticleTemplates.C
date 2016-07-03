/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Will Bainbridge
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

#include "TrackParticle.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class ParticleType>
template<class CloudType>
void Foam::TrackParticle<ParticleType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParticleType::readFields(c);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        TrackParticle<ParticleType>& p = iter();

        p.U_ = U[i];

        i++;
    }
}


template<class ParticleType>
template<class CloudType>
void Foam::TrackParticle<ParticleType>::writeFields(const CloudType& c)
{
    ParticleType::writeFields(c);

    label np = c.size();

    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const TrackParticle<ParticleType>& p = iter();

        U[i] = p.U();

        i++;
    }

    U.write();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
template<class TrackData>
bool Foam::TrackParticle<ParticleType>::move
(
    TrackData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;

    scalar tEnd = (1.0 - this->stepFraction())*trackTime;
    scalar dtMax = tEnd;

    if (tEnd <= SMALL && this->onBoundary())
    {
        // This is a hack to handle particles reaching their endpoint
        // on a processor boundary. If the endpoint is on a processor face
        // it currently gets transferred backwards and forwards infinitely.

        // Remove the particle
        td.keepParticle = false;
    }
    else
    {
        td.keepParticle = true;

        while (td.keepParticle && !td.switchProcessor && tEnd > SMALL)
        {
            scalar dt = min(dtMax, tEnd);

            dt *= this->trackToFace(this->position() + dt*U_, td);

            tEnd -= dt;

            this->stepFraction() = 1.0 - tEnd/trackTime;
        }
    }

    return td.keepParticle;
}


template<class ParticleType>
template<class TrackData>
bool Foam::TrackParticle<ParticleType>::hitPatch
(
    const polyPatch&,
    TrackData& td,
    const label patchi,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    return false;
}


template<class ParticleType>
template<class TrackData>
void Foam::TrackParticle<ParticleType>::hitWedgePatch
(
    const wedgePolyPatch&,
    TrackData& td
)
{
    td.keepParticle = false;
}


template<class ParticleType>
template<class TrackData>
void Foam::TrackParticle<ParticleType>::hitSymmetryPlanePatch
(
    const symmetryPlanePolyPatch&,
    TrackData& td
)
{
    td.keepParticle = false;
}


template<class ParticleType>
template<class TrackData>
void Foam::TrackParticle<ParticleType>::hitSymmetryPatch
(
    const symmetryPolyPatch&,
    TrackData& td
)
{
    td.keepParticle = false;
}


template<class ParticleType>
template<class TrackData>
void Foam::TrackParticle<ParticleType>::hitCyclicPatch
(
    const cyclicPolyPatch&,
    TrackData& td
)
{
    td.keepParticle = false;
}


template<class ParticleType>
template<class TrackData>
void Foam::TrackParticle<ParticleType>::hitProcessorPatch
(
    const processorPolyPatch&,
    TrackData& td
)
{
    td.switchProcessor = true;
}


template<class ParticleType>
template<class TrackData>
void Foam::TrackParticle<ParticleType>::hitWallPatch
(
    const wallPolyPatch& wpp,
    TrackData& td,
    const tetIndices&
)
{
    td.keepParticle = false;
}


template<class ParticleType>
template<class TrackData>
void Foam::TrackParticle<ParticleType>::hitPatch
(
    const polyPatch& wpp,
    TrackData& td
)
{
    td.keepParticle = false;
}


template<class ParticleType>
template<class TrackData>
void Foam::TrackParticle<ParticleType>::correctAfterParallelTransfer
(
    const label patchi,
    TrackData& td
)
{
    ParticleType::correctAfterParallelTransfer(patchi, td);
}



// ************************************************************************* //
