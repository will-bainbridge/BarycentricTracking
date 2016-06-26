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

Class
    Foam::TrackParticle

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef TrackParticle_H
#define TrackParticle_H

#include "particle.H"
#include "autoPtr.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template <class ParticleType>
class TrackParticle;

template <class ParticleType>
Ostream& operator<<(Ostream&, const TrackParticle<ParticleType>&);


/*---------------------------------------------------------------------------*\
                     Class TrackParticle Declaration
\*---------------------------------------------------------------------------*/

template <class ParticleType>
class TrackParticle
:
    public ParticleType
{
    // Private data

        //- Velocity
        vector U_;


public:

    //- Class used to pass tracking data to the trackToFace function
    template <class CloudType>
    class TrackingData
    :
        public ParticleType::template TrackingData<CloudType>
    {
    public:


        // Constructors

            TrackingData(CloudType& cloud)
            :
                ParticleType::template TrackingData<CloudType>(cloud)
            {}

    };


    // Static data members

        //- Runtime type information
        TypeName("TrackParticle");


    // Constructors

        //- Construct from components
        TrackParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label tetFacei,
            const label tetPtI,
            const vector& U
        );

        //- Construct from Istream
        TrackParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new TrackParticle<ParticleType>(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<TrackParticle<ParticleType> > operator()(Istream& is) const
            {
                return autoPtr<TrackParticle<ParticleType> >
                (
                    new TrackParticle<ParticleType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        //- Velocity
        vector& U()
        {
            return U_;
        }

        //- Velocity
        const vector& U() const
        {
            return U_;
        }


        // Tracking

            //- Track all particles for the given amount of time
            template <class TrackData>
            bool move(TrackData& td, const scalar trackTime);

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            template <class TrackData>
            bool hitPatch
            (
                const polyPatch&,
                TrackData& td,
                const label patchi,
                const scalar trackFraction,
                const tetIndices& tetIs
            );

            //- Overridable function to handle the particle hitting a wedge
            template <class TrackData>
            void hitWedgePatch
            (
                const wedgePolyPatch&,
                TrackData& td
            );

            //- Overridable function to handle the particle hitting a
            //  symmetry plane
            template <class TrackData>
            void hitSymmetryPlanePatch
            (
                const symmetryPlanePolyPatch&,
                TrackData& td
            );

            //- Overridable function to handle the particle hitting a
            //  symmetry patch
            template <class TrackData>
            void hitSymmetryPatch
            (
                const symmetryPolyPatch&,
                TrackData& td
            );

            //- Overridable function to handle the particle hitting a cyclic
            template <class TrackData>
            void hitCyclicPatch
            (
                const cyclicPolyPatch&,
                TrackData& td
            );

            //- Overridable function to handle the particle hitting a
            //- processorPatch
            template <class TrackData>
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                TrackData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            template <class TrackData>
            void hitWallPatch
            (
                const wallPolyPatch&,
                TrackData& td,
                const tetIndices&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            template <class TrackData>
            void hitPatch
            (
                const polyPatch&,
                TrackData& td
            );

            //- Convert processor patch addressing to the global equivalents
            //  and set the celli to the face-neighbour
            template <class TrackData>
            void correctAfterParallelTransfer(const label, TrackData&);


        // I-O

            //- Read
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const TrackParticle<ParticleType> &
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "TrackParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TrackParticle.C"
    #include "TrackParticleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
