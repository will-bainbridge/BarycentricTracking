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

Class
    Foam::barycentricParticle

Description
    Base barycentric particle class

\*---------------------------------------------------------------------------*/

#ifndef barycentricParticle_H
#define barycentricParticle_H

#include "vector.H"
#include "Barycentric.H"
#include "Cloud.H"
#include "IDLList.H"
#include "pointField.H"
#include "faceList.H"
#include "OFstream.H"
#include "tetrahedron.H"
#include "FixedList.H"
#include "polyMeshTetDecomposition.H"
#include "barycentricParticleMacros.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class barycentricParticle;

class polyPatch;

class cyclicPolyPatch;
class processorPolyPatch;
class symmetryPlanePolyPatch;
class symmetryPolyPatch;
class wallPolyPatch;
class wedgePolyPatch;

// Forward declaration of friend functions and operators

Ostream& operator<<
(
    Ostream&,
    const barycentricParticle&
);

bool operator==(const barycentricParticle&, const barycentricParticle&);

bool operator!=(const barycentricParticle&, const barycentricParticle&);

/*---------------------------------------------------------------------------*\
                          Class Particle Declaration
\*---------------------------------------------------------------------------*/

class barycentricParticle
:
    public IDLList<barycentricParticle>::link
{
    // Private member data

        //- Size in bytes of the position data
        static const std::size_t sizeofPosition_;

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    template<class CloudType>
    class TrackingData
    {
        // Private data

            //- Reference to the cloud containing (this) particle
            CloudType& cloud_;


    public:

        // Public data

            typedef CloudType cloudType;

            //- Flag to switch processor
            bool switchProcessor;

            //- Flag to indicate whether to keep particle (false = delete)
            bool keepParticle;

            //- The index of the triangle of the tet which has been hit. -1 if
            //  the track has stopped inside the tet.
            label tetTriI;


        // Constructor
        TrackingData(CloudType& cloud)
        :
            cloud_(cloud)
        {}


        // Member functions

            //- Return a reference to the cloud
            CloudType& cloud()
            {
                return cloud_;
            }
    };


protected:

    // Protected data

        //- Reference to the polyMesh database
        const polyMesh& mesh_;

        //- Barycentric coordinates of the particle within the current tet
        Barycentric<scalar> barycentric_;

        //- Index of the cell it is in
        label celli_;

        //- Face index if the particle is on a face otherwise -1
        label facei_;

        //- Fraction of time-step completed
        scalar stepFraction_;

        //- Index of the face that owns the decomposed tet that the
        //  particle is in
        label tetFacei_;

        //- Index of the point on the face that defines the decomposed
        //  tet that the particle is in.  Relative to the face base
        //  point.
        label tetPtI_;

        //- Originating processor id
        label origProc_;

        //- Local particle id on originating processor
        label origId_;


public: // private: // <-- !!! Make private

    // Private Member Functions

        //- Position as a noun. Calculate and return the cartesian position from
        //  the current barycentric coordinates.
        vector position() const;

        //- Position as a verb. Set the current barycentric coordinates from the
        //  given cartesian position.
        void position(const vector& x);


        // Tetrahedra functions

            //- Get indices into the current face for the face-bound vertices of
            //  the current tet.
            void tetFaceIndices
            (
                label& baseI,
                label& vertex1I,
                label& vertex2I
            ) const;

            //- Get indices into the mesh points for the face-bound vertices of
            //  the current tet.
            void tetMeshIndices
            (
                label& baseI,
                label& vertex1I,
                label& vertex2I
            ) const;

            //- Get the vertices of the current tet
            void tetGeometry
            (
                vector& centre,
                vector& base,
                vector& vertex1,
                vector& vertex2
            ) const;

            //- Get the transformation associated with the current tet. This
            //  will convert a barycentric position within the tet to a
            //  cartesian position in the global coordinate system. The
            //  conversion is x = centre + A&y, where x is the cartesian
            //  position, y is the barycentric position, centre is the origin,
            //  and A is the transformation tensor.
            void tetTransform
            (
                vector& centre,
                tensor& A
            ) const;

            //- Get the reverse transform associated with the current tet. The
            //  conversion is detA*y = T&(x - centre). The variables x, y and
            //  centre have the same meaning as for the forward transform. T is
            //  the inverse of the forward transform tensor, A, multiplied by
            //  its determinant, detA. This separation allows the barycentric
            //  tracking algorithm to function on inverted or degenerate
            //  tetrahedra.
            void tetReverseTransform
            (
                vector& centre,
                scalar& detA,
                tensor& T
            ) const;

            //- Get the vertices of the current moving tet. Two values are
            //  returned for each vertex. The first is a constant, and the
            //  second is a linear coefficient of the track fraction.
            void movingTetGeometry
            (
                Pair<vector>& centre,
                Pair<vector>& base,
                Pair<vector>& vertex1,
                Pair<vector>& vertex2
            ) const;

            //- Get the transformation associated with the current, moving, tet.
            //  This is of the same form as for the static case. As with the
            //  moving geometry, a linear function of the tracking fraction is
            //  returned for each component.
            void movingTetTransform
            (
                Pair<vector>& centre,
                Pair<tensor>& A
            ) const;

            //- Get the reverse transformation associated with the current,
            //  moving, tet. This is of the same form as for the static case. As
            //  with the moving geometry, a function of the tracking fraction is
            //  returned for each component. The functions are higher order than
            //  for the forward transform; the determinant is cubic, and the
            //  tensor is quadratic.
            void movingTetReverseTransform
            (
                Pair<vector>& centre,
                FixedList<scalar, 4>& detA,
                FixedList<tensor, 3>& T
            ) const;


        // Transformations

            //- Reflection transform. Corrects the coordinates when the particle
            //  moves between two tets which share a base vertex, but for which
            //  the other two non cell-centre vertices are reversed. All hits
            //  which retain the same face behave this way, as do face hits.
            void reflect();

            //- Rotation transform. Corrects the coordinates when the particle
            //  moves between two tets with different base vertices, but are
            //  otherwise similarly oriented. Hits which change the face within
            //  the cell make use of both this and the reflect transform.
            void rotate(const bool direction);


        // Track

            //- Change tet within a cell. Called after a triangle is hit.
            template<class TrackData>
            void changeTet(const TrackData& td);

            //- Change tet face within a cell. Called by changeTet.
            template<class TrackData>
            void changeFace(const TrackData& td);

            //- Change cell. Called when the particle hits an internal face.
            template<class TrackData>
            void changeCell(const TrackData& td);


protected:

    // Patch interactions

        //- Overridable function to handle the particle hitting a face
        template<class TrackData>
        void hitFace(TrackData& td);

        //- Overridable function to handle the particle hitting a
        //  patch.  Executed before other patch-hitting functions.
        //  trackFraction is passed in to allow mesh motion to
        //  interpolate in time to the correct face state.
        template<class TrackData>
        bool hitPatch
        (
            const polyPatch&,
            TrackData& td,
            const label patchi,
            const scalar trackFraction,
            const tetIndices& tetIs
        );

        //- Overridable function to handle the particle hitting a wedgePatch
        template<class TrackData>
        void hitWedgePatch(const wedgePolyPatch&, TrackData& td);

        //- Overridable function to handle the particle hitting a
        //  symmetryPlanePatch
        template<class TrackData>
        void hitSymmetryPlanePatch
        (
            const symmetryPlanePolyPatch&,
            TrackData& td
        );

        //- Overridable function to handle the particle hitting a
        //  symmetryPatch
        template<class TrackData>
        void hitSymmetryPatch(const symmetryPolyPatch&, TrackData& td);

        //- Overridable function to handle the particle hitting a cyclicPatch
        template<class TrackData>
        void hitCyclicPatch(const cyclicPolyPatch&, TrackData& td);

        //- Overridable function to handle the particle hitting a cyclicAMIPatch
        template<class TrackData>
        void hitCyclicAMIPatch
        (
            const cyclicAMIPolyPatch&,
            TrackData& td,
            const vector& direction
        );

        //- Overridable function to handle the particle hitting a
        //  processorPatch
        template<class TrackData>
        void hitProcessorPatch(const processorPolyPatch&, TrackData& td);

        //- Overridable function to handle the particle hitting a wallPatch
        template<class TrackData>
        void hitWallPatch
        (
            const wallPolyPatch&,
            TrackData& td,
            const tetIndices& tetIs
        );

        //- Overridable function to handle the particle hitting a
        //  general patch
        template<class TrackData>
        void hitPatch(const polyPatch&, TrackData& td);


public:

    // Static data members

        //- Runtime type information
        TypeName("barycentricParticle");

        //- String representation of properties
        DefinePropertyList
        (
            "(Px Py Pz) celli facei stepFraction "
            "tetFacei tetPtI origProc origId"
        );

        //- Cumulative particle counter - used to provode unique ID
        static label particleCount_;


    // Constructors

        //- Construct from components
        barycentricParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label tetFacei,
            const label tetPtI
        );

        //- Construct from components, tetFacei_ and tetPtI_ are not
        //  supplied so they will be deduced by a search
        barycentricParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            bool doCellFacePt = true
        );

        //- Construct from Istream
        barycentricParticle
        (
            const polyMesh& mesh,
            Istream&,
            bool readFields = true
        );

        //- Construct as a copy
        barycentricParticle(const barycentricParticle& p);

        //- Construct as a copy with refernce to a new mesh
        barycentricParticle
        (
            const barycentricParticle& p,
            const polyMesh& mesh
        );

        //- Construct a clone
        virtual autoPtr<barycentricParticle> clone() const
        {
            return autoPtr<barycentricParticle>
            (
                new barycentricParticle(*this)
            );
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

            autoPtr<barycentricParticle> operator()(Istream& is) const
            {
                return autoPtr<barycentricParticle>
                (
                    new barycentricParticle(mesh_, is, true)
                );
            }
        };


    //- Destructor
    virtual ~barycentricParticle()
    {}


    // Member Functions

        // Access

            //- Get unique particle creation id
            inline label getNewParticleID() const;

            //- Return the mesh database
            inline const polyMesh& mesh() const;

            //- Return current particle barycentric coordinates
            inline const Barycentric<scalar>& barycentric() const;

            //- Return current particle barycentric coordinates
            inline Barycentric<scalar>& barycentric();

            //- Return current cell particle is in
            inline label& cell();

            //- Return current cell particle is in
            inline label cell() const;

            //- Return current tet face particle is in
            inline label& tetFace();

            //- Return current tet face particle is in
            inline label tetFace() const;

            //- Return current tet face particle is in
            inline label& tetPt();

            //- Return current tet face particle is in
            inline label tetPt() const;

            //- Return the indices of the current tet that the
            //  particle occupies.
            inline tetIndices currentTetIndices() const;

            //- Return the geometry of the current tet that the
            //  particle occupies.
            inline tetPointRef currentTet() const;

            //- Return the normal of the tri on tetFacei_ for the
            //  current tet.
            inline vector normal() const;

            //- Return the normal of the tri on tetFacei_ for the
            //  current tet at the start of the timestep, i.e. based
            //  on oldPoints
            inline vector oldNormal() const;

            //- Return current face particle is on otherwise -1
            inline label& face();

            //- Return current face particle is on otherwise -1
            inline label face() const;

            //- Return the impact model to be used, soft or hard (default).
            inline bool softImpact() const;

            //- Return the particle current time
            inline scalar currentTime() const;


        // Check

            //- Check the stored cell value (setting if necessary) and
            //  initialise the tetFace and tetPt values
            inline void initCellFacePt();

            //- Have the cell, face, tet points been initialised?
            inline bool hasCellFacePt() const;

            //- Is the particle on the boundary/(or outside the domain)?
            inline bool onBoundary() const;

            //- Is this global face an internal face?
            inline bool internalFace(const label facei) const;

            //- Is this global face a boundary face?
            inline bool boundaryFace(const label facei) const;

            //- Which patch is particle on
            inline label patch(const label facei) const;

            //- Which face of this patch is this particle on
            inline label patchFace
            (
                const label patchi,
                const label facei
            ) const;

            //- Return the fraction of time-step completed
            inline scalar& stepFraction();

            //-  Return the fraction of time-step completed
            inline scalar stepFraction() const;

            //- Return const access to the originating processor id
            inline label origProc() const;

            //- Return the originating processor id for manipulation
            inline label& origProc();

            //- Return const access to  the particle id on originating processor
            inline label origId() const;

            //- Return the particle id on originating processor for manipulation
            inline label& origId();


        // Track

            //- Track particle to end of trajectory
            //  or until it hits the boundary.
            //  On entry 'stepFraction()' should be set to the fraction of the
            //  time-step at which the tracking starts and on exit it contains
            //  the fraction of the time-step completed.
            //  Returns the boundary face index if the track stops at the
            //  boundary, -1 otherwise.
            template<class TrackData>
            label track(const vector& endPosition, TrackData& td);

            //- Track particle to a given position and returns 1.0 if the
            //  trajectory is completed without hitting a face otherwise
            //  stops at the face and returns the fraction of the trajectory
            //  completed.
            //  on entry 'stepFraction()' should be set to the fraction of the
            //  time-step at which the tracking starts.
            template<class TrackData>
            scalar trackToFace(const vector& endPosition, TrackData& td);

            //- ...
            template<class TrackData>
            scalar trackToTri(const vector& endPosition, TrackData& td);

            //- Return the index of the face to be used in the interpolation
            //  routine
            inline label faceInterpolation() const;


    // Transformations

        //- Transform the physical properties of the particle
        //  according to the given transformation tensor
        virtual void transformProperties(const tensor& T);

        //- Transform the physical properties of the particle
        //  according to the given separation vector
        virtual void transformProperties(const vector& separation);

        //- The nearest distance to a wall that
        //  the particle can be in the n direction
        virtual scalar wallImpactDistance(const vector& n) const;


    // Parallel transfer

        //- Convert global addressing to the processor patch
        //  local equivalents
        template<class TrackData>
        void prepareForParallelTransfer(const label patchi, TrackData& td);

        //- Convert processor patch addressing to the global equivalents
        //  and set the celli to the face-neighbour
        template<class TrackData>
        void correctAfterParallelTransfer(const label patchi, TrackData& td);


    // I-O

        //- Read the fields associated with the owner cloud
        template<class CloudType>
        static void readFields(CloudType& c);

        //- Write the fields associated with the owner cloud
        template<class CloudType>
        static void writeFields(const CloudType& c);

        //- Write the particle position and cell
        void writePosition(Ostream&) const;


    // Friend Operators

        friend Ostream& operator<<(Ostream&, const barycentricParticle&);

        friend bool operator==
        (
            const barycentricParticle& pA,
            const barycentricParticle& pB
        );

        friend bool operator!=
        (
            const barycentricParticle& pA,
            const barycentricParticle& pB
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "barycentricParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "barycentricParticleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
