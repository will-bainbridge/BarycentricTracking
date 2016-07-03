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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::TrackParticle<ParticleType>::TrackParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPtI,
    const vector& U
)
:
    ParticleType(mesh, position, celli, tetFacei, tetPtI),
    U_(U)
{}


template<class ParticleType>
Foam::TrackParticle<ParticleType>::TrackParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParticleType(mesh, is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> U_;
        }
        else
        {
            is.read(reinterpret_cast<char&>(U_));
        }
    }

    // Check state of Istream
    is.check
    (
        "TrackParticle<ParticleType>::TrackParticle"
        "(const polyMesh& mesh, Istream& is, bool readFields)"
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParticleType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const TrackParticle<ParticleType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParticleType&>(p) << token::SPACE << p.U_;
    }
    else
    {
        os  << static_cast<const ParticleType&>(p);
        os.write(reinterpret_cast<const char*>(&p.U_));
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const TrackParticle<ParticleType>&)"
    );

    return os;
}


// ************************************************************************* //
