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

#include "barycentricParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::barycentricParticle::propertyList_ =
    Foam::barycentricParticle::propertyList();

const std::size_t Foam::barycentricParticle::sizeofPosition_
(
    offsetof(barycentricParticle, facei_)
  - offsetof(barycentricParticle, position_)
);

const std::size_t Foam::barycentricParticle::sizeofFields_
(
    sizeof(barycentricParticle) - offsetof(barycentricParticle, position_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::barycentricParticle::barycentricParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    mesh_(mesh),
    position_(),
    celli_(-1),
    facei_(-1),
    stepFraction_(0.0),
    tetFacei_(-1),
    tetPtI_(-1),
    origProc_(Pstream::myProcNo()),
    origId_(-1)
{
    if (is.format() == IOstream::ASCII)
    {
        is  >> position_ >> celli_;

        if (readFields)
        {
            is  >> facei_
                >> stepFraction_
                >> tetFacei_
                >> tetPtI_
                >> origProc_
                >> origId_;
        }
    }
    else
    {
        if (readFields)
        {
            is.read(reinterpret_cast<char*>(&position_), sizeofFields_);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&position_), sizeofPosition_);
        }
    }

    // Check state of Istream
    is.check("barycentricParticle::barycentricParticle(Istream&, bool)");
}


void Foam::barycentricParticle::writePosition(Ostream& os) const
{
    if (os.format() == IOstream::ASCII)
    {
        os  << position_ << token::SPACE << celli_;
    }
    else
    {
        os.write(reinterpret_cast<const char*>(&position_), sizeofPosition_);
    }

    // Check state of Ostream
    os.check("barycentricParticle::writePosition(Ostream& os, bool) const");
}


Foam::Ostream& Foam::operator<<(Ostream& os, const barycentricParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << p.position_
            << token::SPACE << p.celli_
            << token::SPACE << p.facei_
            << token::SPACE << p.stepFraction_
            << token::SPACE << p.tetFacei_
            << token::SPACE << p.tetPtI_
            << token::SPACE << p.origProc_
            << token::SPACE << p.origId_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&p.position_),
            barycentricParticle::sizeofFields_
        );
    }

    return os;
}


// ************************************************************************* //
