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

const std::size_t Foam::barycentricParticle::sizeofFields_
(
    sizeof(barycentricParticle) - offsetof(barycentricParticle, barycentric_)
);

struct PositionAndCellI
{
    Foam::vector position_;
    Foam::label celli_;
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::barycentricParticle::barycentricParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    mesh_(mesh),
    barycentric_(),
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
        PositionAndCellI temp;
        is  >> temp.position_ >> temp.celli_;

        celli_ = temp.celli_;
        position(temp.position_);

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
            is.read(reinterpret_cast<char*>(&barycentric_), sizeofFields_);
        }
        else
        {
            PositionAndCellI temp;
            is.read(reinterpret_cast<char*>(&temp), sizeof(PositionAndCellI));

            celli_ = temp.celli_;
            position(temp.position_);
        }
    }

    // Check state of Istream
    is.check("barycentricParticle::barycentricParticle(Istream&, bool)");
}


void Foam::barycentricParticle::writePosition(Ostream& os) const
{
    if (os.format() == IOstream::ASCII)
    {
        os  << position() << token::SPACE << celli_;
    }
    else
    {
        PositionAndCellI temp = {position(), celli_};
        os.write(reinterpret_cast<const char*>(&temp), sizeof(PositionAndCellI));
    }

    // Check state of Ostream
    os.check("barycentricParticle::writePosition(Ostream& os, bool) const");
}


Foam::Ostream& Foam::operator<<(Ostream& os, const barycentricParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << p.barycentric_
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
            reinterpret_cast<const char*>(&p.barycentric_),
            barycentricParticle::sizeofFields_
        );
    }

    return os;
}


// ************************************************************************* //
