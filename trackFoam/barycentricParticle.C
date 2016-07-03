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
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::barycentricParticle::particleCount_ = 0;

const Foam::scalar Foam::barycentricParticle::trackingCorrectionTol = 1e-5;

const Foam::scalar Foam::barycentricParticle::lambdaDistanceToleranceCoeff =
    1e3*SMALL;

const Foam::scalar Foam::barycentricParticle::minStepFractionTol = 1e5*SMALL;

namespace Foam
{
    defineTypeNameAndDebug(barycentricParticle, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::barycentricParticle::barycentricParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPtI
)
:
    mesh_(mesh),
    position_(position),
    celli_(celli),
    facei_(-1),
    stepFraction_(0.0),
    tetFacei_(tetFacei),
    tetPtI_(tetPtI),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID())
{}


Foam::barycentricParticle::barycentricParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    bool doCellFacePt
)
:
    mesh_(mesh),
    position_(position),
    celli_(celli),
    facei_(-1),
    stepFraction_(0.0),
    tetFacei_(-1),
    tetPtI_(-1),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID())
{
    if (doCellFacePt)
    {
        initCellFacePt();
    }
}


Foam::barycentricParticle::barycentricParticle(const barycentricParticle& p)
:
    mesh_(p.mesh_),
    position_(p.position_),
    celli_(p.celli_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    tetFacei_(p.tetFacei_),
    tetPtI_(p.tetPtI_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


Foam::barycentricParticle::barycentricParticle
(
    const barycentricParticle& p,
    const polyMesh& mesh
)
:
    mesh_(mesh),
    position_(p.position_),
    celli_(p.celli_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    tetFacei_(p.tetFacei_),
    tetPtI_(p.tetPtI_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::barycentricParticle::transformProperties(const tensor&)
{}


void Foam::barycentricParticle::transformProperties(const vector&)
{}


Foam::scalar Foam::barycentricParticle::wallImpactDistance(const vector&) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

bool Foam::operator==
(
    const barycentricParticle& pA,
    const barycentricParticle& pB
)
{
    return (pA.origProc() == pB.origProc() && pA.origId() == pB.origId());
}


bool Foam::operator!=
(
    const barycentricParticle& pA,
    const barycentricParticle& pB
)
{
    return !(pA == pB);
}


// ************************************************************************* //
