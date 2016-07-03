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

Application
    trackFoam

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Cloud.H"
#include "dynamicFvMesh.H"
#include "TrackParticle.H"

#include "Barycentric.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef TrackParticle<particle> trackParticle; // <-- drop replacement for particle in here

defineTemplateTypeNameAndDebug(trackParticle, 0);
defineTemplateTypeNameAndDebug(Cloud<trackParticle>, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    #ifdef DyM
        #include "createNamedDynamicFvMesh.H"
    #else
        #include "createMesh.H"
    #endif

    Cloud<trackParticle> cloud
    (
        mesh,
        word("cloud"),
        false
    );

    trackParticle::readFields(cloud);

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        #ifdef DyM
            mesh.update();
            mesh.checkMesh(true);
        #endif

        trackParticle::TrackingData<Cloud<trackParticle>> trackingData(cloud);

        cloud.move(trackingData, runTime.deltaTValue());

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
