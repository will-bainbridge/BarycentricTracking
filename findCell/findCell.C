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
    findCell

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("position");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const vector position = args.argRead<vector>(1);
    const label cellI = mesh.findCell(position);

    if (cellI != -1)
    {
        Info<< "Position " << position << " was found in cell #" << cellI
            << endl;
    }
    else
    {
        Info<< "A cell enclosing position " << position
            << " could not be found." << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
