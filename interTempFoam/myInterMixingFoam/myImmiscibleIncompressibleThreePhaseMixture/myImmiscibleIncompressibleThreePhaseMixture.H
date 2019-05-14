/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::myImmiscibleIncompressibleThreePhaseMixture

Description
    An immiscible incompressible two-phase mixture transport model

SourceFiles
    myImmiscibleIncompressibleThreePhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef myImmiscibleIncompressibleThreePhaseMixture_H
#define myImmiscibleIncompressibleThreePhaseMixture_H

#include "myIncompressibleThreePhaseMixture.H"
#include "myThreePhaseInterfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class myImmiscibleIncompressibleThreePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class myImmiscibleIncompressibleThreePhaseMixture
:
    public myIncompressibleThreePhaseMixture,
    public myThreePhaseInterfaceProperties
{

public:

    // Constructors

        //- Construct from components
        myImmiscibleIncompressibleThreePhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
//ADDED------------------------>
		const volScalarField& T
//<--------------------------ADDED
        );


    //- Destructor
    virtual ~myImmiscibleIncompressibleThreePhaseMixture()
    {}


    // Member Functions

        //- Correct the transport and interface properties
        virtual void correct()
        {
            myIncompressibleThreePhaseMixture::correct();
            myThreePhaseInterfaceProperties::correct();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
