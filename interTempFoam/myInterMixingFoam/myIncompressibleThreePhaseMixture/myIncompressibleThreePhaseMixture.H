/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::myIncompressibleThreePhaseMixture

Description

SourceFiles
    myIncompressibleThreePhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef myIncompressibleThreePhaseMixture_H
#define myIncompressibleThreePhaseMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "IOdictionary.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class myIncompressibleThreePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class myIncompressibleThreePhaseMixture
:
    public IOdictionary,
    public transportModel
{
    // Private data

        word phase1Name_;
        word phase2Name_;
        word phase3Name_;

        volScalarField alpha1_;
        volScalarField alpha2_;
        volScalarField alpha3_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;
//ADDED------------------------------>
        const volScalarField& T_;
//<------------------------------ADDED

        volScalarField nu_;
//ADDED----------------------------->
//	volScalarField kappa_;
//<------------------------------ADDED
        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;
        autoPtr<viscosityModel> nuModel3_;

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;
        dimensionedScalar rho3_;

//ADDED--------------------------->

        dimensionedScalar cp1_;
        dimensionedScalar cp2_;
        dimensionedScalar cp3_;

        dimensionedScalar Pr1_;
        dimensionedScalar Pr2_;
        dimensionedScalar Pr3_;
//<-------------------------ADDED


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();


public:

    // Constructors

        //- Construct from components
        myIncompressibleThreePhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
//ADDED--------------------------->
		const volScalarField& T
//<-------------------------ADDED
        );


    //- Destructor
    ~myIncompressibleThreePhaseMixture()
    {}


    // Member Functions

        const word phase1Name() const
        {
            return phase1Name_;
        }

        const word phase2Name() const
        {
            return phase2Name_;
        }

        const word phase3Name() const
        {
            return phase3Name_;
        }

        const volScalarField& alpha1() const
        {
            return alpha1_;
        }

        volScalarField& alpha1()
        {
            return alpha1_;
        }

        const volScalarField& alpha2() const
        {
            return alpha2_;
        }

        volScalarField& alpha2()
        {
            return alpha2_;
        }

        const volScalarField& alpha3() const
        {
            return alpha3_;
        }

        volScalarField& alpha3()
        {
            return alpha3_;
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //- Return const-access to phase3 density
        const dimensionedScalar& rho3() const
        {
            return rho3_;
        };

//ADDED--------------------------------------->
        const dimensionedScalar& cp1() const
        {
            return cp1_;
        }

        const dimensionedScalar& cp2() const
        {
            return cp2_;
        };

        const dimensionedScalar& cp3() const
        {
            return cp3_;
        };

        const dimensionedScalar& Pr1() const
        {
            return Pr1_;
        }

        const dimensionedScalar& Pr2() const
        {
            return Pr2_;
        };

        const dimensionedScalar& Pr3() const
        {
            return Pr3_;
        };
//<------------------------------------------ADDED

        //- Return the velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the flux
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

//ADDED-------------------------------------->
	const volScalarField& T() const
	{
		return T_;
	}
//<------------------------------------ADDED

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase3 viscosityModel
        const viscosityModel& nuModel3() const
        {
            return nuModel3_();
        }

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

//ADDED------------------------------------->
/*	tmp<volScalarField> kappa() const
	{
		return kappa_;
	}
	tmp<scalarField> kappa(const label patchi) const
	{
		return kappa_.boundaryField()[patchi];
	}
*/
	tmp<surfaceScalarField> kappaf() const;
//<--------------------------------------ADDED

        //- Correct the laminar viscosity
        void correct()
        {
            calcNu();
        }

        //- Read base transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
