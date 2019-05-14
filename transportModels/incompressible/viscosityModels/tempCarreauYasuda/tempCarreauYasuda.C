/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "tempCarreauYasuda.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(tempCarreauYasuda, 0);
    addToRunTimeSelectionTable
    (
        viscosityModel,
        tempCarreauYasuda,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::tempCarreauYasuda::calcNu() const
{
	const volScalarField& T= U_.mesh().lookupObj<volScalarField>("T");
	volScalarField test = Foam::exp(beta_ * (1.0/T -1.0/Tref) );
/*
    return
        nuInf_
      + (nu0_ - nuInf_)
       *pow(scalar(1) + pow(k_*strainRate(), a_), (n_ - 1.0)/a_);
*/
    return max
	(
	nuInf,
	nu0_ * Foam::exp(beta_ * (1.0/T - 1.0/Tref) )
	* pow(lambda_ * test * strainRate() + scalar(1), (n_ - 1.0)/a_)
	);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::tempCarreauYasuda::tempCarreauYasuda
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    tempCarreauYasudaCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    nu0_("nu0", dimViscosity, tempCarreauYasudaCoeffs_),
    nuInf_("nuInf", dimViscosity, tempCarreauYasudaCoeffs_),
    k_("k", dimTime, tempCarreauYasudaCoeffs_),
    n_("n", dimless, tempCarreauYasudaCoeffs_),
	beta_(tempCarreauYasudaCoeffs_.lookup("beta")),
	Tref_(tempCarreauYasudaCoeffs_.lookup("Tref"))
    a_
    (
        tempCarreauYasudaCoeffs_.lookupOrDefault
        (
            "a",
            dimensionedScalar("a", dimless, 2)
        )
    ),
    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::tempCarreauYasuda::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    tempCarreauYasudaCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    tempCarreauYasudaCoeffs_.lookup("nu0") >> nu0_;
    tempCarreauYasudaCoeffs_.lookup("nuInf") >> nuInf_;
    tempCarreauYasudaCoeffs_.lookup("k") >> k_;
    tempCarreauYasudaCoeffs_.lookup("n") >> n_;
	tempCarreauYasudaCoeffs_.lookup("beta") >> beta_;
	tempCarreauYasudaCoeffs_.lookup("Tref") >> Tref_;
    a_ = tempCarreauYasudaCoeffs_.lookupOrDefault
    (
        "a",
        dimensionedScalar("a", dimless, 2)
    );

    return true;
}


// ************************************************************************* //
