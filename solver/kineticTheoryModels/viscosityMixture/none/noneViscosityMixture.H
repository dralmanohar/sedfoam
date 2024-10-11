/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::noneViscosity

Description
    This file is in the \dir kineticTheoryModels/viscosityMixture/none
    subdirectory inside the subdirectory
    \dir kineticTheoryModels/viscosityMixture inside the directory
    \dir kineticTheoryModels

SourceFiles
    noneViscosityMixture.C

\*---------------------------------------------------------------------------*/

#ifndef noneViscosityMixture_H
#define noneViscosityMixture_H

#include "viscosityMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class noneViscosity Declaration
\*---------------------------------------------------------------------------*/

class noneViscosityMixture
:
    public viscosityMixture
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

    //- Construct from components
    noneViscosityMixture(const dictionary& dict);


    //- Destructor
    virtual ~noneViscosityMixture();


    // Member functions


        tmp<volScalarField> mix_nu
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax
        ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
