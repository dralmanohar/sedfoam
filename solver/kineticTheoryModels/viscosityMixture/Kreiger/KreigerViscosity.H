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
    Foam::GidaspowViscosity

Description
    \verbatim
    "A Bubbling Fluidization Model Using Kinetic Theory of Granular Flow",
    Ding, J., Gidaspow, D., AIChE Journal, 36(4):523–538
    \endverbatim

    This file is in the \dir kineticTheoryModels/viscosityMixure/Kreiger
    subdirectory inside the subdirectory
    \dir kineticTheoryModels/viscosityMixture inside the directory
    \dir kineticTheoryModels

SourceFiles
    GidaspowViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef KreigerViscosity_H
#define KreigerViscosity_H

#include "viscosityMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class Kreigerviscosity Declaration
\*---------------------------------------------------------------------------*/

class KreigerViscosity
:
    public viscosityMixture
{

public:

    //- Runtime type information
    TypeName("Kreiger");

    // Constructors

        //- Construct from components
        KreigerViscosity(const dictionary& dict);


    //- Destructor
    virtual ~KreigerViscosity();


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
