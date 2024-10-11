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
    Foam::viscosityMixture

Description
    This file is in the \dir kineticTheoryModels/viscosityMixture/viscosityMixture
    subdirectory inside the subdirectory
    \dir kineticTheoryModels/viscosityMixture inside the directory
    \dir kineticTheoryModels

SourceFiles
    viscosityMixture.C

\*---------------------------------------------------------------------------*/

#ifndef viscosityMixture_H
#define viscosityMixture_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class viscosityMixture Declaration
\*---------------------------------------------------------------------------*/

class viscosityMixture
{
    // Private member functions

        //- Disallow default bitwise copy construct
        viscosityMixture(const viscosityMixture&);

        //- Disallow default bitwise assignment
        void operator=(const viscosityMixture&);


protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("viscosityMixture");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        viscosityMixture,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

    //- Construct from components
    viscosityMixture(const dictionary& dict);


    // Selectors

    static autoPtr<viscosityMixture> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~viscosityMixture();


        // Member Functions


    //Info<<"viscosity mixture = \t"<<viscotyMixture<<endl;
	
        virtual tmp<volScalarField> mix_nu
        (
            const volScalarField& alpha,
            const dimensionedScalar& alphaMax
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
