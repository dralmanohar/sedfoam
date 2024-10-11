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

\*---------------------------------------------------------------------------*/

#include "GidaspowSchillerNaumann_veg.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GidaspowSchillerNaumann_veg, 0);

    addToRunTimeSelectionTable
    (
        vegDragModel,
        GidaspowSchillerNaumann_veg,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GidaspowSchillerNaumann_veg::GidaspowSchillerNaumann_veg
(
    const dictionary& interfaceDict,
    const phaseModel& phasea,
    const phaseModel& phaseb
)
:
    vegDragModel(interfaceDict, phasea, phaseb)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GidaspowSchillerNaumann_veg::~GidaspowSchillerNaumann_veg()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GidaspowSchillerNaumann_veg::K
(
    const volScalarField& Ur
) const
{
   // volScalarField beta(max(scalar(1) - alpha_, scalar(1e-6)));
    volScalarField beta(max(scalar(1) - phaseb_.alpha(), scalar(1e-6)));
    
    //    volScalarField bp(pow(beta, -2.65));

    volScalarField bp(pow(beta, -phaseb_.hExp()));

    //volScalarField Re(max(beta*Ur*phasea_.d()/phaseb_.nu(), scalar(1.0e-3)));
   
    const dimensionedScalar dragConst
    (
        interfaceDict_.getOrDefault
        (
            "dragConst",
            dimensionedScalar("dragConst",
                          dimensionSet(0, 0, 0, 0, 0, 0, 0),
                          1)
        )
    );

    const dimensionedScalar Coe = dragConst;
   
	   Info<<"dragConst in GidaspowSchiller = \t"<<dragConst<<endl; 

    volScalarField Re
    (
        max(beta*Ur*phaseb_.d()*phaseb_.sF()/phasea_.nu(), scalar(1.0e-9))
    );

    volScalarField Cds
    (
        neg(Re - 1000)*(24.0*(1.0 + 0.15*pow(Re, 0.687))/Re)
      + pos(Re - 1000)*0.44
    );

    return 0.75*Coe*Cds*phasea_.rho()*Ur*bp/(phaseb_.d()*phaseb_.sF());
}


/*
	//Info<<"Alphac Max = \t"<<gMax(phaseb_.alpha())<<"\t Min = \t"<<gMin(phaseb_.alpha())<<endl;//<<"\t beta Max = \t"<<max(1 - phaseb_.alpha())<<"\t Min = \t"<<min(1 - phaseb_.alpha())<<endl;
	//Info<<"Alpha Max = \t"<<gMax(alpha_)<<"\t Min = \t"<<gMin(alpha_)<<"\t beta Max = \t"<<max(1 - alpha_ - phaseb_.alpha())<<"\t Min = \t"<<min(1 -alpha_ - phaseb_.alpha())<<endl;

	Info<<" Manohar drag Gida =\t"<<endl;
	Info<<"Manohar Gida spow hExp = \t"<<phaseb_.hExp()<<"\t beta = max =\t"<<gMax(beta)<<"\t beta min = \t"<<gMin(beta)<<"\t nu pahse a = \t"<<phaseb_.nu()<<endl;
	Info<<"Manohar in Manohar file for Cds Max = \t"<<gMax(Cds)<<"\t Min = \t"<<gMin(Cds)<<"\t alpha_ max = \t"<<gMax(phaseb_.alpha())<<"\t alpha min = \t"<<gMin(phaseb_.alpha())<<endl;
	Info<<"Manohar in Manohar file for Alpha  Max = \t"<<gMax(alpha_)<<"\t Min = \t"<<gMin(alpha_)<<"\t alpha_ max = \t"<<gMax(phaseb_.alpha())<<"\t alpha min = \t"<<gMin(phaseb_.alpha())<<endl;
	Info<<"Manohar in Manohar file for Re Max = \t"<<gMax(Re)<<"\t Min = \t"<<gMin(Re)<<"\t Phase b = \t"<<phaseb_.d()<<"phase a rho = \t"<<phasea_.rho()<<endl;
	Info<<"Manohar in Remi Ub Max = \t"<<gMax(Ur)<<"\t Min = \t"<<gMin(Ur)<<endl;

	volScalarField Kmano = 0.75*Cds*phaseb_.rho()*Ur*bp/(phasea_.d()*phasea_.sF());

	//Info<<"Manohar in GidaspowSchillerNeauman recent case max= \t"<<gMax(Kmano)<<"\t min = \t"<<gMin(Kmano)<<"\t max bp = \t"<<gMax(bp)<<"\t min bp = \t"<<gMin(bp)<<"\t hExp = \t"<<phasea_.hExp()<<endl;
    	Info<<" Manohar drag Gida =\t"<<endl;
	Info<<"Manohar Gida spow hExp = \t"<<phaseb_.hExp()<<"\t beta = max =\t"<<gMax(beta)<<"\t beta min = \t"<<gMin(beta)<<"\t nu pahse a = \t"<<phasea_.nu()<<endl;
	Info<<"Manohar in Manohar file for Alpha from phase = \t"<<gMax(phaseb_.alpha())<<"\t Min = \t"<<gMin(phaseb_.alpha())<<endl;
	Info<<"Manohar in Manohar file for Cds Max = \t"<<gMax(Cds)<<"\t Min = \t"<<gMin(Cds)<<endl;
	Info<<"Manohar in Manohar file for Re Max = \t"<<gMax(Re)<<"\t Min = \t"<<gMin(Re)<<"\t Phase b = \t"<<phaseb_.d()<<"phase a rho = \t"<<phasea_.rho()<<endl;
	Info<<"Manohar in Remi Ub Max = \t"<<gMax(Ur)<<"\t Min = \t"<<gMin(Ur)<<"\t phasea_.rho() = \t"<<phasea_.rho()<<endl;
    	
	volScalarField Kmano = 0.75*Cds*phasea_.rho()*Ur*bp/(phaseb_.d()*phaseb_.sF());

	Info<<"Manohar in GidaspowSchillerNeauman recent case max coefficient= \t"<<gMax(Kmano)<<"\t min = \t"<<gMin(Kmano)<<"\t diameter = \t"<<phaseb_.d()<<"\t bp max = \t"<< gMax(bp)<<"\t Min = \t"<<gMin(bp)<<endl;
*/

// ************************************************************************* //
