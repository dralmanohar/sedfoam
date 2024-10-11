/*---------------------------------------------------------------------------*\
Copyright (C) 2015 Cyrille Bonamy, Julien Chauchat, Tian-Jian Hsu
                   and contributors

License
    This file is part of SedFOAM.

    SedFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SedFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SedFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/



Switch MomentumBudget
(
    runTime.controlDict().lookupOrDefault("MomentumBudget", false)
);



Info<<"Manohar Momentum Budget in create File = \t"<<MomentumBudget<<endl;
/*
volTensorField Ei
(
    IOobject
    (
        "Ei",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "Ei",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        tensor::zero
    )
);
*/

/*
volTensorField Ei
(
    IOobject
    (
         "Ei",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    mesh
    
);
*/

volScalarField pMech
(
    IOobject
    (
        "pMech",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pMech",
        dimPressure,
        0
    )
);

volScalarField pA
(
    IOobject
    (
        "pA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pA",
        dimPressure,
        0
    )
);


volScalarField pFF
(
    IOobject
    (
        "pFF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pFF",
        dimPressure,
        0
    )
);

volVectorField rho_mixd
(
    IOobject
    (
        "rho_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "rho_mixd",
        dimensionSet(1, -2, -1, 0, 0, 0, 0),
        vector (0, 0, 0)
    )
);

volTensorField nonl_mixd
(
    IOobject
    (
        "nonl_mixd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "nonl_mixd",
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        tensor::zero
    )
);

volVectorField bouy_mixd
    (
        IOobject
        (
            "bouy_mixd",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "bouy_mixd",
            dimensionSet(1, -2, -2, 0, 0, 0, 0),
            vector (0, 0, 0)
        )
    );


/*
volVectorField rho_mix
(
    IOobject
    (
        "rho_mixd",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),

     mesh,
     dimensionedVector
        (
         "rho_mixd",
          dimensionSet(0, 1, -1, 0, 0, 0, 0),
            0.0
        )
);
*/



volVectorField SolidPressureContr
(
    IOobject
    (
        "SolidPressureContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(p_rbgh)
);

volVectorField divTauaContr
(
    IOobject
    (
        "divTauaContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);


volVectorField ViscStressContr
(
    IOobject
    (
        "ViscStressContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField FrictContr
(
    IOobject
    (
        "FrictContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField phiRaContr
(
    IOobject
    (
        "phiRaContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField phiRa2Contr
(
    IOobject
    (
        "phiRa2Contr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField divRcaContr
(
    IOobject
    (
        "divRcaContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField RcaContr
(
    IOobject
    (
        "RcaContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField BulkContr
(
    IOobject
    (
        "BulkContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField gravityContr
(
    IOobject
    (
        "gravityContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField DragContr
(
    IOobject
    (
        "DragContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField SuspContr
(
    IOobject
    (
        "SuspContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volVectorField ExtContr
(
    IOobject
    (
        "ExtContr",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    fvc::grad(pff)
);

volScalarField Errosion
   		(
       		IOobject
				(
            		"Errosion",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"Errosion",
            			dimensionSet(0, 1, -1, 0, 0, 0, 0),
            			0.0
        		)
 );

volScalarField Deposition
   		(
       		IOobject
				(
            		"Deposition",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"Deposition",
            			dimensionSet(0, 1, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);

/*
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rbgh,
    mesh.solutionDict().subDict("PIMPLE"),
    pRefCell,
    pRefValue
);
if (p_rbgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rbgh = p - rhob*gh;
}
*/

