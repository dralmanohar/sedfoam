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
dimensionedScalar zeroo
(
"zeroo",
dimensionSet(0, 2, -2, 0, 0, 0, 0),
scalar(0.0)
);



Switch favreAveraging
(
    runTime.controlDict().lookupOrDefault("favreAveraging", false)
);


Switch Reset
(
 runTime.controlDict().lookupOrDefault("Reset", false)
 );

Switch TKEBudgetLES
(
    runTime.controlDict().lookupOrDefault("TKEBudgetLES", false)
);

if (Reset)
{
Info<<"Reset options = \t"<<Reset<<endl;
}


Info<<"Reset Mode = \t"<<Reset<<endl;

scalar timeaverage = 0;
scalar runtime     = 0;
scalar deltaT = 0;


volScalarField WcubeMean
(
    IOobject
    (
        "WcubeMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "WcubeMean",
        dimVelocity*dimVelocity*dimVelocity,
        0.0
    )
);

volScalarField UsquarWMean
(
    IOobject
    (
        "UsquarWMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "UsquarWMean",
        dimVelocity*dimVelocity*dimVelocity,
        0.0
    )
);


volScalarField turbDiffMeanDavid
(
    IOobject
    (
        "turbDiffMeanDavid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "turbDiffMeanDavid",
        dimVelocity*dimVelocity*dimVelocity,
        0.0
    )
);

volVectorField turbDiffusionMeanDavid
(
    IOobject
    (
        "turbDiffusionMeanDavid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "turbDiffusionMeanDavid",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        vector::zero
    )
);


volScalarField testVegMean
   		(
       		IOobject
				(
            		"testVegMean",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"testVegMean",
            			dimensionSet(1, -3, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);


volVectorField UbMean
   		(
       		IOobject
       			(
       				"UbMean",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedVector
       			(
       				"UbMean",
       				dimVelocity,
       				vector (0, 0, 0)
       			)
   		);

volScalarField divalphaUaMean
   		(
       		IOobject
				(
            		"divalphaUaMean",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"divalphaUaMean",
            			dimensionSet(0, 0, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);


volScalarField divalphaUafluact
   		(
       		IOobject
				(
            		"divalphaUafluact",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"divalphaUafluact",
            			dimensionSet(0, 0, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);

volScalarField divalphaUa
   		(
       		IOobject
				(
            		"divalphaUa",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"divalphaUa",
            			dimensionSet(0, 0, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);

volScalarField ddtalpha
   		(
       		IOobject
				(
            		"ddtalpha",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"ddtalpha",
            			dimensionSet(0, 0, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);

volScalarField ddtalphaMean
   		(
       		IOobject
				(
            		"ddtalphaMean",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"ddtalphaMean",
            			dimensionSet(0, 0, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);

volScalarField ddtalphaPrimeF
   		(
       		IOobject
				(
            		"ddtalphaPrimeF",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"ddtalphaPrimeF",
            			dimensionSet(0, 0, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);

volScalarField alphaPrimeF
   		(
       		IOobject
				(
            		"alphaPrimeF",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"alphaPrimeF",
            			dimensionSet(0, 0, 0, 0, 0, 0, 0),
            			0.0
        		)
    	);

		
volVectorField UaMean
   		(
       		IOobject
       			(
       				"UaMean",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedVector
       			(
       				"UaMean",
       				dimVelocity,
       				vector (0, 0, 0)
       			)
   		);

volVectorField alphaUaMeandiv
   		(
       		IOobject
       			(
       				"alphaUaMeandiv",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedVector
       			(
       				"alphaUaMeandiv",
       				dimVelocity,
       				vector (0, 0, 0)
       			)
   		);

surfaceScalarField fluxUaalpha
   		(
       		IOobject
       			(
       				"fluxUaalpha",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedScalar
       			(
       				"fluxUaalpha",
       				dimensionSet(0, 3, -1, 0, 0, 0, 0),
       				0.0
       			)
   		);

surfaceScalarField fluxUaMeanalphaMean
   		(
       		IOobject
       			(
       				"fluxUaMeanalphaMean",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedScalar
       			(
       				"fluxUaMeanalphaMean",
       				dimensionSet(0, 3, -1, 0, 0, 0, 0),
       				0.0
       			)
   		);
surfaceScalarField fluxUaPrimealphaPrime
   		(
       		IOobject
       			(
       				"fluxUaPrimealphaPrime",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedScalar
       			(
       				"fluxUaPrimealphaPrime",
       				dimensionSet(0, 3, -1, 0, 0, 0, 0),
       				0.0
       			)
   		);

volScalarField divfluxUaalpha
   		(
       		IOobject
				(
            		"divUaalpha",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"divUaalpha",
            			dimensionSet(0, 1, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);



volVectorField alphaUafluacdiv
   		(
       		IOobject
       			(
       				"alphaUafluacdiv",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedVector
       			(
       				"alphaUafluacdiv",
       				dimVelocity,
       				vector (0, 0, 0)
       			)
   		);

volVectorField alphaUa
   		(
       		IOobject
       			(
       				"alphaUa",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedVector
       			(
       				"alphaUa",
       				dimVelocity,
       				vector (0, 0, 0)
       			)
   		);

///
//
surfaceScalarField fluxUbbeta
   		(
       		IOobject
       			(
       				"fluxUbbeta",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedScalar
       			(
       				"fluxUbbeta",
       				dimensionSet(0, 3, -1, 0, 0, 0, 0),
       				0.0
       			)
   		);


volScalarField beta_bMean
    (
        IOobject
        (
            "beta_bMean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "beta_bMean",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField divbetaUbMean
   		(
       		IOobject
				(
            		"divbetaUbMean",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"divbetaUbMean",
            			dimensionSet(0, 0, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);


volScalarField divbetaUbfluact
   		(
       		IOobject
				(
            		"divbetaUbfluact",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"divbetaUbfluact",
            			dimensionSet(0, 0, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);

volScalarField divbetaUb
   		(
       		IOobject
				(
            		"divbetaUb",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"divbetaUb",
            			dimensionSet(0, 0, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);

volScalarField ddtbeta
   		(
       		IOobject
				(
            		"ddtbeta",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"ddtbeta",
            			dimensionSet(0, 0, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);

volScalarField ddtbetaMean
   		(
       		IOobject
				(
            		"ddtbetaMean",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"ddtbetaMean",
            			dimensionSet(0, 0, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);

volScalarField ddtbetaPrimeF
   		(
       		IOobject
				(
            		"ddtbetaPrimeF",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"ddtbetaPrimeF",
            			dimensionSet(0, 0, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);

volScalarField betaPrimeF
   		(
       		IOobject
				(
            		"betuaPrimeF",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"betaPrimeF",
            			dimensionSet(0, 0, 0, 0, 0, 0, 0),
            			0.0
        		)
    	);

		




volVectorField betaUbMeandiv
   		(
       		IOobject
       			(
       				"betaUbMeandiv",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedVector
       			(
       				"betaUbMeandiv",
       				dimVelocity,
       				vector (0, 0, 0)
       			)
   		);




volScalarField divfluxUbbeta
   		(
       		IOobject
				(
            		"divUbbeta",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
			IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"divUbbeta",
            			dimensionSet(0, 1, -1, 0, 0, 0, 0),
            			0.0
        		)
    	);



volVectorField betaUbfluacdiv
   		(
       		IOobject
       			(
       				"betaUbfluacdiv",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedVector
       			(
       				"betaUbfluacdiv",
       				dimVelocity,
       				vector (0, 0, 0)
       			)
   		);

volVectorField betaUb
   		(
       		IOobject
       			(
       				"betaUb",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedVector
       			(
       				"betaUb",
       				dimVelocity,
       				vector (0, 0, 0)
       			)
   		);

///		
volScalarField TKEMean
(
    IOobject
    (
        "TKEMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TKEMean",
        dimVelocity*dimVelocity,
        0
    )
);


volScalarField TKEMeanProd
(
    IOobject
    (
        "TKEMeanProd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "TKEMeanProd",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0
    )
);

volTensorField strainTensorF
	(
    		IOobject
    			(
        			"strainTensorF",
        			runTime.timeName(),
        			mesh,
        			IOobject::NO_READ,
        			IOobject::NO_WRITE
    			),
    		mesh,
    		dimensionedTensor
    			(
        			"strainTensorF",
        			dimVelocity/dimLength,
        			tensor::zero
    			)
	);

volScalarField viscDissMeanI
(
    IOobject
    (
        "viscDissMeanI",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDissMeanI",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volVectorField betaFluct
(
    IOobject
    (
        "betaFluct",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "betaFluct",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volVectorField betaFluctMean
(
    IOobject
    (
        "betaFluctMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "betaFluctMean",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volScalarField RS2
    (
        IOobject
        (
            "RS2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS2",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField RS4
    (
        IOobject
        (
            "RS4",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS4",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField pMeanF
(
    IOobject
    (
        "pMeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pMeanF",
        dimPressure,
        0
    )
);

volTensorField UaPrime2MeanF
	(
   			IOobject
   				(
					"UaPrime2MeanF",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
   				),
   			mesh,
   			dimensionedTensor
   				(
       				"UaPrime2MeanF",
       				dimVelocity*dimVelocity,
       				tensor::zero
   				)
	);
	
volTensorField UbPrime2MeanF
	(
   			IOobject
   				(
       				"UbPrime2MeanF",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
   				),
   			mesh,
   			dimensionedTensor
   				(
       				"UbPrime2MeanF",
       				dimVelocity*dimVelocity,
       				tensor::zero
  				 )
	);

volTensorField betaSqrUbMean
(
    IOobject
    (
        "betaSqrUbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "betaSqrUbMean",
        dimVelocity*dimVelocity,
        tensor::zero
    )
);


volTensorField alphaSqrUaMean
(
    IOobject
    (
        "alphaSqrUaMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "alphaSqrUaMean",
        dimVelocity*dimVelocity,
        tensor::zero
    )
);

volScalarField RS2Mean
    (
        IOobject
        (
            "RS2Mean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS2Mean",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField RS4Mean
    (
        IOobject
        (
            "RS4Mean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "RS4Mean",
            dimensionSet(0, 2, -2, 0, 0, 0, 0),
            0.0
        )
    );

volScalarField SijFSijF
	(
    		IOobject
    			(
        			"SijFSijF",
        			runTime.timeName(),
        			mesh,
        			IOobject::NO_READ,
        			IOobject::NO_WRITE
    			),
    		mesh,
    		dimensionedScalar
    			(
        			"SijFSijF",
        			dimVelocity/dimLength*dimVelocity/dimLength,
        			0
   			 )
	);

volScalarField SGSDissMean
(
    IOobject
    (
        "SGSDissMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDissMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField SGSDissMeanJ
(
    IOobject
    (
        "SGSDissMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDissMeanJ",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField SGSDissMeanI
(
    IOobject
    (
        "SGSDissMeanI",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDissMeanI",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);


volSymmTensorField fluctstrainTensor
(
    IOobject
    (
        "fluctstrainTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "fluctstrainTensor",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);

volSymmTensorField strainTensor
(
    IOobject
    (
        "strainTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "strainTensor",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);

volSymmTensorField strainTensorM
(
    IOobject
    (
        "strainTensorF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "strainTensorF",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);

volScalarField turbDiffusionMean
(
    IOobject
    (
        "turbDiffusionMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "turbDiffusionMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField turbDiffusionMeanJ
(
    IOobject
    (
        "turbDiffusionMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "turbDiffusionMeanJ",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField pressDiffusionMean
(
    IOobject
    (
        "pressDiffusionMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pressDiffusionMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField pressDiffusionMeanJ
(
    IOobject
    (
        "pressDiffusionMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pressDiffusionMeanJ",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField SGSDiffusionMean
(
    IOobject
    (
        "SGSDiffusionMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDiffusionMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField SGSDiffusionMeanJ
(
    IOobject
    (
        "SGSDiffusionMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "SGSDiffusionMeanJ",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volSymmTensorField fluctstrainTensorF
(
    IOobject
    (
        "fluctstrainTensorF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "fluctstrainTensorF",
        dimVelocity/dimLength,
        symmTensor::zero
    )
);

volVectorField turbDiffMean
(
    IOobject
    (
        "turbDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "turbDiffMean",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volVectorField turbDiffMeanJ
(
    IOobject
    (
        "turbDiffMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "turbDiffMeanJ",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volScalarField SijSijM
	(
    		IOobject
    			(
        			"SijSijM",
        			runTime.timeName(),
        			mesh,
        			IOobject::NO_READ,
        			IOobject::NO_WRITE
    			),
    		mesh,
    		dimensionedScalar
    			(
        			"SijSijM",
        			dimVelocity/dimLength*dimVelocity/dimLength,
        			0
   			 )
	);

volVectorField pressDiffMean
(
    IOobject
    (
        "pressDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "pressDiffMean",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);


volVectorField pressDiffMeanJ
(
    IOobject
    (
        "pressDiffMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "pressDiffMeanJ",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volScalarField viscDiffusionMean
(
    IOobject
    (
        "viscDiffusionMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDiffusionMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volScalarField viscDiffusionMeanJ
(
    IOobject
    (
        "viscDiffusionMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "viscDiffusionMeanJ",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);

volVectorField SGSDiffMean
(
    IOobject
    (
        "SGSDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "SGSDiffMean",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volVectorField SGSDiffMeanJ
(
    IOobject
    (
        "SGSDiffMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "SGSDiffMeanJ",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volVectorField viscDiffMean
(
    IOobject
    (
        "viscDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "viscDiffMean",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volVectorField viscDiffMeanJ
(
    IOobject
    (
        "viscDiffMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "viscDiffMeanJ",
        dimVelocity*dimVelocity*dimVelocity,
        vector::zero
    )
);

volScalarField dragConstbMean
(
    IOobject
    (
        "dragConstbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "dragConstbMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);


volScalarField dragVegMean
(
    IOobject
    (
        "dragVegMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "dragVegMean",
        dimVelocity*dimVelocity*dimVelocity/dimLength,
        0.0
    )
);


volVectorField UbPrimeF
	(
    		IOobject
    			(
        			"UbPrimeF",
        			runTime.timeName(),
        			mesh,
        			IOobject::NO_READ,
        			IOobject::NO_WRITE
    			),
    		mesh,
    		dimensionedVector
    			(
        			"UbPrimeF",
        			dimVelocity,
        			vector (0, 0, 0)
    			)
	);
	
volVectorField UaPrimeF
	(
    		IOobject
    			(
        			"UaPrimeF",
        			runTime.timeName(),
        			mesh,
        			IOobject::NO_READ,
        			IOobject::NO_WRITE
    			),
    		mesh,
    		dimensionedVector
    			(
        			"UaPrimeF",
        			dimVelocity,
        			vector (0, 0, 0)
    			)
	);



volScalarField alpha_aMean
   		(
       		IOobject
				(
            		"alpha_aMean",
            		runTime.timeName(),
            		mesh,
            		IOobject::NO_READ,
					IOobject::NO_WRITE
        		),
        	mesh,
        	dimensionedScalar
        		(
            		"alpha_aMean",
            			dimensionSet(0, 0, 0, 0, 0, 0, 0),
            			0.0
        		)
    	);


volVectorField UaMeanF
	(
   			IOobject
   				(
       				"UaMeanF",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
   				),
   			mesh,
   			dimensionedVector
   				(
					"UaMeanF",
       				dimVelocity,
       				vector (0, 0, 0)
   				)
	);


volVectorField UbMeanF
	(
   			IOobject
   				(
       				"UbMeanF",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
   				),
   			mesh,
   			dimensionedVector
   				(
       				"UbMeanF",
       				dimVelocity,
       				vector (0, 0, 0)
   				)
		);

volVectorField betaUbMean
   		(
       		IOobject
       			(
       				"betaUbMean",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
   			mesh,
   			dimensionedVector
   				(
   					"betaUbMean",
   					dimVelocity,
   					vector (0, 0, 0)
				)
   		);


volVectorField alphaUaMean
   		(
       		IOobject
       			(
       				"alphaUaMean",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
       			),
       		mesh,
       		dimensionedVector
       			(
       				"alphaUaMean",
       				dimVelocity,
       				vector (0, 0, 0)
       			)
   		);

	
volScalarField viscDissMean
	(
   			IOobject
   				(
       				"viscDissMean",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
   				),
   			mesh,
   			dimensionedScalar
   				(
       				"viscDissMean",
       				dimVelocity*dimVelocity*dimVelocity/dimLength,
       				0.0
   				)
	);

volScalarField viscDissMeanJ
	(
   			IOobject
   				(
       				"viscDissMeanJ",
       				runTime.timeName(),
       				mesh,
       				IOobject::NO_READ,
       				IOobject::NO_WRITE
   				),
   			mesh,
   			dimensionedScalar
   				(
       				"viscDissMeanJ",
       				dimVelocity*dimVelocity*dimVelocity/dimLength,
       				0.0
   				)
	);

volScalarField pPrimeF
(
    IOobject
    (
        "pPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "pPrimeF",
        dimPressure,
        0
    )
);



if (Reset)
{

Info<<"Manohar in Reset mode = \t"<<endl;

runtime = 0;

}

else 

{


	runtime = runTime.value();// 0;// this->runTime.timeName();

	Info<<"runtime in create file =\t"<<runtime<<endl;
	Info<<"Manohar is not in Reset Mode = \t"<<runTime.value()<<endl;

IOobject betaFluctMeanHeader
  		 (
        		"betaFluctMean",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ
    		);

if (betaFluctMeanHeader.typeHeaderOk<volVectorField>(true))
		{
    			volVectorField betaFluctMean2
        			(
            				IOobject
            					(
                					"betaFluctMean",
                					runTime.timeName(),
                					mesh,
                					IOobject::MUST_READ,
                					IOobject::NO_WRITE
            					),
            				mesh
        			);
    			betaFluctMean = betaFluctMean2;
		}


IOobject alpha_aMeanHeader
    		(
        		"alpha_aMean",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ
    		);


if (alpha_aMeanHeader.typeHeaderOk<volScalarField>(true))
{
	volScalarField alpha_aMean2
        			(
            				IOobject
            					(
                					"alpha_aMean",
                					runTime.timeName(),
                					mesh,
                					IOobject::MUST_READ,
                					IOobject::NO_WRITE
            					),
            				mesh
        			);
    				alpha_aMean = alpha_aMean2;


          IOobject RS2MeanHeader
    (
        "RS2Mean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (RS2MeanHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField RS2Mean2
        (
            IOobject
            (
                "RS2Mean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        RS2Mean = RS2Mean2;
    }
    // Checks if field RS4Mean is present
    IOobject RS4MeanHeader
    (
        "RS4Mean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (RS4MeanHeader.typeHeaderOk<volScalarField>(true))
    {
        volScalarField RS4Mean2
        (
            IOobject
            (
                "RS4Mean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        RS4Mean = RS4Mean2;
    }

    IOobject viscDiffMeanHeader
    (
        "viscDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (viscDiffMeanHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField viscDiffMean2
(
    IOobject
    (
        "viscDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
viscDiffMean = viscDiffMean2;
}


    IOobject viscDiffMeanJHeader
    (
        "viscDiffMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (viscDiffMeanJHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField viscDiffMeanJ2
(
    IOobject
    (
        "viscDiffMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
viscDiffMeanJ = viscDiffMeanJ2;
}


    IOobject SGSDiffMeanHeader
    (
        "SGSDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDiffMeanHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField SGSDiffMean2
(
    IOobject
    (
        "SGSDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDiffMean = SGSDiffMean2;
}


    IOobject SGSDiffMeanJHeader
    (
        "SGSDiffMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDiffMeanJHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField SGSDiffMeanJ2
(
    IOobject
    (
        "SGSDiffMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDiffMeanJ = SGSDiffMeanJ2;
}


   IOobject WcubeMeanHeader
    (
        "WcubeMeanMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (WcubeMeanHeader.typeHeaderOk<volTensorField>(true))
    {
volScalarField WcubeMean2
(
    IOobject
    (
        "WcubeMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
WcubeMean = WcubeMean2;
}


   IOobject UsquarWMeanHeader
    (
        "UsquarWMeanMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (UsquarWMeanHeader.typeHeaderOk<volTensorField>(true))
    {
volScalarField UsquarWMean2
(
    IOobject
    (
        "UsquarWMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
UsquarWMean = UsquarWMean2;
}



     IOobject turbDiffMeanHeader
    (
        "turbDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffMeanHeader.typeHeaderOk<volTensorField>(true))
    {
 
 volVectorField turbDiffMean2
(
    IOobject
    (
        "turbDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
turbDiffMean = turbDiffMean2;
}


IOobject turbDiffMeanJHeader
    (
        "turbDiffMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffMeanJHeader.typeHeaderOk<volTensorField>(true))
    {
 
 volVectorField turbDiffMeanJ2
(
    IOobject
    (
        "turbDiffMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

turbDiffMeanJ = turbDiffMeanJ2;
}



 IOobject turbDiffMeanDavidHeader
    (
        "turbDiffMeanDavid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (turbDiffMeanDavidHeader.typeHeaderOk<volTensorField>(true))
    {
 
 volVectorField turbDiffMean2David
(
    IOobject
    (
        "turbDiffMeanDavid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
turbDiffMean = turbDiffMean2David;
}




    IOobject dragConstbMeanHeader
    (
        "dragConstbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (dragConstbMeanHeader.typeHeaderOk<volTensorField>(true))
    {
volScalarField dragConstbMean2
(
    IOobject
    (
        "dragConstbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

dragConstbMean = dragConstbMean2;

}


    IOobject dragVegMeanHeader
    (
        "dragVegMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (dragVegMeanHeader.typeHeaderOk<volTensorField>(true))
    {
volScalarField dragVegMean2
(
    IOobject
    (
        "dragVegMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

dragVegMean = dragVegMean2;

}

     IOobject pressDiffMeanHeader
    (
        "pressDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pressDiffMeanHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField pressDiffMean2
(
    IOobject
    (
        "pressDiffMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pressDiffMean = pressDiffMean2;
}


     IOobject pressDiffMeanJHeader
    (
        "pressDiffMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pressDiffMeanJHeader.typeHeaderOk<volTensorField>(true))
    {
volVectorField pressDiffMeanJ2
(
    IOobject
    (
        "pressDiffMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pressDiffMeanJ = pressDiffMeanJ2;
}


	// Checks if field UaMeanF is present
    	IOobject UaMeanFHeader
    				(
        				"UaMeanF",
        				runTime.timeName(),
        				mesh,
        				IOobject::NO_READ
    				);
	
	if (UaMeanFHeader.typeHeaderOk<volVectorField>(true))
    				{
        				volVectorField UaMeanF2
        					(
            						IOobject
            							(
                							"UaMeanF",
                							runTime.timeName(),
                							mesh,
                							IOobject::MUST_READ,
                							IOobject::NO_WRITE
            							),
            						mesh
        					);
        				UaMeanF = UaMeanF2;
        				alphaUaMean = UaMeanF*(alpha_aMean+1e-20);
    				}
	
			// Checks if field UbMeanF is present
	
	IOobject UbMeanFHeader
    				(
        				"UbMeanF",
        				runTime.timeName(),
        				mesh,
        				IOobject::NO_READ
    				);
	
	if (UbMeanFHeader.typeHeaderOk<volVectorField>(true))
    				{
        				volVectorField UbMeanF2
        					(
            						IOobject
            							(
                							"UbMeanF",
                							runTime.timeName(),
                							mesh,
                							IOobject::MUST_READ,
                							IOobject::NO_WRITE
            							),
            						mesh
        					);
        				UbMeanF = UbMeanF2;
        				betaUbMean = UbMeanF*(1-alpha_aMean);
    				}

		
	    
	
			// Checks if field UbPrime2MeanF is present
	
	IOobject UbPrime2MeanFHeader
    				(
        				"UbPrime2MeanF",
        				runTime.timeName(),
        				mesh,
        				IOobject::NO_READ
    				);

 	if (UbPrime2MeanFHeader.typeHeaderOk<volTensorField>(true))
    				{
        				volTensorField UbPrime2MeanF2
        					(
            						IOobject
            							(
                							"UbPrime2MeanF",
                							runTime.timeName(),
                							mesh,
                							IOobject::MUST_READ,
                							IOobject::NO_WRITE
            							),
            						mesh
        					);
        				UbPrime2MeanF = UbPrime2MeanF2;
				        betaSqrUbMean = (UbPrime2MeanF+UbMeanF*UbMeanF)*(1-alpha_aMean);


    				}

	   // Checks if field UaPrime2MeanF is present
    IOobject UaPrime2MeanFHeader
    (
        "UaPrime2MeanF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (UaPrime2MeanFHeader.typeHeaderOk<volTensorField>(true))
    {
        volTensorField UaPrime2MeanF2
        (
            IOobject
            (
                "UaPrime2MeanF",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        UaPrime2MeanF = UaPrime2MeanF2;
        alphaSqrUaMean = (UaPrime2MeanF+UaMeanF*UaMeanF)
                        *(alpha_aMean+1e-20);
    }

    				
  	IOobject viscDissMeanHeader
    		(
        		"viscDissMean",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ
    		);
    	
	if (viscDissMeanHeader.typeHeaderOk<volScalarField>(true))
    		{
        		volScalarField viscDissMean2
        			(
            				IOobject
            					(
                					"viscDissMean",
                					runTime.timeName(),
                					mesh,
                					IOobject::MUST_READ,
                					IOobject::NO_WRITE
            					),
            				mesh
        			);
        		viscDissMean = viscDissMean2;
    		}
    		
    		
  	IOobject viscDissMeanJHeader
    		(
        		"viscDissMeanJ",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ
    		);
    	
	if (viscDissMeanJHeader.typeHeaderOk<volScalarField>(true))
    		{
        		volScalarField viscDissMeanJ2
        			(
            				IOobject
            					(
                					"viscDissMeanJ",
                					runTime.timeName(),
                					mesh,
                					IOobject::MUST_READ,
                					IOobject::NO_WRITE
            					),
            				mesh
        			);
        		viscDissMeanJ = viscDissMeanJ2;
    		}



IOobject pMeanHeader
    (
        "pMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

 if (pMeanHeader.typeHeaderOk<volTensorField>(true))
    {

    volScalarField pMean2
(
    IOobject
    (
        "pMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
  pMeanF = pMean2;
}


     IOobject pPrimeFHeader
    (
        "pPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (pPrimeFHeader.typeHeaderOk<volTensorField>(true))
    {

  volScalarField pPrimeF2
(
    IOobject
    (
        "pPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
pPrimeF = pPrimeF2;
}


     IOobject SGSDissMeanHeader
    (
        "SGSDissMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDissMeanHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField SGSDissMean2
(
    IOobject
    (
        "SGSDissMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDissMean = SGSDissMean2;
}


     IOobject SGSDissMeanJHeader
    (
        "SGSDissMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (SGSDissMeanJHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField SGSDissMeanJ2
(
    IOobject
    (
        "SGSDissMeanJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
SGSDissMeanJ = SGSDissMeanJ2;
}

///

     IOobject fluxUbbetaHeader
    (
        "fluxUbbeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (fluxUbbetaHeader.typeHeaderOk<volTensorField>(true))
    {

surfaceScalarField fluxUbbeta2
(
    IOobject
    (
        "fluxUbbeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
fluxUbbeta = fluxUbbeta2;
}


IOobject beta_bMeanHeader
    (
        "beta_bMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (beta_bMeanHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField beta_bMean2
(
    IOobject
    (
        "beta_bMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
beta_bMean = beta_bMean2;
}


IOobject divbetaUbMeanHeader
    (
        "divbetaUbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divbetaUbMeanHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField divbetaUbMean2
(
    IOobject
    (
        "divbetaUbMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divbetaUbMean = divbetaUbMean2;
}


IOobject divbetaUbfluactHeader
    (
        "divbetaUbfluact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divbetaUbfluactHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField divbetaUbfluact2
(
    IOobject
    (
        "divbetaUbfluact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divbetaUbfluact = divbetaUbfluact2;
}



IOobject divbetaUbHeader
    (
        "divbetaUb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divbetaUbHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField divbetaUb2
(
    IOobject
    (
        "divbetaUb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divbetaUb = divbetaUb2;
}


IOobject ddtbetaHeader
    (
        "ddtbeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (ddtbetaHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField ddtbeta2
(
    IOobject
    (
        "ddtbeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
ddtbeta = ddtbeta2;
}


IOobject ddtbetaMeanHeader
    (
        "ddtbetaMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (ddtbetaMeanHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField ddtbetaMean2
(
    IOobject
    (
        "ddtbetaMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
ddtbetaMean = ddtbetaMean2;
}


IOobject ddtbetaPrimeFHeader
    (
        "ddtbetaPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (ddtbetaPrimeFHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField ddtbetaPrimeF2
(
    IOobject
    (
        "ddtbetaPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
ddtbetaPrimeF = ddtbetaPrimeF2;
}

IOobject betaPrimeFHeader
    (
        "betaPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (betaPrimeFHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField ddtbetaPrimeF2
(
    IOobject
    (
        "betaPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
ddtbetaPrimeF = ddtbetaPrimeF2;
}


IOobject divfluxUbbetaHeader
    (
        "divfluxUbbeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divfluxUbbetaHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField divfluxUbbeta2
(
    IOobject
    (
        "divfluxUbbeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divfluxUbbeta = divfluxUbbeta2;
}
/*
IOobject fluxUbbetaHeader
    (
        "fluxUbbeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (fluxUbbetaHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField fluxUbbeta2
(
    IOobject
    (
        "fluxUbbeta",
        runTime.timeName(),
        mesh,translate
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
fluxUbbeta = fluxUbbeta2;
}
*/

IOobject betaUbMeandivHeader
    (
        "betaUbMeandiv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (betaUbMeandivHeader.typeHeaderOk<volTensorField>(true))
    {

volVectorField betaUbMeandiv2
(
    IOobject
    (
        "betaUbMeandiv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
betaUbMeandiv = betaUbMeandiv2;
}

IOobject betaUbfluacdivHeader
    (
        "betaUbfluacdiv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (betaUbfluacdivHeader.typeHeaderOk<volTensorField>(true))
    {

volVectorField betaUbfluacdiv2
(
    IOobject
    (
        "betaUbfluacdiv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
betaUbfluacdiv = betaUbfluacdiv2;
}


IOobject betaUbHeader
    (
        "betaUb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (betaUbHeader.typeHeaderOk<volTensorField>(true))
    {

volVectorField betaUb2
(
    IOobject
    (
        "betaUb",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
betaUb = betaUb2;
}

///

     IOobject fluxUaalphaHeader
    (
        "fluxUaalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (fluxUaalphaHeader.typeHeaderOk<volTensorField>(true))
    {

surfaceScalarField fluxUaalpha2
(
    IOobject
    (
        "fluxUaalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
fluxUaalpha = fluxUaalpha2;
}


IOobject alpha_aMeanHeader
    (
        "alpha_aMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (alpha_aMeanHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField alpha_aMean2
(
    IOobject
    (
        "alpha_aMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
alpha_aMean = alpha_aMean2;
}


IOobject divalphaUaMeanHeader
    (
        "divalphaUaMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divalphaUaMeanHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField divalphaUaMean2
(
    IOobject
    (
        "divalphaUaMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divalphaUaMean = divalphaUaMean2;
}



IOobject divalphaUafluactHeader
    (
        "divalphaUafluact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divalphaUafluactHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField divalphaUafluact2
(
    IOobject
    (
        "divalphaUafluact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divalphaUafluact = divalphaUafluact2;
}



IOobject divalphaUaHeader
    (
        "divalphaUa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divalphaUaHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField divalphaUa2
(
    IOobject
    (
        "divalphaUa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divalphaUa = divalphaUa2;
}


IOobject ddtalphaHeader
    (
        "ddtalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (ddtalphaHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField ddtalpha2
(
    IOobject
    (
        "ddtalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
ddtalpha = ddtalpha2;
}


IOobject ddtalphaMeanHeader
    (
        "ddtalphaMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (ddtalphaMeanHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField ddtalphaMean2
(
    IOobject
    (
        "ddtalphaMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
ddtalphaMean = ddtalphaMean2;
}


IOobject ddtalphaPrimeFHeader
    (
        "ddtalphaPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (ddtalphaPrimeFHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField ddtalphaPrimeF2
(
    IOobject
    (
        "ddtalphaPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
ddtalphaPrimeF = ddtalphaPrimeF2;
}

IOobject alphaPrimeFHeader
    (
        "alphaPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (alphaPrimeFHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField alphaPrimeF2
(
    IOobject
    (
        "alphaPrimeF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
alphaPrimeF = alphaPrimeF2;
}


IOobject divfluxUaalphaHeader
    (
        "divfluxUaalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divfluxUaalphaHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField divfluxUaalpha2
(
    IOobject
    (
        "divfluxUaalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divfluxUaalpha = divfluxUaalpha2;
}

/*
IOobject divfluxUaalphaMeanHeader
    (
        "divfluxUaalphaMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (divfluxUaalphaMeanHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField divfluxUaalphaMean2
(
    IOobject
    (
        "divfluxUaalphaMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
divfluxUaalphaMean = divfluxUaalphaMean2;
}



IOobject fluxUaalphaHeader
    (
        "fluxUaalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (fluxUaalphaHeader.typeHeaderOk<volTensorField>(true))
    {

volScalarField fluxUaalpha2
(
    IOobject
    (
        "fluxUaalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
fluxUaalpha = fluxUaalpha2;
}
*/

IOobject alphaUaMeandivHeader
    (
        "alphaUaMeandiv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (alphaUaMeandivHeader.typeHeaderOk<volTensorField>(true))
    {

volVectorField alphaUaMeandiv2
(
    IOobject
    (
        "alphaUaMeandiv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
alphaUaMeandiv = alphaUaMeandiv2;
}

/*
IOobject alphaUafluactHeader
    (
        "alphaUafluact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (alphaUafluactHeader.typeHeaderOk<volTensorField>(true))
    {

volVectorField alphaUafluact2
(
    IOobject
    (
        "alphaUafluact",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
alphaUafluact = alphaUafluact2;
}

*/
IOobject alphaUaHeader
    (
        "alphaUa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );
    if (alphaUaHeader.typeHeaderOk<volTensorField>(true))
    {

volVectorField alphaUa2
(
    IOobject
    (
        "alphaUa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
alphaUa = alphaUa2;
}


			
}


	

}

//Info<<"Manohar run Time = \t"<<timeaverage<<"\t delta t = \t"<<runTime.value()<<endl;

//runtime += runTime.deltaTValue();
