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


/*    This terms compute rho U mean */

Info<<"runTime before the addition = \t"<<runtime<<endl;

//runtime = runTime.value();
deltaT = runTime.deltaTValue();

Info<<"timeaverage in Mass budget = \t"<<endl;

timeaverage = runtime  + deltaT;

betaUbMean  = (runtime* betaUbMean  + beta* Ub* deltaT) / (timeaverage);
alphaUaMean = (runtime* alphaUaMean + alpha* Ua* deltaT) / (timeaverage);

UbMean = (runtime* UbMean + Ub*deltaT) / (timeaverage);

/* alpha Mean */
alpha_aMean = (runtime* alpha_aMean + alpha* deltaT) / (timeaverage);
beta_bMean  = (runtime* beta_bMean  + beta* deltaT) / (timeaverage);
/* *****************  Pressure */

/* Favre avrage fields */

UbMeanF = betaUbMean  / (1 - alpha_aMean);

/* Favre fluactuations in velocities and Pressure  */

UbPrimeF = Ub - UbMeanF;
alphaPrimeF = alpha - alpha_aMean;
betaPrimeF = beta - beta_bMean;
pPrimeF  = p_rbgh - pMeanF; 

/* *******************  for UbPrimeMean2F or UaPrimeMean2F     ***********  */

betaSqrUbMean  = (runtime* betaSqrUbMean +  beta* Ub* Ub* deltaT ) / (timeaverage);

UbPrime2MeanF = (betaSqrUbMean / (1 - alpha_aMean)) - UbMeanF* UbMeanF;

/* *************** Q and R factor ************ */

volScalarField UbPrimeFx = UbPrimeF.component(vector::X);
volScalarField UbPrimeFy = UbPrimeF.component(vector::Y);

/* Mass flux computation ****               */

alphaUaMeandiv = alpha_aMean* UaMean;
alphaUafluacdiv = alphaPrimeF* UaPrimeF;
alphaUa         = alpha* Ua;

divalphaUaMean   = fvc::div(alphaUaMeandiv);
divalphaUafluact = fvc::div(alphaUafluacdiv);
divalphaUa        = fvc::div(alphaUa);


ddtalpha       = fvc::ddt(alpha);
ddtalphaMean   = fvc::ddt(alpha_aMean);
ddtalphaPrimeF = fvc::ddt(alphaPrimeF);

fluxUaalpha            = fvc::flux(alphaUa);
//fluxUaMeanalphaMean    = fvc::flux(alpha_aMean, UaMean);
//fluxUaPrimealphaPrime  = fvm::flux(alphaPrimeF, UaPrimeF);

//divfluxUaalpha = fvc::div(fluxUaalpha);

//fluxUaMeanalphaMean    = fvc::flux(alpha_aMean, UaMean);

/* To write some fundamental variables for TKE budget */

if (runTime.outputTime())
{
	ddtalpha.write();
	ddtalphaMean.write();
	ddtalphaPrimeF.write();
	fluxUaalpha.write();
	divalphaUaMean.write();
	divalphaUafluact.write();
	divalphaUa.write();
	alphaUaMeandiv.write();
	alphaUafluacdiv.write();
	alphaUa.write();
}


runtime += deltaT; 

