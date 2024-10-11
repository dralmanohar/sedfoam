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

Info<<"timeaverage = \t"<<runtime + deltaT<<"run time delta value = \t"<<runTime.deltaTValue()<<endl;

timeaverage = runtime  + deltaT;

Info<<"timeaverage = \t"<<timeaverage<<endl;
Info<<"deltaT = \t"<<deltaT<<endl;
Info<<"run Time =\t"<<runTime.value()<<endl;
Info<<"run time from the code = \t"<<runtime<<endl;
Info<<"viscDiss for time = \t"<<sum(viscDissMean)<<endl;
Info<<"viscDissJ for time = \t"<<sum(viscDissMeanJ)<<endl;
//Info<<"alpha_aMean = \t"<<sum(alpha_aMean)<<endl;

betaUbMean  = (runtime* betaUbMean  + beta* Ub* deltaT) / (timeaverage);
alphaUaMean = (runtime* alphaUaMean + alpha* Ua* deltaT) / (timeaverage);

UaMean = (runtime* UaMean + Ua*deltaT) / (timeaverage);
UbMean = (runtime* UbMean + Ub*deltaT) / (timeaverage);

/* alpha Mean */
alpha_aMean = (runtime* alpha_aMean + alpha* deltaT) / (timeaverage);
beta_bMean  = (runtime* beta_bMean  + beta* deltaT) / (timeaverage);
/* *****************  Pressure */

pMeanF = (runtime* pMeanF + p_rbgh* deltaT ) / (timeaverage);

/* Favre avrage fields */

UbMeanF = betaUbMean  / (1 - alpha_aMean);
UaMeanF = alphaUaMean / (alpha_aMean + 1e-20);

/* Favre fluactuations in velocities and Pressure  */

UbPrimeF = Ub - UbMeanF;
UaPrimeF = Ua - UaMeanF;
alphaPrimeF = alpha - alpha_aMean;
betaPrimeF = beta - beta_bMean;
pPrimeF  = p_rbgh - pMeanF; 




/* *******************  for UbPrimeMean2F or UaPrimeMean2F     ***********  */

betaSqrUbMean  = (runtime* betaSqrUbMean +  beta* Ub* Ub* deltaT ) / (timeaverage);
alphaSqrUaMean = (runtime* alphaSqrUaMean + alpha* Ua* Ua* deltaT) / (timeaverage);


UbPrime2MeanF = (betaSqrUbMean / (1 - alpha_aMean)) - UbMeanF* UbMeanF;
UaPrime2MeanF = (alphaSqrUaMean / (alpha_aMean + 1e-20)) - UaMeanF* UaMeanF;


/* *************** Q and R factor ************ */

volScalarField UbPrimeFx = UbPrimeF.component(vector::X);
volScalarField UbPrimeFy = UbPrimeF.component(vector::Y);

UsquarWMean = (runtime* UsquarWMean + (UbPrimeFx*UbPrimeFx*UbPrimeFy)* deltaT) / (timeaverage);
WcubeMean   = (runtime* WcubeMean   + (UbPrimeFy*UbPrimeFy*UbPrimeFy)* deltaT) / (timeaverage);

RS2  = min(UbPrimeFx* UbPrimeFy* neg(UbPrimeFx), zeroo);
RS4  = min(UbPrimeFx* UbPrimeFy* pos(UbPrimeFx), zeroo);


RS2Mean = (runtime* RS2Mean + RS2* deltaT) / (timeaverage);
RS4Mean = (runtime* RS4Mean + RS4* deltaT) / (timeaverage);

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

betaUbMeandiv = beta_bMean* UbMean;
betaUbfluacdiv = betaPrimeF* UbPrimeF;
betaUb         = beta* Ub;

divbetaUbMean    = fvc::div(betaUbMeandiv);
divbetaUbfluact  = fvc::div(betaUbfluacdiv);
divbetaUb        = fvc::div(betaUb);


ddtbeta       = fvc::ddt(beta);
ddtbetaMean   = fvc::ddt(beta_bMean);
ddtbetaPrimeF = fvc::ddt(betaPrimeF);

fluxUbbeta            = fvc::flux(betaUb);
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

if (runTime.outputTime())
{

	alpha_aMean.write();
	UbMeanF.write();
	UaMeanF.write();
	RS2.write();
	RS4.write();
	UbPrime2MeanF.write();
	UaPrime2MeanF.write();
}

/* ******************* Mean Sij or Tij *********** */

//Info<<"CalcFavre beta = \t"<<beta<<"\t viscosity = \t"<<phaseb->nu()<<endl;

strainTensor        = symm(fvc::grad(Ub));
strainTensorM = symm(fvc::grad(UbMeanF));
fluctstrainTensorF  = symm(fvc::grad(UbPrimeF));

SijFSijF = fluctstrainTensorF && fluctstrainTensorF;

SijSijM   = strainTensorM &&  strainTensorM;

viscDissMeanI = 2* beta* phaseb->nu()* (strainTensor && fvc::grad(UbPrimeF));
SGSDissMeanI  = 2* beta* turbulenceb->nut()* (strainTensor && fvc::grad(UbPrimeF));


/*           ************** Manohar Formula ************* */
if (TKEBudgetLES)
{


	Info<< "Computing TKE budget for LES "<<endl;

	TKEMean = 0.5* tr(UbPrime2MeanF);

	TKEMeanProd =  - UbPrime2MeanF && fvc::grad(UbMeanF); 

	turbDiffMean = (runtime* turbDiffMean - 0.5* beta* (UbPrimeF* magSqr(UbPrimeF))* deltaT ) / (timeaverage);
	
	turbDiffMeanDavid = (runtime* turbDiffMeanDavid - 0.5* beta* (UsquarWMean + WcubeMean)* deltaT ) / (timeaverage);

	pressDiffMean = (runtime* pressDiffMean - (beta/phaseb->rho())* (UbPrimeF * p_rbgh)* deltaT ) / (timeaverage); 

	SGSDiffMean   = (runtime* SGSDiffMean - 2* beta* turbulenceb->nut()* (UbPrimeF & strainTensor )* deltaT ) / (timeaverage);
	
	viscDiffMean   = (runtime* viscDiffMean - 2* beta* phaseb->nu()* (UbPrimeF & strainTensor )* deltaT ) / (timeaverage);

	viscDissMean   = (runtime* viscDissMean - (2* beta* phaseb->nu()* (strainTensor && (fvc::grad( UbPrimeF))))* deltaT) / (timeaverage);
//	viscDissMeanM  = (runtime* viscDissMean - (2* beta* phaseb->nu()* SijFSijF)* deltaT) / (timeaverage);

	SGSDissMean   = (runtime* SGSDissMean  + (2* beta* turbulenceb->nut()* (strainTensor && (fvc::grad(UbPrimeF)))* deltaT ) )/ (timeaverage);
//	SGSDissMeanM   = (runtime* SGSDissMean  - (2* beta* turbulenceb->nut()* SijFSijF)* deltaT ) / (timeaverage);

	dragConstbMean = (runtime* dragConstbMean - (UbPrimeF & (alpha* beta* draga->K(mag(Ub-Ua))* (Ub-Ua))/rhob)* deltaT ) / (timeaverage);
	
	//dragVegMean = (runtime* dragVegMean - (UbPrimeF & (alpha* gamma* Kveg* (Ub - Uc))/rhob)* deltaT ) / (timeaverage);
	
	dragVegMean = (runtime* dragVegMean - (UbPrimeF & (beta* gamma* Kveg* (Ub - Uc))/rhob)* deltaT ) / (timeaverage);
}

//Info<<"testVegMean = \t"<<testVegMean.value()<<endl;
//Info<<"dragVegMean = \t"<<dragVegMean<<endl;

/* Writing the TKE budget variables */

//Info<<"Kveg in favre average in tke budget = \t"<<dragVegMean<<endl;

if (runTime.outputTime())
{

	TKEMean.write();
	TKEMeanProd.write();
	viscDissMeanI.write();
	SGSDissMeanI.write();
	turbDiffMean.write();
	viscDissMean.write();
	SGSDissMean.write();
	dragConstbMean.write();
	dragVegMean.write();
	
	turbDiffusionMean = fvc::div(turbDiffMean);
	turbDiffusionMeanDavid = fvc::grad(turbDiffMeanDavid);
	pressDiffusionMean = fvc::div(pressDiffMean);
	SGSDiffusionMean = fvc::div(SGSDiffMean);
	viscDiffusionMean = fvc::div(viscDiffMean);

        turbDiffusionMean.write();
        pressDiffusionMean.write();
        SGSDiffusionMean.write();
	viscDiffusionMean.write();
}

/*     *********** Julien formulation  ***************  */


if (TKEBudgetLES)
{


	Info<< "Computing TKE budget for LES Julien Formula "<<endl;

	pressDiffMeanJ = (runtime* pressDiffMeanJ - (beta/phaseb->rho())* (UbPrimeF * pPrimeF)* deltaT ) / (timeaverage); 

	SGSDiffMeanJ   = (runtime* SGSDiffMeanJ - 2* beta* turbulenceb->nut()* (Ub & strainTensor )* deltaT ) / (timeaverage);
	
	viscDiffusionMeanJ   = (runtime* viscDiffusionMeanJ +  beta* phaseb->nu()* (fvc::laplacian(TKEMean) )* deltaT ) / (timeaverage);
	
	viscDissMeanJ  = (runtime* viscDissMeanJ - (2* beta* phaseb->nu()* SijFSijF)* deltaT) / (timeaverage);

	SGSDissMeanJ   = (runtime* SGSDissMeanJ  - (2* beta* turbulenceb->nut()* SijFSijF)* deltaT ) / (timeaverage);
}


//Info<<"kinematic viscosity = \t"<<phaseb->nu()<<"\t LES viscosity = \t"<<turbulenceb->nut()<<endl;

/* Writing the TKE budget variables */

if (runTime.outputTime())
{

	viscDissMeanJ.write();
	SGSDissMeanJ.write();
		
	pressDiffusionMeanJ = fvc::div(pressDiffMeanJ);
	SGSDiffusionMeanJ = fvc::div(SGSDiffMeanJ);
	//viscDiffusionMeanJ = fvc::div(viscDiffMeanJ);

        pressDiffusionMeanJ.write();
        SGSDiffusionMeanJ.write();
	viscDiffusionMeanJ.write();
}


runtime += deltaT; 


