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
Info<<"viscDiss for time = \t"<<sum(viscDissMean_b)<<endl;
Info<<"viscDissJ for time = \t"<<sum(viscDissMeanJ_b)<<endl;
//Info<<"alpha_aMean = \t"<<sum(alpha_aMean)<<endl;

betaUbMean  = (runtime* betaUbMean  + beta* Ub* deltaT) / (timeaverage);

UbMean = (runtime* UbMean + Ub*deltaT) / (timeaverage);

/* alpha Mean */
alpha_aMean = (runtime* alpha_aMean + alpha* deltaT) / (timeaverage);
beta_bMean  = (runtime* beta_bMean  + beta* deltaT) / (timeaverage);
/* *****************  Pressure */

pMean_b = (runtime* pMean_b + p_rbgh* deltaT ) / (timeaverage);

/* Favre avrage fields */

UbMeanF = betaUbMean  / (1 - alpha_aMean);

/* Favre fluactuations in velocities and Pressure  */

UbPrimeF = Ub - UbMeanF;
pPrimeF_b  = p_rbgh - pMean_b; 

/* *******************  for UbPrimeMean2F or UaPrimeMean2F     ***********  */

betaSqrUbMean  = (runtime* betaSqrUbMean +  beta* Ub* Ub* deltaT ) / (timeaverage);

UbPrime2MeanF = (betaSqrUbMean / (1 - alpha_aMean)) - UbMeanF* UbMeanF;

/* *************** Q and R factor ************ */

volScalarField UbPrimeFx = UbPrimeF.component(vector::X);
volScalarField UbPrimeFy = UbPrimeF.component(vector::Y);

UbsquarWbMean = (runtime* UbsquarWbMean + (UbPrimeFx*UbPrimeFx*UbPrimeFy)* deltaT) / (timeaverage);
WbcubeMean   = (runtime* WbcubeMean   + (UbPrimeFy*UbPrimeFy*UbPrimeFy)* deltaT) / (timeaverage);

RS2_b  = min(UbPrimeFx* UbPrimeFy* neg(UbPrimeFx), zeroo);
RS4_b  = min(UbPrimeFx* UbPrimeFy* pos(UbPrimeFx), zeroo);


RS2Mean_b = (runtime* RS2Mean_b + RS2_b* deltaT) / (timeaverage);
RS4Mean_b = (runtime* RS4Mean_b + RS4_b* deltaT) / (timeaverage);

/* To write some fundamental variables for TKE budget */


if (runTime.outputTime())
{

	alpha_aMean.write();
	UbMeanF.write();
	RS2_b.write();
	RS4_b.write();
	UbPrime2MeanF.write();
}

/* ******************* Mean Sij or Tij *********** */

//Info<<"CalcFavre beta = \t"<<beta<<"\t viscosity = \t"<<phaseb->nu()<<endl;

strainTensor_b        = symm(fvc::grad(Ub));
strainTensorM_b = symm(fvc::grad(UbMeanF));
fluctstrainTensor_b  = symm(fvc::grad(UbPrimeF));

SijbSijb = fluctstrainTensor_b && fluctstrainTensor_b;

SijSij_b   = strainTensorM_b &&  strainTensorM_b;

viscDissMeanI_b = 2* beta* phaseb->nu()* (strainTensor_b && fvc::grad(UbPrimeF));
SGSDissMeanI_b  = 2* beta* turbulenceb->nut()* (strainTensor_b && fvc::grad(UbPrimeF));


/*           ************** Manohar Formula ************* */
if (favreAveraging_fluid)
{


	Info<< "Computing TKE budget for LES "<<endl;

	TKEMean_b = 0.5* tr(UbPrime2MeanF);

	TKEMeanProd_b =  - UbPrime2MeanF && fvc::grad(UbMeanF); 

	turbDiffMean_b = (runtime* turbDiffMean_b - 0.5* beta* (UbPrimeF* magSqr(UbPrimeF))* deltaT ) / (timeaverage);
	
	turbDiffMeanDavid = (runtime* turbDiffMeanDavid - 0.5* beta* (UbsquarWbMean + WbcubeMean)* deltaT ) / (timeaverage);

	pressDiffMean_b = (runtime* pressDiffMean_b - (beta/phaseb->rho())* (UbPrimeF * p_rbgh)* deltaT ) / (timeaverage); 

	SGSDiffMean_b   = (runtime* SGSDiffMean_b - 2* beta* turbulenceb->nut()* (UbPrimeF & strainTensor_b )* deltaT ) / (timeaverage);
	
	viscDiffMean_b   = (runtime* viscDiffMean_b - 2* beta* phaseb->nu()* (UbPrimeF & strainTensor_b )* deltaT ) / (timeaverage);

	viscDissMean_b   = (runtime* viscDissMean_b - (2* beta* phaseb->nu()* (strainTensor_b && (fvc::grad( UbPrimeF))))* deltaT) / (timeaverage);
//	viscDissMeanM  = (runtime* viscDissMean - (2* beta* phaseb->nu()* SijFSijF)* deltaT) / (timeaverage);

	SGSDissMean_b   = (runtime* SGSDissMean_b  + (2* beta* turbulenceb->nut()* (strainTensor_b && (fvc::grad(UbPrimeF)))* deltaT ) )/ (timeaverage);
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

	TKEMean_b.write();
	TKEMeanProd_b.write();
	viscDissMeanI_b.write();
	SGSDissMeanI_b.write();
	turbDiffMean_b.write();
	viscDissMean_b.write();
	SGSDissMean_b.write();
	dragConstbMean.write();
	dragVegMean.write();
	
	turbDiffusionMean_b = fvc::div(turbDiffMean_b);
	turbDiffusionMeanDavid = fvc::grad(turbDiffMeanDavid);
	pressDiffusionMean_b = fvc::div(pressDiffMean_b);
	SGSDiffusionMean_b = fvc::div(SGSDiffMean_b);
	viscDiffusionMean_b = fvc::div(viscDiffMean_b);

        turbDiffusionMean_b.write();
        pressDiffusionMean_b.write();
        SGSDiffusionMean_b.write();
	viscDiffusionMean_b.write();
}

/*     *********** Julien formulation  ***************  */


if (favreAveraging_fluid)
{


	Info<< "Computing TKE budget for LES Julien Formula "<<endl;

	pressDiffMeanJ_b = (runtime* pressDiffMeanJ_b - (beta/phaseb->rho())* (UbPrimeF * pPrimeF_b)* deltaT ) / (timeaverage); 

	SGSDiffMeanJ_b   = (runtime* SGSDiffMeanJ_b - 2* beta* turbulenceb->nut()* (Ub & strainTensor_b )* deltaT ) / (timeaverage);
	
	viscDiffusionMeanJ_b   = (runtime* viscDiffusionMeanJ_b +  beta* phaseb->nu()* (fvc::laplacian(TKEMean_b) )* deltaT ) / (timeaverage);
	
	viscDissMeanJ_b  = (runtime* viscDissMeanJ_b - (2* beta* phaseb->nu()* SijbSijb)* deltaT) / (timeaverage);

	SGSDissMeanJ_b   = (runtime* SGSDissMeanJ_b  - (2* beta* turbulenceb->nut()* SijbSijb)* deltaT ) / (timeaverage);
}


//Info<<"kinematic viscosity = \t"<<phaseb->nu()<<"\t LES viscosity = \t"<<turbulenceb->nut()<<endl;

/* Writing the TKE budget variables */

if (runTime.outputTime())
{

	viscDissMeanJ_b.write();
	SGSDissMeanJ_b.write();
		
	pressDiffusionMeanJ_b = fvc::div(pressDiffMeanJ_b);
	SGSDiffusionMeanJ_b = fvc::div(SGSDiffMeanJ_b);
	//viscDiffusionMeanJ = fvc::div(viscDiffMeanJ);

        pressDiffusionMeanJ_b.write();
        SGSDiffusionMeanJ_b.write();
	viscDiffusionMeanJ_b.write();
}


runtime += deltaT; 


