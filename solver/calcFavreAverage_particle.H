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
//deltaT = runTime.deltaTValue();

Info<<"timeaverage = \t"<<runtime + deltaT<<"run time delta value = \t"<<runTime.deltaTValue()<<endl;

//timeaverage = runtime  + deltaT;

Info<<"timeaverage = \t"<<timeaverage<<endl;
Info<<"deltaT = \t"<<deltaT<<endl;
Info<<"run Time =\t"<<runTime.value()<<endl;
Info<<"run time from the code = \t"<<runtime<<endl;
Info<<"viscDiss for time = \t"<<sum(viscDissMean_b)<<endl;
Info<<"viscDissJ for time = \t"<<sum(viscDissMeanJ_b)<<endl;
//Info<<"alpha_aMean = \t"<<sum(alpha_aMean)<<endl;

alphaUaMean  = (runtime* alphaUaMean  + alpha* Ua* deltaT) / (timeaverage);

UaMean = (runtime* UaMean + Ua*deltaT) / (timeaverage);

/* alpha Mean */
alpha_aMean = (runtime* alpha_aMean + alpha* deltaT) / (timeaverage);
/* *****************  Pressure */

pMean_a = (runtime* pMean_a + pa* deltaT ) / (timeaverage);

/* Favre avrage fields */

UaMeanF = alphaUaMean  / (alpha_aMean + 1e-20);

/* Favre fluactuations in velocities and Pressure  */

UaPrimeF = Ua - UaMeanF;
pPrimeF_a  = pa - pMean_a; 

/* *******************  for UbPrimeMean2F or UaPrimeMean2F     ***********  */

alphaSqrUaMean  = (runtime* alphaSqrUaMean +  alpha* Ua* Ua* deltaT ) / (timeaverage);

UaPrime2MeanF = (alphaSqrUaMean / (alpha_aMean + 1e-20)) - UaMeanF* UaMeanF;

/* *************** Q and R factor ************ */

volScalarField UaPrimeFx = UaPrimeF.component(vector::X);
volScalarField UaPrimeFy = UaPrimeF.component(vector::Y);

RS2_a  = min(UaPrimeFx* UaPrimeFy* neg(UaPrimeFx), zeroo);
RS4_a  = min(UaPrimeFx* UaPrimeFy* pos(UaPrimeFx), zeroo);


RS2Mean_a = (runtime* RS2Mean_a + RS2_a* deltaT) / (timeaverage);
RS4Mean_a = (runtime* RS4Mean_a + RS4_a* deltaT) / (timeaverage);

/* To write some fundamental variables for TKE budget */


if (runTime.outputTime())
{

	alpha_aMean.write();
	UaMeanF.write();
	RS2_a.write();
	RS4_a.write();
	UaPrime2MeanF.write();
	UaMean.write();
	UaMeanF.write();
}

/* ******************* Mean Sij or Tij *********** */

//Info<<"CalcFavre beta = \t"<<beta<<"\t viscosity = \t"<<phaseb->nu()<<endl;

strainTensor_a        = symm(fvc::grad(Ua));
strainTensorM_a = symm(fvc::grad(UaMeanF));
fluctstrainTensor_a  = symm(fvc::grad(UaPrimeF));

SijaSija = fluctstrainTensor_a && fluctstrainTensor_a;

SijSij_a   = strainTensorM_a &&  strainTensorM_a;

viscDissMeanI_a = 2* alpha* phasea->nu()* (strainTensor_a && fvc::grad(UaPrimeF));
SGSDissMeanI_a  = 2* alpha* turbulencea->nut()* (strainTensor_a && fvc::grad(UaPrimeF));


/*           ************** Manohar Formula ************* */
if (favreAveraging_particle)
{


	Info<< "Computing TKE budget for LES "<<endl;

	TKEMean_a = 0.5* tr(UaPrime2MeanF);

	TKEMeanProd_a =  - UaPrime2MeanF && fvc::grad(UaMeanF); 

	turbDiffMean_a = (runtime* turbDiffMean_a - 0.5* alpha* (UaPrimeF* magSqr(UaPrimeF))* deltaT ) / (timeaverage);
	
	pressDiffMean_a = (runtime* pressDiffMean_a - (alpha/phasea->rho())* (UaPrimeF * p_rbgh)* deltaT ) / (timeaverage); 
	
	pressDiffMeankin_a = (runtime* pressDiffMeankin_a - (1/phasea->rho())* (UaPrimeF * pa)* deltaT ) / (timeaverage); 
	
	pressDiffMeancon_a = (runtime* pressDiffMeancon_a - (1/phasea->rho())* (UaPrimeF * pff)* deltaT ) / (timeaverage); 
	
	viscDiffMean_a   = (runtime* viscDiffMean_a - 2* alpha* (phasea->nu()* (UaPrimeF & strainTensor_a ))* deltaT ) / (timeaverage);
	
	viscDissMean_a  = (runtime* viscDissMean_a - (2* alpha* (phasea->nu() + nua_kin + nuFra)* SijaSija)* deltaT) / (timeaverage);
	
	viscDissMeankin_a  = (runtime* viscDissMeankin_a - (2* alpha* (nua_kin)* SijaSija)* deltaT) / (timeaverage);
	
	viscDissMeancon_a  = (runtime* viscDissMeancon_a - (2* alpha* nuFra* SijaSija)* deltaT) / (timeaverage);
	
	SGSDissMeanJ_a   = (runtime* SGSDissMeanJ_a  - (2* alpha* turbulencea->nut()* SijaSija)* deltaT ) / (timeaverage);

	SGSDiffMean_a   = (runtime* SGSDiffMean_a - 2* alpha* turbulencea->nut()* (UaPrimeF & strainTensor_a )* deltaT ) / (timeaverage);

	SGSDissMean_a   = (runtime* SGSDissMean_a  + (2* alpha* turbulencea->nut()* (strainTensor_a && (fvc::grad(UaPrimeF)))* deltaT ) )/ (timeaverage);
	
	
//	pressDiffMeanJ_a = (runtime* pressDiffMeanJ_a - (alpha/phasea->rho())* (UaPrimeF * pPrimeF_a)* deltaT ) / (timeaverage); 

//	SGSDiffMeanJ_a   = (runtime* SGSDiffMeanJ_a - 2* alpha* turbulencea->nut()* (Ua & strainTensor_a )* deltaT ) / (timeaverage);
	
//	viscDiffusionMeanJ_b   = (runtime* viscDiffusionMeanJ_b +  beta* phaseb->nu()* (fvc::laplacian(TKEMean_b) )* deltaT ) / (timeaverage);
	
	viscDissMeanJ_a  = (runtime* viscDissMeanJ_a - (2* alpha* phasea->nu()* SijaSija)* deltaT) / (timeaverage);

//	SGSDissMeanJ_b   = (runtime* SGSDissMeanJ_b  - (2* beta* turbulenceb->nut()* SijbSijb)* deltaT ) / (timeaverage);

	
}

//Info<<"testVegMean = \t"<<testVegMean.value()<<endl;
//Info<<"dragVegMean = \t"<<dragVegMean<<endl;

/* Writing the TKE budget variables */

//Info<<"Kveg in favre average in tke budget = \t"<<dragVegMean<<endl;


if (runTime.outputTime())
{

	TKEMean_a.write();
	TKEMeanProd_a.write();
	viscDissMeanI_a.write();
	SGSDissMeanI_a.write();
	turbDiffMean_a.write();
	viscDissMean_a.write();
	viscDissMeanJ_a.write();
	viscDissMeankin_a.write();
	viscDissMeancon_a.write();
	SGSDissMean_a.write();
	
	turbDiffusionMean_a = fvc::div(turbDiffMean_a);
	pressDiffusionMean_a = fvc::div(pressDiffMean_a);
	pressDiffusionMeankin_a = fvc::div(pressDiffMeankin_a);
	pressDiffusionMeancon_a = fvc::div(pressDiffMeancon_a);
	SGSDiffusionMean_a = fvc::div(SGSDiffMean_a);
	viscDiffusionMean_a = fvc::div(viscDiffMean_a);

        turbDiffusionMean_a.write();
        pressDiffusionMean_a.write();
        SGSDiffusionMean_a.write();
	viscDiffusionMean_a.write();
}


 
/*     *********** Julien formulation  ***************  */

/*
if (TKEBudgetLES)
{


	Info<< "Computing TKE budget for LES Julien Formula "<<endl;

	pressDiffMeanJ_b = (runtime* pressDiffMeanJ_b - (beta/phaseb->rho())* (UbPrimeF * pPrimeF_b)* deltaT ) / (timeaverage); 

	SGSDiffMeanJ_b   = (runtime* SGSDiffMeanJ_b - 2* beta* turbulenceb->nut()* (Ub & strainTensor_b )* deltaT ) / (timeaverage);
	
	viscDiffusionMeanJ_b   = (runtime* viscDiffusionMeanJ_b +  beta* phaseb->nu()* (fvc::laplacian(TKEMean_b) )* deltaT ) / (timeaverage);
	
	viscDissMeanJ_b  = (runtime* viscDissMeanJ_b - (2* beta* phaseb->nu()* SijbSijb)* deltaT) / (timeaverage);

	SGSDissMeanJ_b   = (runtime* SGSDissMeanJ_b  - (2* beta* turbulenceb->nut()* SijbSijb)* deltaT ) / (timeaverage);
}


//Info<<"kinematic viscosity = \t"<<phaseb->nu()<<"\t LES viscosity = \t"<<turbulenceb->nut()<<endl;
*/
/* Writing the TKE budget variables */
/*
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
*/

//runtime += deltaT; 


