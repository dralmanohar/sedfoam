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

alphaF = alpha - alpha_aMean;

Mean_alphaUa     = (runtime* Mean_alphaUa  + alpha* Ua* deltaT) / (timeaverage);
Mean_UaFalphaF   = (runtime* Mean_UaFalphaF  + alphaF* UaPrimeF* deltaT) / (timeaverage);
alphaFMean       = (runtime* alphaFMean  + alphaF* deltaT) / (timeaverage);
MeanUaF          = (runtime* MeanUaF  + UaPrimeF* deltaT) / (timeaverage);


if (runTime.outputTime())
{
	alphaF.write();
	alphaFMean.write();
	Mean_alphaUa.write();
	Mean_UaFalphaF.write();
	MeanUaF.write();
}


//runtime += deltaT; 

