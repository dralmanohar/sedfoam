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

if (favreAveraging_fluid)
{
    #include "calcFavreAverage_fluid.H"
    #include "Random.H"
}

if (favreAveraging_particle)
{
    #include "calcFavreAverage_particle.H"
    #include "Random.H"
}


if (favreAveraging_mass)
{
    #include "calcFavreAverage_mass.H"
    #include "Random.H"
}

if (MomentumBudget)
{       
	Info<<"Momentum Budget in Write IO file "<<endl;
	#include "calcmomentumbudget.H"
}

if (runTime.outputTime())
{
    if (writeTau == 1){
        K.write();
        Kveg.write();
        Taua = fvc::reconstruct(fvc::interpolate(rhoa*(alpha*nuEffa+nuFra))*mesh.magSf()*fvc::snGrad(Ua)) +
            rhoa*(alpha*nuEffa+nuFra)*dev(gradUaT);
        Taub = fvc::reconstruct(fvc::interpolate(rhob*beta*nuEffb)*mesh.magSf()*fvc::snGrad(Ub)) +
            rhob*(scalar(1.0)-alpha)*nuEffb*dev(gradUbT);
        Taua.write();
        Taub.write();
    }


    if (granularRheology.on())
    {
        pa.write();
        granularRheology.muI().write();
        if (granularRheology.dilatancy())
        {
            alphaPlastic.write();
            delta.write();
        }
    }
    if (kineticTheory.on())
    {
        pa.write();
        kineticTheory.Theta().write();
    }
    if (SUSlocal)
    {
        SUS.write();
    }
    pff.write();
    runTime.write();
}
