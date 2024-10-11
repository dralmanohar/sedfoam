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


Info<<"Manohar in calc Momentum budget files = \t"<<endl;

Info<<"Ufall velocity from momentum calculation manohar  = \t"<<max(UfallField.internalField()).value()<<endl;

tensor Ei(1.,0,0,0,1.,0,0,0,1.);

if (runTime.outputTime())
{
Errosion = alpha* U.component(vector::Y);
Deposition = alpha* UfallField;

Errosion.write();
Deposition.write();


//Taua = Ei; //+ spherSigmaSGSa); // - (nuEffa+nuFra/(alpha+alphaSmall))*gradUaT
 

//Taua = rhoa*(alpha*nuEffa + nuFra)*(fvc::grad(Ua) + gradUaT - ((2./3.)*Ei)*tr(fvc::grad(Ua)))  + lambdaUa*tr(fvc::grad(Ua)); 

//Taua =  fvc::reconstruct(fvc::interpolate(rhoa*(alpha*nuEffa+nuFra)*fvc::grad(Ua))*mesh.magSf())+ rhoa*(alpha*nuEffa+nuFra)*dev(gradUaT);

//Taub = rhob*(scalar(1.0) - alpha)*nuEffb*(fvc::grad(Ub) + gradUbT - ((2./3.)*Ei)*tr(fvc::grad(Ub)));  

//fvc::reconstruct(fvc::interpolate(rhob*(scalar(1.0)-alpha)*nuEffb)*mesh.magSf()*fvc::snGrad(Ub)); // + fvc::reconstruct(fvc::interpolate(rhob*(scalar(1.0) - alpha)*nuEffb)*mesh.magSf()*fvc::snGrad(Ub)().T()); //rhob*(scalar(1.0)-alpha)*nuEffb*dev(gradUbT);


//Info<<"nuFra = \t"<<nuFra<<endl;

//Taua.write();
//Taub.write();
//

Taua = fvc::reconstruct(fvc::interpolate(rhoa*(alpha*nuEffa+nuFra)) *mesh.magSf()*fvc::snGrad(Ua)) + rhoa*(alpha*nuEffa+nuFra)*dev(gradUaT);

Taub = fvc::reconstruct(fvc::interpolate(rhob*(scalar(1.0)-alpha)*nuEffb) *mesh.magSf()*fvc::snGrad(Ub)) + rhob*(scalar(1.0)-alpha)*nuEffb*dev(gradUbT);

Taua.write();
Taub.write();

volTensorField Rca(
       "Rca",
      ((2.0/3.0)*I)*((nuEffa+nuFra/(alpha+alphaSmall))*tr(gradUaT) + spherSigmaSGSa) - (nuEffa+nuFra/(alpha+alphaSmall))*gradUaT
      );

volTensorField Rca2(
      "Rca2",
       ((2.0/3.0)*I)*((alpha*nuEffa+nuFra)*tr(gradUaT) + alpha*spherSigmaSGSa) - (alpha*nuEffa+nuFra)*gradUaT
       );
        
volVectorField Ra(
       "Ra",
        -nuEffa*fvc::grad(alpha)/(alpha+alphaSmall)
       );
        
surfaceScalarField phiRa(
      -fvc::interpolate(nuEffa)*mesh.magSf()*(fvc::snGrad(alpha))
       /fvc::interpolate(alpha+ alphaSmall)
       );

Rca.write();
Rca2.write();
Ra.write();


pMech = p;
pA = pa;
pFF = pff;
SolidPressureContr = -fvc::grad(pa+pff);
ViscStressContr    =  rhoa*alpha*fvc::laplacian(nuEffa, Ua);
FrictContr         =  rhoa*fvc::laplacian(nuFra, Ua);
phiRaContr         = -rhoa*alpha*fvc::div(phiRa, Ua, "div(phiRa,Ua)");
phiRa2Contr        =  rhoa*alpha*fvc::Sp(fvc::div(phiRa), Ua);
divRcaContr        = -rhoa*alpha*fvc::div(Rca);
RcaContr           = -rhoa*alpha/(alpha+alphaSmall)*(fvc::grad(alpha) & Rca);
BulkContr          =  fvc::grad(lambdaUa*tr(gradUaT));
divTauaContr       =  ViscStressContr + FrictContr + phiRaContr + phiRa2Contr + divRcaContr + RcaContr + BulkContr;
        
gravityContr       =  alpha/(alpha+alphaSmall)*(rhoa-rhob)*g;
DragContr          =  alpha*beta*K*(Ub-Ua);
ExtContr           =  rhoa*alpha*gradPOSC/(tilt*rhob +(1.0-tilt)*rhoa);
SuspContr          = -SUS*K*beta*turbulenceb->nut()*((SUS_I*iso-SUS_A*aniso) & fvc::grad(alpha));

pMech.write();
pA.write();
pFF.write();
SolidPressureContr.write();
ViscStressContr.write();
FrictContr.write();
phiRaContr.write();
phiRa2Contr.write();
divRcaContr.write();
RcaContr.write();
BulkContr.write();
divTauaContr.write();
gravityContr.write();
DragContr.write();
ExtContr.write();
SuspContr.write();


/* mixture properties for diagnostic*/

rho_mixd  = rhoa* alpha* Ua + rhob* beta* Ub;
nonl_mixd  = rhoa* alpha* Ua * Ua  + rhob* beta* Ub * Ub;
bouy_mixd = (rhoa* alpha  + rhob* beta)* g;


rho_mixd.write();
nonl_mixd.write();
bouy_mixd.write();
rho.write();

}
//divTauaContr.write();
//gravityContr.write();
//DragContr.write();
//ExtContr.write();
//SuspContr.write();

