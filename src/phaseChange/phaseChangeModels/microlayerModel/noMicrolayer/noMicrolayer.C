/*---------------------------------------------------------------------------*\
            Copyright (c) 2017-2019, German Aerospace Center (DLR)
-------------------------------------------------------------------------------
License
    This file is part of the VoFLibrary source code library, which is an
	unofficial extension to OpenFOAM.
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


#include "noMicrolayer.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{ 
    defineTypeNameAndDebug(noMicrolayer, 0);
    addToRunTimeSelectionTable(microlayerModel,noMicrolayer, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noMicrolayer::noMicrolayer
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const solidThermo& solid,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    reconstructionSchemes& surf,
    const dictionary& dict
)
:
    microlayerModel
    (
        typeName,
        phase1,
        phase2,
        solid,
        p,
        satModel,
        surf,
        dict
    )
{
}
// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * *  //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //

Foam::tmp<Foam::fvScalarMatrix> Foam::noMicrolayer::hSourceML()
{
    const volScalarField& hsolid = solid_.he();  

    volScalarField sourceCoeff
    (
        IOobject
        (
            "sourceCoeff",
            hsolid.mesh().time().timeName(),
            hsolid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        hsolid.mesh(),
        dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0),
        "zeroGradient"
    );


    tmp<fvScalarMatrix> hSource(fvm::Sp(sourceCoeff, hsolid));

    return hSource;

}


void Foam::noMicrolayer::initialiseML()
{   
   
}

void Foam::noMicrolayer::updateML()
{       
    
}

Foam::tmp<Foam::volScalarField> Foam::noMicrolayer::energySourceML()
{

    const volScalarField& TSat = satModel_.TSat(); 
    const volScalarField& k1 = phase1_.kappa(); 
    const dimensionedScalar kmax ("kFluidMax",k1.dimensions(), gMax(k1.internalField()));
    const dimensionedScalar dummyLength ("dummyLength",dimLength, 1.0);

    tmp<volScalarField> sourceCoeff(TSat*0.0/(pow(dummyLength,2)/kmax));    
    volScalarField& sourceCoeffRef = sourceCoeff.ref();    

    sourceCoeffRef.correctBoundaryConditions();

    return sourceCoeff;

}



Foam::tmp<Foam::volScalarField> 
Foam::noMicrolayer::massSourceML( volScalarField& rhoSource)
{
    
    tmp<volScalarField> massSourceML(rhoSource * 0.0);
    rhoSource *= 0.0;
    return massSourceML;
}


Foam::tmp<Foam::volScalarField>
Foam::noMicrolayer::alphaSourceML( volScalarField& rhoSource)
{
    dimensionedScalar scale("scale", dimless/dimDensity,0.0);
    tmp<volScalarField> alphaSourceML(rhoSource * scale);

   return alphaSourceML;
}
