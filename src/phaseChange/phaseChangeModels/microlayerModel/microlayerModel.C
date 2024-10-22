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

#include "microlayerModel.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(microlayerModel, 0);
    defineRunTimeSelectionTable(microlayerModel, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::microlayerModel::microlayerModel
(
    const word& type,
    const phaseModel& phase1,
    const phaseModel& phase2,
    const solidThermo& solid,
    const volScalarField& p,
    singleComponentSatProp& satModel,
    reconstructionSchemes& surf,
    const dictionary& dict
)
:
    dictionary(dict),
    microlayerModelCoeffs_(optionalSubDict(type + "Coeffs")),
    phase1_(phase1),
    phase2_(phase2),
    solid_(solid),
    p_(p),
    satModel_(satModel),
    surf_(surf),
    dml_
    (
        IOobject
        (
            "dml",
            phase1_.mesh().time().timeName(),
            phase1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase1_.mesh(),
        dimLength
    )
    
{}

// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //


const Foam::dictionary&
Foam::microlayerModel::modelDict() const
{
    return microlayerModelCoeffs_;
}

Foam::dictionary&
Foam::microlayerModel::modelDict()
{
    return microlayerModelCoeffs_;
}

Foam::volScalarField& Foam::microlayerModel::dml()
{
    return dml_;
}

Foam::tmp<Foam::volScalarField> Foam::microlayerModel::dml() const 
{
    return dml_;
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
