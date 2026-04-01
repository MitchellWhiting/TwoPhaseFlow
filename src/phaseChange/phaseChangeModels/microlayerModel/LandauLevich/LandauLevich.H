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

Class
    Foam::LandauLevich

Description
    interface heat resistance model based on the LandauLevich equation

SourceFiles
    LandauLevich.C

\*---------------------------------------------------------------------------*/

#ifndef LandauLevich_H
#define LandauLevich_H

#include "typeInfo.H"
#include "microlayerModel.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "scalar.H"

#include "autoPtr.H"

#include "fvScalarMatrix.H"
#include "fvm.H"

#include "reconstructionSchemes.H"
#include "singleComponentSatProp.H"

#include "mixedFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class LandauLevich Declaration
\*---------------------------------------------------------------------------*/

class LandauLevich : public microlayerModel
{

private:

    // volScalarField dml_;

    scalar evapCoeff_;
    dimensionedScalar Rgas_;
    scalar sigma_;
    string fluidPatch_;
    string solidPatch_;
    scalar Jakob_;
    scalar prevT_;
    scalar prevR_;
    scalar prevRbase_;
    scalar prevdRdt_;
    label prevCell_;
    vector origin_;


    // Private Member Functions

        //- Disallow copy construct
        LandauLevich(const LandauLevich&);

        //- Disallow default bitwise assignment
        void operator=(const LandauLevich&);

public:

    //- Runtime type information
    TypeName("LandauLevich");

    // Constructors

        //- Construct from components
        LandauLevich
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const solidThermo& solid,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );

    //- Destructor
    virtual ~LandauLevich()
    {}

    // Member Functions
        virtual tmp<fvScalarMatrix> hSourceML();

        virtual void initialiseML();

        virtual void updateML();
        
        virtual tmp<volScalarField> energySourceML();
        
        virtual tmp<volScalarField> massSourceML(volScalarField& rhoSource);

        virtual tmp<volScalarField> alphaSourceML(volScalarField& rhoSource);
 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
