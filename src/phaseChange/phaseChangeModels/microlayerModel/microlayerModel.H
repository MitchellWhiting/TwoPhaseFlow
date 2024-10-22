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
    Foam::microlayerEnergySourceTermModel

Description
    base class of the energy source terms models

SourceFiles
    microlayerEnergySourceTermModel.C

\*---------------------------------------------------------------------------*/

#ifndef microlayerModel_H
#define microlayerModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"
#include "scalar.H"

#include "phaseModel.H"
#include "solidThermo.H"
#include "fvMatricesFwd.H"
#include "fvc.H"
#include "fvMesh.H"
#include "fvCFD.H"


#include "reconstructionSchemes.H"
#include "singleComponentSatProp.H"

// #include "dictionary.H"
// #include "surfaceFields.H"
// #include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class microlayerModel Declaration
\*---------------------------------------------------------------------------*/

class microlayerModel
    : public dictionary //IOobject //IOdictionary //dictionary
{
private:
    //-
        dictionary microlayerModelCoeffs_;

protected:

        //- Thermo-package of phase 1
        const phaseModel& phase1_;

        //- Thermo-package of phase 2
        const phaseModel& phase2_;

        const solidThermo& solid_;

        //- reference to pressure
        const volScalarField& p_;

        //- reference to saturation model
        singleComponentSatProp& satModel_;

        //- reference to reconstruction scheme
        reconstructionSchemes& surf_;

        volScalarField dml_;


    // Private Member Functions

        //- Disallow copy construct
        microlayerModel(const microlayerModel&);

        //- Disallow default bitwise assignment
        void operator=(const microlayerModel&);

public:

    //-
    dictionary& modelDict();

    //-
    const dictionary& modelDict() const;

public:

    //- Runtime type information
    TypeName("microlayerModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            microlayerModel,
            components,
            (
                const phaseModel& phase1,
                const phaseModel& phase2,
                const solidThermo& solid,
                const volScalarField& p,
                singleComponentSatProp& satModel,
                reconstructionSchemes& surf,
                const dictionary& dict
            ),
            (phase1,phase2,solid,p,satModel,surf,dict)
        );

    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<microlayerModel> New
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const solidThermo& solid,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );

    // Constructors

        //- Construct from components
        microlayerModel
        (
            const word& type,
            const phaseModel& phase1,
            const phaseModel& phase2,
            const solidThermo& solid,
            const volScalarField& p,
            singleComponentSatProp& satModel,
            reconstructionSchemes& surf,
            const dictionary& dict
        );

    //- Destructor
    virtual ~microlayerModel()
    {}

    // Member Functions

        virtual tmp<fvScalarMatrix> hSourceML() = 0;

        virtual void initialiseML() = 0;
        
        virtual void updateML() = 0;

        virtual tmp<volScalarField> energySourceML() = 0;
        
        virtual tmp<volScalarField> massSourceML(volScalarField& rhoSource) = 0;

        virtual tmp<volScalarField> alphaSourceML(volScalarField& rhoSource) = 0;

        virtual volScalarField& dml();

        virtual tmp<volScalarField> dml() const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
