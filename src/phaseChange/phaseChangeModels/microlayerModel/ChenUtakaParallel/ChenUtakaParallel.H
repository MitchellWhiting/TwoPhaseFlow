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
    Foam::ChenUtakaParallel

Description
    interface heat resistance model based on the ChenUtakaParallel equation

SourceFiles
    ChenUtakaParallel.C

\*---------------------------------------------------------------------------*/

#ifndef ChenUtakaParallel_H
#define ChenUtakaParallel_H

#include "typeInfo.H"
#include "microlayerModel.H"
#include "bubbleGrowthRate.H"

#include "volFields.H"
#include "dimensionedScalar.H"
#include "scalar.H"

#include "autoPtr.H"

#include "fvScalarMatrix.H"
#include "fvm.H"

#include "reconstructionSchemes.H"
#include "singleComponentSatProp.H"

#include "mixedFvPatchFields.H"
#include "mappedPatchBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class ChenUtakaParallel Declaration
\*---------------------------------------------------------------------------*/

class ChenUtakaParallel : public microlayerModel
{

private:

    // volScalarField dml_;
    
    // Stored Mass Source Term (kg/m3/s)
    volScalarField MSource_;

    scalar evapCoeff_;
    dimensionedScalar Rgas_;
    string fluidPatch_;
    string solidPatch_;
    string method_;
    scalar gradient_;
    scalar coefficient_;
    vector origin_;
    bool debug_;

    // MappedPatchBase for data exchange
    autoPtr<mappedPatchBase> mapper_;
    bool mapperInitialized_;
    label fluidPatchSize_; // For AMR detection

    // Private Member Functions
        template <typename T>
        void distributeField(const labelList& neiProcs,Field<T>& field);
        
        labelList getNeibourProcs(const boolList& nextToInterface);

        //- Disallow copy construct
        ChenUtakaParallel(const ChenUtakaParallel&);

        //- Disallow default bitwise assignment
        void operator=(const ChenUtakaParallel&);
        
        void updateMapper();

public:

    //- Runtime type information
    TypeName("ChenUtakaParallel");

    // Constructors

        //- Construct from components
        ChenUtakaParallel
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
    virtual ~ChenUtakaParallel();

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
