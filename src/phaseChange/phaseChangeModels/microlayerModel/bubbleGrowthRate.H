/*---------------------------------------------------------------------------*\
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
    Foam:: bubbleGrowthRate

Description
    Computes and outputs bubble growth rate and acceleration along a patch for
    calculations of microalyer thickenss
SourceFiles
    bubbleGrowthRate.C

\*---------------------------------------------------------------------------*/

#ifndef bubbleGrowthRate_H
#define bubbleGrowthRate_H

#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"
#include "fvMesh.H"
#include "IOobject.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"
#include "singleComponentSatProp.H"
#include "surfaceFields.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                Class ChenUtaka Declaration
\*---------------------------------------------------------------------------*/
    class BubbleGrowthRate : public dictionary
    {
    private:
        scalar prevTime_;
        scalar prevRad_;
        scalar prevRate_;
        scalar currTime_;

        // Private Member Functions

            //- Disallow copy construct
            BubbleGrowthRate(const BubbleGrowthRate&);

            //- Disallow default bitwise assignment
            void operator=(const BubbleGrowthRate&);

    public:
        //- Runtime type information
        TypeName("BubbleGrowthRate");

        // Constructors
        //- Construct from components
        BubbleGrowthRate
        (
            const string& fluidPatch,
            const phaseModel& phase1,
            const Time& runTime,
            const scalar& startTime,
            const scalar& startRad,
            const scalar& startRate,
            const dictionary& dict
        );

        //- Destructor
        virtual ~BubbleGrowthRate()
        {}

        // Member Functions

        virtual scalar computeBubbleGrowthRate
        (    
            const fvMesh& mesh,
            const volScalarField& alpha,
            const word& patchName,
            const point& bubbleCenter,
            scalar& R,
            scalar& dRdt,
            scalar& d2Rdt2,
            label& triggeredCell
        );

        virtual scalar computeBubbleGrowthAcceleration() = 0;


    };
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
