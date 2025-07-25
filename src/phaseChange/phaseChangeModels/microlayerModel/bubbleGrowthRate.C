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

\*---------------------------------------------------------------------------*/
#include "bubbleGrowthRate.H"
#include "fvMesh.H"
#include "vector.H"
#include "Time.H"
#include "volFields.H"
#include "fvc.H"
#include "OFstream.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// namespace Foam
// { 
//     defineTypeNameAndDebug(BubbleGrowthRate, 0);
//     defineRunTimeSelectionTable(BubbleGrowthRate, components);
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::BubbleGrowthRate::BubbleGrowthRate
(
    const string& fluidPatch,
    const phaseModel& phase1,
    const Time& runTime,
    const scalar& startTime,
    const scalar& startRad,
    const scalar& startRate,
    const dictionary& dict
)
:
    dictionary(dict),
    prevRad_(startRad),
    prevRate_(startRate),
    prevTime_(startTime),
    currTime_(runTime.timeOutputValue())
{
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Foam::BubbleGrowthRate::~BubbleGrowthRate() 
// {}
// * * * * * * * * * * * * * * Public Access Member Functions  * * * * * * * //

Foam::scalar
Foam::BubbleGrowthRate::computeBubbleGrowthRate(    
    const fvMesh& mesh,
    const volScalarField& alpha,
    const word& patchName,
    const point& bubbleCenter,
    scalar& R,
    scalar& dRdt,
    scalar& d2Rdt2,
    label& triggeredCell
)
{

    label patchID = mesh.boundaryMesh().findPatchID(patchName);

    if (patchID < 0)
    {
        FatalErrorInFunction
            << "Patch " << patchName << " not found!" << nl << exit(FatalError);
    }

    const polyPatch& pp = mesh.boundaryMesh()[patchID];
    const labelList& faceCells = pp.faceCells();
    const vectorField& Cf = pp.faceCentres();

    scalar maxRadius = 0.0;
    triggeredCell = -1;

    forAll(faceCells, i)
    {
        label cI = faceCells[i];
        if (alpha[cI] < 0.001)
        {
            scalar r = mag(Cf[i] - bubbleCenter);
            if (r > maxRadius)
            {
                maxRadius = r;
                triggeredCell = cI;
            }
        }
    }

    scalar currentTime = mesh.time().timeOutputValue();
    scalar dt = currentTime - prevTime_;

    dRdt = (maxRadius - prevRad_) / dt;
    d2Rdt2 = (dRdt - prevRate_) / dt;
    R = maxRadius;

    prevRad_ = R;
    prevTime_ = currentTime;
    prevRate_ = dRdt;

    Info << "\n[calculateBubbleRadius] R = " << R
         << ", dRdt = " << dRdt
         << ", d2Rdt2 = " << d2Rdt2 << endl;

    if (triggeredCell >= 0)
    {
        Info << "[calculateBubbleRadius] Cell " << triggeredCell
             << " just transitioned to gas (alpha < 0.001)" << endl;
    }
    return R, dRdt, currentTime;
}
