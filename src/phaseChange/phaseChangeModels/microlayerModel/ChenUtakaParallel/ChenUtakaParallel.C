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


#include "ChenUtakaParallel.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "Pstream.H"
#include "PstreamBuffers.H"
#include <algorithm>
#include <vector>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ChenUtakaParallel, 0);
    addToRunTimeSelectionTable(microlayerModel,ChenUtakaParallel, components);

    // Removed old performParallelCoupling which used custom structs

    void Foam::ChenUtakaParallel::updateMapper()
    {
        // Find fluid patch
        const label fluidPatchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);
        if (fluidPatchID == -1) return;
        
        const polyPatch& pp = phase1_.mesh().boundaryMesh()[fluidPatchID];
        
        // AMR Check: If initialized but size changed, reset
        if (mapperInitialized_)
        {
            if (pp.size() != fluidPatchSize_)
            {
                if (debug_) 
                {
                    Info << "ChenUtakaParallel: Mesh topology change detected (AMR). " 
                         << "Patch size changed from " << fluidPatchSize_ << " to " << pp.size() 
                         << ". Recreating mapper." << endl;
                }
                mapperInitialized_ = false;
                // fluidPatchSize_ will be updated below
            }
            else
            {
                return; // Valid mapper exists
            }
        }

        // Check if it is already a mapped patch (e.g. user defined BC)
        if (isA<mappedPatchBase>(pp))
        {
             // If user already set it as mappedPatch, we might conflict if we create another mapper.
             // But usually for microlayer, it is a standard wall patch.
        }

        // Create auxiliary mapper
        // We map FROM solid TO fluid.
        // Sample region: solid mesh name
        // Sample patch: solidPatch_
        
        const word& sampleRegion = solid_.T().mesh().name();
        const word& samplePatch = solidPatch_;
        
        if (debug_)
        {
            Info << "ChenUtakaParallel: Creating mappedPatchBase from " 
                 << sampleRegion << "::" << samplePatch 
                 << " to " << phase1_.mesh().name() << "::" << fluidPatch_ 
                 << " (size: " << pp.size() << ")" << endl;
        }

        mapper_.reset
        (
            new mappedPatchBase
            (
                pp,
                sampleRegion,
                mappedPatchBase::NEARESTPATCHFACE, // sampleMode
                samplePatch,
                0.0 // offset
            )
        );
        
        fluidPatchSize_ = pp.size();
        mapperInitialized_ = true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ChenUtakaParallel::ChenUtakaParallel
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
    ),
    MSource_
    (
        IOobject
        (
            "MSource_ChenUtaka",
            phase1.mesh().time().timeName(),
            phase1.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        phase1.mesh(),
        dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    evapCoeff_(modelDict().lookupOrDefault<scalar>("evapCoeff",1)),
    Rgas_(modelDict().lookupOrDefault<scalar>("Rgas",1)),
    fluidPatch_(modelDict().lookupOrDefault<string>("fluidPatch","fluid_to_solid")), //Only working for 1 patch at the moment
    solidPatch_(modelDict().lookupOrDefault<string>("solidPatch","solid_to_fluid")),
    method_(modelDict().lookupOrDefault<string>("initialisationMethod","ChenUtakaParallel")), //Only working for 1 patch at the moment
    gradient_(modelDict().lookupOrDefault<scalar>("gradient",1)),
    coefficient_(modelDict().lookupOrDefault<scalar>("coefficient_",1)),
    origin_(modelDict().lookupOrDefault<vector>("origin",vector(0, 0, 0))),
    mapper_(nullptr),
    mapperInitialized_(false),
    fluidPatchSize_(-1),
    debug_(modelDict().lookupOrDefault<bool>("debug", false))
{
    if (phase2_.thermo().incompressible())  //Is this needed?
    {
        Rgas_.value() = modelDict().get<scalar>("Rgas");
    }
}

Foam::ChenUtakaParallel::~ChenUtakaParallel()
{
    if (mapper_.valid())
    {
        mapper_.clear();
    }
}

// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * *  //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
// Foam::tmp<Foam::fvScalarMatrix> Foam::ChenUtakaParallel::dmlInitial()
// {
//
// }


Foam::tmp<Foam::fvScalarMatrix> Foam::ChenUtakaParallel::hSourceML()
{
    if (debug_) Info<< "I'm inside hSourceML"<< endl;
    
    updateMapper();

    //Fluid Patch
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();
    const label fluidPatchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);

    // Solid Patch
    const fvMesh& solidFvMesh = solid_.T().mesh();
    const label solidPatchID = solidFvMesh.boundaryMesh().findPatchID(solidPatch_);
    
    if (fluidPatchID == -1 || solidPatchID == -1 || !mapper_.valid())
    {
        if (debug_) Info << "ChenUtakaParallel: Missing patches or mapper invalid in hSourceML" << endl;
        // Return zero source if patches missing
        volScalarField sourceCoeff
        (
            IOobject("sourceCoeff", solidFvMesh.time().timeName(), solidFvMesh, IOobject::NO_READ, IOobject::NO_WRITE),
            solidFvMesh,
            dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0),
            "zeroGradient"
        );
        return tmp<fvScalarMatrix>(fvm::Sp(sourceCoeff, solid_.he()));
    }

    // 1. Map Solid Data to Fluid Patch
    // We need Tsolid, hsolid, alpha, Cp on the fluid patch
    const fvPatch& solidPatch = solidFvMesh.boundary()[solidPatchID];
    const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
    
    // Prepare fields on SOLID patch
    scalarField TsolidFld = solid_.T().boundaryField()[solidPatchID];
    scalarField hsolidFld = solid_.he().boundaryField()[solidPatchID];
    
    // Corrected physics: Use thermal conductivity ksolid instead of alpha*Cp
    // alpha = k / (rho * Cp)  => k = alpha * rho * Cp
    // If we use solid_.kappa(), it's direct.
    // solidThermo usually has kappa() returning volScalarField.
    scalarField kSolidFld = solid_.kappa()().boundaryField()[solidPatchID];

    // Distribute to Fluid Patch
    // mapper_->distribute(fld) pulls data FROM sample (solid) TO local (fluid)
    
    scalarField Tsolid_on_Fluid = TsolidFld;
    mapper_->distribute(Tsolid_on_Fluid);
    
    scalarField hsolid_on_Fluid = hsolidFld;
    mapper_->distribute(hsolid_on_Fluid);
    
    scalarField kSolid_on_Fluid = kSolidFld;
    mapper_->distribute(kSolid_on_Fluid);
    
    // 2. Calculate qml and sourceCoeff on Fluid Patch
    scalarField sourceCoeff_on_Fluid(fluidPatch.size(), 0.0);
    
    const volScalarField& TSat = satModel_.TSat();
    const volScalarField& k1 = phase1_.kappa();
    const volScalarField& rho2 = phase2_.thermo().rho();
    
    scalar kMax = gMax(k1.internalField());
    
    const dimensionedScalar RgasDim ("Rgas", dimGasConstant, Rgas_.value());
    volScalarField RintField
    (
        (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*RgasDim,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
    );

    const labelUList& fluidFaceCells = fluidPatch.faceCells();
    
    // We need a second pass of mapping for geometry if we want exact match
    // ySolid is distance from solid cell center to face
    scalarField ySolidFld(solidPatch.size());
    const labelUList& solidFaceCells = solidPatch.faceCells();
    forAll(solidPatch, i)
    {
        ySolidFld[i] = mag(solidFvMesh.C()[solidFaceCells[i]] - solidPatch.Cf()[i]);
    }
    scalarField ySolid_on_Fluid = ySolidFld;
    mapper_->distribute(ySolid_on_Fluid);
    
    label processedCount = 0;
    label activeCount = 0;
    label skippedPhaseCount = 0;
    label skippedDmlCount = 0;
    scalar maxQ = 0.0;

    forAll(fluidPatch, faceI)
    {
        label cellI = fluidFaceCells[faceI];
        scalar phase1Val = phase1_[cellI];
        scalar dmlVal = dml_[cellI];
        
        if (phase1Val >= 0.001) 
        {
            skippedPhaseCount++;
            continue;
        }
        if (dmlVal <= 1e-10) 
        {
            skippedDmlCount++;
            continue;
        }
        
        processedCount++;

        scalar T_s = Tsolid_on_Fluid[faceI];
        scalar h_s = hsolid_on_Fluid[faceI];
        scalar y_s = ySolid_on_Fluid[faceI];
        scalar k_s = kSolid_on_Fluid[faceI];
        
        // Corrected kdsolid (Conductance W/m2K)
        scalar kdsolid = k_s / (y_s + SMALL);
        
        scalar RintVal = RintField.boundaryField()[fluidPatchID][faceI];
        
        // kdfluid (Conductance W/m2K)
        scalar kdfluid = 1.0 / ((dmlVal / kMax) + RintVal);
        
        scalar TSatVal = TSat.boundaryField()[fluidPatchID][faceI];
        
        scalar Twall = (kdsolid * T_s + kdfluid * TSatVal) / (kdsolid + kdfluid);
        scalar qml = kdfluid * (Twall - TSatVal);
        
        scalar yFluid = 2.0 * mag(mesh.C()[cellI] - fluidPatch.Cf()[faceI]);
        
        sourceCoeff_on_Fluid[faceI] = qml / (yFluid * h_s + SMALL);
        
        if (mag(qml) > maxQ) maxQ = mag(qml);
        if (mag(qml) > SMALL) activeCount++;
    }
    
    if (debug_)
    {
        label totalProcessed = processedCount;
        label totalActive = activeCount;
        scalar globalMaxQ = maxQ;
        
        reduce(totalProcessed, sumOp<label>());
        reduce(totalActive, sumOp<label>());
        reduce(globalMaxQ, maxOp<scalar>());
        
        if (Pstream::master())
        {
            Info << "ChenUtakaParallel hSourceML: Total Processed: " << totalProcessed 
                 << ", Total Active Q > 0: " << totalActive 
                 << ", Global Max Qml: " << globalMaxQ << endl;
        }
    }
    
    // 3. Map sourceCoeff back to Solid
    scalarField sourceCoeff_on_Solid = sourceCoeff_on_Fluid;
    mapper_->reverseDistribute(sourceCoeff_on_Solid);
    
    // 5. Apply to solid field
    volScalarField sourceCoeff
    (
        IOobject("sourceCoeff", solidFvMesh.time().timeName(), solidFvMesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        solidFvMesh,
        dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0),
        "zeroGradient"
    );
    
    forAll(solidPatch, i)
    {
        label cellI = solidFaceCells[i];
        sourceCoeff[cellI] = sourceCoeff_on_Solid[i];
    }
    
    sourceCoeff.correctBoundaryConditions();
    tmp<fvScalarMatrix> hSource(fvm::Sp(sourceCoeff, solid_.he()));

    return hSource;
}

void Foam::ChenUtakaParallel::initialiseML()
{
    // Support both naming conventions for compatibility
    if (method_ == "ChenUtakaParallel")
    {
        Info<< "ChenUtakaParallel: Executing initialiseML with method " << method_ << endl;

        const label patchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);
        if (patchID == -1) 
        {
            WarningInFunction << "Fluid patch " << fluidPatch_ << " not found!" << endl;
            return;
        }

        const polyPatch& cPatch = phase1_.mesh().boundaryMesh()[patchID];
        const labelUList& faceCells = cPatch.faceCells();

        vector point (0, 0, 0);
        scalar distance (0.0);

        scalar maxInitDml = 0.0;

        forAll(faceCells, cellI)
        { 
            const label globalCellI = faceCells[cellI];
            point = phase1_.mesh().C()[globalCellI];
            distance = mag(point - origin_);

            scalar addedDml = gradient_*distance;
            dml_[globalCellI] += addedDml;
            
            if (addedDml > maxInitDml) maxInitDml = addedDml;
        }
        
        reduce(maxInitDml, maxOp<scalar>());
        Info << "ChenUtakaParallel: Initialised dml. Global max added = " << maxInitDml << endl;
    }
    else
    {
        Info << "ChenUtakaParallel: Skipping initialiseML. Method '" << method_ << "' does not match 'ChenUtakaParallel'." << endl;
    }

}

void Foam::ChenUtakaParallel::updateML()
{
    if (debug_) Info<< "I'm inside updateML"<< endl;
    
    updateMapper();

    //Fluid Patch
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();
    const label fluidPatchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);
    
    // Reset MSource_
    MSource_ = dimensionedScalar("0", MSource_.dimensions(), 0.0);

    if (fluidPatchID != -1 && mapper_.valid())
    {
        const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
        const labelUList& fluidFaceCells = fluidPatch.faceCells();
        
        // Similar to hSourceML but now we update dml and MSource
        
        // We need Tsolid, ySolid from solid
        const fvMesh& solidFvMesh = solid_.T().mesh();
        const label solidPatchID = solidFvMesh.boundaryMesh().findPatchID(solidPatch_);
        const fvPatch& solidPatch = solidFvMesh.boundary()[solidPatchID];
        
        scalarField TsolidFld = solid_.T().boundaryField()[solidPatchID];
        scalarField kSolidFld = solid_.kappa()().boundaryField()[solidPatchID];
        
        scalarField ySolidFld(solidPatch.size());
        const labelUList& solidFaceCells = solidPatch.faceCells();
        forAll(solidPatch, i) ySolidFld[i] = mag(solidFvMesh.C()[solidFaceCells[i]] - solidPatch.Cf()[i]);
        
        scalarField Tsolid_on_Fluid = TsolidFld;
        mapper_->distribute(Tsolid_on_Fluid);
        
        scalarField kSolid_on_Fluid = kSolidFld;
        mapper_->distribute(kSolid_on_Fluid);
        
        scalarField ySolid_on_Fluid = ySolidFld;
        mapper_->distribute(ySolid_on_Fluid);
        
        const volScalarField& TSat = satModel_.TSat();
        const volScalarField& k1 = phase1_.kappa();
        const volScalarField& rho1 = phase1_.thermo().rho();
        const volScalarField& rho2 = phase2_.thermo().rho();
        
        scalar rho1Max = gMax(rho1.internalField());
        scalar kMax = gMax(k1.internalField());
        scalar deltaT = mesh.time().deltaTValue();
        
        const dimensionedScalar RgasDim ("Rgas", dimGasConstant, Rgas_.value());
        volScalarField RintField
        (
            (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*RgasDim,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
        );
        
        label updatedCount = 0;
        label skippedPhaseCount = 0;
        label skippedDmlCount = 0;
        
        // Debugging first active cell
        bool debugPrinted = false;

        forAll(fluidPatch, faceI)
        {
            label cellI = fluidFaceCells[faceI];
            scalar phase1Val = phase1_[cellI];
            scalar dmlVal = dml_[cellI];
            
            if (phase1Val < 0.5) // Check phase
            {
                if (dmlVal > 1e-10)
                {
                    scalar T_s = Tsolid_on_Fluid[faceI];
                    scalar y_s = ySolid_on_Fluid[faceI];
                    scalar k_s = kSolid_on_Fluid[faceI];
                    
                    // Corrected physics: kdsolid = k / y
                    scalar kdsolid = k_s / (y_s + SMALL);
                    
                    scalar RintVal = RintField.boundaryField()[fluidPatchID][faceI];
                    scalar kdfluid = 1.0 / ((dmlVal / kMax) + RintVal);
                    scalar TSatVal = TSat.boundaryField()[fluidPatchID][faceI];
                    
                    scalar Twall = (kdsolid * T_s + kdfluid * TSatVal) / (kdsolid + kdfluid);
                    scalar qml = kdfluid * (Twall - TSatVal);
                    
                    scalar LVal = satModel_.L().boundaryField()[fluidPatchID][faceI];
                    scalar mSourceFlux = (LVal > SMALL) ? (qml / LVal) : 0.0;
                    
                    scalar dmlNew = max(0.0, dmlVal - mSourceFlux * deltaT / rho1Max);
                    
                    dml_[cellI] = dmlNew;
                    
                    scalar faceArea = fluidPatch.magSf()[faceI];
                    scalar cellVol = mesh.V()[cellI];
                    MSource_[cellI] = mSourceFlux * faceArea / cellVol;
                    
                    updatedCount++;
                    
                    if (debug_ && !debugPrinted && updatedCount == 1)
                    {
                        Info << "ChenUtakaParallel DEBUG Single Face:" << endl
                             << "  FaceI: " << faceI << " dml: " << dmlVal << " -> " << dmlNew << endl
                             << "  Tsolid: " << T_s << " TSat: " << TSatVal << " Twall: " << Twall << endl
                             << "  ySolid: " << y_s << " kSolid: " << k_s << " kdsolid: " << kdsolid << endl
                             << "  Rint: " << RintVal << " kdfluid: " << kdfluid << endl
                             << "  qml: " << qml << " mSourceFlux: " << mSourceFlux << endl;
                        debugPrinted = true;
                    }
                }
                else
                {
                    // dml too small
                    skippedDmlCount++;
                    dml_[cellI] = 0.0;
                    MSource_[cellI] = 0.0;
                }
            }
            else
            {
                // Not liquid
                skippedPhaseCount++;
                dml_[cellI] = 0.0;
                MSource_[cellI] = 0.0;
            }
        }
        
        if (debug_)
        {
            label totalUpdated = updatedCount;
            reduce(totalUpdated, sumOp<label>());
            
            if (Pstream::master())
            {
                Info << "ChenUtakaParallel updateML: Total Updated dml/MSource on " << totalUpdated << " cells." << endl;
            }
        }
    }

    // Continue with smearing logic (Linear Extrapolation)
    // This part requires parallel reductions for rMax1, etc.
    // to be correct globally.

    scalar rMax1(0.0); // Initialize to 0.0 to match Serial behavior and avoid huge dml if not found
    scalar yMax1(0.0); // local yMax1

    // We only loop over local fluid patch to find local min/max
    if (fluidPatchID != -1)
    {
        const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
        const labelUList& fluidFaceCells = fluidPatch.faceCells();

        forAll(fluidPatch, faceI)
        {
            label fluidCellI = fluidFaceCells[faceI];
            if (phase1_[fluidCellI] < 0.5)
            {
                 scalar rPos = mag(mesh.C()[fluidCellI] - origin_);
                 if (rPos > rMax1) 
                 {
                     rMax1 = rPos;
                     yMax1 = dml_[fluidCellI];
                 }
            }
        }
    }

    // Parallel Reduce
    // Reduce rMax1 and corresponding yMax1
    // We need to find which proc has rMax1 and what its yMax1 is.
    // Standard reduce(rMax1, maxOp) only gives us the global max R.
    
    // 1. Find global rMax1
    scalar globalRMax1 = rMax1;
    reduce(globalRMax1, maxOp<scalar>());

    // 2. Determine who has it.
    // If local rMax1 == globalRMax1 (within tolerance), we broadcast our yMax1
    
    scalar globalYMax1 = 0.0;
    
    // If I hold the max, I propose my yMax. Otherwise I propose -GREAT.
    
    scalar myProposedYMax = -GREAT;
    // Only propose if I actually found a valid max (rMax1 > -GREAT check implicit if initialized to 0.0 or valid)
    if (mag(rMax1 - globalRMax1) < SMALL)
    {
        myProposedYMax = yMax1; // Assuming yMax1 was set locally
    }
    
    reduce(myProposedYMax, maxOp<scalar>());
    globalYMax1 = myProposedYMax;
    
    // If globalYMax1 is still -GREAT (e.g. no fluid patch active?), set to 0
    if (globalYMax1 <= -GREAT/2) globalYMax1 = 0.0;

    rMax1 = globalRMax1;
    yMax1 = globalYMax1;

    // Apply linear extrapolation to all fluid cells (liquid phase)
    if (fluidPatchID != -1)
    {
        const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
        const labelUList& fluidFaceCells = fluidPatch.faceCells();

        forAll(fluidPatch, faceI)
        {
            label fluidCellI = fluidFaceCells[faceI];
            if (phase1_[fluidCellI] > 0.5)
            {
                 scalar rPos = mag(mesh.C()[fluidCellI] - origin_);
                 scalar distance = mag(rPos - rMax1); // Use mag to be safe
                 
                 // Apply extrapolation logic from Serial code
                 dml_[fluidCellI] = gradient_*distance + yMax1; 
            }

            if (method_ == "cooperLloyd" && phase1_[fluidCellI] >= 0.5)
            {
                const volScalarField& rho1 = phase1_.thermo().rho();
                const volScalarField& mu1 = phase1_.thermo().mu();
                dml_[fluidCellI] = coefficient_*pow((mu1[fluidCellI]/rho1[fluidCellI])*phase1_.mesh().time().value() ,0.5);
            }
        }
    }

    if (debug_)
    {
        scalar globalMaxMSource = gMax(MSource_.internalField());
        scalar globalMinMSource = gMin(MSource_.internalField());
        scalar globalMaxDml = gMax(dml_);
        scalar globalMinDml = gMin(dml_);
        
        // gMax/gMin already perform parallel reduction, so we don't need manual reduce()
        // Wait, gMax on volScalarField returns dimensionedScalar which is already reduced?
        // Let's verify. gMax(GeometricField) usually does return a global max.
        // Yes, gMax/gMin in OpenFOAM perform parallel reduction automatically.
        
        if (Pstream::master())
        {
            Info << "ChenUtakaParallel Debug: MSource_ max = " << globalMaxMSource 
                 << ", min = " << globalMinMSource << endl;
            Info << "ChenUtakaParallel Debug: dml_ max = " << globalMaxDml 
                 << ", min = " << globalMinDml << endl;
        }
    }
}

Foam::tmp<Foam::volScalarField> Foam::ChenUtakaParallel::energySourceML()
{
    if (debug_) Info<< "I'm inside energySourceML"<< endl;
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();
    const label fluidPatchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);

    // Prepare Energy Source Field on FLUID MESH
    tmp<volScalarField> tEnergySource
    (
         new volScalarField
         (
            IOobject
            (
                "energySourceML",
                mesh.time().timeName(),
                mesh, // Fluid Mesh
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0) // W/m3
         )
    );
    volScalarField& energySource = tEnergySource.ref();

    // Use stored MSource_ for energy source
    // energySource = MSource_ * L

    if (fluidPatchID != -1)
    {
        const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
        const labelUList& fluidFaceCells = fluidPatch.faceCells();
        
        // MSource_ has been updated in updateML(). 
        const volScalarField& L = satModel_.L();
        
        forAll(fluidFaceCells, i)
        {
            label cellI = fluidFaceCells[i];
            // Check if cellI is valid (it should be, but just in case)
            if (cellI >= 0 && cellI < mesh.nCells())
            {
                energySource[cellI] = MSource_[cellI] * L[cellI];
            }
        }
    }

    return tEnergySource;

}



Foam::tmp<Foam::volScalarField>
Foam::ChenUtakaParallel::massSourceML( volScalarField& rhoSource)
{
    if (debug_) Info<< "I'm inside massSourceML"<< endl;

    // Use stored MSource_ for Smearing
    // The input rhoSource is often empty or dummy, so we overwrite it with MSource_ 
    // or assume psiEqn RHS is MSource_.

    tmp<volScalarField> massSourceML(rhoSource * 0.0);
    volScalarField& massSourceMLRef = massSourceML.ref();

    const fvMesh& mesh = phase1_.mesh();

    dimensionedScalar DPsi
    (
        "DPsi",
        dimensionSet(0,2,0,0,0,0,0),
        3/sqr(gAverage(mesh.nonOrthDeltaCoeffs()))
    );

    // Use MSource_ as the source term for smearing
    dimensionedScalar intPsi0 = fvc::domainIntegrate(MSource_);

    volScalarField psi
    (
        IOobject
        (
            "psi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimDensity/dimTime, 0),
        "zeroGradient"
    );


    //- Smearing of source term field
    fvScalarMatrix psiEqn
    (
        fvm::Sp(scalar(1),psi) - fvm::laplacian(DPsi,psi) == MSource_
    );

    psiEqn.solve();

    // Cut cells with cutoff < alpha1 < 1-cutoff
    // and rescale remaining source term field
    dimensionedScalar intPsiVapor
    (
        "intPsiVapor",
        dimensionSet(1,0,-1,0,0,0,0),
        0.0
    );

    forAll(mesh.C(),celli)
    {
        if (phase1_[celli] < 1e-3)
        {
            intPsiVapor.value() +=
                (1.0-phase1_[celli])*psi[celli]*mesh.V()[celli];
        }
    }

    //- Calculate Nl and Nv
    dimensionedScalar Nv ("Nv", dimless, 2.0);

    reduce(intPsiVapor.value(),sumOp<scalar>());

    if (intPsiVapor.value() > 1e-99)
    {
        Nv = intPsi0/intPsiVapor;
    }

    //- Set source terms in cells with alpha1 < cutoff or alpha1 > 1-cutoff
    forAll(mesh.C(),celli)
    {
        if (phase1_[celli] < 1e-3)
        {
            massSourceMLRef[celli] = Nv.value()*(1.0-phase1_[celli])*psi[celli];
        }
        else
        {
            massSourceMLRef[celli] = 0.0;
        }
    }

    return massSourceML;
}


Foam::tmp<Foam::volScalarField>
Foam::ChenUtakaParallel::alphaSourceML( volScalarField& rhoSource)
{
    // Use MSource_ directly since rhoSource is likely zero
    if (debug_) Info<< "I'm inside alphaSourceML"<< endl;

    tmp<volScalarField> alphaSourceML
    (
        MSource_ / phase1_.thermo().rho()
    );

   return alphaSourceML;
}
