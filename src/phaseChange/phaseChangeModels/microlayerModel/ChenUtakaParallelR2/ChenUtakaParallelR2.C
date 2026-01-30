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


#include "ChenUtakaParallelR2.H"
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
    // Local structs for parallel communication
    namespace {
    struct FluidPacket
    {
        label procID;
        label cellID;
        vector faceCentre;
        vector cellCentre;
        scalar TSat;
        scalar dml;
        scalar Rint;
        scalar kMax;
        scalar rho1Max;
        scalar mu1;
        scalar phase1Val;
        scalar L; // Latent heat

        friend Istream& operator>>(Istream& is, FluidPacket& d)
        {
            return is >> d.procID >> d.cellID >> d.faceCentre >> d.cellCentre
                      >> d.TSat >> d.dml >> d.Rint >> d.kMax
                      >> d.rho1Max >> d.mu1 >> d.phase1Val >> d.L;
        }
        friend Ostream& operator<<(Ostream& os, const FluidPacket& d)
        {
            return os << d.procID << token::SPACE << d.cellID << token::SPACE
                      << d.faceCentre << token::SPACE << d.cellCentre << token::SPACE
                      << d.TSat << token::SPACE << d.dml << token::SPACE
                      << d.Rint << token::SPACE << d.kMax << token::SPACE
                      << d.rho1Max << token::SPACE << d.mu1 << token::SPACE
                      << d.phase1Val << token::SPACE << d.L;
        }

        bool operator==(const FluidPacket& rhs) const
        {
            return procID == rhs.procID &&
                   cellID == rhs.cellID &&
                   faceCentre == rhs.faceCentre &&
                   cellCentre == rhs.cellCentre &&
                   TSat == rhs.TSat &&
                   dml == rhs.dml &&
                   Rint == rhs.Rint &&
                   kMax == rhs.kMax &&
                   rho1Max == rhs.rho1Max &&
                   mu1 == rhs.mu1 &&
                   phase1Val == rhs.phase1Val &&
                   L == rhs.L;
        }

        bool operator!=(const FluidPacket& rhs) const
        {
            return !(*this == rhs);
        }
    };

    struct SolidPacket
    {
        label procID;
        label cellID;
        vector faceCentre;
        vector cellCentre;
        scalar Tsolid;
        scalar hsolid;
        scalar ksolid; 
        scalar Cpsolid;
        // scalar rhosolid; unused in whole context
        scalar alphasolid;

        friend Istream& operator>>(Istream& is, SolidPacket& d)
        {
            return is >> d.procID >> d.cellID >> d.faceCentre >> d.cellCentre
                      >> d.Tsolid >> d.hsolid >> d.ksolid >> d.Cpsolid >> d.alphasolid;
        }
        friend Ostream& operator<<(Ostream& os, const SolidPacket& d)
        {
            return os << d.procID << token::SPACE << d.cellID << token::SPACE
                      << d.faceCentre << token::SPACE << d.cellCentre << token::SPACE
                      << d.Tsolid << token::SPACE << d.hsolid << token::SPACE
                      << d.ksolid << token::SPACE << d.Cpsolid << token::SPACE << d.alphasolid;
        }

        bool operator==(const SolidPacket& rhs) const
        {
            return procID == rhs.procID &&
                   cellID == rhs.cellID &&
                   faceCentre == rhs.faceCentre &&
                   cellCentre == rhs.cellCentre &&
                   Tsolid == rhs.Tsolid &&
                   hsolid == rhs.hsolid &&
                   ksolid == rhs.ksolid &&
                   Cpsolid == rhs.Cpsolid &&
                   alphasolid == rhs.alphasolid;
        }

        bool operator!=(const SolidPacket& rhs) const
        {
            return !(*this == rhs);
        }
    };

    struct FluidResult
    {
        label cellID;
        scalar dmlNew;
        scalar MSource;

        friend Istream& operator>>(Istream& is, FluidResult& d)
        {
            return is >> d.cellID >> d.dmlNew >> d.MSource;
        }
        friend Ostream& operator<<(Ostream& os, const FluidResult& d)
        {
            return os << d.cellID << token::SPACE << d.dmlNew << token::SPACE << d.MSource;
        }

        bool operator==(const FluidResult& rhs) const
        {
            return cellID == rhs.cellID &&
                   dmlNew == rhs.dmlNew &&
                   MSource == rhs.MSource;
        }

        bool operator!=(const FluidResult& rhs) const
        {
            return !(*this == rhs);
        }
    };

    struct SolidResult
    {
        label cellID;
        scalar sourceCoeff;

        friend Istream& operator>>(Istream& is, SolidResult& d)
        {
            return is >> d.cellID >> d.sourceCoeff;
        }
        friend Ostream& operator<<(Ostream& os, const SolidResult& d)
        {
            return os << d.cellID << token::SPACE << d.sourceCoeff;
        }

        bool operator==(const SolidResult& rhs) const
        {
            return cellID == rhs.cellID &&
                   sourceCoeff == rhs.sourceCoeff;
        }

        bool operator!=(const SolidResult& rhs) const
        {
            return !(*this == rhs);
        }
    };
    } // End anonymous namespace

    defineTypeNameAndDebug(ChenUtakaParallelR2, 0);
    addToRunTimeSelectionTable(microlayerModel,ChenUtakaParallelR2, components);

    // Helper function for parallel coupling
    namespace {
    void performParallelCoupling(
        List<FluidPacket>& localFluidPackets,
        List<SolidPacket>& localSolidPackets,
        List<FluidResult>& fluidResults,
        List<SolidResult>& solidResults,
        scalar deltaT,
        const Foam::word& mode,
        bool debug
    )
    {
        // 1. Gather all data to master
        List<List<FluidPacket>> allFluidPackets;
        List<List<SolidPacket>> allSolidPackets;

        if (Pstream::parRun())
        {
            // Manual gather implementation
            if (Pstream::master())
            {
                allFluidPackets.setSize(Pstream::nProcs());
                allFluidPackets[Pstream::myProcNo()] = localFluidPackets;

                allSolidPackets.setSize(Pstream::nProcs());
                allSolidPackets[Pstream::myProcNo()] = localSolidPackets;

                for (int proci = 0; proci < Pstream::nProcs(); ++proci)
                {
                    if (proci != Pstream::masterNo())
                    {
                        IPstream fromProc(Pstream::commsTypes::scheduled, proci);
                        fromProc >> allFluidPackets[proci];
                        fromProc >> allSolidPackets[proci];
                    }
                }
            }
            else
            {
                OPstream toMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
                toMaster << localFluidPackets;
                toMaster << localSolidPackets;
            }
        }
        else
        {
            allFluidPackets.setSize(1);
            allFluidPackets[0] = localFluidPackets;
            allSolidPackets.setSize(1);
            allSolidPackets[0] = localSolidPackets;
        }

        // 2. Master (or single proc) computes interactions
        List<List<FluidResult>> allFluidResults(Pstream::nProcs());
        List<List<SolidResult>> allSolidResults(Pstream::nProcs());

        if (Pstream::master())
        {
            // Collect all solid faces into a searchable structure
            List<SolidPacket> globalSolidFaces;
            forAll(allSolidPackets, proci)
            {
                forAll(allSolidPackets[proci], i)
                {
                    globalSolidFaces.append(allSolidPackets[proci][i]);
                }
            }

            // Iterate over all fluid faces
            label totalFluidFaces = 0;
            label skippedPhase = 0;
            label skippedDml = 0;
            label processed = 0;

            forAll(allFluidPackets, fluidProcI)
            {
                const List<FluidPacket>& fluidFaces = allFluidPackets[fluidProcI];
                List<FluidResult>& procFluidResults = allFluidResults[fluidProcI];
                procFluidResults.setSize(fluidFaces.size());
                
                totalFluidFaces += fluidFaces.size();

                forAll(fluidFaces, i)
                {
                    const FluidPacket& fData = fluidFaces[i];
                    FluidResult& fRes = procFluidResults[i];
                    fRes.cellID = fData.cellID;
                    fRes.dmlNew = fData.dml;
                    fRes.MSource = 0.0;
                    
                    // Restore filtering depending on "mode" word:
                    // 1. Check if we are in gas or interface (phase1 < 0.5)
                    // 2. Check if dml exists (primary condition for this model) - Restored as per serial code
                    if (mode == "hSourceML")
                    {
                        if (fData.phase1Val >= 0.001)
                        {
                            skippedPhase++;
                            continue;
                        }
                        
                        if (fData.dml <= 1e-10)
                        {
                            skippedDml++;
                            continue;
                        }
                        
                        processed++;

                        // Find nearest solid face
                        scalar minDist = GREAT;
                        label bestSolidIdx = -1;

                        forAll(globalSolidFaces, sIdx)
                        {
                            scalar dist = mag(globalSolidFaces[sIdx].faceCentre - fData.faceCentre);
                            
                            if (dist < minDist)
                            {
                                minDist = dist;
                                bestSolidIdx = sIdx;
                            }
                        }

                        // Reasonable tolerance (increased to 5e-3 = 5mm for testing)
                        if (bestSolidIdx != -1 && minDist < 5e-3) 
                        {
                            const SolidPacket& sData = globalSolidFaces[bestSolidIdx];

                            scalar yFluid = 2.0 * mag(fData.cellCentre - fData.faceCentre);
                            scalar ySolid = mag(sData.cellCentre - sData.faceCentre);

                            // scalar kdsolid = sData.ksolid / ySolid;//💡reCheck
                            scalar kdsolid = sData.alphasolid * sData.Cpsolid / ySolid;
                            // Use raw dml as per serial code, since we filtered dml <= 1e-10 above.
                            scalar localDml = fData.dml; 
                            scalar kdfluid = 1.0 / ((localDml / fData.kMax) + fData.Rint);

                            scalar Twall = (kdsolid * sData.Tsolid + kdfluid * fData.TSat)
                                        / (kdsolid + kdfluid);

                            scalar qml = kdfluid * (Twall - fData.TSat);

                            // Fluid side results
                            fRes.MSource = (qml / fData.L); //💡Keep the variable for return structure, but we do not use it here.

                            // Update dml, but prevent it from going negative
                            fRes.dmlNew = max(0.0, fData.dml - fRes.MSource * deltaT / fData.rho1Max); //💡Keep the variable for return structure, but we do not use it here.

                            // Solid side results
                            SolidResult sRes;
                            sRes.cellID = sData.cellID;
                            sRes.sourceCoeff = qml / (yFluid * sData.hsolid);

                            allSolidResults[sData.procID].append(sRes);
                        }
                    }

                    else if (mode == "updateML")
                    {
                        if (fData.phase1Val < 0.5)
                        {
                            skippedPhase++;
                            continue;
                        }
                        
                        if (fData.dml > 1e-10)
                        {
                            skippedDml++;
                            continue;
                        }
                        
                        processed++;

                        // Find nearest solid face
                        scalar minDist = GREAT;
                        label bestSolidIdx = -1;

                        forAll(globalSolidFaces, sIdx)
                        {
                            scalar dist = mag(globalSolidFaces[sIdx].faceCentre - fData.faceCentre);
                            
                            if (dist < minDist)
                            {
                                minDist = dist;
                                bestSolidIdx = sIdx;
                            }
                        }

                        // Reasonable tolerance (increased to 5e-3 = 5mm for testing)
                        if (bestSolidIdx != -1 && minDist < 5e-3) 
                        {
                            const SolidPacket& sData = globalSolidFaces[bestSolidIdx];

                            scalar yFluid = 2.0 * mag(fData.cellCentre - fData.faceCentre);
                            scalar ySolid = mag(sData.cellCentre - sData.faceCentre);

                            // scalar kdsolid = sData.ksolid / ySolid;//💡reCheck
                            scalar kdsolid = sData.alphasolid * sData.Cpsolid / ySolid;
                            // Use raw dml as per serial code, since we filtered dml <= 1e-10 above.
                            scalar localDml = fData.dml; 
                            scalar kdfluid = 1.0 / ((localDml / fData.kMax) + fData.Rint);

                            scalar Twall = (kdsolid * sData.Tsolid + kdfluid * fData.TSat)
                                        / (kdsolid + kdfluid);

                            scalar qml = kdfluid * (Twall - fData.TSat);

                            // Fluid side results
                            fRes.MSource = (qml / fData.L);
                            // Update dml, but prevent it from going negative
                            fRes.dmlNew = max(0.0, fData.dml - fRes.MSource * deltaT / fData.rho1Max);

                            // Solid side results //💡Keep the variable for return structure, but we do not use it here.
                            SolidResult sRes; 
                            sRes.cellID = sData.cellID;
                            sRes.sourceCoeff = qml / (yFluid * sData.hsolid);

                            allSolidResults[sData.procID].append(sRes);
                        }
                    }

                    else if (mode == "energySourceML")
                    {
                        if (fData.phase1Val < 0.001)
                        {
                            skippedPhase++;
                            continue;
                        }
                        
                        if (fData.dml > 1e-10)
                        {
                            skippedDml++;
                            continue;
                        }
                        
                        processed++;

                        // Find nearest solid face
                        scalar minDist = GREAT;
                        label bestSolidIdx = -1;

                        forAll(globalSolidFaces, sIdx)
                        {
                            scalar dist = mag(globalSolidFaces[sIdx].faceCentre - fData.faceCentre);
                            
                            if (dist < minDist)
                            {
                                minDist = dist;
                                bestSolidIdx = sIdx;
                            }
                        }

                        // Reasonable tolerance (increased to 5e-3 = 5mm for testing)
                        if (bestSolidIdx != -1 && minDist < 5e-3) 
                        {
                            const SolidPacket& sData = globalSolidFaces[bestSolidIdx];

                            scalar yFluid = 2.0 * mag(fData.cellCentre - fData.faceCentre);
                            scalar ySolid = mag(sData.cellCentre - sData.faceCentre);

                            // scalar kdsolid = sData.ksolid / ySolid;//💡reCheck
                            scalar kdsolid = sData.alphasolid * sData.Cpsolid / ySolid;
                            // Use raw dml as per serial code, since we filtered dml <= 1e-10 above.
                            scalar localDml = fData.dml; 
                            scalar kdfluid = 1.0 / ((localDml / fData.kMax) + fData.Rint);

                            scalar Twall = (kdsolid * sData.Tsolid + kdfluid * fData.TSat)
                                        / (kdsolid + kdfluid);

                            scalar qml = kdfluid * (Twall - fData.TSat);

                            // Fluid side results
                            fRes.MSource = (qml / yFluid); //💡here MSource means sourceCoeffRef[fluidCellI]

                            // Solid side results //💡Keep the variable for return structure, but we do not use it here.
                            SolidResult sRes;
                            sRes.cellID = sData.cellID;
                            sRes.sourceCoeff = qml / (yFluid * sData.hsolid);

                            allSolidResults[sData.procID].append(sRes);
                        }
                    }
                }
            }
            
            if (debug) // Debug switch controlled by ChenUtakaCoeffs
            {
                Info << "ChenUtakaParallelR2 Coupling Stats: Total=" << totalFluidFaces 
                     << ", Skipped(Phase)=" << skippedPhase 
                     << ", Skipped(dml)=" << skippedDml 
                     << ", Processed=" << processed << endl;
            }
        }

        // 3. Scatter results back
        if (Pstream::parRun())
        {
            // Manual scatter implementation
            if (Pstream::master())
            {
                // Send to slaves
                for (int proci = 0; proci < Pstream::nProcs(); ++proci)
                {
                    if (proci != Pstream::masterNo())
                    {
                        OPstream toProc(Pstream::commsTypes::scheduled, proci);
                        toProc << allFluidResults[proci];
                        toProc << allSolidResults[proci];
                    }
                }

                // Copy to self
                fluidResults = allFluidResults[Pstream::masterNo()];
                solidResults = allSolidResults[Pstream::masterNo()];
            }
            else
            {
                // Receive from master
                IPstream fromMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
                fromMaster >> fluidResults;
                fromMaster >> solidResults;
            }
        }
        else
        {
            fluidResults = allFluidResults[0];
            solidResults = allSolidResults[0];
            return;
        }

    }
    } // End anonymous namespace
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ChenUtakaParallelR2::ChenUtakaParallelR2
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
    ChenUtakaCoeffs_(modelDict().subDict("ChenUtakaCoeffs")),
    evapCoeff_(ChenUtakaCoeffs_.lookupOrDefault<scalar>("evapCoeff",1)),
    Rgas_(ChenUtakaCoeffs_.lookupOrDefault<scalar>("Rgas",1)),
    fluidPatch_(ChenUtakaCoeffs_.lookupOrDefault<string>("fluidPatch","fluid_to_solid")), //Only working for 1 patch at the moment
    solidPatch_(ChenUtakaCoeffs_.lookupOrDefault<string>("solidPatch","solid_to_fluid")),
    method_(ChenUtakaCoeffs_.lookupOrDefault<string>("initialisationMethod","ChenUtakaParallel")), //Only working for 1 patch at the moment
    gradient_(ChenUtakaCoeffs_.lookupOrDefault<scalar>("gradient",1)),
    coefficient_(ChenUtakaCoeffs_.lookupOrDefault<scalar>("coefficient_",1)),
    origin_(ChenUtakaCoeffs_.lookupOrDefault<vector>("origin",vector(0, 0, 0))),
    debug_(ChenUtakaCoeffs_.lookupOrDefault<bool>("debug", false)),
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
    )

{
    if (phase2_.thermo().incompressible())  //Is this needed?
    {
        Rgas_.value() = ChenUtakaCoeffs_.get<scalar>("Rgas");
    }
}
// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * *  //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
// Foam::tmp<Foam::fvScalarMatrix> Foam::ChenUtakaParallelR2::dmlInitial()
// {
//
// }


Foam::tmp<Foam::fvScalarMatrix> Foam::ChenUtakaParallelR2::hSourceML()
{
    // if (debug_) Info<< "I'm inside hSourceML"<< endl;
    Info<< "I'm inside hSourceML"<< endl;
    //Fluid Patch
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();
    const label fluidPatchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);

    // Solid Patch
    const fvMesh& solidFvMesh = solid_.T().mesh();
    const label solidPatchID = solidFvMesh.boundaryMesh().findPatchID(solidPatch_);

    // Prepare Source Coeff Field
    volScalarField sourceCoeff
    (
        IOobject
        (
            "sourceCoeff",
            solidFvMesh.time().timeName(),
            solidFvMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidFvMesh,
        dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0),
        "zeroGradient"
    );

    // prepare Collect Local Data
    List<FluidPacket> localFluidPackets;
    List<SolidPacket> localSolidPackets;

    // Fluid Data
    if (fluidPatchID != -1) // local Fluid Patch
    {
        const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
        const labelUList& fluidFaceCells = fluidPatch.faceCells();

        const volScalarField& TSat = satModel_.TSat();
        const volScalarField& rho1 = phase1_.thermo().rho();
        const volScalarField& rho2 = phase2_.thermo().rho();
        const volScalarField& k1 = phase1_.kappa();

        const dimensionedScalar rho1Max ("rho1Max",rho1.dimensions(), gMax(rho1.internalField()));
        const dimensionedScalar Rgas ("Rgas", dimGasConstant, Rgas_.value());
        const dimensionedScalar kmax ("kFluidMax",k1.dimensions(), gMax(k1.internalField()));

        volScalarField Rint
        (
            (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*Rgas,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
        );

        localFluidPackets.setSize(fluidPatch.size());
        scalar maxLocalDml = 0.0;
        label countDmlPos = 0;
        
        forAll(fluidPatch, faceI) //prepare Fluid Packets
        {
            const label cellI = fluidFaceCells[faceI];
            localFluidPackets[faceI].procID = Pstream::myProcNo();
            localFluidPackets[faceI].cellID = cellI;
            localFluidPackets[faceI].faceCentre = fluidPatch.Cf()[faceI];
            localFluidPackets[faceI].cellCentre = mesh.C()[cellI];
            localFluidPackets[faceI].TSat = TSat[cellI];
            localFluidPackets[faceI].dml = dml_[cellI];
            localFluidPackets[faceI].Rint = Rint[cellI];
            localFluidPackets[faceI].kMax = kmax.value();
            localFluidPackets[faceI].rho1Max = rho1Max.value();
            localFluidPackets[faceI].phase1Val = phase1_[cellI];
            localFluidPackets[faceI].L = satModel_.L()[cellI];
            // mu1 not needed here but struct has it, fill 0 or whatever
            localFluidPackets[faceI].mu1 = 0.0;
            
            if (dml_[cellI] > 1e-10) {
                if (dml_[cellI] > maxLocalDml) maxLocalDml = dml_[cellI];
                countDmlPos++;
            }
        }
        
        if (debug_)
        {        
            if (maxLocalDml > 1e-10)
            {
            Pout << "ChenUtakaParallelR2::hSourceML: Local dml max=" << maxLocalDml 
                 << ", count > 1e-10 =" << countDmlPos << "/" << fluidPatch.size() << endl;
            }
        }
    }

    // Solid Data
    if (solidPatchID != -1)  // local Solid Patch
    {
        const fvPatch& solidPatch = solidFvMesh.boundary()[solidPatchID];
        const labelUList& solidFaceCells = solidPatch.faceCells();

        const volScalarField& Tsolid = solid_.T();
        const volScalarField& hsolid = solid_.he();
        const volScalarField& ksolid = solid_.kappa();
        const volScalarField& Cpsolid = solid_.Cp();
        const volScalarField& alphasolid = solid_.alpha();

        const dimensionedScalar ksolidMax ("kSolidMax",ksolid.dimensions(), gMax(ksolid.internalField()));
        const dimensionedScalar CpsolidMax ("CpSolidMax",Cpsolid.dimensions(), gMax(Cpsolid.internalField()));
        const dimensionedScalar alphasolidMax ("alphasolidMax",alphasolid.dimensions(), gMax(alphasolid.internalField()));

        localSolidPackets.setSize(solidPatch.size());
        forAll(solidPatch, faceI)
        {
            const label cellI = solidFaceCells[faceI];
            localSolidPackets[faceI].procID = Pstream::myProcNo();
            localSolidPackets[faceI].cellID = cellI;
            localSolidPackets[faceI].faceCentre = solidPatch.Cf()[faceI];
            localSolidPackets[faceI].cellCentre = solidFvMesh.C()[cellI];
            localSolidPackets[faceI].Tsolid = Tsolid[cellI];
            localSolidPackets[faceI].hsolid = hsolid[cellI];
            localSolidPackets[faceI].ksolid = ksolidMax.value();
            localSolidPackets[faceI].Cpsolid = CpsolidMax.value();
            localSolidPackets[faceI].alphasolid = alphasolidMax.value();
        }
    }
    else
    {
        if (Pstream::master())
        {
            WarningInFunction
                << "Solid patch '" << solidPatch_ << "' not found in solid mesh! "
                << "Parallel coupling will FAIL. Please check 'solidPatch' entry in dictionary."
                << endl;
        }
    }
    // prepare get results
    List<FluidResult> fluidResults;
    List<SolidResult> solidResults;

    performParallelCoupling(
        localFluidPackets,
        localSolidPackets,
        fluidResults,
        solidResults,
        phase1_.mesh().time().deltaTValue(),
        "hSourceML",
        debug_
    );

    // Apply Results
    // Only apply solid results in hSourceML
    forAll(solidResults, i)
    {
        sourceCoeff[solidResults[i].cellID] = solidResults[i].sourceCoeff;
    }

    sourceCoeff.correctBoundaryConditions();
    tmp<fvScalarMatrix> hSource(fvm::Sp(sourceCoeff, solid_.he()));

    return hSource;
}

void Foam::ChenUtakaParallelR2::initialiseML()
{
    // Support both naming conventions for compatibility
    if (method_ == "ChenUtakaParallelR2" || method_ == "chenUtaka")
    {
        Info<< "ChenUtakaParallelR2: Executing initialiseML with method " << method_ << endl;

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
        scalar dx (0.0);
        scalar dy (0.0);
        scalar dz (0.0);

        scalar maxInitDml = 0.0;

        forAll(faceCells, cellI)
        { 
            const label globalCellI = faceCells[cellI];
            point = phase1_.mesh().C()[globalCellI];
            dx =  pow(origin_[0] - point[0] , 2);
            dy =  pow(origin_[1] - point[1] , 2);
            dz =  pow(origin_[2] - point[2] , 2);

            distance = pow(dx + dy + dz, 0.5);

            scalar addedDml = gradient_*distance;
            dml_[globalCellI] += addedDml;
            
            if (addedDml > maxInitDml) maxInitDml = addedDml;
        }
        
        Info << "ChenUtakaParallelR2: Initialised dml. Local max added = " << maxInitDml * 1e6 << "um" << endl;
    }
    else
    {
        Info << "ChenUtakaParallelR2: Skipping initialiseML. Method '" << method_ << "' does not match 'ChenUtakaParallelR2' or 'chenUtaka'." << endl;
    }

}

void Foam::ChenUtakaParallelR2::updateML()
{
    // if (debug_) Info<< "I'm inside updateML"<< endl;
    Info<< "I'm inside updateML"<< endl;
    //Fluid Patch
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();
    const label fluidPatchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);

    // Solid Patch
    const fvMesh& solidFvMesh = solid_.T().mesh();
    const label solidPatchID = solidFvMesh.boundaryMesh().findPatchID(solidPatch_);

    // prepare Collect Data 
    List<FluidPacket> localFluidPackets;
    List<SolidPacket> localSolidPackets;

    // Fluid Data
    if (fluidPatchID != -1)
    {
        const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
        const labelUList& fluidFaceCells = fluidPatch.faceCells();

        const volScalarField& TSat = satModel_.TSat();
        const volScalarField& rho1 = phase1_.thermo().rho();
        const volScalarField& rho2 = phase2_.thermo().rho();
        const volScalarField& k1 = phase1_.kappa();
        const volScalarField& mu1 = phase1_.thermo().mu();

        const dimensionedScalar rho1Max ("rho1Max",rho1.dimensions(), gMax(rho1.internalField()));
        const dimensionedScalar Rgas ("Rgas", dimGasConstant, Rgas_.value());
        const dimensionedScalar kmax ("kFluidMax",k1.dimensions(), gMax(k1.internalField()));

        volScalarField Rint
        (
            (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*Rgas,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
        );

        localFluidPackets.setSize(fluidPatch.size());
        forAll(fluidPatch, faceI)
        {
            const label cellI = fluidFaceCells[faceI];
            localFluidPackets[faceI].procID = Pstream::myProcNo();
            localFluidPackets[faceI].cellID = cellI;
            localFluidPackets[faceI].faceCentre = fluidPatch.Cf()[faceI];
            localFluidPackets[faceI].cellCentre = mesh.C()[cellI];
            localFluidPackets[faceI].TSat = TSat[cellI];
            localFluidPackets[faceI].dml = dml_[cellI];
            localFluidPackets[faceI].Rint = Rint[cellI];
            localFluidPackets[faceI].kMax = kmax.value();
            localFluidPackets[faceI].rho1Max = rho1Max.value();
            localFluidPackets[faceI].phase1Val = phase1_[cellI];
            localFluidPackets[faceI].L = satModel_.L()[cellI];
            localFluidPackets[faceI].mu1 = mu1[cellI];
        }
    }

    // Solid Data
    if (solidPatchID != -1)
    {
        const fvPatch& solidPatch = solidFvMesh.boundary()[solidPatchID];
        const labelUList& solidFaceCells = solidPatch.faceCells();

        const volScalarField& Tsolid = solid_.T();
        const volScalarField& hsolid = solid_.he();
        const volScalarField& ksolid = solid_.kappa();
        const volScalarField& Cpsolid = solid_.Cp();
        const volScalarField& alphasolid = solid_.alpha();

        const dimensionedScalar ksolidMax ("kSolidMax",ksolid.dimensions(), gMax(ksolid.internalField()));
        const dimensionedScalar CpsolidMax ("CpSolidMax",Cpsolid.dimensions(), gMax(Cpsolid.internalField()));
        const dimensionedScalar alphasolidMax ("alphasolidMax",alphasolid.dimensions(), gMax(alphasolid.internalField()));

        localSolidPackets.setSize(solidPatch.size());
        forAll(solidPatch, faceI)
        {
            const label cellI = solidFaceCells[faceI];
            localSolidPackets[faceI].procID = Pstream::myProcNo();
            localSolidPackets[faceI].cellID = cellI;
            localSolidPackets[faceI].faceCentre = solidPatch.Cf()[faceI];
            localSolidPackets[faceI].cellCentre = solidFvMesh.C()[cellI];
            localSolidPackets[faceI].Tsolid = Tsolid[cellI];
            localSolidPackets[faceI].hsolid = hsolid[cellI];
            localSolidPackets[faceI].ksolid = ksolidMax.value();
            localSolidPackets[faceI].Cpsolid = CpsolidMax.value();
        }
    }
    // prepare get results
    List<FluidResult> fluidResults;
    List<SolidResult> solidResults;

    performParallelCoupling(
        localFluidPackets,
        localSolidPackets,
        fluidResults,
        solidResults,
        phase1_.mesh().time().deltaTValue(),
        "updateML",
        debug_
    );

    // Apply Fluid Results (dml update and MSource storage)
    // Reset MSource_ first
    MSource_ = dimensionedScalar("0", MSource_.dimensions(), 0.0);

    // volScalarField MSource(satModel_.TSat()*0.0); // Dummy for MSource -> REMOVE THIS

    if (fluidPatchID != -1)
    {
        const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
        
        // Ensure sizes match
        if (fluidResults.size() == fluidPatch.size())
        {
            forAll(fluidResults, i)
            {
                label faceI = i; 
                label cellI = fluidResults[i].cellID; 

                if (phase1_[cellI] < 0.5)
                {
                    if (dml_[cellI] >  1e-10)
                    {
                     dml_[cellI] = fluidResults[i].dmlNew;
                    //  MSource_[cellI] =  fluidResults[i].MSource; 
                     // Calculate Volumetric Mass Source (kg/m3/s)
                     // fluidResults[i].MSource is Flux (kg/m2/s)
                    //  scalar mSourceFlux = fluidResults[i].MSource; 
                    //  scalar faceArea = fluidPatch.magSf()[faceI];
                    //  scalar cellVol = mesh.V()[cellI];
                     
                    //  MSource_[cellI] = mSourceFlux * faceArea / cellVol;
                    }
                }
                else
                {
                     dml_[cellI] = 0.0;
                    //  MSource_[cellI] = 0.0;
                }
            }
        }
    }

    // Continue with smearing logic (Linear Extrapolation)
    // This part requires parallel reductions for xMin1, xMax1, etc.
    // to be correct globally.

    scalar xMax1(0.0); // Initialize to 0.0 to match Serial behavior and avoid huge dml if not found
    scalar xMin1(GREAT);
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
                 scalar xPos = mesh.C()[fluidCellI].x();
                 if (xPos < xMin1) xMin1 = xPos;
                 if (xPos > xMax1) 
                 {
                     xMax1 = xPos;
                     yMax1 = dml_[fluidCellI];
                 }
            }
        }
    }

    // Parallel Reduce
    reduce(xMin1, minOp<scalar>());
    // Reduce xMax1 and corresponding yMax1
    // We need to find which proc has xMax1 and what its yMax1 is.
    // Standard reduce(xMax1, maxOp) only gives us the global max X.
    
    // 1. Find global xMax1
    scalar globalXMax1 = xMax1;
    reduce(globalXMax1, maxOp<scalar>());

    // 2. Determine who has it.
    // If local xMax1 == globalXMax1 (within tolerance), we broadcast our yMax1
    
    scalar globalYMax1 = 0.0;
    
    // If I hold the max, I propose my yMax. Otherwise I propose -GREAT.
    
    scalar myProposedYMax = -GREAT;
    // Only propose if I actually found a valid max (xMax1 > -GREAT check implicit if initialized to 0.0 or valid)
    if (mag(xMax1 - globalXMax1) < SMALL)
    {
        myProposedYMax = yMax1; // Assuming yMax1 was set locally
    }
    
    reduce(myProposedYMax, maxOp<scalar>());
    globalYMax1 = myProposedYMax;
    
    // If globalYMax1 is still -GREAT (e.g. no fluid patch active?), set to 0
    if (globalYMax1 <= -GREAT/2) globalYMax1 = 0.0;

    xMax1 = globalXMax1;
    yMax1 = globalYMax1;
    scalar globalmidpoint = 0.0;
    globalmidpoint = (xMin1 + xMax1)/2.0;

    if (debug_)
    {
        Info << "xMin = " << xMin1 << ", xMid = " << globalmidpoint << ", xMax = " << xMax1 << endl;
    }

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
                 scalar point = phase1_.mesh().C()[fluidCellI].x();
                 scalar distance = mag(point - xMax1); // Use mag to be safe
                 
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
        Info << "ChenUtakaParallelR2 Debug: MSource_ max = " << gMax(MSource_.internalField()) 
             << ", min = " << gMin(MSource_.internalField()) << endl;
        Info << "ChenUtakaParallelR2 Debug: dml_ max = " << gMax(dml_) 
             << ", min = " << gMin(dml_) << endl;
    }
}

Foam::tmp<Foam::volScalarField> Foam::ChenUtakaParallelR2::energySourceML()
{
    // if (debug_) Info<< "I'm inside energySourceML"<< endl;
    Info<< "I'm inside energySourceML"<< endl;
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();
    const label fluidPatchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);

    const fvMesh& solidFvMesh = solid_.T().mesh();
    const label solidPatchID = solidFvMesh.boundaryMesh().findPatchID(solidPatch_);

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
            dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0) 
         )
    );
    volScalarField& energySource = tEnergySource.ref();

    // Calculate MSource_ for energy source

    // prepare Collect Local Data
    List<FluidPacket> localFluidPackets;
    List<SolidPacket> localSolidPackets;

    // Fluid Data
    if (fluidPatchID != -1) // local Fluid Patch
    {
        const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
        const labelUList& fluidFaceCells = fluidPatch.faceCells();

        const volScalarField& TSat = satModel_.TSat();
        const volScalarField& rho1 = phase1_.thermo().rho();
        const volScalarField& rho2 = phase2_.thermo().rho();
        const volScalarField& k1 = phase1_.kappa();

        const dimensionedScalar rho1Max ("rho1Max",rho1.dimensions(), gMax(rho1.internalField()));
        const dimensionedScalar Rgas ("Rgas", dimGasConstant, Rgas_.value());
        const dimensionedScalar kmax ("kFluidMax",k1.dimensions(), gMax(k1.internalField()));

        volScalarField Rint
        (
            (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*Rgas,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
        );

        localFluidPackets.setSize(fluidPatch.size());
        
        forAll(fluidPatch, faceI) //prepare Fluid Packets
        {
            const label cellI = fluidFaceCells[faceI];
            localFluidPackets[faceI].procID = Pstream::myProcNo();
            localFluidPackets[faceI].cellID = cellI;
            localFluidPackets[faceI].faceCentre = fluidPatch.Cf()[faceI];
            localFluidPackets[faceI].cellCentre = mesh.C()[cellI];
            localFluidPackets[faceI].TSat = TSat[cellI];
            localFluidPackets[faceI].dml = dml_[cellI];
            localFluidPackets[faceI].Rint = Rint[cellI];
            localFluidPackets[faceI].kMax = kmax.value();
            localFluidPackets[faceI].rho1Max = rho1Max.value();
            localFluidPackets[faceI].phase1Val = phase1_[cellI];
            localFluidPackets[faceI].L = satModel_.L()[cellI];
            // mu1 not needed here but struct has it, fill 0 or whatever
            localFluidPackets[faceI].mu1 = 0.0;
        }
    }

    // Solid Data
    if (solidPatchID != -1)  // local Solid Patch
    {
        const fvPatch& solidPatch = solidFvMesh.boundary()[solidPatchID];
        const labelUList& solidFaceCells = solidPatch.faceCells();

        const volScalarField& Tsolid = solid_.T();
        const volScalarField& hsolid = solid_.he();
        const volScalarField& ksolid = solid_.kappa();
        const volScalarField& Cpsolid = solid_.Cp();
        const volScalarField& alphasolid = solid_.alpha();

        const dimensionedScalar ksolidMax ("kSolidMax",ksolid.dimensions(), gMax(ksolid.internalField()));
        const dimensionedScalar CpsolidMax ("CpSolidMax",Cpsolid.dimensions(), gMax(Cpsolid.internalField()));
        const dimensionedScalar alphasolidMax ("alphasolidMax",alphasolid.dimensions(), gMax(alphasolid.internalField()));

        localSolidPackets.setSize(solidPatch.size());
        forAll(solidPatch, faceI)
        {
            const label cellI = solidFaceCells[faceI];
            localSolidPackets[faceI].procID = Pstream::myProcNo();
            localSolidPackets[faceI].cellID = cellI;
            localSolidPackets[faceI].faceCentre = solidPatch.Cf()[faceI];
            localSolidPackets[faceI].cellCentre = solidFvMesh.C()[cellI];
            localSolidPackets[faceI].Tsolid = Tsolid[cellI];
            localSolidPackets[faceI].hsolid = hsolid[cellI];
            localSolidPackets[faceI].ksolid = ksolidMax.value();
            localSolidPackets[faceI].Cpsolid = CpsolidMax.value();
            localSolidPackets[faceI].alphasolid = alphasolidMax.value();
        }
    }
    else
    {
        if (Pstream::master())
        {
            WarningInFunction
                << "Solid patch '" << solidPatch_ << "' not found in solid mesh! "
                << "Parallel coupling will FAIL. Please check 'solidPatch' entry in dictionary."
                << endl;
        }
    }
    // prepare get results
    List<FluidResult> fluidResults;
    List<SolidResult> solidResults;

    performParallelCoupling(
        localFluidPackets,
        localSolidPackets,
        fluidResults,
        solidResults,
        phase1_.mesh().time().deltaTValue(),
        "energySourceML",
        debug_
    );

    if (fluidPatchID != -1)
    {
        const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
        const labelUList& fluidFaceCells = fluidPatch.faceCells();
        const volScalarField& L = satModel_.L();
        
        forAll(fluidFaceCells, i)
        {
            label cellI = fluidFaceCells[i];
            // Check if cellI is valid (it should be, but just in case)
            if (cellI >= 0 && cellI < mesh.nCells())
            {
                energySource[cellI] = MSource_[cellI];
            }
        }
    }

    return tEnergySource;
}



Foam::tmp<Foam::volScalarField>
Foam::ChenUtakaParallelR2::massSourceML( volScalarField& rhoSource)
{
    // if (debug_) Info<< "I'm inside massSourceML"<< endl;
    Info<< "I'm inside massSourceML"<< endl;

    tmp<volScalarField> massSourceML(rhoSource * 0.0);
    volScalarField& massSourceMLRef = massSourceML.ref();

    const fvMesh& mesh = phase1_.mesh();

    dimensionedScalar DPsi
    (
        "DPsi",
        dimensionSet(0,2,0,0,0,0,0),
        3/sqr(gAverage(mesh.nonOrthDeltaCoeffs()))
    );

    dimensionedScalar intPsi0 = fvc::domainIntegrate(rhoSource);

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
        fvm::Sp(scalar(1),psi) - fvm::laplacian(DPsi,psi) == rhoSource
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
Foam::ChenUtakaParallelR2::alphaSourceML( volScalarField& rhoSource)
{
    Info<< "I'm inside alphaSourceML"<< endl;

    tmp<volScalarField> alphaSourceML
    (
        rhoSource / phase1_.thermo().rho()
    );

   return alphaSourceML;
}
