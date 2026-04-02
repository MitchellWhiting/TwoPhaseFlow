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


#include "LandauLevichParallel.H"
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

#include "cubicEqn.H"

namespace Foam
{
    defineTypeNameAndDebug(LandauLevichParallel, 0);
    addToRunTimeSelectionTable(microlayerModel,LandauLevichParallel, components);

    void Foam::LandauLevichParallel::updateMapper()
    {
        const label fluidPatchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);
        if (fluidPatchID == -1)
        {
            return;
        }

        const polyPatch& pp = phase1_.mesh().boundaryMesh()[fluidPatchID];

        if (mapperInitialized_)
        {
            if (pp.size() == fluidPatchSize_)
            {
                return;
            }
            mapperInitialized_ = false;
        }

        const word& sampleRegion = solid_.T().mesh().name();
        const word& samplePatch = solidPatch_;

        mapper_.reset
        (
            new mappedPatchBase
            (
                pp,
                sampleRegion,
                mappedPatchBase::NEARESTPATCHFACE,
                samplePatch,
                0.0
            )
        );

        fluidPatchSize_ = pp.size();
        mapperInitialized_ = true;
    }
}

Foam::LandauLevichParallel::LandauLevichParallel
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
    evapCoeff_(modelDict().lookupOrDefault<scalar>("evapCoeff",1)),
    Rgas_(modelDict().lookupOrDefault<scalar>("Rgas",1)),
    sigma_(modelDict().lookupOrDefault<scalar>("sigma",1)),
    fluidPatch_(modelDict().lookupOrDefault<string>("fluidPatch","fluid_to_solid")),
    solidPatch_(modelDict().lookupOrDefault<string>("solidPatch","solid_to_fluid")),
    Jakob_(modelDict().lookupOrDefault<scalar>("Jakob",1)),
    prevT_(modelDict().lookupOrDefault<scalar>("prevT",0)),
    prevR_(modelDict().lookupOrDefault<scalar>("Rinit",0)),
    prevRbase_(modelDict().lookupOrDefault<scalar>("RbaseInit",0)),
    prevdRdt_(modelDict().lookupOrDefault<scalar>("prevdRdt",0)),
    prevCell_(modelDict().lookupOrDefault<label>("prevCell",0)),
    origin_(modelDict().lookupOrDefault<vector>("origin",vector(0, 0, 0))),
    debug_(modelDict().lookupOrDefault<bool>("debug", false)),
    mapper_(nullptr),
    mapperInitialized_(false),
    fluidPatchSize_(-1)
{
    if (phase2_.thermo().incompressible())
    {
        Rgas_.value() = modelDict().get<scalar>("Rgas");
    }
}

Foam::LandauLevichParallel::~LandauLevichParallel()
{
    if (mapper_.valid())
    {
        mapper_.clear();
    }
}

Foam::tmp<Foam::fvScalarMatrix> Foam::LandauLevichParallel::hSourceML()
{
    if (debug_) Info<< "LandauLevichParallel: Entering hSourceML" << endl;
    updateMapper();

    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();
    const label fluidPatchID = mesh.boundaryMesh().findPatchID(fluidPatch_);
    const fvMesh& solidFvMesh = solid_.T().mesh();
    const label solidPatchID = solidFvMesh.boundaryMesh().findPatchID(solidPatch_);

    if (fluidPatchID == -1 || solidPatchID == -1 || !mapper_.valid())
    {
        volScalarField sourceCoeff
        (
            IOobject
            (
                "sourceCoeff",
                solidFvMesh.time().timeName(),
                solidFvMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidFvMesh,
            dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0),
            "zeroGradient"
        );

        return tmp<fvScalarMatrix>(fvm::Sp(sourceCoeff, solid_.he()));
    }

    const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
    const fvPatch& solidPatch = solidFvMesh.boundary()[solidPatchID];
    const labelUList& fluidFaceCells = fluidPatch.faceCells();
    const labelUList& solidFaceCells = solidPatch.faceCells();

    scalarField TsolidFld = solid_.T().boundaryField()[solidPatchID];
    scalarField hsolidFld = solid_.he().boundaryField()[solidPatchID];
    scalarField kSolidFld = solid_.kappa()().boundaryField()[solidPatchID];

    scalarField ySolidFld(solidPatch.size());
    forAll(solidPatch, i)
    {
        ySolidFld[i] = mag(solidFvMesh.C()[solidFaceCells[i]] - solidPatch.Cf()[i]);
    }

    scalarField Tsolid_on_Fluid = TsolidFld;
    mapper_->distribute(Tsolid_on_Fluid);
    scalarField hsolid_on_Fluid = hsolidFld;
    mapper_->distribute(hsolid_on_Fluid);
    scalarField kSolid_on_Fluid = kSolidFld;
    mapper_->distribute(kSolid_on_Fluid);
    scalarField ySolid_on_Fluid = ySolidFld;
    mapper_->distribute(ySolid_on_Fluid);

    const volScalarField& TSat = satModel_.TSat();
    const volScalarField& rho2 = phase2_.thermo().rho();
    const volScalarField& k1 = phase1_.kappa();
    scalar kmax = gMax(k1.internalField());
    const dimensionedScalar Rgas ("Rgas", dimGasConstant, Rgas_.value());

    volScalarField Rint
    (
        (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*Rgas,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
    );

    scalarField sourceCoeff_on_Fluid(fluidPatch.size(), 0.0);

    forAll(fluidPatch, faceI)
    {
        const label fluidCellI = fluidFaceCells[faceI];

        if (phase1_[fluidCellI] >= 0.001 || dml_[fluidCellI] <= 1e-10)
        {
            continue;
        }

        scalar kdsolid = kSolid_on_Fluid[faceI]/(ySolid_on_Fluid[faceI] + SMALL);
        scalar kdfluid = 1.0/((dml_[fluidCellI]/kmax) + Rint.boundaryField()[fluidPatchID][faceI]);
        scalar TSatVal = TSat.boundaryField()[fluidPatchID][faceI];
        scalar Twall = (kdsolid*Tsolid_on_Fluid[faceI] + kdfluid*TSatVal)/(kdsolid + kdfluid);
        scalar qml = kdfluid*(Twall - TSatVal);
        scalar yFluid = 2*mag(mesh.C()[fluidCellI] - fluidPatch.Cf()[faceI]);

        sourceCoeff_on_Fluid[faceI] = qml/(yFluid*(hsolid_on_Fluid[faceI] + SMALL));
    }

    scalarField sourceCoeff_on_Solid = sourceCoeff_on_Fluid;
    mapper_->reverseDistribute(sourceCoeff_on_Solid);

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

    forAll(solidPatch, i)
    {
        const label cellI = solidFaceCells[i];
        sourceCoeff[cellI] = sourceCoeff_on_Solid[i];
    }

    sourceCoeff.correctBoundaryConditions();

    return tmp<fvScalarMatrix>(fvm::Sp(sourceCoeff, solid_.he()));
}

void Foam::LandauLevichParallel::initialiseML()
{
    if (debug_) Info<< "LandauLevichParallel: Entering initialiseML" << endl;

    const label patchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);
    if (patchID == -1)
    {
        return;
    }

    const polyPatch& cPatch = phase1_.mesh().boundaryMesh()[patchID];
    const labelUList& faceCells = cPatch.faceCells();

    const volScalarField& rho1 = phase1_.thermo().rho();
    const dimensionedScalar rho1Max ("rho1Max",rho1.dimensions(), gMax(rho1.internalField()));
    const volScalarField& k1 = phase1_.kappa();
    const dimensionedScalar kMax ("kFluidMax",k1.dimensions(), gMax(k1.internalField()));
    const volScalarField& mu1 = phase1_.thermo().mu();
    const dimensionedScalar mu1Max ("mu1Max",mu1.dimensions(), gMax(mu1.internalField()));
    const volScalarField& Cp1 = phase1_.thermo().Cp();
    const dimensionedScalar Cp1Max ("Cp1Max",Cp1.dimensions(), gMax(Cp1.internalField()));
    const dimensionedScalar alphaDiff = kMax/(Cp1Max*rho1Max);

    scalar localMaxR = -GREAT;
    scalar localMaxdRdt = 0.0;
    label localMaxRCell = -1;

    forAll(faceCells, i)
    {
        const label cellI = faceCells[i];
        scalar R = mag(phase1_.mesh().C()[cellI] - origin_);

        if (R <= prevR_ && phase1_[cellI] < 0.1)
        {
            scalar t = pow(R/(2*Jakob_*pow((3*alphaDiff.value()/constant::mathematical::pi),0.5)), 2);
            scalar dRdt = pow((3*alphaDiff.value()/constant::mathematical::pi),0.5)*Jakob_*pow(t,-0.5);
            scalar d2Rdt2 = -0.5*pow((3*alphaDiff.value()/constant::mathematical::pi),0.5)*Jakob_*pow(t,-1.5);

            if (R > localMaxR)
            {
                localMaxR = R;
                localMaxRCell = cellI;
                localMaxdRdt = dRdt;
            }

            const scalar a = (rho1Max.value()/sigma_)*(pow(dRdt/R,2) - (d2Rdt2/(3*R)));
            const scalar b = -((rho1Max.value()*d2Rdt2)/(2*sigma_));
            const scalar c = 1/(R);

            Roots<3> roots(cubicEqn(a, b, c, -1).roots());
            scalar x_bar = 0;
            for (label ir = 0; ir < 3; ++ir)
            {
                scalar root = roots[ir];
                if (root > 0)
                {
                    x_bar = root;
                    break;
                }
            }

            scalar Rm = 1/(3*a*pow(x_bar,2) + 2*b*x_bar + c);
            dml_[cellI] = 1.34*Rm*pow((mu1Max.value()*dRdt/sigma_),(2.0/3.0));
        }
    }

    scalar globalMaxR = localMaxR;
    reduce(globalMaxR, maxOp<scalar>());

    scalar globalMaxdRdt = -GREAT;
    label globalMaxRCell = -1;
    if (mag(localMaxR - globalMaxR) < SMALL)
    {
        globalMaxdRdt = localMaxdRdt;
        globalMaxRCell = localMaxRCell;
    }
    reduce(globalMaxdRdt, maxOp<scalar>());
    reduce(globalMaxRCell, maxOp<label>());

    if (globalMaxR > -GREAT/2)
    {
        prevR_ = globalMaxR;
        prevRbase_ = globalMaxR;
        prevdRdt_ = globalMaxdRdt;
        prevCell_ = globalMaxRCell;
    }

    if (debug_ && Pstream::master())
    {
        Info << "LandauLevichParallel initialiseML: prevR=" << prevR_
             << " prevRbase=" << prevRbase_
             << " prevdRdt=" << prevdRdt_ << endl;
    }
}

void Foam::LandauLevichParallel::updateML()
{
    if (debug_) Info<< "LandauLevichParallel: Entering updateML" << endl;
    updateMapper();

    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();
    const label fluidPatchID = mesh.boundaryMesh().findPatchID(fluidPatch_);
    const fvMesh& solidFvMesh = solid_.T().mesh();
    const label solidPatchID = solidFvMesh.boundaryMesh().findPatchID(solidPatch_);

    if (fluidPatchID == -1 || solidPatchID == -1 || !mapper_.valid())
    {
        return;
    }

    const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
    const fvPatch& solidPatch = solidFvMesh.boundary()[solidPatchID];
    const labelUList& fluidFaceCells = fluidPatch.faceCells();
    const labelUList& solidFaceCells = solidPatch.faceCells();

    scalarField TsolidFld = solid_.T().boundaryField()[solidPatchID];
    scalarField kSolidFld = solid_.kappa()().boundaryField()[solidPatchID];
    scalarField ySolidFld(solidPatch.size());
    forAll(solidPatch, i)
    {
        ySolidFld[i] = mag(solidFvMesh.C()[solidFaceCells[i]] - solidPatch.Cf()[i]);
    }

    scalarField Tsolid_on_Fluid = TsolidFld;
    mapper_->distribute(Tsolid_on_Fluid);
    scalarField kSolid_on_Fluid = kSolidFld;
    mapper_->distribute(kSolid_on_Fluid);
    scalarField ySolid_on_Fluid = ySolidFld;
    mapper_->distribute(ySolid_on_Fluid);

    const volScalarField& TSat = satModel_.TSat();
    const volScalarField& rho1 = phase1_.thermo().rho();
    const volScalarField& rho2 = phase2_.thermo().rho();
    const volScalarField& k1 = phase1_.kappa();
    const volScalarField& mu1 = phase1_.thermo().mu();

    const dimensionedScalar rho1Max ("rho1Max",rho1.dimensions(), gMax(rho1.internalField()));
    const dimensionedScalar kmax ("kFluidMax",k1.dimensions(), gMax(k1.internalField()));
    const dimensionedScalar mu1max ("mu1Max",mu1.dimensions(), gMax(mu1.internalField()));
    const dimensionedScalar Rgas ("Rgas", dimGasConstant, Rgas_.value());

    volScalarField Rint
    (
        (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*Rgas,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
    );

    scalar localMaxRbase = -GREAT;
    forAll(fluidPatch, faceI)
    {
        const label fluidCellI = fluidFaceCells[faceI];
        if (phase1_[fluidCellI] < 0.1)
        {
            const scalar Rbase = mag(phase1_.mesh().C()[fluidCellI] - origin_);
            if (Rbase > localMaxRbase)
            {
                localMaxRbase = Rbase;
            }
        }
    }

    scalar globalMaxRbase = localMaxRbase;
    reduce(globalMaxRbase, maxOp<scalar>());

    bool updatedRadius = false;
    if (globalMaxRbase > prevRbase_ + SMALL)
    {
        const dimensionedScalar bubbleV(72*fvc::domainIntegrate(1.0 - phase1_));
        const scalar R = (pow((6.0*bubbleV.value())/(constant::mathematical::pi), 1.0/3.0))/2;
        const scalar currentTime = mesh.time().timeOutputValue();
        const scalar dt = currentTime - prevT_;

        if (dt > SMALL && R > SMALL)
        {
            const scalar dRdt = (R - prevR_)/dt;
            const scalar d2Rdt2 = (dRdt - prevdRdt_)/dt;

            const scalar a = (rho1Max.value()/sigma_)*(pow(dRdt/R,2) - (d2Rdt2/(3*R)));
            const scalar b = -((rho1Max.value()*d2Rdt2)/(2*sigma_));
            const scalar c = 1/R;

            Roots<3> roots(cubicEqn(a, b, c, -1).roots());
            scalar x_bar = 0;
            for (label ir = 0; ir < 3; ++ir)
            {
                const scalar root = roots[ir];
                if (root > 0)
                {
                    x_bar = root;
                    break;
                }
            }

            if (x_bar > SMALL)
            {
                const scalar Rm = 1/(3*a*pow(x_bar,2) + 2*b*x_bar + c);
                const scalar dmlNew = max
                (
                    scalar(0.0),
                    1.34*Rm*pow((mu1max.value()*dRdt/sigma_),(2.0/3.0))
                );

                forAll(fluidPatch, faceI)
                {
                    const label fluidCellI = fluidFaceCells[faceI];
                    if (phase1_[fluidCellI] < 0.1)
                    {
                        const scalar Rbase = mag(phase1_.mesh().C()[fluidCellI] - origin_);
                        if (mag(Rbase - globalMaxRbase) < 1e-10)
                        {
                            dml_[fluidCellI] = dmlNew;
                        }
                    }
                }

                prevRbase_ = globalMaxRbase;
                prevR_ = R;
                prevT_ = currentTime;
                prevdRdt_ = dRdt;
                prevCell_ = -1;
                updatedRadius = true;
            }
        }
    }

    volScalarField MSource(TSat*0.0/(Rint*satModel_.L()));

    label updatedCount = 0;
    label skippedPhaseCount = 0;
    label skippedDmlCount = 0;
    bool debugPrinted = false;

    forAll(fluidPatch, faceI)
    {
        const label fluidCellI = fluidFaceCells[faceI];

        if (phase1_[fluidCellI] < 0.1)
        {
            if (dml_[fluidCellI] > 1e-10)
            {
                scalar kdsolid = kSolid_on_Fluid[faceI]/(ySolid_on_Fluid[faceI] + SMALL);
                scalar kdfluid = 1/((dml_[fluidCellI]/(kmax.value())) + Rint.boundaryField()[fluidPatchID][faceI]);
                scalar TSatVal = TSat.boundaryField()[fluidPatchID][faceI];
                scalar Twall = (kdsolid*Tsolid_on_Fluid[faceI] + kdfluid*TSatVal)/(kdsolid + kdfluid);

                MSource[fluidCellI] = ((Twall - TSatVal)*kdfluid)/satModel_.L()[fluidCellI];
                dml_[fluidCellI] = max
                (
                    scalar(0.0),
                    dml_[fluidCellI] - MSource[fluidCellI]*phase1_.mesh().time().deltaTValue()/rho1Max.value()
                );
                updatedCount++;

                if (debug_ && !debugPrinted && updatedCount == 1)
                {
                    Info << "LandauLevichParallel DEBUG Single Face:" << endl
                         << "  FaceI: " << faceI << " dml: " << dml_[fluidCellI] << endl
                         << "  Tsolid: " << Tsolid_on_Fluid[faceI]
                         << " TSat: " << TSatVal
                         << " Twall: " << Twall << endl
                         << "  kdsolid: " << kdsolid
                         << " kdfluid: " << kdfluid
                         << " MSource: " << MSource[fluidCellI] << endl;
                    debugPrinted = true;
                }
            }
            else
            {
                dml_[fluidCellI] = 0.0;
                MSource[fluidCellI] = 0.0;
                skippedDmlCount++;
            }
        }
        else
        {
            skippedPhaseCount++;
            dml_[fluidCellI] = 0.0;
            MSource[fluidCellI] = 0.0;
        }
    }

    if (debug_)
    {
        label totalUpdated = updatedCount;
        reduce(totalUpdated, sumOp<label>());
        scalar globalMaxDml = gMax(dml_);
        scalar globalMinDml = gMin(dml_);
        scalar globalMaxMSource = gMax(MSource.internalField());
        scalar globalMinMSource = gMin(MSource.internalField());

        if (Pstream::master())
        {
            Info << "LandauLevichParallel updateML: Total Updated dml/MSource on "
                 << totalUpdated << " cells." << endl;
            Info << "LandauLevichParallel updateML: skippedPhase=" << skippedPhaseCount
                 << " skippedDml=" << skippedDmlCount
                 << " radiusEvent=" << (updatedRadius ? 1 : 0) << endl;
            Info << "LandauLevichParallel Debug: MSource max = " << globalMaxMSource
                 << ", min = " << globalMinMSource << endl;
            Info << "LandauLevichParallel Debug: dml max = " << globalMaxDml
                 << ", min = " << globalMinDml << endl;
        }
    }
}

Foam::tmp<Foam::volScalarField> Foam::LandauLevichParallel::energySourceML()
{
    if (debug_) Info<< "LandauLevichParallel: Entering energySourceML" << endl;
    updateMapper();

    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh();
    const label fluidPatchID = mesh.boundaryMesh().findPatchID(fluidPatch_);
    const fvMesh& solidFvMesh = solid_.T().mesh();
    const label solidPatchID = solidFvMesh.boundaryMesh().findPatchID(solidPatch_);

    if (fluidPatchID == -1 || solidPatchID == -1 || !mapper_.valid())
    {
        dimensionedScalar dummyLength ("dummyLength",dimLength,1.0);
        const volScalarField& TSat = satModel_.TSat();
        const volScalarField& rho2 = phase2_.thermo().rho();
        const dimensionedScalar Rgas ("Rgas", dimGasConstant, Rgas_.value());
        volScalarField Rint
        (
            (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*Rgas,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
        );
        return tmp<volScalarField>(TSat*0.0/(Rint*dummyLength));
    }

    const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];
    const fvPatch& solidPatch = solidFvMesh.boundary()[solidPatchID];
    const labelUList& fluidFaceCells = fluidPatch.faceCells();
    const labelUList& solidFaceCells = solidPatch.faceCells();

    scalarField TsolidFld = solid_.T().boundaryField()[solidPatchID];
    scalarField kSolidFld = solid_.kappa()().boundaryField()[solidPatchID];
    scalarField ySolidFld(solidPatch.size());
    forAll(solidPatch, i)
    {
        ySolidFld[i] = mag(solidFvMesh.C()[solidFaceCells[i]] - solidPatch.Cf()[i]);
    }

    scalarField Tsolid_on_Fluid = TsolidFld;
    mapper_->distribute(Tsolid_on_Fluid);
    scalarField kSolid_on_Fluid = kSolidFld;
    mapper_->distribute(kSolid_on_Fluid);
    scalarField ySolid_on_Fluid = ySolidFld;
    mapper_->distribute(ySolid_on_Fluid);

    const volScalarField& TSat = satModel_.TSat();
    const volScalarField& rho2 = phase2_.thermo().rho();
    const volScalarField& k1 = phase1_.kappa();
    const dimensionedScalar kmax ("kFluidMax",k1.dimensions(), gMax(k1.internalField()));
    const dimensionedScalar Rgas ("Rgas", dimGasConstant, Rgas_.value());

    volScalarField Rint
    (
        (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*Rgas,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
    );

    const dimensionedScalar dummyLength ("dummyLength",dimLength, 1.0);
    tmp<volScalarField> sourceCoeff(TSat*0.0/(Rint*dummyLength));
    volScalarField& sourceCoeffRef = sourceCoeff.ref();

    forAll(fluidPatch, faceI)
    {
        const label fluidCellI = fluidFaceCells[faceI];

        if (phase1_[fluidCellI] < 0.001 && dml_[fluidCellI] > 1e-10)
        {
            scalar kdsolid = kSolid_on_Fluid[faceI]/(ySolid_on_Fluid[faceI] + SMALL);
            scalar kdfluid = 1/((dml_[fluidCellI]/(kmax.value())) + Rint.boundaryField()[fluidPatchID][faceI]);
            scalar TSatVal = TSat.boundaryField()[fluidPatchID][faceI];
            scalar Twall = (kdsolid*Tsolid_on_Fluid[faceI] + kdfluid*TSatVal)/(kdsolid + kdfluid);
            scalar qml = kdfluid*(Twall - TSatVal);
            scalar yDimFluid = 2*mag(mesh.C()[fluidCellI] - fluidPatch.Cf()[faceI]);
            sourceCoeffRef[fluidCellI] = qml/yDimFluid;
        }
        else
        {
            sourceCoeffRef[fluidCellI] = 0.0;
        }
    }

    sourceCoeffRef.correctBoundaryConditions();
    return sourceCoeff;
}

Foam::tmp<Foam::volScalarField>
Foam::LandauLevichParallel::massSourceML(volScalarField& rhoSource)
{
    if (debug_) Info<< "LandauLevichParallel: Entering massSourceML" << endl;

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

    fvScalarMatrix psiEqn
    (
        fvm::Sp(scalar(1),psi) - fvm::laplacian(DPsi,psi) == rhoSource
    );

    psiEqn.solve();

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
            intPsiVapor.value() += (1.0-phase1_[celli])*psi[celli]*mesh.V()[celli];
        }
    }

    dimensionedScalar Nv ("Nv", dimless, 2.0);
    reduce(intPsiVapor.value(),sumOp<scalar>());

    if (intPsiVapor.value() > 1e-99)
    {
        Nv = intPsi0/intPsiVapor;
    }

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
Foam::LandauLevichParallel::alphaSourceML(volScalarField& rhoSource)
{
    if (debug_) Info<< "LandauLevichParallel: Entering alphaSourceML" << endl;

    tmp<volScalarField> alphaSourceML
    (
        rhoSource / phase1_.thermo().rho()
    );

    return alphaSourceML;
}
