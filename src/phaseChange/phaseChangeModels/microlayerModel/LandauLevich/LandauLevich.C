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


#include "LandauLevich.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"

#include "cubicEqn.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{ 
    defineTypeNameAndDebug(LandauLevich, 0);
    addToRunTimeSelectionTable(microlayerModel,LandauLevich, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LandauLevich::LandauLevich
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
    fluidPatch_(modelDict().lookupOrDefault<string>("fluidPatch","fluid_to_solid")), //Only working for 1 patch at the moment
    solidPatch_(modelDict().lookupOrDefault<string>("solidPatch","solid_to_fluid")),
    Jakob_(modelDict().lookupOrDefault<scalar>("Jakob",1)),
    prevT_(modelDict().lookupOrDefault<scalar>("prevT",0)),
    prevR_(modelDict().lookupOrDefault<scalar>("Rinit",0)),
    prevRbase_(modelDict().lookupOrDefault<scalar>("RbaseInit",0)),
    prevdRdt_(modelDict().lookupOrDefault<scalar>("prevdRdt",0)),
    prevCell_(modelDict().lookupOrDefault<label>("prevCell",0)),
    origin_(modelDict().lookupOrDefault<vector>("origin",vector(0, 0, 0)))

{
    if (phase2_.thermo().incompressible())  //Is this needed?
    {
        Rgas_.value() = modelDict().get<scalar>("Rgas");
    }

}
// * * * * * * * * * * * * * * Protected Access Member Functions  * * * * *  //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //

Foam::tmp<Foam::fvScalarMatrix> Foam::LandauLevich::hSourceML()
{
    Info<< "I'm inside hSourceML"<< endl;
    //Fluid Patch
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh(); 
    const label fluidPatchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);
    const fvPatch& fluidPatch = mesh.boundary()[fluidPatchID];

    //Solid Patch
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>(fluidPatch.patch());
    const polyMesh& solidMesh = mpp.sampleMesh();
    const fvMesh& solidFvMesh = refCast<const fvMesh>(solidMesh);
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& solidPatch = solidFvMesh.boundary()[samplePatchi];

    // const polyPatch& fluidPatch = phase1_.mesh().boundaryMesh()[fluidPatchID];
    // const labelUList& fluidFaceCells = fluidPatch.faceCells();
    // const labelUList& solidFaceCells = solidPatch.faceCells();



    //Fluid Fields
    const volScalarField& TSat = satModel_.TSat(); 
    const volScalarField& rho1 = phase1_.thermo().rho();
    const volScalarField& rho2 = phase2_.thermo().rho(); 
    const volScalarField& k1 = phase1_.kappa(); //Can change to kappaEff later (need alphat from turbulence model)

    const dimensionedScalar rho1Max ("rho1Max",rho1.dimensions(), gMax(rho1.internalField()));
    const dimensionedScalar Rgas ("Rgas", dimGasConstant, Rgas_.value());
    const dimensionedScalar kmax ("kFluidMax",k1.dimensions(), gMax(k1.internalField()));

    //Solid Fields
    const volScalarField& Tsolid = solid_.T(); 
    const volScalarField& hsolid = solid_.he();  
    const volScalarField& ksolid = solid_.kappa(); 
    const volScalarField& Cpsolid = solid_.Cp(); 
    const volScalarField& rhosolid = solid_.rho(); 
    const volScalarField& alphasolid = solid_.alpha(); 

    const dimensionedScalar ksolidMax ("kSolidMax",ksolid.dimensions(), gMax(ksolid.internalField()));
    const dimensionedScalar CpsolidMax ("CpSolidMax",Cpsolid.dimensions(), gMax(Cpsolid.internalField()));
    const dimensionedScalar rhosolidMax ("rhosolidMax",rhosolid.dimensions(), gMax(rhosolid.internalField()));
    const dimensionedScalar alphasolidMax ("alphasolidMax",alphasolid.dimensions(), gMax(alphasolid.internalField()));

    volScalarField Rint
    (
        (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*Rgas,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
    );

    volScalarField sourceCoeff
    (
        IOobject
        (
            "sourceCoeff",
            solidMesh.time().timeName(),
            solidFvMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidFvMesh,
        dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0),
        "zeroGradient"
    );

    // volScalarField MSource(TSat*0.0/(Rint*satModel_.L()));

    const labelList& fluidFaceCells = fluidPatch.faceCells();
    const labelList& solidFaceCells = solidPatch.faceCells();
    const vectorField& fluidFaceCentres = fluidPatch.Cf();

    forAll(fluidPatch, faceI)
    {
        const scalar fluidCellI = fluidFaceCells[faceI];
        const scalar solidCellI = solidFaceCells[faceI];

        // Ensure both cells are valid local cells (avoid ghost cells)
        if (fluidCellI < 0 || solidCellI < 0) continue;

        const dimensionedScalar dmlCell ("dmlCell", dimLength, dml_[fluidCellI]);
        const scalar dmlVal = dml_[fluidCellI];

        // Skip inactive cells or zero phase region
        if (phase1_[fluidCellI] >= 0.001 || dmlVal <= 1e-10)
            continue;


        const vector& faceC = fluidFaceCentres[faceI];
        const vector& fluidCellC = mesh.C()[fluidCellI];
        const vector& solidCellC = solidFvMesh.C()[solidCellI];
    
        const scalar yFluid = 2 * mag(fluidCellC - faceC);
        const scalar ySolid = mag(solidCellC - faceC);
    
        const scalar kdSolid = alphasolidMax.value() * CpsolidMax.value() / ySolid;
        const scalar kdFluid = 1.0 / ((dmlVal / kmax.value()) + Rint[fluidCellI]);
    
        const scalar Twall = (kdSolid * Tsolid[solidCellI] + kdFluid * TSat[fluidCellI])
                             / (kdSolid + kdFluid);
    
        const scalar qml = kdFluid * (Twall - TSat[fluidCellI]);
    
        sourceCoeff[solidCellI] = qml / (yFluid * hsolid[solidCellI]);

    }
    
    sourceCoeff.correctBoundaryConditions();

    tmp<fvScalarMatrix> hSource(fvm::Sp(sourceCoeff, hsolid));

    return hSource;

}

void Foam::LandauLevich::initialiseML()
{   
    Info<< "I'm inside initialiseML"<< endl;

    const label patchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);
    const polyPatch& cPatch = phase1_.mesh().boundaryMesh()[patchID];
    const labelUList& faceCells = cPatch.faceCells();

    //Fluid Fields
    const volScalarField& rho1 = phase1_.thermo().rho();
    const dimensionedScalar rho1Max ("rho1Max",rho1.dimensions(), gMax(rho1.internalField()));
    const volScalarField& k1 = phase1_.kappa(); //Can change to kappaEff later (need alphat from turbulence model)
    const dimensionedScalar kMax ("kFluidMax",k1.dimensions(), gMax(k1.internalField()));
    const volScalarField& mu1 = phase1_.thermo().mu();
    const dimensionedScalar mu1Max ("mu1Max",mu1.dimensions(), gMax(mu1.internalField()));
    const volScalarField& Cp1 = phase1_.thermo().Cp();
    const dimensionedScalar Cp1Max ("Cp1Max",Cp1.dimensions(), gMax(Cp1.internalField()));

    const dimensionedScalar alphaDiff = kMax/(Cp1Max*rho1Max);

    scalar maxR = 0.0;
    scalar maxRCell = 0;
    scalar maxdRdt = 0.0;
    forAll(faceCells, cellI)
    {
        scalar R = mag(phase1_.mesh().C()[cellI] - origin_);

        if (R <= prevR_ && phase1_[cellI] < 0.1)
        {
            scalar t = pow(R/(2 * Jakob_ * pow((3*alphaDiff.value()/constant::mathematical::pi),0.5)), 2);
            scalar dRdt = pow((3*alphaDiff.value()/constant::mathematical::pi),0.5) * Jakob_ * pow(t,-0.5);
            scalar d2Rdt2 = -0.5*pow((3*alphaDiff.value()/constant::mathematical::pi),0.5) * Jakob_ * pow(t,-1.5);

            if (R > maxR)
            {
                maxR = R;
                maxRCell = cellI;
                maxdRdt = dRdt;
            }

            Info << "\n[initialisedBubbleRadius] R = " << R
            << ", dRdt = " << dRdt 
            << ", d2Rdt2 = " << d2Rdt2 
            << ", t = " << t << endl;

            const scalar a = (rho1Max.value()/sigma_)*(pow(dRdt/R,2) - (d2Rdt2/(3*R))); 
            const scalar b = -((rho1Max.value()*d2Rdt2)/(2*sigma_)) ;
            const scalar c = 1/(R);

            Roots<3> roots(cubicEqn(a, b, c, -1).roots());
            scalar x_bar = 0;
            for (label i = 0; i < 3; ++i) 
            {
                scalar root = roots[i];
                if (root > 0)
                {
                    x_bar = root;
                    break;
                }  
            }

            scalar Rm = 1/(3*a*pow(x_bar,2) + 2*b*x_bar + c); // dimensioned scalar could be more appropriate here
            dml_[cellI] = 1.34*Rm*pow((mu1Max.value()*dRdt/sigma_),(2.0/3.0)); 
        }

    }
    prevR_ = maxR; 
    prevRbase_ = maxR; 
    prevdRdt_ = maxdRdt;
    prevCell_ = maxRCell;
}

void Foam::LandauLevich::updateML()
{       
    Info<< "I'm inside updateML"<< endl;

    //Fluid Patch
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh(); 
    const label fluidPatchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);
    const polyPatch& fluidPatch = phase1_.mesh().boundaryMesh()[fluidPatchID];
    const labelUList& fluidFaceCells = fluidPatch.faceCells();

    //Solid Patch
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>(fluidPatch);
    const polyMesh& solidMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& solidPatch = refCast<const fvMesh>(solidMesh).boundary()[samplePatchi];
    const labelUList& solidFaceCells = solidPatch.faceCells();

    //Solid Fields
    const volScalarField& Tsolid = solid_.T(); 
    const volScalarField& ksolid = solid_.kappa(); 
    const volScalarField& Cpsolid = solid_.Cp(); 

    const volScalarField& alphasolid = solid_.alpha(); 

    const dimensionedScalar ksolidMax ("kSolidMax",ksolid.dimensions(), gMax(ksolid.internalField()));
    const dimensionedScalar CpsolidMax ("CpSolidMax",Cpsolid.dimensions(), gMax(Cpsolid.internalField()));
    const dimensionedScalar alphasolidMax ("alphasolidMax",alphasolid.dimensions(), gMax(alphasolid.internalField()));

    //Fluid Fields
    const volScalarField& TSat = satModel_.TSat(); 
    const volScalarField& rho1 = phase1_.thermo().rho();
    const dimensionedScalar rho1Max ("rho1Max",rho1.dimensions(), gMax(rho1.internalField()));
    const volScalarField& rho2 = phase2_.thermo().rho(); 
    // const volScalarField& T2 = phase2_.thermo().T(); 
    const dimensionedScalar Rgas ("Rgas", dimGasConstant, Rgas_.value());
    const volScalarField& k1 = phase1_.kappa(); //Can change to kappaEff later (need alphat from turbulence model)
    const dimensionedScalar kmax ("kFluidMax",k1.dimensions(), gMax(k1.internalField()));
    const volScalarField& mu1 = phase1_.thermo().mu();
    const dimensionedScalar mu1max ("mu1Max",mu1.dimensions(), gMax(mu1.internalField()));

    volScalarField Rint
    (
        (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*Rgas,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
    );

    forAll(fluidPatch, faceI)
    {
        const scalar fluidCellI = fluidFaceCells[faceI];
        if (phase1_[fluidCellI] < 0.1)
        {
            scalar Rbase = mag(phase1_.mesh().C()[fluidCellI] - origin_);  
            // Check if the radius is greater than the previous base radius
            if (Rbase > prevRbase_)
            {
                // Check if the cell is the same as the previous one
                label triggeredCell = fluidCellI;
                if (triggeredCell != prevCell_)
                {
                    dimensionedScalar bubbleV (72*fvc::domainIntegrate(1.0 - phase1_)); //wedge of 5 degrees, 360/5 = 72
                    Info << "\n bubbleV = " << bubbleV.value() << endl;
                    scalar R = (pow((6.0 * bubbleV.value()) / (constant::mathematical::pi), 1.0 / 3.0))/2;
                    scalar currentTime = mesh.time().timeOutputValue();
                    scalar dt = currentTime - prevT_;
                    scalar dRdt = (R - prevR_) / dt;
                    scalar d2Rdt2 = (dRdt - prevdRdt_) / dt;

                    Info << "\n[calculateBubbleRadius] R = " << R
                        << ", dRdt = " << dRdt
                        << ", d2Rdt2 = " << d2Rdt2 
                        << ", dt = " << dt << endl;

                    // update old rates for future timestep
                    prevRbase_ = Rbase;
                    prevR_ = R;
                    prevT_ = currentTime;
                    prevdRdt_ = dRdt;
                    prevCell_ = triggeredCell;

                    // Info << "updated old values" << endl;

                    const scalar a = (rho1Max.value()/sigma_)*(pow(dRdt/R,2) - (d2Rdt2/(3*R)));   
                    const scalar b = -((rho1Max.value()*d2Rdt2)/(2*sigma_)) ;
                    const scalar c = 1/(R);

                    Info << "coefficients a: " << a << "  : " << b << "  c: " << c << endl;

                    // Info << "found coefficients" << endl;

                    Roots<3> roots(cubicEqn(a, b, c, -1).roots());
                    Info << "found roots" << roots << endl;
                    scalar x_bar = 0;
                    for (label i = 0; i < 3; ++i) 
                    {
                        scalar root = roots[i];
                        if (root > 0)
                        {
                            x_bar = root;
                            break;
                        }  
                    }
                    Info << "found x_bar " << x_bar << endl;

                    scalar Rm = 1/(3*a*pow(x_bar,2) + 2*b*x_bar + c); // dimensioned scalar could be more appropriate here
                    Info << "found Rm " << Rm*1000 << " [mm] " << endl;

                    dml_[fluidCellI] = 1.34*Rm*pow((mu1max.value()*dRdt/sigma_),(2.0/3.0)); 
                    Info << "updated dml " << dml_[fluidCellI] << endl;
                }
                }
            }
        }
    




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    volScalarField MSource(TSat*0.0/(Rint*satModel_.L()));

    forAll(fluidPatch, faceI)
    {
        const scalar fluidCellI = fluidFaceCells[faceI];
        const scalar solidCellI = solidFaceCells[faceI];

        //for linear extrapolation of microlayer


        if (phase1_[fluidCellI] < 0.1)
        {
            if (dml_[fluidCellI] >  1e-10)
            {
                const vector solidCellC =   Tsolid.mesh().C()[solidCellI];
                const vector faceC = fluidPatch.faceCentres()[faceI];
                scalar yDimSolid = Foam::mag(solidCellC - faceC);
                // Info<< "yDimSolid  " << yDimSolid << endl;
                scalar kdsolid = alphasolidMax.value()*CpsolidMax.value()/yDimSolid;
                // Info<< "kdsolid  " << kdsolid << endl;

                const vector fluidCellC =   mesh.C()[fluidCellI];
                scalar yDimFluid = 2*Foam::mag(fluidCellC - faceC); 
                // Info<< "yDimFluid  " << yDimFluid << endl;

                // scalar kdfluid = kmax.value()/(dml_[fluidCellI]); 
                scalar kdfluid =   1/((dml_[fluidCellI]/(kmax.value())) + Rint[fluidCellI]);
                // Info<< "kmax  " << kmax.value() << endl;
                // Info<< "kdfluid  " << kdfluid << endl;

                scalar Twall = (kdsolid*Tsolid[solidCellI]+kdfluid*TSat[fluidCellI])/(kdsolid + kdfluid);

                MSource[fluidCellI] = ((Twall - TSat[fluidCellI])*kdfluid)/satModel_.L()[fluidCellI]; 
                // MSource[fluidCellI] = ((Tsolid.boundaryField()[samplePatchi][faceI] - TSat[fluidCellI])*kdfluid.value())/satModel_.L()[fluidCellI]; 

                dml_[fluidCellI] = dml_[fluidCellI] - MSource[fluidCellI]*phase1_.mesh().time().deltaTValue()/rho1Max.value();
            }
            else
            {
                // If dml is too small, set it to zero to avoid negative values
                dml_[fluidCellI] = 0.0;
                MSource[fluidCellI] = 0.0;
            }
        }
    }
}

Foam::tmp<Foam::volScalarField> Foam::LandauLevich::energySourceML()
{
    Info<< "I'm inside energySourceML"<< endl;
    //Fluid Patch
    surf_.reconstruct(false);
    const fvMesh& mesh = phase1_.mesh(); 
    const label fluidPatchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);
    const polyPatch& fluidPatch = phase1_.mesh().boundaryMesh()[fluidPatchID];

    // const List< polyPatch> allPatchList(Pstream::nProcs());
    // // if (Pstream::master()){ v = something(); } // <- must do on master

    // // Pstream::gather(fluidPatch, polyPatch); // <- root process gathers

    // allPatchList[Pstream::myProcNo()] =& phase1_.mesh().boundaryMesh()[fluidPatchID];
    // Pstream::gatherList(allPatchList);

    const labelUList& fluidFaceCells = fluidPatch.faceCells();
    // const List<label> fluidFaceCells = fluidPatch.faceCells();
    // Pstream::gather(fluidPatch);
    //Solid Patch
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(fluidPatch);
    const polyMesh& solidMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& solidPatch = refCast<const fvMesh>(solidMesh).boundary()[samplePatchi];
    const labelUList& solidFaceCells = solidPatch.faceCells();
    
    // Pstream::gatherList(solidFaceCells);
    //Solid Fields
    const volScalarField& Tsolid = solid_.T(); 
    const volScalarField& ksolid = solid_.kappa(); 
    const volScalarField& Cpsolid = solid_.Cp(); 
    const volScalarField& rhosolid = solid_.rho(); 
    const volScalarField& alphasolid = solid_.alpha(); 

    // const dimensionedScalar ksolidMax ("kSolidMax",ksolid.dimensions(), gMax(ksolid.internalField()));
    // const dimensionedScalar CpsolidMax ("CpSolidMax",Cpsolid.dimensions(), gMax(Cpsolid.internalField()));
    // const dimensionedScalar rhosolidMax ("rhosolidMax",rhosolid.dimensions(), gMax(rhosolid.internalField()));
    // const dimensionedScalar alphasolidMax ("alphasolidMax",alphasolid.dimensions(), gMax(alphasolid.internalField()));
    // returnReduce(max(ksolid.internalField()).value(), maxOp<scalar>());

    const dimensionedScalar ksolidMax ("kSolidMax",ksolid.dimensions(), returnReduce(gMax(ksolid.internalField()), maxOp<scalar>()));
    const dimensionedScalar CpsolidMax ("CpSolidMax",Cpsolid.dimensions(), returnReduce(gMax(Cpsolid.internalField()), maxOp<scalar>()));
    const dimensionedScalar rhosolidMax ("rhosolidMax",rhosolid.dimensions(), returnReduce(gMax(rhosolid.internalField()), maxOp<scalar>()));
    const dimensionedScalar alphasolidMax ("alphasolidMax",alphasolid.dimensions(), returnReduce(gMax(alphasolid.internalField()), maxOp<scalar>()));


    //Fluid Fields
    const volScalarField& TSat = satModel_.TSat(); 
    const volScalarField& rho1 = phase1_.thermo().rho();
    const dimensionedScalar rho1Max ("rho1Max",rho1.dimensions(), gMax(rho1.internalField()));
    const volScalarField& rho2 = phase2_.thermo().rho(); 
    const dimensionedScalar Rgas ("Rgas", dimGasConstant, Rgas_.value());
    const volScalarField& k1 = phase1_.kappa(); 
    const dimensionedScalar kmax ("kFluidMax",k1.dimensions(), gMax(k1.internalField()));

    // Info<< "Got all fields"<< endl;

    volScalarField Rint
    (
        (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*Rgas,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
    );

    const dimensionedScalar dummyLength ("dummyLength",dimLength, 1.0);
    // dimensionedScalar yDimFluid = ("yDimFluid", dimLength, 1.0); 
    // volScalarField kdfluid(phase1_.kappa()*0.0/dummyLength);

    // volScalarField kdh(kdfluid/dummyLength);    

    tmp<volScalarField> sourceCoeff(TSat*0.0/(Rint*dummyLength));    
    volScalarField& sourceCoeffRef = sourceCoeff.ref();    
    
    volScalarField MSource(TSat*0.0/(Rint*satModel_.L()));    


    forAll(fluidPatch, faceI)
    {
        const scalar fluidCellI = fluidFaceCells[faceI];
        const scalar solidCellI = solidFaceCells[faceI];

        if (phase1_[fluidCellI] < 0.001 and dml_[fluidCellI] > 1e-10)
        {
            // Pout<< "if"<< endl;
            const vector solidCellC =   Tsolid.mesh().C()[solidCellI];
            const vector faceC = fluidPatch.faceCentres()[faceI];
            scalar yDimSolid = Foam::mag(solidCellC - faceC);
            // Info<< "yDimSolid  " << yDimSolid << endl;
            scalar kdsolid = alphasolidMax.value()*CpsolidMax.value()/yDimSolid;
            // Info<< "kdsolid  " << kdsolid << endl;

            const vector fluidCellC =   mesh.C()[fluidCellI];
            scalar yDimFluid = 2*Foam::mag(fluidCellC - faceC); 
            // Info<< "yDimFluid  " << yDimFluid << endl;

            // scalar kdfluid = kmax.value()/(dml_[fluidCellI]); 
            scalar kdfluid =   1/((dml_[fluidCellI]/(kmax.value())) + Rint[fluidCellI]);
            // Info<< "Rint  " <<  Rint[fluidCellI] << endl;
            // Info<< "kmax  " << kmax.value() << endl;
            // Info<< "kdfluid  " << kdfluid << endl;

            scalar Twall = (kdsolid*Tsolid[solidCellI]+kdfluid*TSat[fluidCellI])/(kdsolid + kdfluid);
            // Info<< "TSat  " << TSat[fluidCellI] << endl;
            // Info<< "Twall  " << Twall << endl;
            // Info<< "Tsolid  " << Tsolid[solidCellI] << endl;

            scalar qml =  kdfluid*(Twall - TSat[fluidCellI]);
            // Info<< "qml  " << qml << endl;

            sourceCoeffRef[fluidCellI] =  qml/(yDimFluid);
            // Info<< "sourceCoeffRef[fluidCellI]  " << sourceCoeffRef[fluidCellI] << endl;
            // sourceCoeffRef[fluidCellI] =  qml;

        }

        else
        {
            // Pout<< "else"<< endl;
            sourceCoeffRef[fluidCellI] = 0.0;
            // dmlCell[fluidCellI] = 0.0;
        }
    } 
    

    // Pout<< "exited for loop"<< endl;
    // reduce(sourceCoeff, sumOp<volScalarField>());

    // sourceCoeffRef.ref() *= mag(surf_.normal().internalField())/TSat.mesh().V();
    sourceCoeffRef.correctBoundaryConditions();
    // Pout<< "just before return"<< endl;
    // Info<< "minSourceCoeff " <<gMin(sourceCoeffRef.internalField())<< endl;
    return sourceCoeff;

}



Foam::tmp<Foam::volScalarField> 
Foam::LandauLevich::massSourceML( volScalarField& rhoSource)
{
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
Foam::LandauLevich::alphaSourceML( volScalarField& rhoSource)
{
    Info<< "I'm inside alphaSourceML"<< endl;

    tmp<volScalarField> alphaSourceML
    (
        rhoSource / phase1_.thermo().rho() 
    );

   return alphaSourceML;
}
