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


#include "ChenUtaka.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{ 
    defineTypeNameAndDebug(ChenUtaka, 0);
    addToRunTimeSelectionTable(microlayerModel,ChenUtaka, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ChenUtaka::ChenUtaka
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
    fluidPatch_(modelDict().lookupOrDefault<string>("fluidPatch","fluid_to_solid")), //Only working for 1 patch at the moment
    solidPatch_(modelDict().lookupOrDefault<string>("solidPatch","solid_to_fluid")),
    method_(modelDict().lookupOrDefault<string>("initialisationMethod","chenUtaka")), //Only working for 1 patch at the moment
    gradient_(modelDict().lookupOrDefault<scalar>("gradient",1)),
    coefficient_(modelDict().lookupOrDefault<scalar>("coefficient_",1)),
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
// Foam::tmp<Foam::fvScalarMatrix> Foam::ChenUtaka::dmlInitial()
// {

// }


Foam::tmp<Foam::fvScalarMatrix> Foam::ChenUtaka::hSourceML()
{
    Info<< "I'm inside hSourceML"<< endl;
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
    const fvMesh& solidFvMesh = refCast<const fvMesh>(solidMesh);
    const labelUList& solidFaceCells = solidPatch.faceCells();

    //Fluid Fields
    const volScalarField& TSat = satModel_.TSat(); 
    const volScalarField& rho1 = phase1_.thermo().rho();
    const dimensionedScalar rho1Max ("rho1Max",rho1.dimensions(), gMax(rho1.internalField()));
    const volScalarField& rho2 = phase2_.thermo().rho(); 
    const dimensionedScalar Rgas ("Rgas", dimGasConstant, Rgas_.value());
    const volScalarField& k1 = phase1_.kappa(); //Can change to kappaEff later (need alphat from turbulence model)
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

    // volScalarField sourceCoeff
    // (
    //     IOobject
    //     (
    //         "sourceCoeff",
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar("0", dimensionSet(1,-3,-1,0,0,0,0), 0.0),
    //     "zeroGradient"
    // );

    volScalarField MSource(TSat*0.0/(Rint*satModel_.L()));


    forAll(fluidPatch, faceI)
    {
        const scalar fluidCellI = fluidFaceCells[faceI];
        const scalar solidCellI = solidFaceCells[faceI];

        const dimensionedScalar dmlCell ("dmlCell", dimLength, dml_[fluidCellI]);

        if (phase1_[fluidCellI] < 0.001 and dmlCell.value() > 1e-10)
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

            // scalar kdh = kdfluid/yDimFluid;
            // Info<< " d/k " << (dml_[fluidCellI]/(kmax.value()))  <<endl;
            // Info<< " Rint " << Rint[fluidCellI]  <<endl;

            // sourceCoeff[solidCellI] = -(kdh.value()*kdsolid.value()*(1/CpsolidMax.value()))/(kdSum.value());
            // source1[solidCellI] = -kdh.value()*TSat[fluidCellI]*((kdfluid.value())/((kdSum.value())) - 1);

            // sourceCoeff[solidCellI] = (kdh.value())/((CpsolidMax.value()));
            // source1[solidCellI] = TSat[fluidCellI]*kdh.value();
            // scalar Twall = (kdsolid.value()*Tsolid[solidCellI]+kdfluid*TSat[fluidCellI])/(kdsolid.value()+kdfluid);
            scalar Twall = (kdsolid*Tsolid[solidCellI]+kdfluid*TSat[fluidCellI])/(kdsolid + kdfluid);

            // scalar Twall = Tsolid.boundaryField()[samplePatchi][faceI];
            scalar qml =  kdfluid*(Twall - TSat[fluidCellI]);
            sourceCoeff[solidCellI] =  qml/(yDimFluid*hsolid[solidCellI]);

        }
            
        else
        {
            sourceCoeff[solidCellI] = 0.0;
        }

    }
    
    // reduce(TSourceML, sumOp<scalar>());

    // tmp<volScalarField> energySource(TSource*(Twall-TSat));
    // volScalarField& energySourceRef = energySource.ref();
    // energySourceRef.ref() *= mag(surf_.normal().internalField())/TSat.mesh().V();

    sourceCoeff.correctBoundaryConditions();
    // sourceCoeff *= 1/solidFvMesh.V();

    tmp<fvScalarMatrix> hSource(fvm::Sp(sourceCoeff, hsolid));

    return hSource;

}


void Foam::ChenUtaka::initialiseML()
{   
    if (method_ == "chenUtaka")
    {
        Info<< "I'm inside initialiseML"<< endl;

        const label patchID = phase1_.mesh().boundaryMesh().findPatchID(fluidPatch_);
        const polyPatch& cPatch = phase1_.mesh().boundaryMesh()[patchID];
        const labelUList& faceCells = cPatch.faceCells();

        vector point (0, 0, 0);
        scalar distance (0.0);
        scalar dx (0.0);
        scalar dy (0.0);
        scalar dz (0.0);

        
        forAll(faceCells, cellI)
        { //Add condition for dml>0 - v

            point = phase1_.mesh().C()[cellI];
            dx =  pow(origin_[0] - point[0] , 2);
            dy =  pow(origin_[1] - point[1] , 2);
            dz =  pow(origin_[2] - point[2] , 2);

            distance = pow(dx + dy + dz, 0.5);

            dml_[cellI] += gradient_*distance;
        }
    }

}

void Foam::ChenUtaka::updateML()
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

    volScalarField Rint
    (
        (2-evapCoeff_)/(2*evapCoeff_)*(pow(2*constant::mathematical::pi*Rgas,0.5))*pow(TSat,1.5)/(pow(satModel_.L(),2)*rho2)
    );


    volScalarField MSource(TSat*0.0/(Rint*satModel_.L()));

    scalar xMax1(0.0); // Cell next to TPL
    scalar xMin1(1.0);

    scalar xMax2(0.0);
    scalar xMin2(0.0);
    scalar x1(0.0);
    scalar x2(0.0);
    scalar x3(0.0);
    scalar x4(0.0);

    scalar midpoint(0.0);

    scalar yMax1(0.0); // dml of xMax1
    scalar yMax2(0.0); // dml of xMax2

    scalar yMin1(0.0);
    scalar yMin2(0.0);

    scalar dx1 (30e-6);
    scalar dx2 (60e-6);

    forAll(fluidPatch, faceI)
    {
        const scalar fluidCellI = fluidFaceCells[faceI];
        const scalar solidCellI = solidFaceCells[faceI];

        //for linear extrapolation of microlayer


        if (phase1_[fluidCellI] < 0.5)
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
                dml_[fluidCellI] = 0.0;
            }
           
            if (mesh.C()[fluidCellI].x() < xMin1)
            {
                xMin1 = mesh.C()[fluidCellI].x();
                yMin1 = dml_[fluidCellI];
            }

            if ( mesh.C()[fluidCellI].x() > xMax1)
            {
                xMax1 = mesh.C()[fluidCellI].x();
                yMax1 = dml_[fluidCellI];
            }

            midpoint = (xMax1 - xMin1)/2;


            // if (mesh.C()[fluidCellI].x() < xMin2 and mesh.C()[fluidCellI].x() > xMin1)
            // {
            //     xMin2 = mesh.C()[fluidCellI].x();
            //     yMin2 = dml_[fluidCellI];
            // }

            // if (dml_[fluidCellI] > yMax1)
            // xMin1 = xMin2;
            // yMin1 = yMin2;

            // xMin2 = xMax2;
            // yMin2 = yMax2;

            // xMax2 = xMax1;
            // yMax2 = yMax1;




            // else if (mesh.C()[fluidCellI].x() > xMax1 - dx2 && mesh.C()[fluidCellI].x() <= xMax1 - dx1)
            // {
            //     xMin1 = max(mesh.C()[fluidCellI].x(),0);
            //     yMin1 = max(dml_[fluidCellI],0);
            // }

            // // Info << "xMin2 : " << xMin2 << endl;
            // Info << "xMin1 : " << xMin1 << endl;


        }

            
    }

    Info << "Min" << xMin1 << endl;
    Info << "Mid" << midpoint << endl;
    Info << "Max" << xMax1 << endl;

    forAll(fluidPatch, faceI)
    {
        const scalar fluidCellI = fluidFaceCells[faceI];

        // // if (mesh.C()[fluidCellI].x() > xMax1 - dx2 && mesh.C()[fluidCellI].x() <= xMax1 - dx1)
        // if (mesh.C()[fluidCellI].x() > x1 && mesh.C()[fluidCellI].x() <= xMax1 - midpoint)
        // {
        //     x1 = mesh.C()[fluidCellI].x();
        //     y1 = max(dml_[fluidCellI],0);
        // }
        // Info << "xMin1 : " << xMin1 << endl;


        if (phase1_[fluidCellI] > 0.5)
        {
            scalar point = phase1_.mesh().C()[fluidCellI].x();

            // dml_[fluidCellI] = yMax1 + (yMax1 - yMax2)*(point - xMax1)/(xMax1 - xMax2);
            // Info<< "xMax1 - xMin1 : "<< xMax1 - xMin1 << endl;
            // dml_[fluidCellI] = max(yMax1 + (yMax1 - yMin1)*(point - xMax1)/(xMax1 - xMin1), yMax1);
            // dml_[fluidCellI] = yMin2 + (yMin1 - yMin2)*(point - xMin2)/(xMin1 - xMin2);

            scalar distance = point - xMax1;
            Info << "distance : " << distance << endl;

            dml_[fluidCellI] = gradient_*distance + yMax1;
        }
        // Info << "xMin1 : " << xMin1 << endl;

        if (method_ == "cooperLloyd" && phase1_[fluidCellI] >= 0.5)
        {
            dml_[fluidCellI] = coefficient_*pow((mu1[fluidCellI]/rho1[fluidCellI])*phase1_.mesh().time().value() ,0.5);  //
        }
    }


    
}

Foam::tmp<Foam::volScalarField> Foam::ChenUtaka::energySourceML()
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

    // const dimensionedVector solidCellC0 = Tsolid.mesh().C()[solidFaceCells[0]];
    // const dimensionedVector solidCellC1 = Tsolid.mesh().C()[solidFaceCells[1]];
    // // const dimensionedVector faceC = fluidPatch.faceCentres()[0];    
    // dimensionedScalar d_solid("dSolid", dimLength, 0.5*Foam::mag(solidCellC0 - solidCellC1).value());
    // Info<< "d_solid"<<d_solid.dimensions() <<endl;
    // Info<< "alphasolidMax"<<alphasolidMax.dimensions() <<alphasolidMax.value()<< endl;

    // volScalarField kdsolid (alphasolidMax.value()*ksolid/d_solid);
    // volScalarField kdsolid
    // (
    //     IOobject
    //     (
    //         "kdsolid",
    //         solidMesh.time().timeName(),
    //         solidMesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     solidMesh,
    //     dimensionedScalar("0", dimensionSet(1,0,-3,-1,0,0,0), alphasolidMax.value()*ksolidMax.value()/d_solid.value()),
    //     "zeroGradient"
    // );

    // Info<< "Made fields"<< endl;
    // Pout << "Hello from processor " << Pstream::myProcNo() << "! I am working on "<< fluidPatch.size() << " cells" << endl;

    // forAll(fluidPatch, faceI)
    // {
    //     const scalar fluidCellI = fluidFaceCells[faceI];

    //     if (phase1_[fluidCellI] < 0.001 and dml_[fluidCellI] > 1e-10)
    //     {
    //         const dimensionedVector faceC = fluidPatch.faceCentres()[faceI];
    //         const dimensionedVector fluidCellC = mesh.C()[fluidCellI];
    //         yDimFluid.value() = 2*Foam::mag(fluidCellC - faceC).value(); 
    //         kdfluid[fluidCellI] = kmax.value()/(dml_[fluidCellI]);
    //         kdh[fluidCellI] = kdfluid[fluidCellI]/yDimFluid.value();
    //     }
    // } 
    // Info<< "minkdfluid " <<gMin(kdfluid.internalField())<< endl;
    // Info<< "minkdh " <<gMin(kdh.internalField())<< endl;
    // Info<< "minkdsolid " <<gMin(kdsolid.internalField())<< endl;
    // Info<< "minsum " <<gMin(kdsolid.internalField()) + gMin(kdfluid.internalField())<< endl;

    // volScalarField Twall((kdsolid*Tsolid + kdfluid*TSat)/(kdsolid + kdfluid));
    // Info<< "minTwall " << gMin(Twall.internalField()) << endl;
    // sourceCoeffRef =  kdh*(Twall - TSat);


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
Foam::ChenUtaka::massSourceML( volScalarField& rhoSource)
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
Foam::ChenUtaka::alphaSourceML( volScalarField& rhoSource)
{
    Info<< "I'm inside alphaSourceML"<< endl;

    tmp<volScalarField> alphaSourceML
    (
        rhoSource / phase1_.thermo().rho() 
    );

   return alphaSourceML;
}
