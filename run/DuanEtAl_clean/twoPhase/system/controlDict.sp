/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  plus                                  |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     myBastard;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         0.7;

deltaT         1e-10;

writeControl    adjustableRunTime;

writeInterval   0.01;

purgeWrite      2;

writeFormat     binary;

writePrecision  10;

writeCompression no;

timeFormat      general;

timePrecision   10;

runTimeModifiable yes;

adjustTimeStep  yes;

maxCo           0.1;
maxDi           1e8;
maxCapillaryNum 1e8;
maxAlphaCo      0.1;
maxDeltaT       1.0;

libs ("libpostProcess.so"
      "libdynamicLoadBalanceFvMesh.so"	

);

// ************************************************************************* //
functions
{

//----------------------------------------------------------------------------//
//- Function to sample the interface.
 /*       freeSurface
   {
       type            surfaces;
       libs            ("libsampling.so");
       surfaceFormat  raw;
       executeControl  timeStep;
       executeInterval 500;
       writeControl    timeStep;
       writeInterval 500;
       sampleOnExecute true;
       region fluid;
       fields
       (
           alpha.water
       );
       surfaces
       (
           freeSurface
           {
               type        isoSurfaceCell;
               isoField    alpha.water;
               isoValue    0.5;
               interpolate true;
               regularise  false;
           }
       );
       interpolationScheme cellPointFace;
   }
    */
    interfaceEnergyFluxes1
    {
    type        interfaceEnergyFluxes;
    libs        ("libfieldFunctionObjects.so");
    writeControl    timeStep;
    writeInterval   100;
    region fluid;
    patches     ("fluid_to_solid");
    };
    interfaceRegion1
    {
    type        interfaceRegion;
    libs        ("libfieldFunctionObjects.so");
    writeControl    timeStep;
    writeInterval   100;
    region fluid;
    nLayers 3;
    patches     ("fluid_to_solid");
    };
    // wallHeatFluxSolid1
    // {
    // type        wallHeatFluxSolid;
    // libs        ("libfieldFunctionObjects.so");
    // writeControl    timeStep;
    // writeInterval   100;
    // region fluid;
    // patches     ("fluid_to_solid");
    // };
    // superheated1
    // {
    // type        superHeated;
    // libs        ("libfieldFunctionObjects.so");
    // executeControl  timeStep;
    // executeInterval 100;
    // writeControl    timeStep;
    // writeInterval   100;
    // region fluid;
    // TSat    uniform 298.15;
    // patches     ("fluid_to_solid");
    // };

    bubbleIntegral
   {
       libs        ("libutilityFunctionObjects.so");
        type coded;
       // Name of on-the-fly generated functionObject
        name bubbleIntegral;
        executeControl  timeStep;
        executeInterval 100;
        writeControl    timeStep;
        writeInterval   100;
        region fluid;
       
       codeOptions
       #{
           -I$(LIB_SRC)/finiteVolume/lnInclude \
           -I$(LIB_SRC)/OpenFOAM/lnInclude
       #};
        codeInclude
       #{
           #include "volFieldsFwd.H"
           #include "OFstream.H"
           #include <iostream>
       #};
       
       codeData
       #{
             autoPtr<OFstream> outputFilePtr;
       #};
           
       codeRead
       #{
          
            outputFilePtr.reset(new OFstream("bubbleIntegral.dat"));
          
       #};
       
       codeExecute
       #{
       
            const volScalarField& alpha1 = mesh().lookupObject<volScalarField>("alpha.water");
            
            scalar x (0.0);
            scalar y (0.0);
            scalar V (0.0);
            scalar A (0.0);
            scalar intTot (0.0);
            scalar pi (3.1415926535897);

            scalar Deq (0.0);

            const faceList & ff = alpha1.mesh().faces();
            const pointField & pp = alpha1.mesh().points();
            
            forAll(alpha1,cellI)
            {
                y = alpha1.mesh().C()[cellI].y();
                if (y>0){        

                    V = alpha1.mesh().V()[cellI];
                    x = alpha1.mesh().C()[cellI].x();

                    const cell & cc = alpha1.mesh().cells()[cellI];
                    labelList pLabels(cc.labels(ff));
                    pointField pLocal(pLabels.size(), vector::zero);

                    forAll (pLabels, pointI)
                        pLocal[pointI] = pp[pLabels[pointI]];

                    scalar zDim = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));

                    A = (V/(zDim));
                    intTot += (1.0 - alpha1[cellI])*2*pi*x*A;
                        
                }
            }

            reduce(intTot, sumOp<scalar>());

            Deq =  std::cbrt((6*intTot)/(pi));
            //reduce(Deq, sumOp<scalar>());

                       
            if (Pstream::myProcNo() == 0)
            {
                outputFilePtr() << mesh().time().timeName() << " " << intTot << " " << Deq << endl;
            }
            
            Info << "intTot = "  << mesh().time().timeName() << " " << intTot << " " << Deq << endl;
            
       #};
       
       codeWrite
       #{   #};
   }

}
