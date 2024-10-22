if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if (oCorr == 0)
    {
        scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

        mesh.update();

        if (mesh.changing())
        {
            MRF.update();

            Info<< "Execution time for mesh.update() = "
                << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                << " s" << endl;
        }
    }



if(overwrite)
{
    continue;
}
{
    if (i == 0)
    {
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
        {
            fvScalarMatrix hEqn
            (
                fvm::ddt(betav*rho, h)
            - (
                thermo.isotropic()
                ? fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
                : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
                )
            ==
             fvOptions(rho, h)
            - phaseChangeSolid.hSourceML()

            );

            hEqn.relax();

            fvOptions.constrain(hEqn);

            hEqn.solve(mesh.solver(h.select(finalIter)));

            fvOptions.correct(h);
        }

        phaseChangeSolid.updateML();

        thermo.correct();

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;
    }

    else
    {    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
        {
            fvScalarMatrix hEqn
            (
                fvm::ddt(betav*rho, h)
            - (
                thermo.isotropic()
                ? fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
                : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
                )
            ==
                fvOptions(rho, h)
            );

            hEqn.relax();

            fvOptions.constrain(hEqn);

            hEqn.solve(mesh.solver(h.select(finalIter)));

            fvOptions.correct(h);
        }

        thermo.correct();

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;
    }
}


if (finalIter)
{
    mesh.data::remove("finalIteration");
}
