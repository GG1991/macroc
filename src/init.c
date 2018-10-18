/*
 *  This source code is part of MacroC: a finite element code
 *  to solve macrostructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#include "macroc.h"


PetscErrorCode init()
{
    PetscErrorCode ierr;
    char mess[64];

    int rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    sprintf(mess, "Problem size %d\n", nproc);
    print0(mess);

    final_time = FINAL_TIME;
    ts = TIME_STEPS;
    dt = DT;

    PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL);
    PetscOptionsGetInt(NULL, NULL, "-ts", &ts, NULL);

    DMBoundaryType bx = DM_BOUNDARY_NONE, by = DM_BOUNDARY_NONE,
                   bz = DM_BOUNDARY_NONE;
    ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_BOX,
                        NX, NY, NZ,
                        PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                        DIM, 1, NULL, NULL, NULL, &DA);
    ierr = DMDASetUniformCoordinates(DA, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    CHKERRQ(ierr);
    ierr = DMSetMatType(DA, MATAIJ); CHKERRQ(ierr);
    ierr = DMSetFromOptions(DA); CHKERRQ(ierr);
    ierr = DMSetUp(DA); CHKERRQ(ierr);
    ierr = DMCreateMatrix(DA, &A); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(DA, &u); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(DA, &b); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(DA, &du); CHKERRQ(ierr);

    ierr = VecZeroEntries(u); CHKERRQ(ierr);
    ierr = VecZeroEntries(b); CHKERRQ(ierr);
    ierr = VecZeroEntries(du); CHKERRQ(ierr);

    dx = 1.;
    dy = 1.;
    dz = 1.;
    wg = dx * dy * dz / NPE;

    // Initializes <materials> declared in <micropp_c_wrapper.h>
    micropp_C_material_create();
    micropp_C_material_set(0, 1.0e7, 0.25, 1.0e4, 1.0e7, 0);
    micropp_C_material_set(1, 1.0e7, 0.25, 1.0e4, 1.0e7, 1);
    print0("Material Values:\n");

    if(!rank) {
        micropp_C_material_print(0);
        micropp_C_material_print(1);
    }

    PetscInt nex, ney, nez;
    ierr = DMDAGetElementsSizes(DA, &nex, &ney, &nez); CHKERRQ(ierr);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,
                            "rank:%d\tne:%d\tnex:%d\tney:%d\tnez:%d\n",
                            rank, (int) (nex * ney * nez), (int) nex,
                            (int) ney, (int) nez);
    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

    // Initializes <micro> declared in <micropp_c_wrapper.h>
    int ngpl = nex * ney * nez * NGP;
    int size[3] = { 5, 5, 5 };
    int type = 1;
    double params[4] = { 1., 1., 1., .5 };
    micropp_C_create3(ngpl, size, type, params);

    return ierr;
}


PetscErrorCode finish()
{
    PetscErrorCode ierr;
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = VecDestroy(&u); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = VecDestroy(&du); CHKERRQ(ierr);
    ierr = PetscFinalize();
    return ierr;
}
