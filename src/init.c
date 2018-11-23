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
	PetscInt micro_n = 6;
	PetscInt micro_type = 1;
	PetscReal micro_mat_1[4] = { 1.0e7, 0.25, 1.0e4, 1.0e7 };
	PetscReal micro_mat_2[4] = { 1.0e9, 0.25, 1.0e4, 1.0e7 };

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	/*
	 * MIC_SPHERE = 0
	 * MIC_LAYER_Y = 1
	 * MIC_CILI_FIB_Z = 3
	 * MIC_CILI_FIB_XZ = 4
	 * MIC_QUAD_FIB_XYZ = 5
	 * MIC_QUAD_FIB_XZ = 6
	 * MIC_QUAD_FIB_XZ_BROKEN_X = 7
	 */

	final_time = FINAL_TIME;
	ts = TIME_STEPS;
	dt = DT;

	lx = LX;
	ly = LY;
	lz = LZ;
	NX = NX_CONST;
	NY = NY_CONST;
	NZ = NZ_CONST;

	/* Read command line options */

	vtu_freq = VTU_FREQ;
	newton_max_its = NEWTON_MAX_ITS;
	newton_min_tol = NEWTON_MIN_TOL;
	bc_type = BC_BENDING;

	PetscOptionsGetReal(NULL, NULL, "-dt", &dt, NULL);
	PetscOptionsGetReal(NULL, NULL, "-lx", &lx, NULL);
	PetscOptionsGetReal(NULL, NULL, "-ly", &ly, NULL);
	PetscOptionsGetReal(NULL, NULL, "-lz", &lz, NULL);
	PetscOptionsGetReal(NULL, NULL, "-new_tol", &newton_min_tol, NULL);
	PetscOptionsGetInt(NULL, NULL, "-ts", &ts, NULL);
	PetscOptionsGetInt(NULL, NULL, "-vtu_freq", &vtu_freq, NULL);
	PetscOptionsGetInt(NULL, NULL, "-new_its", &newton_max_its, NULL);
	PetscOptionsGetInt(NULL, NULL, "-bc_type", &bc_type, NULL);

	PetscInt nmax = 4;
	PetscOptionsGetInt(NULL, NULL, "-micro_n", &micro_n, NULL);
	PetscOptionsGetInt(NULL, NULL, "-micro_type", &micro_type, NULL);
	PetscOptionsGetRealArray(NULL, NULL, "-micro_mat_1", micro_mat_1, &nmax, NULL);
	PetscOptionsGetRealArray(NULL, NULL, "-micro_mat_2", micro_mat_2, &nmax, NULL);

	ierr = DMDACreate3d(PETSC_COMM_WORLD,
			    DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
			    DMDA_STENCIL_BOX,
			    NX, NY, NZ,
			    PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			    DIM, 1, NULL, NULL, NULL, &da);

	ierr = DMSetMatType(da, MATAIJ); CHKERRQ(ierr);
	ierr = DMSetFromOptions(da); CHKERRQ(ierr);
	ierr = DMSetUp(da); CHKERRQ(ierr);
	ierr = DMCreateMatrix(da, &A); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da, &u); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da, &b); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da, &du); CHKERRQ(ierr);

	ierr = VecSetOption(u, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
	ierr = VecSetOption(b, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

	ierr = VecZeroEntries(u); CHKERRQ(ierr);
	ierr = VecZeroEntries(b); CHKERRQ(ierr);
	ierr = VecZeroEntries(du); CHKERRQ(ierr);

	ierr = DMDAGetInfo(da, 0, &NX, &NY, &NZ,
			   &NP_X, &NP_Y, &NP_Y, 0, 0, 0, 0, 0, 0); CHKERRQ(ierr);

	PetscPrintf(PETSC_COMM_WORLD, "Boundary Condition : ");
	switch (bc_type) {
		case (BC_BENDING):
			PetscPrintf(PETSC_COMM_WORLD, "BC_BENDING\n");
			break;
		case (BC_CIRCLE):
			PetscPrintf(PETSC_COMM_WORLD, "BC_CIRCLE\n");
			break;
		default:
			break;
	}

	PetscPrintf(PETSC_COMM_WORLD,
		    "Number of CPUs     : %d\n", nproc);
	PetscPrintf(PETSC_COMM_WORLD,
		    "Number of Elements : %ld\n", (NX - 1) * (NY - 1) * (NZ- 1));
	PetscPrintf(PETSC_COMM_WORLD,
		    "Number of Nodes    : %ld\n", NX * NY * NZ);
	PetscPrintf(PETSC_COMM_WORLD,
		    "Number of DOFs     : %ld\n\n", (NX * NY * NZ) * DIM);
	PetscPrintf(PETSC_COMM_WORLD,
		    "NP_X : %ld\tNP_Y : %ld\tNP_Z : %ld\n", NP_X, NP_Y, NP_Z);
	PetscPrintf(PETSC_COMM_WORLD,
		    "NX   : %ld\tNY   : %ld\tNZ   : %ld\n\n", NX, NY, NZ);


	dx = lx / (NX - 1);
	dy = ly / (NY - 1);
	dz = lz / (NZ - 1);
	wg = dx * dy * dz / NPE;
	rad = 1.;


	/* Prepares the solver */

	KSPType ksptype;
	PetscReal rtol, abstol, dtol;
	PetscInt maxits;

	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
	ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
	ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
	ierr = KSPSetUp(ksp); CHKERRQ(ierr);

	ierr = KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &maxits);
	ierr = KSPGetType(ksp, &ksptype);
	PetscPrintf(PETSC_COMM_WORLD,
		    "KSP Info: type = %s\trtol = %e\t"
		    "abstol = %e\tdtol = %e\tmaxits = %d\n\n",
		    ksptype, rtol, abstol, dtol, maxits);


	ierr = DMDAGetElementsSizes(da, &nex, &ney, &nez); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da, &si_local, &sj_local, &sk_local,
				   &nx_local, &ny_local, &nz_local); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da, &si_ghost, &sj_ghost, &sk_ghost,
				   &nx_ghost, &ny_ghost, &nz_ghost); CHKERRQ(ierr);


	PetscSynchronizedPrintf(PETSC_COMM_WORLD,
				"rank:%d\tne:%d\tnex:%d\tney:%d\tnez:%d\n",
				rank, (int) (nex * ney * nez),
				(int) nex, (int) ney, (int) nez);
	PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);


	/* Measure the balance between processes */

	int min, max;
	minmax_elems_across_mpis(da, &min, &max);
	PetscPrintf(PETSC_COMM_WORLD,
		    "Min : %d Max : %d Unbalance (Max - Min) / Max = %3.1lf %\n",
		    min, max, 1. * (max - min) / (1. * max) * 100.);

	ierr = bc_init(da, &index_dirichlet, &nbcs,
		       &index_dirichlet_positive, &nbcs_positive);


	/* Prepares micropp (Constructor calls) */

	// Initializes <materials> declared in <micropp_c_wrapper.h>
	micropp_C_material_set(0,
			       micro_mat_1[0], micro_mat_1[1],
			       micro_mat_1[2], micro_mat_1[3], 1);
	micropp_C_material_set(1,
			       micro_mat_2[0], micro_mat_2[1],
			       micro_mat_2[2], micro_mat_2[3], 0);
	PetscPrintf(PETSC_COMM_WORLD, "Material Values : \n");

	if(!rank) {
		micropp_C_material_print(0);
		micropp_C_material_print(1);
	}

	// Initializes <micro> declared in <micropp_c_wrapper.h>
	int ngpl = nex * ney * nez * NGP;
	int size[3] = { micro_n, micro_n, micro_n };
	double params[4] = { 1., 1., 1., .5 };
	micropp_C_create3(ngpl, size, micro_type, params);

	if(!rank)
		micropp_C_print_info();

	return ierr;
}


PetscErrorCode finish()
{
	PetscErrorCode ierr;
	ierr = MatDestroy(&A); CHKERRQ(ierr);
	ierr = VecDestroy(&u); CHKERRQ(ierr);
	ierr = VecDestroy(&b); CHKERRQ(ierr);
	ierr = VecDestroy(&du); CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

	ierr = bc_finish(index_dirichlet); CHKERRQ(ierr);

	ierr = PetscFinalize();
	return ierr;
}
