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

#ifndef MACROC_H
#define MACROC_H

#include <stdio.h>

#include "micropp_c_wrapper.h"

#include "petscksp.h"
#include "petscdm.h"
#include "petscdmda.h"

#define NGP            8
#define NPE            8
#define NVOI           6
#define DIM            3
#define NEWTON_MIN_TOL 1.0e-1
#define NEWTON_MAX_ITS 2

#define FINAL_TIME     10.0
#define TIME_STEPS     1
#define VTU_FREQ       -1
#define DT             0.001
#define NX_CONST       5
#define NY_CONST       5
#define NZ_CONST       5
#define LX             10.0
#define LY             1.0
#define LZ             1.0

#define U_MAX          -0.8

static char help[] = "FE code to solve macroscopic problems with PETSc.\n";

#define CONSTXG        0.577350269189626

enum {BC_BENDING, BC_CIRCLE};

static double xg[8][3] = {
	{ -CONSTXG, -CONSTXG, -CONSTXG },
	{ +CONSTXG, -CONSTXG, -CONSTXG },
	{ +CONSTXG, +CONSTXG, -CONSTXG },
	{ -CONSTXG, +CONSTXG, -CONSTXG },
	{ -CONSTXG, -CONSTXG, +CONSTXG },
	{ +CONSTXG, -CONSTXG, +CONSTXG },
	{ +CONSTXG, +CONSTXG, +CONSTXG },
	{ -CONSTXG, +CONSTXG, +CONSTXG } };

double lx, ly, lz, dx, dy, dz;
double wg;

PetscReal dt, final_time;
PetscReal newton_min_tol;
PetscInt ts;
PetscInt vtu_freq;
PetscInt newton_max_its;
PetscInt bc_type;
PetscInt *index_dirichlet;
PetscInt *index_dirichlet_positive;
PetscInt nbcs;
PetscInt nbcs_positive;
PetscReal rad;


PetscInt NX, NY, NZ;                   /* Total number of nodes in each
					  direction
					  (THE SAME for all processors) */

PetscInt NP_X, NP_Y, NP_Z;             /* Number of processors in each
					  direction
					  (THE SAME for all processors) */

PetscInt nx_ghost, ny_ghost, nz_ghost; /* Number of nodes in each direction
					  including ghost nodes
					  (DIFFERENT for each processors) */

PetscInt nx_local, ny_local, nz_local; /* Number of nodes in each direction
					  (DIFFERENT for each processors) */


PetscInt si_ghost, sj_ghost, sk_ghost; /* Number of nodes in each direction
					  including ghost nodes
					  (DIFFERENT for each processors) */

PetscInt si_local, sj_local, sk_local; /* Number of nodes in each direction
					  (DIFFERENT for each processors) */


PetscInt nex, ney, nez;                /* Number of elements in each direction
					  (DIFFERENT for each processors) */


PetscInt nelem;                        /* Total number of elements
					  (DIFFERENT for each processors) */


DM da;
PC pc;
KSP ksp;
Mat A;
Vec u, du, b;

PetscErrorCode init();
PetscErrorCode finish();
PetscErrorCode set_strains();
PetscErrorCode assembly_jac(Mat A);
PetscErrorCode assembly_res(Vec b);
PetscErrorCode solve_Ax(KSP ksp, Vec b, Vec x);

void calc_B(int gp, double B[6][NPE * DIM]);

PetscErrorCode write_pvtu(const char *filename);
PetscErrorCode minmax_elems_across_mpis(DM da, int *min, int *max);

double get_displacement(int time_s);
PetscErrorCode apply_bc_on_u(double U, Vec u);
PetscErrorCode bc_apply_on_u_bending(double U, Vec u);

PetscErrorCode apply_bc_on_jac(Mat A);
PetscErrorCode bc_apply_on_u_bending(double U, Vec u);
PetscErrorCode bc_apply_on_u_circle(double U, Vec u);

PetscErrorCode apply_bc_on_res(Vec b);

PetscErrorCode bc_init(DM da, PetscInt **_index_dirichlet, PetscInt *_nbcs, PetscInt **_index_dirichlet_positive, PetscInt *_nbcs_positive);
PetscErrorCode bc_init_bending(DM da, PetscInt **_index_dirichlet, PetscInt *_nbcs);
PetscErrorCode bc_init_circle(DM da, PetscInt **_index_dirichlet, PetscInt *_nbcs);
PetscErrorCode bc_finish(PetscInt *index_dirichlet);

PetscInt get_non_linear_gps(void);

double calc_force(DM da);
PetscErrorCode calc_force_bending(DM da, PetscReal *_mpi_force);
PetscErrorCode calc_force_circle(DM da, PetscReal *_mpi_force);

#endif
