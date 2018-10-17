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

#define NGP        8
#define NPE        8
#define NVOI       6
#define DIM        3
#define NEWTON_TOL 1.0e-1
#define NEWTON_ITS 1
#define CONSTXG    0.577350269189626

#define FINAL_TIME 10.0
#define TIME_STEPS 1
#define DT         0.0001
#define NX         2
#define NY         2
#define NZ         2

#define print0(mess) { if(!rank) printf("%s", mess); }
#define gpi(ex, ey, ez, gp) ( (ez * nexl * neyl + ey * nexl + ex) * NGP + gp)

static char help[] = "FE code to solve macroscopic problems with PETSc.\n";

static double xg[8][3] = {
    { -CONSTXG, -CONSTXG, -CONSTXG },
    { +CONSTXG, -CONSTXG, -CONSTXG },
    { +CONSTXG, +CONSTXG, -CONSTXG },
    { -CONSTXG, +CONSTXG, -CONSTXG },
    { -CONSTXG, -CONSTXG, +CONSTXG },
    { +CONSTXG, -CONSTXG, +CONSTXG },
    { +CONSTXG, +CONSTXG, +CONSTXG },
    { -CONSTXG, +CONSTXG, +CONSTXG } };

int rank, nproc, nproc_x, nproc_y, nproc_z;

int nelem, nelem_local;
int ngpl;

double lx, ly, lz, dx, dy, dz;
double wg;

PetscInt ts;
PetscReal dt, final_time;
double UY, U_max;

DM DA;
Mat A;
Vec u, du, b;

int init();
int finish();
int set_bc(int time_step, Vec u);
int set_strains();
int assembly_jac(Mat A);
int assembly_res(Vec b);
int solve_Ax();

void calc_B(int gp, double B[6][NPE * DIM]);

#endif
