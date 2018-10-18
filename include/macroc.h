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
#define NEWTON_TOL -1.0e-1
#define NEWTON_ITS 1
#define CONSTXG    0.577350269189626

#define FINAL_TIME 10.0
#define TIME_STEPS 1
#define DT         0.0001
#define NX         100
#define NY         100
#define NZ         100

#define U_MAX      0.8

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

double lx, ly, lz, dx, dy, dz;
double wg;

PetscInt ts;
PetscReal dt, final_time;

DM DA;
Mat A;
Vec u, du, b;

PetscErrorCode init();
PetscErrorCode finish();
PetscErrorCode set_bc(int time_step, Vec u);
PetscErrorCode set_strains();
PetscErrorCode assembly_jac(Mat A);
PetscErrorCode assembly_res(Vec b);
PetscErrorCode solve_Ax(Mat A, Vec b, Vec x);

void calc_B(int gp, double B[6][NPE * DIM]);

#endif
