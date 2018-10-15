/*
 *  This source code is part of MacroC: a finite element code
 *  to solve macrostructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Guido Giuntoli <gagiuntoli@gmail.com>
 *                         Based on the PETSc example develop by:
 *                         Dave May
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

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>

#define NEWTON_TOL 1.0e-1
#define NEWTON_ITS 4

#define print0(mess) { if(!rank) printf("%s", mess); }

static char help[] = "FE code to solve macroscopic problems with PETSc.\n";
int rank, nproc;

PetscErrorCode solve_elasticity_2d(PetscInt mx,PetscInt my);
PetscErrorCode init(int nx, int ny, int nz);
PetscErrorCode set_bc(int time_step);
PetscErrorCode assembly_jac(void);
PetscErrorCode assembly_res(void);
PetscErrorCode solve_Ax(void);

#endif
