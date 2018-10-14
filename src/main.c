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

static char help[] = "FE code to solve macroscopic problems with PETSc.\n";

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>

PetscErrorCode solve_elasticity_2d(PetscInt mx,PetscInt my);

int main(int argc,char **args)
{
  PetscErrorCode ierr;
  PetscInt       mx,my;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  mx   = my = 10;
  ierr = PetscOptionsGetInt(NULL,NULL,"-mx",&mx,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-my",&my,NULL);CHKERRQ(ierr);
  ierr = solve_elasticity_2d(mx,my);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
