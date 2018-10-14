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

#define NEWTON_TOL 1.0e-1
#define NEWTON_ITS 4

#include <stdio.h>

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>

#define print0(mess) { if(!rank) printf("%s", mess); }

PetscErrorCode solve_elasticity_2d(PetscInt mx,PetscInt my);
PetscErrorCode init(int nx, int ny, int nz);
PetscErrorCode set_bc(int time_step);
PetscErrorCode assembly_jac(void);
PetscErrorCode assembly_res(void);
PetscErrorCode solve_Ax(void);

int main(int argc,char **args)
{
  	PetscErrorCode ierr;
  	PetscInt       nx = 10, ny = 10, nz = 10;
  	int time_s, newton_it, tsteps = 1;
	double norm;
  	int rank, nproc;
  	char mess[64];

  	ierr = PetscInitialize(&argc,&args,(char*)0,help); if(ierr) return ierr;
  	ierr = PetscOptionsGetInt(NULL,NULL,"-nx",&nx,NULL);CHKERRQ(ierr);
  	ierr = PetscOptionsGetInt(NULL,NULL,"-ny",&ny,NULL);CHKERRQ(ierr);
  	ierr = PetscOptionsGetInt(NULL,NULL,"-zy",&nz,NULL);CHKERRQ(ierr);
  	ierr = PetscOptionsGetInt(NULL,NULL,"-ts",&tsteps,NULL);CHKERRQ(ierr);

  	ierr = init(nx, ny, nz);

  	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  	sprintf(mess, "Problem size %d\n", nproc);
  	print0(mess);

  	for(time_s = 0; time_s < tsteps; ++time_s) {

		ierr = set_bc(time_s);

        newton_it = 0;
  		while(newton_it < NEWTON_ITS ) {

	  		ierr = assembly_res();
        	/* norm = |b| */

	  		ierr = assembly_jac();
	  		ierr = solve_Ax();
        	/* u = u + du */

        	newton_it ++;
        }
  	}

  	ierr = solve_elasticity_2d(nx,ny);CHKERRQ(ierr);
  	ierr = PetscFinalize();
  	return ierr;
}
