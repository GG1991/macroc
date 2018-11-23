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


int main(int argc,char **args)
{
	PetscErrorCode ierr;
	double t1, t2;
	double norm;
	int time_s, newton_it;
	FILE *file_out;

	ierr = PetscInitialize(&argc, &args, (char*)0, help);
	if(ierr)
		return ierr;

	PetscFOpen(PETSC_COMM_WORLD, "f_vs_d.dat", "w", &file_out);
	PetscPrintf(PETSC_COMM_WORLD, "\nMacroC : A HPC for FE2 Multi-scale Simulations\n\n");
	PetscPrintf(PETSC_COMM_WORLD, " >>> Starting Initialization : \n\n");

	ierr = init();

	PetscPrintf(PETSC_COMM_WORLD, "\n <<< Finishing Initialization.\n\n\n");
	PetscPrintf(PETSC_COMM_WORLD, " >>> Starting Calculation.\n");

	t1 = MPI_Wtime();

	for(time_s = 0; time_s < ts; ++time_s) {

		PetscPrintf(PETSC_COMM_WORLD, "\n\nTime Step = %d\n", time_s);

		double U = get_displacement(time_s);
		ierr = apply_bc_on_u(U, u);

		newton_it = 0;
		while(newton_it < newton_max_its) {

			PetscPrintf(PETSC_COMM_WORLD, "\nNewton Iteration = %d\n", newton_it);
			PetscPrintf(PETSC_COMM_WORLD, "Homogenizing MicroPP\n");
			ierr = set_strains();
			micropp_C_homogenize();

			PetscPrintf(PETSC_COMM_WORLD, "Assemblying RHS\n");

			ierr = assembly_res(b);
			ierr = VecNorm(b, NORM_2, &norm); CHKERRQ(ierr);
			PetscPrintf(PETSC_COMM_WORLD, "|RES| = %e\n", norm);
			if (norm < newton_min_tol) break;

			ierr = assembly_jac(A);
			ierr = solve_Ax(ksp, b, du);

			ierr = VecAXPY(u, 1., du); CHKERRQ(ierr);

			newton_it ++;
		}
		micropp_C_update_vars();


		/* Output the solution */

		double force = calc_force(da);

		PetscFPrintf(PETSC_COMM_WORLD, file_out, "%d\t%e\t%e\t%e\n",
			     time_s, time_s * dt, U, force);

		int64_t non_linear_micro_gps = get_non_linear_gps();
		PetscPrintf(PETSC_COMM_WORLD, "Non-Linear Gauss points : %ld\n", non_linear_micro_gps);


		if (vtu_freq > 0 && time_s % vtu_freq == 0) {

			char file_prefix[PETSC_MAX_PATH_LEN];
			ierr = PetscSNPrintf(file_prefix,
					     sizeof(file_prefix), "solution_%d",
					     time_s); CHKERRQ(ierr);

			ierr = write_pvtu(file_prefix);
		}
	}

	t2 = MPI_Wtime();

	PetscPrintf(PETSC_COMM_WORLD,
		    "\n\n <<< Calculation Finished OK\n\n"
		    "Elapsed time : %f\n", t2 - t1);

	ierr = PetscFClose(PETSC_COMM_WORLD, file_out);

	ierr = finish();
	return ierr;
}
