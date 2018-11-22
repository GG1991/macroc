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


double calc_force(DM da)
{
	PetscErrorCode ierr;
	PetscReal mpi_force;

	switch (bc_type) {

		case (BC_BENDING):
			ierr = calc_force_bending(da, &mpi_force);
			break;

		case (BC_CIRCLE):
			ierr = calc_force_circle(da, &mpi_force);
			break;

		default:
			ierr = PETSC_ERR_ARG_WRONG;
			break;
	}

	double force = 0.0;

	ierr = MPI_Reduce(&mpi_force, &force, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	return force;
}


/*
   Integrates the average stress in the elements next to the boundary x = Lx

*/

PetscErrorCode calc_force_bending(DM da, PetscReal *_mpi_force)
{

	PetscErrorCode ierr;
	PetscInt gp, gpi;
	PetscInt e, ey, ez;
	PetscReal mpi_force = 0.0;
	PetscReal stress_ave[6];
	PetscReal stress[6];

	PetscInt nex, ney, nez;
	ierr = DMDAGetElementsSizes(da, &nex, &ney, &nez); CHKERRQ(ierr);

	PetscInt si, sj, sk;
	PetscInt nx, ny, nz;
	ierr = DMDAGetCorners(da, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);

	if (si + nx == NX) {
		for(ey = 0; ey < ney; ++ey) {
			for(ez = 0; ez < nez; ++ez) {

				e = (nex - 1) + ey * nex + ez * (nex * ney);

				memset(stress_ave, 0., NVOI * sizeof(PetscReal));
				for(gp = 0; gp < NGP; ++gp) {

					gpi = e * NGP + gp;
					micropp_C_get_stress3(gpi, stress);

					int i;
					for (i = 0; i < NVOI; ++i)
						stress_ave[i] += stress[i];

				}
				mpi_force += stress_ave[3] * dy * dz;

			}
		}
	}

	/* int rank;
	 * ierr  = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	 * PetscSynchronizedPrintf(PETSC_COMM_WORLD, "rank:%d\tforce_mpi:%lf\n", rank, mpi_force);
	 * PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	 */

	*_mpi_force = mpi_force;

	return ierr;
}


/*
   Integrates the average stress in the elements in the boundary y = Ly
   and in the circle x^2 + z^2 < rad

*/

PetscErrorCode calc_force_circle(DM da, PetscReal *_mpi_force)
{

	PetscErrorCode ierr;
	PetscInt gp, gpi;
	PetscInt e, ex, ez;
	PetscReal mpi_force = 0.0;
	PetscReal stress_ave[6];
	PetscReal stress[6];

	PetscInt nex, ney, nez;
	ierr = DMDAGetElementsSizes(da, &nex, &ney, &nez); CHKERRQ(ierr);

	PetscInt si, sj, sk;
	PetscInt nx, ny, nz;
	ierr = DMDAGetGhostCorners(da, &si, &sj, &sk, NULL, NULL, NULL); CHKERRQ(ierr);
	ierr = DMDAGetCorners(da, NULL, NULL, NULL, &nx, &ny, &nz); CHKERRQ(ierr);

	if (sj + ny == NY) {

		for(ex = 0; ex < nex; ++ex) {
			for(ez = 0; ez < nez; ++ez) {

				double x = lx / 2. - ((si + ex) * dx + dx / 2.);
				double z = lz / 2. - ((sk + ez) * dz + dz / 2.);

				if ((x * x + z * z) < 1 * (rad * rad)) {

					e = ex + (ney - 1) * nex + ez * (nex * ney);

					memset(stress_ave, 0., NVOI * sizeof(PetscReal));
					for(gp = 0; gp < NGP; ++gp) {

						gpi = e * NGP + gp;
						micropp_C_get_stress3(gpi, stress);

						int i;
						for (i = 0; i < NVOI; ++i)
							stress_ave[i] += stress[i];

					}
					mpi_force += stress_ave[1] * dx * dz;
				}

			}
		}
	}

	*_mpi_force = mpi_force;

	return ierr;
}
