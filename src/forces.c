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


PetscErrorCode calc_force(DM da, Vec b, PetscReal *force)
{
	PetscErrorCode ierr;
	PetscReal force_per_mpi;

	switch (bc_type) {
		case (BC_BENDING):
			ierr = calc_force_bending(da, &force_per_mpi);
			break;

		case (BC_CIRCLE):
			ierr = calc_force_circle(da, &force_per_mpi);
			break;

		default:
			break;
	}

	*force = 0.0;
	ierr = MPI_Reduce(&force_per_mpi, force, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	return ierr;
}


PetscErrorCode calc_force_bending(DM da, PetscReal *_force_per_mpi)
{
	/* Integrates the average stress in the elements next to the boundary x = Lx */

	PetscErrorCode ierr;
	PetscInt gp, gpi;
	PetscInt e, ey, ez;
	PetscReal force_per_mpi = 0.0;
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
				force_per_mpi += stress_ave[3] * dy * dz;

			}
		}
	}

	/* int rank;
	 * ierr  = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	 * PetscSynchronizedPrintf(PETSC_COMM_WORLD, "rank:%d\tforce_mpi:%lf\n", rank, force_per_mpi);
	 * PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	 */

	*_force_per_mpi = force_per_mpi;

	return ierr;
}

PetscErrorCode calc_force_circle(DM da, PetscReal *_force_per_mpi)
{
	PetscErrorCode ierr;
	PetscInt gp, gpi;
	PetscInt e, ex, ez;
	PetscReal force_per_mpi = 0.0;
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
					force_per_mpi += stress_ave[1] * dx * dz;
				}

			}
		}
	}

	*_force_per_mpi = force_per_mpi;

	return ierr;
}
