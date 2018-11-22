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


PetscErrorCode apply_bc_on_u(int time_step, Vec u)
{
	PetscErrorCode ierr;
	PetscReal time = time_step * dt;

	double U;
	if(time < final_time / 2.)
		U = U_MAX * (time / final_time);
	else
		U = U_MAX;

	if (bc_type == BC_BENDING) {

		ierr = bc_apply_on_u_bending(U, u);

	} else if (bc_type == BC_CIRCLE) {

		ierr = bc_apply_on_u_circle(U, u);

	}


	//VecView(u, PETSC_VIEWER_STDOUT_WORLD);
	return ierr;
}


PetscErrorCode bc_apply_on_u_bending(double U, Vec u)
{
	PetscErrorCode ierr;
	PetscReal *bc_vals;
	PetscInt i, j, k, d;
	PetscInt si, sj, sk;
	PetscInt nx, ny, nz;
	PetscInt M, N, P;

	ISLocalToGlobalMapping ltogm;
	const PetscInt *g_idx;
	ierr = DMGetLocalToGlobalMapping(da, &ltogm); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingGetIndices(ltogm, &g_idx); CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &M, &N, &P, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	ierr = DMDAGetGhostCorners(da, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);

	bc_vals = malloc(nbcs * sizeof(PetscReal));

	PetscInt index = 0;
	if (si == 0) /* X = 0 */
		for (k = 0; k < nz; ++k)
			for (j = 0; j < ny; ++j)
				for (d = 0; d < DIM; ++d)
					bc_vals[index++] =  0.;

	if (si + nx == M) /* X = LX */
		for (k = 0; k < nz; ++k)
			for (j = 0; j < ny; ++j)
				for (d = 0; d < DIM; ++d)
					bc_vals[index++] = (d == 1) ? U : 0.;

	ierr = VecSetValues(u, nbcs, index_dirichlet, bc_vals, INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(u); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingRestoreIndices(ltogm, &g_idx); CHKERRQ(ierr);

	free(bc_vals);
	return ierr;
}


PetscErrorCode bc_apply_on_u_circle(double U, Vec u)
{
	PetscErrorCode ierr;
	PetscReal *bc_vals;
	PetscInt i, j, k, d;
	PetscInt si, sj, sk;
	PetscInt nx, ny, nz;
	PetscInt M, N, P;

	ISLocalToGlobalMapping ltogm;
	const PetscInt *g_idx;
	ierr = DMGetLocalToGlobalMapping(da, &ltogm); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingGetIndices(ltogm, &g_idx); CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &M, &N, &P, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	ierr = DMDAGetGhostCorners(da, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);

	bc_vals = malloc(nbcs * sizeof(PetscReal));

	PetscInt index = 0;

	if (si == 0 && sj == 0) /* X = 0 & Y = 0 (ALONG Z) */
		for (k = 0; k < nz; ++k)
			for (d = 0; d < DIM; ++d)
				bc_vals[index++] =  0.;

	if (si + nx == M && sj == 0) /* X = LX & Y = 0 (ALONG Z) */
		for (k = 0; k < nz; ++k)
			for (d = 0; d < DIM; ++d)
				bc_vals[index++] =  0.;

	if (sk == 0 && sj == 0) /* Z = 0 & Y = 0 (ALONG X) */
		for (i = 1; i < nx - 1; ++i)
			for (d = 0; d < DIM; ++d)
				bc_vals[index++] =  0.;

	if (sk + nz == P && sj == 0) /* Z = LZ & Y = 0 (ALONG X) */
		for (i = 1; i < nx - 1; ++i)
			for (d = 0; d < DIM; ++d)
				bc_vals[index++] =  0.;


	if (sj + ny == N) { /* Y = LY (INSIDE CIRCLE) */
		for (i = 0; i < nx; ++i) {
			for (k = 0; k < nz; ++k) {
				double x = lx / 2. - ((si + i) * dx + dx / 2.);
				double z = lz / 2. - ((sk + k) * dz + dz / 2.);
				if ((x * x + z * z) < (rad * rad))
					bc_vals[index++] = U;
			}
		}
	}

	ierr = VecSetValues(u, nbcs, index_dirichlet, bc_vals, INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(u); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingRestoreIndices(ltogm, &g_idx); CHKERRQ(ierr);

	free(bc_vals);
	return ierr;
}


PetscErrorCode bc_init(DM da, PetscInt **_index_dirichlet, PetscInt *_nbcs,
		       PetscInt **_index_dirichlet_positive,
		       PetscInt *_nbcs_positive)
{
	PetscErrorCode ierr;
	PetscInt i;
	PetscInt nbcs, *index_dirichlet;
	PetscInt *index_dirichlet_positive, nbcs_positive;

	if (bc_type == BC_BENDING) {

		ierr = bc_init_bending(da, _index_dirichlet, _nbcs);

	} else if (bc_type == BC_CIRCLE) {

		ierr = bc_init_circle(da, _index_dirichlet, _nbcs);

	}
	nbcs = *_nbcs;
	index_dirichlet = *_index_dirichlet;

	/* <nbcs_positive> and <index_dirichlet_positive> are used for
	 * MatSetValues that does not acept negative indeces.
	 */
	nbcs_positive = 0;
	for (i = 0; i < nbcs; ++i)
		if (index_dirichlet[i] >= 0)
			nbcs_positive ++;

	index_dirichlet_positive = malloc(nbcs_positive * sizeof(PetscInt));

	for (i = 0; i < nbcs; ++i)
		if (index_dirichlet[i] >= 0)
			index_dirichlet_positive[i] = index_dirichlet[i];

	*_nbcs_positive = nbcs_positive;
	*_index_dirichlet_positive = index_dirichlet_positive;

	return ierr;
}


PetscErrorCode bc_init_bending(DM da, PetscInt **_index_dirichlet,
			       PetscInt *_nbcs)
{
	PetscErrorCode ierr;
	PetscInt si, sj, sk;
	PetscInt nx, ny, nz;
	PetscInt i, j, k, d;
	PetscInt M, N, P;
	PetscInt nbcs;
	PetscInt *ix;

	ISLocalToGlobalMapping ltogm;
	const PetscInt *g_idx;
	ierr = DMGetLocalToGlobalMapping(da, &ltogm); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingGetIndices(ltogm, &g_idx); CHKERRQ(ierr);

	ierr = DMDAGetInfo(da, 0, &M, &N, &P, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	ierr = DMDAGetGhostCorners(da, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);

	nbcs = 2 * ny * nz * DIM;

	ix = malloc(nbcs * sizeof(PetscInt));
	for (i = 0; i < nbcs; ++i)
		ix[i] = -1;

	PetscInt index = 0;

	if (si == 0) {
		i = 0; /* X = 0 */
		for (k = 0; k < nz; ++k) {
			for (j = 0; j < ny; ++j) {
				for (d = 0; d < DIM; ++d) {
					PetscInt local_id = i + j * nx + k * nx * ny;
					ix[index++] = g_idx[local_id * DIM + d];
				}
			}
		}
	}


	if (si + nx == M) {
		i = nx - 1; /* X = LX */
		for (k = 0; k < nz; ++k) {
			for (j = 0; j < ny; ++j) {
				for (d = 0; d < DIM; ++d) {
					PetscInt local_id = i + j * nx + k * nx * ny;
					ix[index++] = g_idx[local_id * DIM + d];
				}
			}
		}
	}

	*_index_dirichlet = ix;
	*_nbcs = nbcs;

	ierr = ISLocalToGlobalMappingRestoreIndices(ltogm, &g_idx); CHKERRQ(ierr);

	return ierr;
}

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
	PetscErrorCode ierr;
	PetscInt gp, gpi;
	PetscInt e, ey, ez;
	PetscReal force_per_mpi = 0.0;
	PetscReal stress_ave[6];
	PetscReal stress[6];
	PetscInt M, N, P;

	PetscInt nex, ney, nez;
	ierr = DMDAGetElementsSizes(da, &nex, &ney, &nez); CHKERRQ(ierr);

	PetscInt si, sj, sk;
	PetscInt nx, ny, nz;
	ierr = DMDAGetInfo(da, 0, &M, &N, &P, 0, 0, 0, 0, 0, 0, 0, 0, 0); CHKERRQ(ierr);
	ierr = DMDAGetCorners(da, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);

	if (si + nx == M) {
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
	PetscInt M, N, P;

	PetscInt nex, ney, nez;
	ierr = DMDAGetElementsSizes(da, &nex, &ney, &nez); CHKERRQ(ierr);

	PetscInt si, sj, sk;
	PetscInt nx, ny, nz;
	ierr = DMDAGetInfo(da, 0, &M, &N, &P, 0, 0, 0, 0, 0, 0, 0, 0, 0); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(da, &si, &sj, &sk, NULL, NULL, NULL); CHKERRQ(ierr);
	ierr = DMDAGetCorners(da, NULL, NULL, NULL, &nx, &ny, &nz); CHKERRQ(ierr);

	int rank;
	ierr  = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	if (sj + ny == N) {

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

	/*
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "rank:%d\tforce_mpi:%lf\n", rank, force_per_mpi);
	PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	*/

	*_force_per_mpi = force_per_mpi;

	return ierr;
}


PetscErrorCode bc_init_circle(DM da, PetscInt **_index_dirichlet,
			      PetscInt *_nbcs)
{
	PetscErrorCode ierr;
	PetscInt si, sj, sk;
	PetscInt nx, ny, nz;
	PetscInt i, j, k, d;
	PetscInt M, N, P;
	PetscInt nbcs;
	PetscInt *ix;

	ISLocalToGlobalMapping ltogm;
	const PetscInt *g_idx;
	ierr = DMGetLocalToGlobalMapping(da, &ltogm); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingGetIndices(ltogm, &g_idx); CHKERRQ(ierr);

	ierr = DMDAGetInfo(da, 0, &M, &N, &P, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	ierr = DMDAGetGhostCorners(da, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);

	nbcs = (2 * nx + 2 * nz) * DIM + nx * nz;

	ix = malloc(nbcs * sizeof(PetscInt));
	for (i = 0; i < nbcs; ++i)
		ix[i] = -1;

	PetscInt index = 0;

	if (si == 0 && sj == 0) { /* X = 0 & Y = 0 (ALONG Z) */
		i = 0;
		j = 0;
		for (k = 0; k < nz; ++k) {
			PetscInt local_id = i + j * nx + k * nx * ny;
			for (d = 0; d < DIM; ++d)
				ix[index++] = g_idx[local_id * DIM + d];
		}
	}

	if (si + nx == M && sj == 0) { /* X = LX & Y = 0 (ALONG Z) */
		i = nx - 1;
		j = 0;
		for (k = 0; k < nz; ++k) {
			PetscInt local_id = i + j * nx + k * nx * ny;
			for (d = 0; d < DIM; ++d)
				ix[index++] = g_idx[local_id * DIM + d];
		}
	}

	if (sk == 0 && sj == 0) { /* Z = 0 & Y = 0 (ALONG X) */
		k = 0;
		j = 0;
		for (i = 1; i < nx - 1; ++i) {
			PetscInt local_id = i + j * nx + k * nx * ny;
			for (d = 0; d < DIM; ++d)
				ix[index++] = g_idx[local_id * DIM + d];
		}
	}

	if (sk + nz == P && sj == 0) { /* Z = LZ & Y = 0 (ALONG X) */
		k = nz - 1;
		j = 0;
		for (i = 1; i < nx - 1; ++i) {
			PetscInt local_id = i + j * nx + k * nx * ny;
			for (d = 0; d < DIM; ++d)
				ix[index++] = g_idx[local_id * DIM + d];
		}
	}


	if (sj + ny == N) { /* Y = LY (INSIDE CIRCLE) */
		j = ny - 1;
		for (i = 0; i < nx; ++i) {
			for (k = 0; k < nz; ++k) {
				double x = lx / 2. - ((si + i) * dx + dx / 2.);
				double z = lz / 2. - ((sk + k) * dz + dz / 2.);
				const int d = 1;
				if ((x * x + z * z) < (rad * rad)) {
					PetscInt local_id = i + j * nx + k * nx * ny;
					ix[index++] = g_idx[local_id * DIM + d];
				}
			}
		}
	}

	*_index_dirichlet = ix;
	*_nbcs = nbcs;
	return ierr;
}


PetscErrorCode apply_bc_on_jac(Mat A)
{
	PetscErrorCode ierr;
	ierr = MatZeroRowsColumns(A, nbcs_positive, index_dirichlet_positive,
				  1., NULL, NULL);
	return ierr;
}


PetscErrorCode apply_bc_on_res(Vec b)
{
	PetscErrorCode ierr;

	PetscReal *zeros = calloc(nbcs, sizeof(PetscReal));

	ierr = VecSetValues(b, nbcs, index_dirichlet, zeros, INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

	free(zeros);
	return ierr;
}


PetscErrorCode bc_finish(PetscInt *index_dirichlet)
{
	PetscErrorCode ierr = 0;

	free(index_dirichlet);
	free(index_dirichlet_positive);

	return ierr;
}
