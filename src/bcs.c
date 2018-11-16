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
	}

	//VecView(u, PETSC_VIEWER_STDOUT_WORLD);
	return ierr;
}


PetscErrorCode bc_apply_on_u_bending(double U, Vec u)
{
	PetscErrorCode ierr;
	PetscInt *ix;
	PetscReal *bc_vals;
	PetscInt i, j, k, d;
	PetscInt si, sj, sk;
	PetscInt nx, ny, nz;
	PetscInt M, N, P;
	PetscInt nbcs;

	ISLocalToGlobalMapping ltogm;
	const PetscInt *g_idx;
	ierr = DMGetLocalToGlobalMapping(da, &ltogm); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingGetIndices(ltogm, &g_idx); CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &M, &N, &P, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	ierr = DMDAGetGhostCorners(da, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);

	nbcs = 2 * ny * nz * DIM;
	ix = malloc(nbcs * sizeof(PetscInt));
	bc_vals = malloc(nbcs * sizeof(PetscReal));
	for (i = 0; i < nbcs; ++i)
		ix[i] = -1;

	PetscInt index = 0;
	if (si == 0) {
		i = 0; /* X = 0 */
		for (k = 0; k < nz; ++k) {
			for (j = 0; j < ny; ++j) {
				for (d = 0; d < DIM; ++d) {

					PetscInt local_id = i + j * nx + k * nx * ny;

					ix[index] = g_idx[local_id * DIM + d];
					bc_vals[index] =  0.;
					index ++;
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

					ix[index] = g_idx[local_id * DIM + d];
					bc_vals[index] = (d == 1) ? U : 0.;
					index ++;
				}
			}
		}
	}

	ierr = VecSetValues(u, nbcs, ix, bc_vals, INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(u); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingRestoreIndices(ltogm, &g_idx); CHKERRQ(ierr);

	free(ix);
	free(bc_vals);
	return ierr;
}


PetscErrorCode apply_bc_on_jac(Mat A)
{
	PetscErrorCode ierr;
	PetscInt si, sj, sk;
	PetscInt nx, ny, nz;
	PetscInt M, N, P;
	PetscInt i, j, k, d;
	PetscInt nbcs = 0;
	const PetscInt *g_idx;
	PetscInt *rows;
	ISLocalToGlobalMapping ltogm;

	ierr = DMGetLocalToGlobalMapping(da, &ltogm); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingGetIndices(ltogm, &g_idx); CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &M, &N, &P, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	ierr = DMDAGetGhostCorners(da, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);

	ierr = PetscMalloc1(ny * nz * DIM, &rows); CHKERRQ(ierr);

	for (i = 0; i < ny * nz * DIM; ++i)
		rows[i] = -1;

	i = 0;
	for (k = 0; k < nz; ++k) {
		for (j = 0; j < ny; ++j) {
			for (d = 0; d < DIM; ++d) {

				PetscInt local_id = i + j * nx + k * nx * ny;
				PetscInt index = (k * ny + j) * DIM + d;
				rows[index] = g_idx[local_id * DIM + d];
			}
		}
	}

	nbcs = 0;
	if (si == 0)
		nbcs = ny * nz * DIM;

	MatZeroRowsColumns(A, nbcs, rows, 1., NULL, NULL);

	for (i = 0; i < ny * nz * DIM; ++i)
		rows[i] = -1;

	i = nx - 1;
	for (k = 0; k < nz; ++k) {
		for (j = 0; j < ny; ++j) {
			for (d = 0; d < DIM; ++d) {

				PetscInt local_id = i + j * nx + k * nx * ny;
				PetscInt index = (k * ny + j) * DIM + d;
				rows[index] = g_idx[local_id * DIM + d];
			}
		}
	}

	nbcs = 0;
	if (si + nx == M)
		nbcs = ny * nz * DIM;

	MatZeroRowsColumns(A, nbcs, rows, 1., NULL, NULL);

	ierr = ISLocalToGlobalMappingRestoreIndices(ltogm, &g_idx); CHKERRQ(ierr);
	PetscFree(rows);

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
}


PetscErrorCode bc_init(DM da, PetscInt **index_dirichlet, PetscInt *nbcs)
{
	PetscErrorCode ierr = 0;
	if (bc_type == BC_BENDING) {
		ierr = bc_init_bending(da, index_dirichlet, nbcs);
	}
	return ierr;
}


PetscErrorCode bc_init_bending(DM da, PetscInt **_index_dirichlet, PetscInt *_nbcs)
{
	PetscErrorCode ierr;
	PetscInt si, sj, sk;
	PetscInt nx, ny, nz;
	PetscInt i, j, k, d;
	PetscInt M, N, P;
	PetscInt nbcs;
	PetscInt index = 0;
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

	if (si == 0) {
		i = 0; /* X = 0 */
		for (k = 0; k < nz; ++k) {
			for (j = 0; j < ny; ++j) {
				for (d = 0; d < DIM; ++d) {
					PetscInt local_id = i + j * nx + k * nx * ny;
					ix[index] = g_idx[local_id * DIM + d];
					index ++;
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
					ix[index] = g_idx[local_id * DIM + d];
					index ++;
				}
			}
		}
	}

	*_index_dirichlet = ix;
	*_nbcs = nbcs;
	return ierr;
}


PetscErrorCode bc_finish(PetscInt *index_dirichlet)
{
	PetscErrorCode ierr = 0;

	free(index_dirichlet);

	return ierr;
}
