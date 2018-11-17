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
	if (si == 0) { /* X = 0 */
		for (k = 0; k < nz; ++k) {
			for (j = 0; j < ny; ++j) {
				for (d = 0; d < DIM; ++d) {

					bc_vals[index] =  0.;
					index ++;

				}
			}
		}
	}


	if (si + nx == M) { /* X = LX */
		for (k = 0; k < nz; ++k) {
			for (j = 0; j < ny; ++j) {
				for (d = 0; d < DIM; ++d) {
					bc_vals[index] = (d == 1) ? U : 0.;
					index ++;
				}
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
	if (si == 0 && sj == 0) { /* Y = 0 and borders */
		for (k = 0; k < nz; ++k) {
			for (j = 0; j < ny; ++j) {
				for (d = 0; d < DIM; ++d) {

					bc_vals[index] =  0.;
					index ++;

				}
			}
		}
	}


	if (si + nx == M) { /* X = LX */
		for (k = 0; k < nz; ++k) {
			for (j = 0; j < ny; ++j) {
				for (d = 0; d < DIM; ++d) {
					bc_vals[index] = (d == 1) ? U : 0.;
					index ++;
				}
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
