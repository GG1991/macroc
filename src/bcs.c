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


/* Sets the displacement values on Dirichlet indeces
 * acording to the time step.
 */

PetscErrorCode apply_bc_on_u(double U, Vec u)
{
	PetscErrorCode ierr;

	if (bc_type == BC_BENDING) {

		ierr = bc_apply_on_u_bending(U, u);

	} else if (bc_type == BC_CIRCLE) {

		ierr = bc_apply_on_u_circle(U, u);

	}

	//VecView(u, PETSC_VIEWER_STDOUT_WORLD);
	return ierr;
}


/*
 * Ramp function to simulate Loading / Unloading
 */

double get_displacement(int time_s)
{
	double time = time_s * dt;
	double U;

	if(time < final_time / 2.)
		U = U_MAX * (time / final_time);
	else
		U = U_MAX;

	return U;
}


PetscErrorCode bc_apply_on_u_bending(double U, Vec u)
{
	PetscErrorCode ierr;
	PetscReal *bc_vals;
	PetscInt i, j, k, d;

	bc_vals = malloc(nbcs * sizeof(PetscReal));

	PetscInt index = 0;

	/* Surface X = 0 */
	if (si_ghost == 0)
		for (k = 0; k < nz_ghost; ++k)
			for (j = 0; j < ny_ghost; ++j)
				for (d = 0; d < DIM; ++d)
					bc_vals[index++] =  0.;

	/* Surface X = LX */
	if (si_ghost + nx_ghost == NX)
		for (k = 0; k < nz_ghost; ++k)
			for (j = 0; j < ny_ghost; ++j)
				for (d = 0; d < DIM; ++d)
					bc_vals[index++] = (d == 1) ? U : 0.;

	ierr = VecSetValues(u, nbcs, index_dirichlet, bc_vals, INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(u); CHKERRQ(ierr);

	free(bc_vals);
	return ierr;
}


PetscErrorCode bc_apply_on_u_circle(double U, Vec u)
{
	PetscErrorCode ierr;
	PetscReal *bc_vals;
	PetscInt i, j, k, d;

	bc_vals = malloc(nbcs * sizeof(PetscReal));

	PetscInt index = 0;

	/* X = 0 & Y = 0 (ALONG Z) */
	if (si_ghost == 0 && sj_ghost == 0)
		for (k = 0; k < nz_ghost; ++k)
			for (d = 0; d < DIM; ++d)
				bc_vals[index++] =  0.;

	/* X = LX & Y = 0 (ALONG Z) */
	if (si_ghost + nx_ghost == NX && sj_ghost == 0)
		for (k = 0; k < nz_ghost; ++k)
			for (d = 0; d < DIM; ++d)
				bc_vals[index++] =  0.;

	/* Z = 0 & Y = 0 (ALONG X) */
	if (sk_ghost == 0 && sj_ghost == 0)
		for (i = 1; i < nx_ghost - 1; ++i)
			for (d = 0; d < DIM; ++d)
				bc_vals[index++] =  0.;

	/* Z = LZ & Y = 0 (ALONG X) */
	if (sk_ghost + nz_ghost == NZ && sj_ghost == 0)
		for (i = 1; i < nx_ghost - 1; ++i)
			for (d = 0; d < DIM; ++d)
				bc_vals[index++] =  0.;

	/* Y = LY (INSIDE CIRCLE) */
	if (sj_ghost + ny_ghost == NY) {
		for (i = 0; i < nx_ghost; ++i) {
			for (k = 0; k < nz_ghost; ++k) {
				double x = lx / 2. - ((si_ghost + i) * dx + dx / 2.);
				double z = lz / 2. - ((sk_ghost + k) * dz + dz / 2.);
				if ((x * x + z * z) < (rad * rad))
					bc_vals[index++] = U;
			}
		}
	}

	ierr = VecSetValues(u, nbcs, index_dirichlet, bc_vals, INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(u); CHKERRQ(ierr);

	free(bc_vals);
	return ierr;
}


/* creates :
 * <nbcs>, <index_dirichlet>, <nbcs_positive> and <index_dirichlet_positive>
 * for being used in the BC settings on <A>, <u> and <b>.
 */

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

	/* filters positive values of <nbcs> and <index_dirichlet>
	 * to create  <nbcs_positive> and <index_dirichlet_positive>
	 * for being used in MatSetValues in BC setting
	 * (does not acept negative indeces).
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
	PetscInt i, j, k, d;
	PetscInt nbcs;
	PetscInt *ix;

	ISLocalToGlobalMapping ltogm;
	const PetscInt *g_idx;
	ierr = DMGetLocalToGlobalMapping(da, &ltogm); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingGetIndices(ltogm, &g_idx); CHKERRQ(ierr);

	nbcs = 2 * ny_ghost * nz_ghost * DIM;

	ix = malloc(nbcs * sizeof(PetscInt));
	for (i = 0; i < nbcs; ++i)
		ix[i] = -1;

	PetscInt index = 0;

	/* Surface X = 0 */
	if (si_ghost == 0) {
		i = 0;
		for (k = 0; k < nz_ghost; ++k) {
			for (j = 0; j < ny_ghost; ++j) {
				for (d = 0; d < DIM; ++d) {
					PetscInt local_id = i + j * nx_ghost + k * nx_ghost * ny_ghost;
					ix[index++] = g_idx[local_id * DIM + d];
				}
			}
		}
	}

	/* Surface X = LX */
	if (si_ghost + nx_ghost == NX) {
		i = nx_ghost - 1;
		for (k = 0; k < nz_ghost; ++k) {
			for (j = 0; j < ny_ghost; ++j) {
				for (d = 0; d < DIM; ++d) {
					PetscInt local_id = i + j * nx_ghost + k * nx_ghost * ny_ghost;
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


PetscErrorCode bc_init_circle(DM da, PetscInt **_index_dirichlet,
			      PetscInt *_nbcs)
{
	PetscErrorCode ierr;
	PetscInt i, j, k, d;
	PetscInt nbcs;
	PetscInt *ix;

	ISLocalToGlobalMapping ltogm;
	const PetscInt *g_idx;
	ierr = DMGetLocalToGlobalMapping(da, &ltogm); CHKERRQ(ierr);
	ierr = ISLocalToGlobalMappingGetIndices(ltogm, &g_idx); CHKERRQ(ierr);

	nbcs = (2 * nx_ghost + 2 * nz_ghost) * DIM + nx_ghost * nz_ghost;

	ix = malloc(nbcs * sizeof(PetscInt));
	for (i = 0; i < nbcs; ++i)
		ix[i] = -1;

	PetscInt index = 0;

	/* X = 0 & Y = 0 (ALONG Z) */
	if (si_ghost == 0 && sj_ghost == 0) {
		i = 0;
		j = 0;
		for (k = 0; k < nz_ghost; ++k) {
			PetscInt local_id = i + j * nx_ghost + k * nx_ghost * ny_ghost;
			for (d = 0; d < DIM; ++d)
				ix[index++] = g_idx[local_id * DIM + d];
		}
	}

	/* X = LX & Y = 0 (ALONG Z) */
	if (si_ghost + nx_ghost == NX && sj_ghost == 0) {
		i = nx_ghost - 1;
		j = 0;
		for (k = 0; k < nz_ghost; ++k) {
			PetscInt local_id = i + j * nx_ghost + k * nx_ghost * ny_ghost;
			for (d = 0; d < DIM; ++d)
				ix[index++] = g_idx[local_id * DIM + d];
		}
	}

	/* Z = 0 & Y = 0 (ALONG X) */
	if (sk_ghost == 0 && sj_ghost == 0) {
		k = 0;
		j = 0;
		for (i = 1; i < nx_ghost - 1; ++i) {
			PetscInt local_id = i + j * nx_ghost + k * nx_ghost * ny_ghost;
			for (d = 0; d < DIM; ++d)
				ix[index++] = g_idx[local_id * DIM + d];
		}
	}

	/* Z = LZ & Y = 0 (ALONG X) */
	if (sk_ghost + nz_ghost == NZ && sj_ghost == 0) {
		k = nz_ghost - 1;
		j = 0;
		for (i = 1; i < nx_ghost - 1; ++i) {
			PetscInt local_id = i + j * nx_ghost + k * nx_ghost * ny_ghost;
			for (d = 0; d < DIM; ++d)
				ix[index++] = g_idx[local_id * DIM + d];
		}
	}

	/* Circle in Surface <Y = LY>  &  <X^2 + Z^2 < R^2> */
	if (sj_ghost + ny_ghost == NY) {
		j = ny_ghost - 1;
		for (i = 0; i < nx_ghost; ++i) {
			for (k = 0; k < nz_ghost; ++k) {
				double x = lx / 2. - ((si_ghost + i) * dx + dx / 2.);
				double z = lz / 2. - ((sk_ghost + k) * dz + dz / 2.);
				const int d = 1;
				if ((x * x + z * z) < (rad * rad)) {
					PetscInt local_id = i + j * nx_ghost + k * nx_ghost * ny_ghost;
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
