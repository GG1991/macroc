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


PetscErrorCode set_strains()
{
	PetscErrorCode ierr;
	PetscInt ie, nelem, npe, gpi;
	PetscInt i, j, gp, n, d;
	double strain[NVOI];
	double u_e[NPE * DIM];
	double B[NVOI][NPE * DIM];

	Vec u_loc;
	PetscScalar *u_arr;
	const PetscInt *eix;
	ierr = DMGetLocalVector(da, &u_loc); CHKERRQ(ierr);
	ierr = VecZeroEntries(u_loc); CHKERRQ(ierr);
	ierr = VecGetArray(u_loc, &u_arr);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da, u, INSERT_VALUES, u_loc); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, u, INSERT_VALUES, u_loc); CHKERRQ(ierr);
	ierr = DMDAGetElements(da, &nelem, &npe, &eix); CHKERRQ(ierr);

	for(ie = 0; ie < nelem; ++ie) {

		for(n = 0; n < NPE; ++n)
			for(d = 0; d < DIM; ++d)
				u_e[n * DIM + d] = u_arr[eix[ie * NPE + n] * DIM + d];

		for(gp = 0; gp < NGP; ++gp) {

			calc_B(gp, B);
			memset(strain, 0., NVOI * sizeof(double));
			for(i = 0; i < NVOI; ++i)
				for(j = 0; j < NPE * DIM; ++j)
					strain[i] += B[i][j] * u_e[j];

			gpi = ie * NGP + gp;
			micropp_C_set_strain3(gpi, strain);
		}
	}

	ierr = VecRestoreArray(u_loc, &u_arr); CHKERRQ(ierr);
	ierr = VecDestroy(&u_loc); CHKERRQ(ierr);
	return ierr;
}


PetscErrorCode assembly_jac(Mat A)
{
	PetscErrorCode ierr;
	PetscInt ie, nelem;
	PetscInt n, d, gp, gpi, npe;
	PetscInt ix[NPE * DIM];
	PetscInt i, j, k, l;
	double ctan[NVOI * NVOI];
	double Ae[NPE * DIM * NPE * DIM];
	double B[NVOI][NPE * DIM];

	ierr = MatZeroEntries(A); CHKERRQ(ierr);

	const PetscInt *eix;
	ierr = DMDAGetElements(da, &nelem, &npe, &eix); CHKERRQ(ierr);

	for(ie = 0; ie < nelem; ++ie) {

		memset(Ae, 0., NPE * DIM * NPE * DIM * sizeof(double));
		for(gp = 0; gp < NGP; ++gp) {

			calc_B(gp, B);
			gpi = ie * NGP + gp;
			micropp_C_get_ctan3(gpi, ctan);

			for (i = 0; i < NPE * DIM; ++i)
				for (j = 0; j < NPE * DIM; ++j)
					for (k = 0; k < NVOI; ++k)
						for (l = 0; l < NVOI; ++l)
							Ae[NPE * DIM * i + j] +=
								B[k][i] * ctan[k * NVOI + l] * B[l][j] * wg;

		}
		for(n = 0; n < NPE; ++n)
			for(d = 0; d < DIM; ++d)
				ix[n * DIM + d] = eix[ie * NPE + n] * DIM + d;

		ierr = MatSetValuesLocal(A, NPE * DIM, ix, NPE * DIM, ix, Ae,
					 ADD_VALUES); CHKERRQ(ierr);
	}
	ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	ierr = apply_bc_on_jac(A);

	//MatView(A, PETSC_VIEWER_DRAW_WORLD);

	return ierr;
}


PetscErrorCode assembly_res(Vec b, PetscReal *force)
{
	PetscErrorCode ierr;
	PetscInt ie, nelem;
	PetscInt i, j, k;
	PetscInt n, d, gp, gpi, npe;
	PetscInt ix[NPE * DIM];
	double stress[NVOI];
	double be[NPE * DIM * NPE * DIM];
	double B[NVOI][NPE * DIM];

	ierr = VecZeroEntries(b); CHKERRQ(ierr);

	Vec b_loc;
	PetscScalar *b_arr;
	ierr = DMGetLocalVector(da, &b_loc); CHKERRQ(ierr);
	ierr = VecZeroEntries(b_loc); CHKERRQ(ierr);
	ierr = VecGetArray(b_loc, &b_arr);CHKERRQ(ierr);

	const PetscInt *eix;
	ierr = DMDAGetElements(da, &nelem, &npe, &eix); CHKERRQ(ierr);

	for(ie = 0; ie < nelem; ++ie) {

		memset(be, 0., NPE * DIM * sizeof(double));
		for(gp = 0; gp < NGP; ++gp) {

			calc_B(gp, B);
			gpi = ie * NGP + gp;
			micropp_C_get_stress3(gpi, stress);

			for (i = 0; i < NPE * DIM; ++i)
				for (j = 0; j < NVOI; ++j)
					be[i] += B[j][i] * stress[j] * wg;
		}

		for(n = 0; n < NPE; ++n)
			for(d = 0; d < DIM; ++d)
				ix[n * DIM + d] = eix[ie * NPE + n] * DIM + d;

		for(n = 0; n < NPE * DIM; ++n)
			b_arr[ix[n]] += be[n];
	}

	ierr = DMLocalToGlobalBegin(da, b_loc, ADD_VALUES, b); CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da, b_loc, ADD_VALUES, b); CHKERRQ(ierr);
	ierr = VecRestoreArray(b_loc, &b_arr); CHKERRQ(ierr);
	ierr = VecDestroy(&b_loc);

	ierr = calc_force(da, b, force);

	//VecView(b, PETSC_VIEWER_STDOUT_WORLD);

	ierr = apply_bc_on_res(b);

	//VecView(b, PETSC_VIEWER_STDOUT_WORLD);

	ierr = VecScale(b, -1.); CHKERRQ(ierr);

	return ierr;
}


PetscErrorCode solve_Ax(KSP ksp, Vec b, Vec x)
{
	PetscErrorCode ierr;
	PetscInt its;
	PetscReal rnorm;

	ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(ksp, &its);
	ierr = KSPGetResidualNorm(ksp, &rnorm);
	PetscPrintf(PETSC_COMM_WORLD,
		    "KSP : |Ax - b|/|Ax| = %e\tIts = %d\n", rnorm, its);

	return ierr;
}


void calc_B(int gp, double B[6][NPE * DIM])
{
	int i;
	double dx = 1., dy = 1., dz = 1.;

	const double dsh[NPE][DIM] = {
		{
			-(1 - xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
			-(1 - xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
			-(1 - xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{
			+(1 - xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
			-(1 + xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
			-(1 + xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{
			+(1 + xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
			+(1 + xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
			-(1 + xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{
			-(1 + xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
			+(1 - xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
			-(1 - xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{
			-(1 - xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
			-(1 - xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
			+(1 - xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{
			+(1 - xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
			-(1 + xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
			+(1 + xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{
			+(1 + xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
			+(1 + xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
			+(1 + xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{
			-(1 + xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
			+(1 - xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
			+(1 - xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz } };

	for (i = 0; i < NPE; ++i) {
		B[0][i * DIM    ] = dsh[i][0];
		B[0][i * DIM + 1] = 0;
		B[0][i * DIM + 2] = 0;
		B[1][i * DIM    ] = 0;
		B[1][i * DIM + 1] = dsh[i][1];
		B[1][i * DIM + 2] = 0;
		B[2][i * DIM    ] = 0;
		B[2][i * DIM + 1] = 0;
		B[2][i * DIM + 2] = dsh[i][2];
		B[3][i * DIM    ] = dsh[i][1];
		B[3][i * DIM + 1] = dsh[i][0];
		B[3][i * DIM + 2] = 0;
		B[4][i * DIM    ] = dsh[i][2];
		B[4][i * DIM + 1] = 0;
		B[4][i * DIM + 2] = dsh[i][0];
		B[5][i * DIM    ] = 0;
		B[5][i * DIM + 1] = dsh[i][2];
		B[5][i * DIM + 2] = dsh[i][1];
	}
}
