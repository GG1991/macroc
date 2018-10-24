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


PetscErrorCode set_bc(int time_step, Vec u)
{
    PetscErrorCode ierr;
    PetscInt *ix;
    PetscReal time = time_step * dt;
    PetscReal *bc_vals;
    PetscInt i, j, k, d;
    PetscInt si, sj, sk;
    PetscInt nx, ny, nz;
    PetscInt M, N, P;
    PetscInt nbcs;

    double UY;
    if(time < final_time / 2.)
        UY = U_MAX * (time / final_time);
    else
        UY = U_MAX;

    ISLocalToGlobalMapping ltogm;
    const PetscInt *g_idx;
    ierr = DMGetLocalToGlobalMapping(da, &ltogm); CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingGetIndices(ltogm, &g_idx); CHKERRQ(ierr);
    ierr = DMDAGetInfo(da, 0, &M, &N, &P, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    ierr = DMDAGetGhostCorners(da, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);

    ix = malloc(ny * nz * DIM * sizeof(PetscInt));
    bc_vals = malloc(ny * nz * DIM * sizeof(PetscReal));

    i = 0; /* X = 0 */
    for (k = 0; k < nz; ++k) {
        for (j = 0; j < ny; ++j) {
            for (d = 0; d < DIM; ++d) {

                PetscInt local_id = i + j * nx + k * nx * ny;
                PetscInt index = (k * ny + j) * DIM + d;

                ix[index] = g_idx[local_id * DIM + d];
                bc_vals[index] =  0.;
            }
        }
    }

    nbcs = 0;
    if (si == 0)
        nbcs = ny * nz * DIM;

    ierr = VecSetValues(u, nbcs, ix, bc_vals, INSERT_VALUES); CHKERRQ(ierr);

    i = nx - 1; /* X = LX */
    for (k = 0; k < nz; ++k) {
        for (j = 0; j < ny; ++j) {
            for (d = 0; d < DIM; ++d) {

                PetscInt local_id = i + j * nx + k * nx * ny;
                PetscInt index = (k * ny + j) * DIM + d;

                ix[index] = g_idx[local_id * DIM + d];
                bc_vals[index] = (d == 1) ? UY : 0.;
            }
        }
    }

    nbcs = 0;
    if (si + nx == M)
        nbcs = ny * nz * DIM;

    ierr = VecSetValues(u, nbcs, ix, bc_vals, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(u); CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingRestoreIndices(ltogm, &g_idx); CHKERRQ(ierr);

    //VecView(u, PETSC_VIEWER_STDOUT_WORLD);

    free(ix);
    free(bc_vals);
    return ierr;
}


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


int assembly_jac(Mat A)
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

    /* Boundary Conditions */
    PetscInt si, sj, sk;
    PetscInt nx, ny, nz;
    PetscInt M, N, P;
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

    //MatView(A, PETSC_VIEWER_DRAW_WORLD);

    ierr = ISLocalToGlobalMappingRestoreIndices(ltogm, &g_idx); CHKERRQ(ierr);
    PetscFree(rows);
    return ierr;
}


PetscErrorCode assembly_res(Vec b)
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

    /* Boundary Conditions */
    PetscInt si, sj, sk;
    PetscInt nx, ny, nz;
    PetscInt M, N, P;
    PetscInt nbcs = 0;
    ISLocalToGlobalMapping ltogm;
    const PetscInt *g_idx;
    ierr = DMGetLocalToGlobalMapping(da, &ltogm); CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingGetIndices(ltogm, &g_idx); CHKERRQ(ierr);

    ierr = DMDAGetInfo(da, 0, &M, &N, &P, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    ierr = DMDAGetGhostCorners(da, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);

    PetscInt *rows = malloc(ny * nz * DIM * sizeof(PetscInt));
    double *zeros = calloc(ny * nz * DIM, sizeof(double));

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

    ierr = VecSetValues(b, nbcs, rows, zeros, INSERT_VALUES); CHKERRQ(ierr);

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

    ierr = VecSetValues(b, nbcs, rows, zeros, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

    //VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    ierr = VecScale(b, -1.); CHKERRQ(ierr);

    free(rows);
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
