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


int set_bc(int time_step, Vec u)
{
    int ierr;
    int *bc_global_ids;
    double UY;
    double time = time_step * dt;
    double *bc_vals;
    const PetscInt *g_idx;
    int i, j, k, d;
    int si, sj, sk;
    int nx, ny, nz;

    ISLocalToGlobalMapping ltogm;

    if(time < final_time / 2.){
        UY = U_max * (time / final_time);
    } else {
        UY = U_max;
    }
    UY = 5.;

    ierr = DMGetLocalToGlobalMapping(DA, &ltogm); CHKERRQ(ierr);
    ierr = ISLocalToGlobalMappingGetIndices(ltogm, &g_idx); CHKERRQ(ierr);
    //ierr = DMDAGetGhostCorners(DA, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);
    ierr = DMDAGetCorners(DA, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);

    bc_global_ids = malloc(2 * ny * nz * DIM * sizeof(int));
    bc_vals = malloc(2 * ny * nz * DIM * sizeof(double));

    /* init the entries to -1 so VecSetValues will ignore them */
    for (i = 0; i < 2 * ny * nz * DIM; ++i){
        bc_global_ids[i] = -1;
    }

    int count = 0;

    i = 0; /* X = 0 */
    for (j = 0; j < ny; j++) {
        for (k = 0; k < nz; k++) {
            for (d = 0; d < DIM; d++) {

                int local_id = i + j * nx + k * nx * ny;

                bc_global_ids[count] = g_idx[DIM * local_id + d];
                bc_vals[count] =  0.;
                count ++;
            }
        }
    }

    i = nx - 1; /* X = LX */
    for (j = 0; j < ny; j++) {
        for (k = 0; k < nz; k++) {
            for (d = 0; d < DIM; d++) {

                int local_id = i + j * nx + k * nx * ny;

                bc_global_ids[count] = g_idx[DIM * local_id + d];
                bc_vals[count] =  (d == 1) ? UY : 0.;
                count ++;
            }
        }
    }

    ierr = ISLocalToGlobalMappingRestoreIndices(ltogm, &g_idx); CHKERRQ(ierr);

    if(u){
        ierr = VecSetValues(u, count, bc_global_ids, bc_vals, INSERT_VALUES);
        CHKERRQ(ierr);
        ierr = VecAssemblyBegin(u); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(u); CHKERRQ(ierr);
    }

    /* Debug */
    /*
    int rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    printf("rank:%d\tsi:%d\tsj:%d\tsk:%d\tnx:%d\tny:%d\tnz:%d\n",
           rank, si, sj, sk, nx, ny, nz);
    ierr = VecView(u, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    */

    free(bc_global_ids);
    free(bc_vals);

    return ierr;
}


int set_strains()
{
    int ex, ey, ez, gp;
    int ierr = 0;
    int gpi_;
    double strain[6];

    for(ex = 0; ex < nexl; ++ex) {
        for(ey = 0; ey < neyl; ++ey) {
            for(ez = 0; ex < nezl; ++ez) {

                for(gp = 0; gp < NGP; ++gp) {

                    strain[0] = 0.0;
//                    gpi_ = gpi(ex, ey, ez, gp);
//                    micropp_C_set_strain3(gpi_, strain);

                }

            }
        }
    }

    return ierr;
}


int assembly_jac(Mat A)
{
    int ierr = 0;
    int i, j, k, l;
    int ie, nelem;
    int n, d, gp, npe;
    int ix[NPE * DIM];
    double ctan[36];
    double Ae[NPE * DIM * NPE * DIM];
    double B[6][NPE * DIM];

    const PetscInt *eix, *eixp;

    ierr = MatZeroEntries(A);CHKERRQ(ierr);
    DMDAGetElements(DA, &nelem, &npe, &eix);
    for(ie = 0; ie < nelem; ++ie) {

        memset(Ae, 0., NPE * DIM * NPE * DIM * sizeof(double));
        for(gp = 0; gp < NGP; ++gp) {

            calc_B(gp, B);
            micropp_C_get_ctan3(gp, ctan);

            for (i = 0; i < NPE; ++i) {
                for (j = 0; j < NPE; ++j) {
                    for (k = 0; k < NVOI; ++k) {
                        for (l = 0; l < NVOI; ++l) {
                            Ae[NPE * i + j] += B[k][i] * ctan[k * NVOI + l] \
                                               * B[l][j] * wg;
                        }
                    }
                }
            }

        }
        eixp = &eix[ie * NPE];
        for(n = 0; n < NPE; ++n) {
            for(d = 0; d < DIM; ++d) {
                ix[n * DIM + d] = eix[n] * DIM + d;
            }
        }

        ierr = MatSetValuesLocal(A, NPE * DIM, ix, NPE * DIM, ix, Ae,
                                 ADD_VALUES);
        CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    return ierr;
}


int assembly_res(Vec b)
{
    int ierr = 0;
    int ie, nelem;
    int i, j;
    int n, d, gp, npe;
    int ix[NPE * DIM];
    double stress[6];
    double be[NPE * DIM * NPE * DIM];
    double B[6][NPE * DIM];

    const PetscInt *eix, *eixp;

    ierr = VecZeroEntries(b);
    ierr = DMDAGetElements(DA, &nelem, &npe, &eix); CHKERRQ(ierr);

    for(ie = 0; ie < nelem; ++ie) {

        memset(be, 0., NPE * DIM * sizeof(double));
        for(gp = 0; gp < NGP; ++gp) {

            calc_B(gp, B);
            micropp_C_get_stress3(gp, stress);

            for (i = 0; i < NPE * DIM; ++i) {
                for (j = 0; j < NVOI; ++j) {
                    be[i] += B[j][i] * stress[j] * wg;
                }
            }

        }

        eixp = &eix[ie * NPE];
        for(n = 0; n < NPE; ++n) {
            for(d = 0; d < DIM; ++d) {
                ix[n * DIM + d] = eix[n] * DIM + d;
            }
        }

        ierr = VecSetValuesLocal(b, NPE * DIM, ix, be, ADD_VALUES);
        CHKERRQ(ierr);
    }
    VecAssemblyBegin(b); CHKERRQ(ierr);
    VecAssemblyEnd(b); CHKERRQ(ierr);

    return ierr;
}


int solve_Ax(Mat A, Vec b, Vec x)
{
    int ierr;
    KSP ksp;
    PC pc;

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCBDDC); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = KSPSetUp(ksp); CHKERRQ(ierr);

    ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

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
