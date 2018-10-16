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


int set_bc(int time_step)
{
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

int assembly_jac(Mat A)
{
    int ierr = 0;
    int ex, ey, ez, gp;
    int sex, sey, sez;
    int i, j, k, l;
    double ctan[36];
    double Ae[NPE * DIM * NPE * DIM];
    double B[6][NPE * DIM];

    MatStencil u_eqn[NPE * DIM];

    ierr = DMDAGetElementsCorners(DA, &sex, &sey, &sez); CHKERRQ(ierr);
    for(ex = 0; ex < nexl; ++ex) {
        for(ey = 0; ey < neyl; ++ey) {
            for(ez = 0; ez < nezl; ++ez) {

                memset(Ae, 0., NPE * DIM * NPE * DIM * sizeof(double));
                for(gp = 0; gp < NGP; ++gp) {

                    calc_B(gp, B);
                    micropp_C_get_ctan3(gp, ctan);

                    for (i = 0; i < NPE; ++i) {
                        for (j = 0; j < NPE; ++j) {
                            for (k = 0; k < NVOI; ++k) {
                                for (l = 0; l < NVOI; ++l) {
                                    Ae[NPE * i + j] += \
                                                       B[k][i] \
                                                       * ctan[k * NVOI + l] \
                                                       * B[l][j] * wg;
                                }
                            }
                        }
                    }

                }

                get_elem_stencil(u_eqn, sex + ex, sey + ey, sez + ez);
                ierr = MatSetValuesStencil(A, NPE * DIM, u_eqn, NPE * DIM,
                                           u_eqn, Ae, ADD_VALUES);
                CHKERRQ(ierr);
            }
        }
    }
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    return ierr;
}

#define LOCIX(ex, ey, ez) ((ez) * (nexl + 1)* (neyl + 1) + (ey) * (nexl + 1) + (ex)) 
void get_elem_ix_local(int ex, int ey, int ez, PetscInt ixloc[NPE * DIM])
{
    int d;

    for(d = 0; d < DIM; ++d){
        ixloc[0 * DIM + d] = LOCIX(ex    , ey    , ez    ) * DIM + d;
        ixloc[1 * DIM + d] = LOCIX(ex + 1, ey    , ez    ) * DIM + d;
        ixloc[2 * DIM + d] = LOCIX(ex + 1, ey + 1, ez    ) * DIM + d;
        ixloc[3 * DIM + d] = LOCIX(ex    , ey + 1, ez    ) * DIM + d;
        ixloc[4 * DIM + d] = LOCIX(ex    , ey    , ez + 1) * DIM + d;
        ixloc[5 * DIM + d] = LOCIX(ex + 1, ey    , ez + 1) * DIM + d;
        ixloc[6 * DIM + d] = LOCIX(ex + 1, ey + 1, ez + 1) * DIM + d;
        ixloc[7 * DIM + d] = LOCIX(ex    , ey + 1, ez + 1) * DIM + d;
    }

}

int assembly_res(Vec b)
{
    int ierr = 0;
    int ex, ey, ez, gp;
    int sex, sey, sez;
    int i, j;
    double stress[6];
    double be[NPE * DIM * NPE * DIM];
    double B[6][NPE * DIM];
    PetscInt ix_loc[NPE * DIM];

    MatStencil u_eqn[NPE * DIM];

    ierr = DMDAGetElementsCorners(DA, &sex, &sey, &sez); CHKERRQ(ierr);
    for(ex = 0; ex < nexl; ++ex) {
        for(ey = 0; ey < neyl; ++ey) {
            for(ez = 0; ez < nezl; ++ez) {

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

                get_elem_ix_local(ex, ey, ez, ix_loc);
                VecSetValuesLocal(b, NPE * DIM, ix_loc, be, ADD_VALUES);
            }
        }
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

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


void get_elem_stencil(MatStencil s_u[NPE * DIM], int ei, int ej, int ek)
{
    /* displacement */
    /* node 0 */
    s_u[0].i = ei ; s_u[0].j = ej; s_u[0].k = ek; s_u[0].c = 0; /* Ux0 */
    s_u[1].i = ei ; s_u[1].j = ej; s_u[1].k = ek; s_u[1].c = 1; /* Uy0 */
    s_u[2].i = ei ; s_u[2].j = ej; s_u[2].k = ek; s_u[2].c = 2; /* Uz0 */

    /* node 1 */
    s_u[3].i = ei + 1; s_u[3].j = ej; s_u[3].k = ek; s_u[3].c = 0; /* Ux0 */
    s_u[4].i = ei + 1; s_u[4].j = ej; s_u[4].k = ek; s_u[4].c = 1; /* Uy0 */
    s_u[5].i = ei + 1; s_u[5].j = ej; s_u[5].k = ek; s_u[5].c = 2; /* Uz0 */

    /* node 2 */
    s_u[6].i = ei + 1; s_u[6].j = ej + 1; s_u[6].k = ek; s_u[6].c = 0; /* Ux0 */
    s_u[7].i = ei + 1; s_u[7].j = ej + 1; s_u[7].k = ek; s_u[7].c = 1; /* Uy0 */
    s_u[8].i = ei + 1; s_u[8].j = ej + 1; s_u[8].k = ek; s_u[8].c = 2; /* Uz0 */

    /* node 3 */
    s_u[9].i = ei; s_u[9].j = ej + 1; s_u[9].k = ek; s_u[9].c = 0; /* Ux0 */
    s_u[10].i = ei; s_u[10].j = ej + 1; s_u[10].k = ek; s_u[10].c = 1; /* Uy0 */
    s_u[11].i = ei; s_u[11].j = ej + 1; s_u[11].k = ek; s_u[11].c = 2; /* Uz0 */

    /* node 4 */
    s_u[12].i = ei ; s_u[12].j = ej; s_u[12].k = ek + 1; s_u[12].c = 0; /* Ux0 */
    s_u[13].i = ei ; s_u[13].j = ej; s_u[13].k = ek + 1; s_u[13].c = 1; /* Uy0 */
    s_u[14].i = ei ; s_u[14].j = ej; s_u[14].k = ek + 1; s_u[14].c = 2; /* Uz0 */

    /* node 5 */
    s_u[15].i = ei + 1; s_u[15].j = ej; s_u[15].k = ek + 1; s_u[15].c = 0; /* Ux0 */
    s_u[16].i = ei + 1; s_u[16].j = ej; s_u[16].k = ek + 1; s_u[16].c = 1; /* Uy0 */
    s_u[17].i = ei + 1; s_u[17].j = ej; s_u[17].k = ek + 1; s_u[17].c = 2; /* Uz0 */

    /* node 6 */
    s_u[18].i = ei + 1; s_u[18].j = ej + 1; s_u[18].k = ek + 1; s_u[18].c = 0; /* Ux0 */
    s_u[19].i = ei + 1; s_u[19].j = ej + 1; s_u[19].k = ek + 1; s_u[19].c = 1; /* Uy0 */
    s_u[20].i = ei + 1; s_u[20].j = ej + 1; s_u[20].k = ek + 1; s_u[20].c = 2; /* Uz0 */

    /* node 7 */
    s_u[21].i = ei; s_u[21].j = ej + 1; s_u[21].k = ek + 1; s_u[21].c = 0; /* Ux0 */
    s_u[22].i = ei; s_u[22].j = ej + 1; s_u[22].k = ek + 1; s_u[22].c = 1; /* Uy0 */
    s_u[23].i = ei; s_u[23].j = ej + 1; s_u[23].k = ek + 1; s_u[23].c = 2; /* Uz0 */
}
