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


int main(int argc,char **args)
{
    int ierr;
    int time_s, newton_it;
    double norm;
    char mess[64];

    ierr = PetscInitialize(&argc,&args,(char*)0,help); if(ierr) return ierr;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    ierr = init();

    sprintf(mess, "Problem size %d\n", nproc);
    print0(mess);

    for(time_s = 0; time_s < tsteps; ++time_s) {

        //ierr = set_bc(time_s);

        newton_it = 0;
        while(newton_it < NEWTON_ITS) {

            //ierr = set_strains();
            //micropp_C_homogenize();

            ierr = assembly_res(b);
            /* norm = |b| */

            ierr = assembly_jac(A);
            //ierr = solve_Ax();
            /* u = u + du */

            newton_it ++;
        }
    }

    ierr = finish();
    return ierr;
}
