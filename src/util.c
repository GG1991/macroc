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


PetscErrorCode minmax_elems_across_mpis(DM da, int *min, int *max) {

	PetscErrorCode ierr;
	int _min, _max;
	PetscInt nex, ney, nez;
	ierr = DMDAGetElementsSizes(da, &nex, &ney, &nez); CHKERRQ(ierr);
	int nelem = nex * ney * nez;

	int *nelems_all = NULL;

	int rank, nproc;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	if (rank == 0)
		nelems_all = malloc(nproc * sizeof(int));

	MPI_Gather(&nelem, 1, MPI_INT, nelems_all, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		_min = nelems_all[0];
		_max = nelems_all[0];
		int i;
		for (i = 1; i < nproc; ++i) {
			if (nelems_all[i] < _min)
				_min = nelems_all[i];

			if (nelems_all[i] > _max)
				_max = nelems_all[i];
		}
	}

	MPI_Bcast(&_min, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&_max, 1, MPI_INT, 0, MPI_COMM_WORLD);

	*min = _min;
	*max = _max;

	return ierr;
}


/*
 * Returns the number of non_linear Gauss points at the micro-scale
 */

PetscInt get_non_linear_gps(void)
{
	PetscErrorCode ierr;
	PetscInt non_linear = 0;

	PetscInt mpi_non_linear = micropp_C_get_non_linear_gps();

	ierr = MPI_Reduce(&mpi_non_linear, &non_linear, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

	return non_linear;
}
