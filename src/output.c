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


PetscErrorCode write_vtk(const char file_prefix[])
{
    int rank;
    char vtk_filename[PETSC_MAX_PATH_LEN];
    FILE *vtk_fp = NULL;
    PetscErrorCode ierr;
    PetscInt si, sj, sk, nx, ny, nz, i;
    PetscInt N;
    PetscInt memory_offset;
    DM cda;
    Vec coords;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    ierr = PetscSNPrintf(vtk_filename, 
                         sizeof(vtk_filename), "subdo-%s-p%1.4d.vts",
                         file_prefix, rank); CHKERRQ(ierr);
    vtk_fp = fopen(vtk_filename,"w");
    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "<?xml version=\"1.0\"?>\n");

    ierr = DMDAGetGhostCorners(DA, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);
    N = nx * ny * nz;

    PetscFPrintf(PETSC_COMM_SELF, vtk_fp,
                 "<VTKFile type=\"StructuredGrid\" version=\"0.1\" "
                 "byte_order=\"LittleEndian\">\n");
    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, 
                 "  <StructuredGrid WholeExtent=\"%D %D %D %D %D %D\">\n",
                 si, si + nx - 1, sj, sj + ny - 1, sk, sk + nz - 1);
    PetscFPrintf(PETSC_COMM_SELF, vtk_fp,
                 "    <Piece Extent=\"%D %D %D %D %D %D\">\n",
                 si, si + nx - 1, sj, sj + ny - 1, sk, sk + nz - 1);

    memory_offset = 0;

    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "      <CellData></CellData>\n");
    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "      <Points>\n");

    PetscFPrintf(PETSC_COMM_SELF, vtk_fp,
                 "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                 "format=\"appended\" offset=\"%d\" />\n", memory_offset);
    memory_offset = memory_offset + sizeof(PetscInt) + \
                    sizeof(PetscScalar) * N * DIM;

    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "      </Points>\n");
//    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "      <PointData Scalars=\" ");
//
//    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "disp");
//
//    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "\">\n");
//    PetscFPrintf(PETSC_COMM_SELF, vtk_fp,
//                 "        <DataArray type=\"Float64\" Name=\"%s\" "
//                 "format=\"appended\" offset=\"%d\"/>\n", "disp", memory_offset);
//
//    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "      </PointData>\n");
    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "    </Piece>\n");
    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "  </StructuredGrid>\n");

    Vec u_loc;
    PetscScalar *u_arr;

    ierr = DMGetLocalVector(DA, &u_loc); CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(DA, u, INSERT_VALUES, u_loc); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(DA, u,INSERT_VALUES, u_loc); CHKERRQ(ierr);
    ierr = VecGetArray(u_loc, &u_arr); CHKERRQ(ierr);

    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "  <AppendedData encoding=\"raw\">\n");
    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "_");

    ierr = DMGetCoordinateDM(DA, &cda); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(DA, &coords); CHKERRQ(ierr);

    /* write coordinates */
    {
        int length = sizeof(PetscScalar) * N * DIM;
        PetscScalar *allcoords;

        fwrite(&length, sizeof(int), 1, vtk_fp);
        ierr = VecGetArray(coords, &allcoords); CHKERRQ(ierr);
        fwrite(allcoords, sizeof(PetscScalar), 3 * N, vtk_fp);
        ierr = VecRestoreArray(coords, &allcoords); CHKERRQ(ierr);
    }

    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "\n  </AppendedData>\n");
    PetscFPrintf(PETSC_COMM_SELF, vtk_fp, "</VTKFile>\n");

    fclose(vtk_fp);
    return ierr;
}
