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


PetscErrorCode write_pvtu(const char file_prefix[])
{
    int rank, nproc;
    char name_pvtu[PETSC_MAX_PATH_LEN];
    char name_vtu[PETSC_MAX_PATH_LEN];
    FILE *fp = NULL;
    PetscErrorCode ierr;
    PetscInt si, sj, sk, nx, ny, nz;
    PetscInt i, j, k, e;
    PetscInt n, N;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);
    ierr = PetscSNPrintf(name_pvtu,
                         sizeof(name_pvtu), "%s.pvtu",
                         file_prefix, rank); CHKERRQ(ierr);

    fp = fopen(name_pvtu,"w");

    PetscFPrintf(PETSC_COMM_SELF, fp,
                 "<?xml version=\"1.0\"?>\n"
                 "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
                 "byte_order=\"LittleEndian\">\n"
                 "<PUnstructuredGrid GhostLevel=\"0\">\n"
                 "<PPoints>\n"
                 "  <PDataArray type=\"Float64\" Name=\"Position\" "
                 "NumberOfComponents=\"3\"/>\n"
                 "</PPoints>\n"
                 "<PCells>\n"
                 "  <PDataArray type=\"Int32\" Name=\"connectivity\" "
                 "NumberOfComponents=\"1\"/>\n"
                 "  <PDataArray type=\"Int32\" Name=\"offsets\"      "
                 "NumberOfComponents=\"1\"/>\n"
                 "  <PDataArray type=\"UInt8\" Name=\"types\"        "
                 "NumberOfComponents=\"1\"/>\n"
                 "</PCells>\n"
                 "<PPointData Vectors=\"displ\">\n"
                 "  <PDataArray type=\"Float64\" Name=\"displ\"    "
                 "NumberOfComponents=\"3\" />\n"
                 "</PPointData>\n"
                 "<PCellData>\n"
                 "  <PDataArray type=\"Int32\"   Name=\"part\"   "
                 "NumberOfComponents=\"1\"/>\n"
//                 "<PDataArray type=\"Float64\" Name=\"strain\" "
//                 "NumberOfComponents=\"6\"/>\n"
//                 "<PDataArray type=\"Float64\" Name=\"stress\" "
//                 "NumberOfComponents=\"6\"/>\n"
                 "</PCellData>\n");

    for (i = 0; i < nproc; ++i) {
      sprintf(name_vtu, "%s-subdo-%d", file_prefix, i);
      PetscFPrintf(PETSC_COMM_SELF, fp,
                   "  <Piece Source=\"%s.vtu\"/>\n", name_vtu);
    }
    PetscFPrintf(PETSC_COMM_SELF, fp,
                 "</PUnstructuredGrid>\n"
                 "</VTKFile>\n");

    fclose(fp);

    //------------------------------------------------------------

    sprintf(name_vtu, "%s-subdo-%d.vtu", file_prefix, rank);
    fp = fopen(name_vtu, "w");

    PetscInt nelem, npe;
    const PetscInt *eix, *eixp;

    ierr = DMDAGetGhostCorners(DA, &si, &sj, &sk, &nx, &ny, &nz); CHKERRQ(ierr);
    ierr = DMDAGetElements(DA, &nelem, &npe, &eix); CHKERRQ(ierr);
    N = nx * ny * nz;

    PetscFPrintf(PETSC_COMM_SELF, fp,
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n"
            "<UnstructuredGrid>\n"
            "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n"
            "<Points>\n",
            N, nelem);

    PetscFPrintf(PETSC_COMM_SELF, fp, "<DataArray type=\"Float64\" "
                 "Name=\"Position\" NumberOfComponents=\"3\" "
                 "format=\"ascii\">\n");

    for (k = sk ; k < sk + nz; ++k) {
        for (j = sj ; j < sj + ny; ++j) {
            for (i = si ; i < si + nx; ++i) {
                PetscFPrintf(PETSC_COMM_SELF, fp, "%01.6e\t%01.6e\t%01.6e\n",
                             i * dx, j * dy, k * dz);
            }
        }
    }
    PetscFPrintf(PETSC_COMM_SELF, fp,
                 "</DataArray>\n"
                 "</Points>\n"
                 "<Cells>\n");

    PetscFPrintf(PETSC_COMM_SELF, fp,
                 "<DataArray type=\"Int32\" Name=\"connectivity\" "
                 "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for (e = 0; e < nelem ; ++e) {
        for (n = 0 ; n < npe; ++n)
            PetscFPrintf(PETSC_COMM_SELF, fp, "%-6d\t", eix[e * npe + n]);
        PetscFPrintf(PETSC_COMM_SELF, fp, "\n");
    }
    PetscFPrintf(PETSC_COMM_SELF, fp, "</DataArray>\n");


    PetscFPrintf(PETSC_COMM_SELF, fp,
                 "<DataArray type=\"Int32\" Name=\"offsets\" "
                 "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for (e = 1; e < nelem + 1; ++e) {
        PetscFPrintf(PETSC_COMM_SELF, fp, "%d\t", e * npe);
    }
    PetscFPrintf(PETSC_COMM_SELF, fp, "\n</DataArray>\n");

    PetscFPrintf(PETSC_COMM_SELF, fp,
                 "<DataArray type=\"UInt8\"  Name=\"types\" "
                 "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for (e = 0; e < nelem ; ++e) {
        PetscFPrintf(PETSC_COMM_SELF, fp, "12\t");
    }
    PetscFPrintf(PETSC_COMM_SELF, fp, "\n</DataArray>\n");

    PetscFPrintf(PETSC_COMM_SELF, fp, "</Cells>\n");
    PetscFPrintf(PETSC_COMM_SELF, fp, "<PointData Vectors=\"displ\">\n");
    PetscFPrintf(PETSC_COMM_SELF, fp,
                 "<DataArray type=\"Float64\" Name=\"displ\" "
                 "NumberOfComponents=\"3\" format=\"ascii\" >\n");

    PetscScalar *u_arr;
    ierr = VecGetArray(u, &u_arr);
    for (n = 0; n < N; ++n) {
        PetscFPrintf(PETSC_COMM_SELF, fp, "%01.6e\t%01.6e\t%01.6e\n",
                     u_arr[n * DIM + 0], u_arr[n * DIM + 1], u_arr[n * DIM + 2]);
    }
    VecRestoreArray(u, &u_arr);
    PetscFPrintf(PETSC_COMM_SELF, fp, "</DataArray>\n");
    PetscFPrintf(PETSC_COMM_SELF, fp, "</PointData>\n");
    PetscFPrintf(PETSC_COMM_SELF, fp, "<CellData>\n");

    PetscFPrintf(PETSC_COMM_SELF, fp,
                 "<DataArray type=\"Int32\" Name=\"part\" "
                 "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for (e = 0; e < nelem; ++e)
        PetscFPrintf(PETSC_COMM_SELF, fp, "%d\t", rank);
    PetscFPrintf(PETSC_COMM_SELF, fp, "\n</DataArray>\n");

    PetscFPrintf(PETSC_COMM_SELF, fp, "</CellData>\n");

    PetscFPrintf(PETSC_COMM_SELF, fp,
                 "</Piece>\n"
                 "</UnstructuredGrid>\n"
                 "</VTKFile>\n");

    fclose(fp);

    return ierr;
}
