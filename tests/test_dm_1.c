
static char help[] = "Tests VecView() contour plotting for 2d DMDAs.\n\n";

/*
  VecView() on DMDA vectors first puts the Vec elements into global natural ordering before printing (or plotting)
them. In 2d 5 by 2 DMDA this means the numbering is

     5  6   7   8   9
     0  1   2   3   4

Now the default split across 2 processors with the DM  is (by rank)

    0  0   0  1   1
    0  0   0  1   1

So the global PETSc ordering is

    3  4  5   8  9
    0  1  2   6  7

Use the options
     -da_grid_x <nx> - number of grid points in x direction, if M < 0
     -da_grid_y <ny> - number of grid points in y direction, if N < 0
     -da_processors_x <MX> number of processors in x directio
     -da_processors_y <MY> number of processors in x direction
*/

#include <petscdm.h>
#include <petscdmda.h>

int main(int argc,char **argv)
{
    PetscMPIInt      rank;
    PetscInt         M = 10,N = 8;
    PetscErrorCode   ierr;
    PetscBool        flg = PETSC_FALSE;
    DM               da;
    PetscViewer      viewer;
    Vec              local,global;
    PetscScalar      value;
    DMBoundaryType   bx    = DM_BOUNDARY_NONE,by = DM_BOUNDARY_NONE;
    DMDAStencilType  stype = DMDA_STENCIL_BOX;

    ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
    ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,0,"",300,0,300,300,&viewer);
    CHKERRQ(ierr);

    ierr = PetscOptionsGetBool(NULL,NULL,"-star_stencil",&flg,NULL);CHKERRQ(ierr);
    if (flg) stype = DMDA_STENCIL_STAR;

    /* Create distributed array and get vectors */
    ierr = DMDACreate2d(PETSC_COMM_WORLD, bx, by, stype, M, N,
                        PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &da);
    CHKERRQ(ierr);
    ierr = DMSetFromOptions(da);CHKERRQ(ierr);
    ierr = DMSetUp(da);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&global);CHKERRQ(ierr);
    ierr = DMCreateLocalVector(da,&local);CHKERRQ(ierr);

    ierr  = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    value = rank + 1;
    if(rank == 0) {
        VecSetValue(global,1,value, ADD_VALUES);
    } else if(rank == 1) {
        VecSetValue(global,8,value, ADD_VALUES);
    }

    VecAssemblyBegin(global);
    VecAssemblyEnd(global);

    int x0, y0, z0, dx, dy, dz;
    DMDAGetCorners(da,&x0,&y0,&z0,&dx,&dy,&dz);
    int nel, npe;
    const PetscInt *e;
    DMDAGetElements(da,&nel,&npe,&e);
    //DMDAGetGhostCorners(da,&x0,&y0,&z0,&dx,&dy,&dz);
    printf("rank:%d x0:%d\ty0:%d\tz0:%d\tdx:%d\tdy:%d\tdz:%d\n",rank,x0,y0,z0,dx,dy,dz);
    printf("rank:%d nel:%d\tnpe:%d\n",rank,nel,npe);
    int i, d;
    for(i = 0; i < nel; ++i){
        printf("e %d : %d %d %d %d\n",i,e[npe*i],e[npe*i+1],e[npe*i+2],e[npe*i+3]);
    }
    printf("\n");

    flg  = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,NULL, "-view_global", &flg,NULL);CHKERRQ(ierr);
    if (flg) { /* view global vector in natural ordering */
        ierr = VecView(global,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
    ierr = DMView(da,viewer);CHKERRQ(ierr);
    ierr = VecView(global,viewer);CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = VecDestroy(&local);CHKERRQ(ierr);
    ierr = VecDestroy(&global);CHKERRQ(ierr);
    ierr = DMDestroy(&da);CHKERRQ(ierr);
    ierr = PetscFinalize();
    return ierr;
}
