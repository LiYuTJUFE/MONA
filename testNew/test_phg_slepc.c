/*
 * =====================================================================================
 *
 *       Filename:  test_mg.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018年12月26日 15时50分13秒
 *
 *         Author:  Li Yu (liyu@tjufe.edu.cn), 
 *   Organization:  
 *
 * =====================================================================================
 */

#include <stdio.h>
#include "gcge.h"
#include "pase.h"
#include "gcge_app_slepc.h"
#include <slepceps.h>

static int num_levels = 3;
/* 这两个函数在get_mat_phg.c中 */
int CreateMatrixPHG (void **matA, void **matB, void **dofU, void **mapM, void **gridG, int argc, char *argv[]);
int DestroyMatrixPHG(void **matA, void **matB, void **dofU, void **mapM, void **gridG, int argc, char *argv[]);

static char help[] = "Test MultiGrid.\n";
void PETSCPrintMat(Mat A, char *name);
void PETSCPrintVec(Vec x);
void PETSCPrintBV (BV  x, char *name);

void Multigrid_LinearSolverCreate(PASE_MULTIGRID *multi_grid, GCGE_OPS *gcge_ops, 
      int num_vecs, void *A, void *B);
void GCGE_SOLVER_SetMultigridLinearSolver(GCGE_SOLVER *solver, PASE_MULTIGRID multi_grid);
void GCGE_MultiGrid_LinearSolver(void *Matrix, void **b, void **x, int *start, int *end, struct GCGE_OPS_ *ops);
void Multigrid_LinearSolverDestroy(PASE_MULTIGRID *multi_grid, int num_vecs);
/* 
 *  Description:  测试PASE_MULTIGRID
 */

int
main ( int argc, char *argv[] )
{
   /* PetscInitialize */
   SlepcInitialize(&argc,&argv,(char*)0,help);
   PetscErrorCode ierr;

   /* 测试矩阵声明 */
   Mat      petsc_mat_A, petsc_mat_B;

   /* 得到一个PHG矩阵, 并将之转换为PETSC矩阵 */
   void *phg_mat_A, *phg_mat_B, *phg_dof_U, *phg_map_M, *phg_grid_G;
   CreateMatrixPHG (&phg_mat_A, &phg_mat_B, &phg_dof_U, &phg_map_M, &phg_grid_G, argc, argv);
   MatrixConvertPHG2PETSC((void **)(&petsc_mat_A), &phg_mat_A);
   MatrixConvertPHG2PETSC((void **)(&petsc_mat_B), &phg_mat_B);
   DestroyMatrixPHG(&phg_mat_A, &phg_mat_B, &phg_dof_U, &phg_map_M, &phg_grid_G, argc, argv);

#if 0
   /* MAT_GLOBAL_SUM MAT_GLOBAL_MAX MAT_LOCAL*/
   PetscInt nrows, ncols;
   MatInfo mat_info;
   MatGetInfo(petsc_mat_A, MAT_GLOBAL_SUM, &mat_info);
   MatGetSize(petsc_mat_A, &nrows, &ncols);
   printf("nz_used of A for global = %g, nrows = %d, ncols = %d\n", 
	 mat_info.nz_used, nrows, ncols);
   MatGetInfo(petsc_mat_A, MAT_LOCAL, &mat_info);
   MatGetLocalSize(petsc_mat_A, &nrows, &ncols);
   printf("nz_used of A for local = %g, nrows = %d, ncols = %d\n", 
	 mat_info.nz_used, nrows, ncols);
   MatGetInfo(petsc_mat_B, MAT_GLOBAL_SUM, &mat_info);
   MatGetSize(petsc_mat_B, &nrows, &ncols);
   printf("nz_used of B for global = %g, nrows = %d, ncols = %d\n", 
	 mat_info.nz_used, nrows, ncols);
   MatGetInfo(petsc_mat_B, MAT_LOCAL, &mat_info);
   MatGetLocalSize(petsc_mat_B, &nrows, &ncols);
   printf("nz_used of B for local = %g, nrows = %d, ncols = %d\n", 
	 mat_info.nz_used, nrows, ncols);
#endif

   /* slepc 求解器 */
   EPS            eps;
   ST             st;          /* spectral transformation context */
   EPSType        type;
   PetscBool      flag,terse;
   PetscInt        nev;
   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Create the eigensolver and set various options
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   /*
      Create eigensolver context
      */
   //printf("EPSCreate\n");
   ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

   /*
      Set operators. In this case, it is a generalized eigenvalue problem
      */
   //printf("EPSSetOperators\n");
   ierr = EPSSetOperators(eps,petsc_mat_A,petsc_mat_B);CHKERRQ(ierr);
   ierr = EPSSetProblemType(eps,EPS_GHEP);CHKERRQ(ierr);

   /*
      Set solver parameters at runtime
      */
   //printf("EPSSetFromOptions\n");
   ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

   ierr = PetscObjectTypeCompareAny((PetscObject)eps,&flag,EPSBLOPEX,EPSLOBPCG,EPSRQCG,"");CHKERRQ(ierr);
   if (flag) {
      ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL);CHKERRQ(ierr);
   } else {
      /*
	 Select portion of spectrum
	 */
      ierr = EPSSetTarget(eps,0.0);CHKERRQ(ierr);
      ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE);CHKERRQ(ierr);
      /*
	 Use shift-and-invert to avoid solving linear systems with a singular B
	 in case nulldim>0
	 */
      ierr = PetscObjectTypeCompareAny((PetscObject)eps,&flag,EPSGD,EPSJD,"");CHKERRQ(ierr);
      if (!flag) {
	 ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
	 ierr = STSetType(st,STSINVERT);CHKERRQ(ierr);
      }
   }

   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Solve the eigensystem
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   //printf("EPSSolve\n");
   ierr = EPSSolve(eps);CHKERRQ(ierr);

   /*
Optional: Get some information from the solver and display it
*/
   ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
   ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);

   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Display solution and clean up
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   /* show detailed info unless -terse option is given by user */
   ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
   if (terse) {
      ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
   } else {
      ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
      ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
   }
   ierr = EPSDestroy(&eps);CHKERRQ(ierr);

   ierr = MatDestroy(&petsc_mat_A);
   ierr = MatDestroy(&petsc_mat_B);

   /* PetscFinalize */
   ierr = SlepcFinalize();
   return 0;
}
void PETSCPrintMat(Mat A, char *name)
{
    PetscErrorCode ierr;
    int nrows = 0, ncols = 0;
    ierr = MatGetSize(A, &nrows, &ncols);
    int row = 0, col = 0;
    const PetscInt *cols;
    const PetscScalar *vals;
    for(row=0; row<nrows; row++) {
        ierr = MatGetRow(A, row, &ncols, &cols, &vals);
        for(col=0; col<ncols; col++) {
            GCGE_Printf("%s(%d, %d) = %18.15e;\n", name, row+1, cols[col]+1, vals[col]);
        }
        ierr = MatRestoreRow(A, row, &ncols, &cols, &vals);
    }
}

void PETSCPrintVec(Vec x)
{
    PetscErrorCode ierr;
    int size = 0;
    int i = 0;
    ierr = VecGetSize(x, &size);
    const PetscScalar *array;
    ierr = VecGetArrayRead(x, &array);
    for(i=0; i<size; i++)
    {
        GCGE_Printf("%18.15e\t", array[i]);
    }
    ierr = VecRestoreArrayRead(x, &array);
    GCGE_Printf("\n");
}

void PETSCPrintBV(BV x, char *name)
{
    PetscErrorCode ierr;
    int n = 0, i = 0;
    ierr = BVGetSizes(x, NULL, NULL, &n);
    Vec xs;
    GCGE_Printf("%s = [\n", name);
    for(i=0; i<n; i++)
    {
        ierr = BVGetColumn(x, i, &xs);
	PETSCPrintVec(xs);
        ierr = BVRestoreColumn(x, i, &xs);
    }
    GCGE_Printf("];\n");
}
