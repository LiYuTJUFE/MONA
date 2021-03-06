/**
 *    @file  ops_eig_sol_gcg.c
 *   @brief  特征值求解器 GCG 
 *
 *  特征值求解器 GCG
 *
 *  @author  Yu Li, liyu@tjufe.edu.cn
 *
 *       Created:  2020/8/18
 *      Revision:  none
 */

#include	<stdio.h>
#include	<stdlib.h>
#include    <math.h>
#include    <memory.h>
#include    <assert.h>
#include    <time.h>
#include    <string.h> 

#include    "ops_eig_sol_gcg.h"

#define     DEBUG 0
#define     TIME_GCG 1

typedef struct TimeGCG_ {
	double initX_time;
	double checkconv_time;
	double compP_time;
	double compRR_time;
	double rr_matW_time;
	double dsyevx_time;
	double compRV_time;
	double compW_time;
	double compX_time;            
	double linsol_time;
    double time_total;
} TimeGCG;

struct TimeGCG_ time_gcg = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};


static int sizeN, startN, endN;
static int sizeP, startP, endP;
static int sizeW, startW, endW;
static int sizeC, sizeX , sizeV, endX;


static void   **mv_ws[3]; 
static double *dbl_ws; 
static int    *int_ws;
static struct OPS_ *ops_gcg;
static struct GCGSolver_ *gcg_solver;

static void InitializeX(void **V, void **ritz_vec, void *B, int nevGiven)
{	
#if TIME_GCG
#if USE_MPI
    time_gcg.initX_time -= MPI_Wtime();
#else
    time_gcg.initX_time -= clock();
#endif
#endif
	int start[2], end[2];	
	start[0] = 0; end[0] = nevGiven;
	start[1] = 0; end[1] = nevGiven;
	ops_gcg->MultiVecAxpby(1.0,ritz_vec,0,V,start,end,ops_gcg);	
#if DEBUG
	ops_gcg->Printf("V\n");	
	ops_gcg->MultiVecView(V,0,sizeV,ops_gcg);
#endif
	ops_gcg->Printf("sizeX = %d, nevGiven = %d, %s\n",
		sizeX,nevGiven,gcg_solver->initX_orth_method);
	if (0 == strcmp("mgs", gcg_solver->initX_orth_method))
			MultiVecOrthSetup_ModifiedGramSchmidt(
				gcg_solver->initX_orth_block_size,
				gcg_solver->initX_orth_max_reorth,
				gcg_solver->initX_orth_zero_tol,
				//mv_ws[0],dbl_ws,ops_gcg);
				ritz_vec,gcg_solver->dbl_ws,ops_gcg);
	else if (0 == strcmp("bgs", gcg_solver->initX_orth_method))
			MultiVecOrthSetup_BinaryGramSchmidt(
				gcg_solver->initX_orth_block_size,
				gcg_solver->initX_orth_max_reorth,
				gcg_solver->initX_orth_zero_tol,
				//mv_ws[0],dbl_ws,ops_gcg);
				ritz_vec,gcg_solver->dbl_ws,ops_gcg);
	else
			MultiVecOrthSetup_ModifiedGramSchmidt(
				gcg_solver->initX_orth_block_size,
				gcg_solver->initX_orth_max_reorth,
				gcg_solver->initX_orth_zero_tol,
				//mv_ws[0],dbl_ws,ops_gcg);
				ritz_vec,gcg_solver->dbl_ws,ops_gcg);

	ops_gcg->MultiVecOrth(V,0,&nevGiven,B,ops_gcg);
	ops_gcg->Printf("sizeX = %d, nevGiven = %d\n",sizeX,nevGiven);
#if DEBUG
	ops_gcg->Printf("before Orth\n");
#endif
	ops_gcg->MultiVecSetRandomValue(V,nevGiven,sizeX,ops_gcg);
	//MultiVecOrthSetup_BinaryGramSchmidt(2,1e-16,mv_ws[0],dbl_ws,ops_gcg);
	//ops_gcg->MultiVecOrth(V,nevGiven,&endX,B,ops_gcg);
	//MultiVecOrthSetup_ModifiedGramSchmidt(2,1e-16,mv_ws[0],dbl_ws,ops_gcg);
	ops_gcg->MultiVecOrth(V,nevGiven,&endX,B,ops_gcg);
#if DEBUG
	ops_gcg->Printf("after Orth\n");
#endif
	assert(endX==sizeX);
	/* 多次正交化, 保证有 sizeX 个正交向量 */
	//int pre_endX;	
	//while (endX < sizeX) {
	//	ops_gcg->MultiVecSetRandomValue(V,endX,sizeX,ops_gcg);
	//	pre_endX = endX; endX = sizeX;
	//	ops_gcg->MultiVecOrth(V,pre_endX,&endX,B,ops_gcg);
	//}
#if DEBUG
	ops_gcg->Printf("nevGiven = %d\n",nevGiven);
	ops_gcg->MultiVecView(V,0,endX,ops_gcg);
#endif
#if TIME_GCG
#if USE_MPI
    time_gcg.initX_time += MPI_Wtime();
#else
    time_gcg.initX_time += clock();
#endif
#endif
	return;
}
static void ComputeRitzVec(void **ritz_vec, void **V, double *ss_evec)
{
#if TIME_GCG
#if USE_MPI
    time_gcg.compRV_time -= MPI_Wtime();
#else
    time_gcg.compRV_time -= clock();
#endif
#endif 
	int start[2], end[2]; double *coef;
	start[0] = startN; end[0] = endW;
	start[1] = startN; end[1] = endX;
	coef     = ss_evec;
#if DEBUG
	ops_gcg->Printf("startN = %d, endW = %d, endX = %d\n",startN,endW,endX);	
	ops_gcg->Printf("coef: (%d * %d)\n",end[0]-start[0],end[1]-start[1]);	
	int row, col;
	for (row = 0; row < end[0]-start[0]; ++row) {
		for (col = 0; col < end[1]-start[1]; ++col) {
			ops_gcg->Printf("%6.4e\t",coef[row+col*(sizeV-sizeC)]);
		}	
		ops_gcg->Printf("\n");
	}
	ops_gcg->Printf("V:\n");
	ops_gcg->MultiVecView(V,start[0],end[0],ops_gcg);
	ops_gcg->Printf("startN = %d, endW = %d, endX = %d\n",startN,endW,endX);	
	ops_gcg->Printf("V = %p, (%d, %d), ritz_vec = %p (%d, %d)\n",V,start[0],end[0],ritz_vec,start[1],end[1]);	
#endif
	ops_gcg->MultiVecLinearComb(V,ritz_vec,0, 
			start,end,coef,sizeV-sizeC,NULL,0,ops_gcg);
			
#if DEBUG	
	ops_gcg->Printf("ritz vec:\n");
	ops_gcg->MultiVecView(ritz_vec,start[1],end[1],ops_gcg);
#endif
#if TIME_GCG
#if USE_MPI
    time_gcg.compRV_time += MPI_Wtime();
#else
    time_gcg.compRV_time += clock();
#endif
#endif
	return;
}
static int CheckConvergence(void *A, void *B, double *ss_eval, void **ritz_vec, 
	int numCheck, double *tol, int *offset)
{
#if TIME_GCG
#if USE_MPI
    time_gcg.checkconv_time -= MPI_Wtime();
#else
    time_gcg.checkconv_time -= clock();
#endif
#endif
#if DEBUG
	ops_gcg->Printf("numCheck = %d\n", numCheck);
#endif
	int start[2], end[2], idx; double *inner_prod;
	int nevConv;
	start[0] = startN; end[0] = start[0]+numCheck;
	start[1] = 0     ; end[1] = numCheck;	
	ops_gcg->MatDotMultiVec(A,ritz_vec,mv_ws[0],start,end,ops_gcg);	
	ops_gcg->MatDotMultiVec(B,ritz_vec,mv_ws[1],start,end,ops_gcg);	
	/* lambda Bx */
	ops_gcg->MultiVecLinearComb(NULL,mv_ws[1],0,start,end,
			NULL,0,ss_eval+startN,1,ops_gcg);
	start[0] = 0     ; end[0] = numCheck;
	start[1] = 0     ; end[1] = numCheck;
	/* Ax - lambda Bx */
	ops_gcg->MultiVecAxpby(-1.0,mv_ws[1],1.0,mv_ws[0],start,end,ops_gcg);
	/* 不使用 ss_evec 部分 */
	inner_prod = dbl_ws+(sizeV-sizeC)*sizeW;
	ops_gcg->MultiVecInnerProd('D',mv_ws[0],mv_ws[0],0,
			start,end,inner_prod,1,ops_gcg);
	for (idx = 0; idx < numCheck; ++idx) {
		inner_prod[idx] = sqrt(inner_prod[idx]);
		ops_gcg->Printf("GCG: [%d] %6.14e (%6.4e, %6.4e)\n",
			startN+idx,ss_eval[startN+idx],
			inner_prod[idx], inner_prod[idx]/fabs(ss_eval[startN+idx]));
	}
	for (idx = 0; idx < numCheck; ++idx) {
		/* 绝对残量 和 相对残量 需分别小于 tol[0] 和 tol[1] */
		if (inner_prod[idx] > tol[0] || 
				inner_prod[idx] > fabs(ss_eval[startN+idx])*tol[1]) break;
	}	
	for ( ; idx > 0; --idx) {
		/* 最后一个收敛的特征值与第一个不收敛的特征值不是重根 */
		if ( fabs((ss_eval[startN+idx-1]-ss_eval[startN+idx])/ss_eval[startN+idx-1]) 
		      > gcg_solver->gapMin) {
			break;
		}
	}
	nevConv = sizeC+idx;
	
	/* offset[0] 为未收敛块的个数, offset[2n-1] <= idx < offset[2n]
	 * idx 是不收敛的标号 1 <= n <= offset[0] */
	int state, num_unconv;
	/* 1 1 0 0 1 1 1 1 0 0 1 0 1 0 0 0 0 0 0 */
	offset[0] = 0; state = 1; num_unconv = 0; 
	for (idx = 0; idx < numCheck; ++idx) {
		/* 这一个是不收敛的 */
		if (inner_prod[idx] > tol[0] || 
				inner_prod[idx] > fabs(ss_eval[startN+idx])*tol[1]) {
			/* 上一个是收敛的 */
			if (state) {
				offset[ offset[0]*2+1 ] = startN+idx;
				state = 0;
			}
			++num_unconv;
			if (num_unconv == sizeN) {
				offset[ offset[0]*2+2 ] = startN+idx+1;
				++offset[0];
				break;
			}
		}
		else {
			/* 上一个是不收敛的 */
			if (!state) {
				offset[ offset[0]*2+2 ] = startN+idx;
				++offset[0];
				state = 1;
			}
		}
	}
	if (num_unconv < sizeN) {
		if (state == 1) {
			offset[ offset[0]*2+1 ] = startN+numCheck;	
		}
		offset[ offset[0]*2+2 ] = startN+numCheck+sizeN-num_unconv;
		offset[ offset[0]*2+2 ] = offset[ offset[0]*2+2 ] < endX?
				offset[ offset[0]*2+2 ]:endX;
		assert(offset[ offset[0]*2+1 ]<offset[ offset[0]*2+2 ]);
		++offset[0];	
	}
	
#if TIME_GCG
#if USE_MPI
    time_gcg.checkconv_time += MPI_Wtime();
#else
    time_gcg.checkconv_time += clock();
#endif
#endif
	for (idx = 0; idx < offset[0]; ++idx) {
		ops_gcg->Printf("offset [%d,%d)\n",
			offset[idx*2+1],offset[idx*2+2]);
	}
	assert(offset[0]>0); 
	return nevConv;
}
static void ComputeP(void **V, double *ss_evec, int *offset)
{
#if TIME_GCG
#if USE_MPI
    time_gcg.compP_time -= MPI_Wtime();
#else
    time_gcg.compP_time -= clock();
#endif
#endif
	int length, incx, incy, ldm, block_size;
	int nrows, idx, col, start[2], end[2]; 
	double *source, *destin, *mat, *coef;
	
	/* 复制 n 部分对应的列 */
	ops_gcg->Printf("offset[0] = %d, sizeP = %d\n", offset[0], sizeP);	
	block_size = 0;
	for (idx = 0; idx < offset[0]; ++idx) {
		length = (sizeV-sizeC)*(offset[idx*2+2]-offset[idx*2+1]);		
		source = ss_evec+(sizeV-sizeC)*(offset[idx*2+1]-sizeC) ; incx = 1;
		destin = ss_evec+(sizeV-sizeC)*(sizeX-sizeC+block_size); incy = 1;
		dcopy(&length,source,&incx,destin,&incy);		
		block_size += offset[idx*2+2]-offset[idx*2+1];
		ops_gcg->Printf("offset [%d, %d)\n", offset[idx*2+1],offset[idx*2+2]);	
	}
	sizeP = block_size;
	/* 置零 np 部分 */
	for (idx = 0; idx < offset[0]; ++idx) {
		length = (offset[idx*2+2]-offset[idx*2+1]);
		destin = ss_evec+(sizeV-sizeC)*(sizeX-sizeC)+(offset[idx*2+1]-sizeC);
		for (col = 0; col < sizeP; ++col) {
			memset(destin,0,length*sizeof(double));
			destin += sizeV-sizeC;
		}			
	}
	
	
	/* 小规模正交化 */
	mat    = ss_evec; 
	nrows  = sizeV-sizeC; ldm  = sizeV-sizeC ;
	startP = sizeX-sizeC; endP = startP+sizeP;
#if DEBUG	
	ops_gcg->Printf("sizeC = %d, sizeN = %d, sizeX = %d, sizeP = %d, sizeW = %d\n",
			sizeC,sizeN,sizeX,sizeP,sizeW);	
	ops_gcg->Printf("startP = %d, endP = %d, startW = %d, endW = %d, sizeV = %d\n",
			startP,endP,startW,endW,sizeV);
	int row, col, ncols;
	for (row = 0; row < nrows; ++row) {
		for (col = 0; col < endP; ++col) {
			ops_gcg->Printf("%6.4e\t",mat[row+ldm*col]);
		}
		ops_gcg->Printf("\n");
	}
#endif

	ops_gcg->Printf("startP = %d, endP = %d, sizeP = %d, startW = %d, endW = %d, sizeW = %d, sizeV = %d\n",
			startP,endP,sizeP,startW,endW,sizeW,sizeV);
	double *orth_dbl_ws = ss_evec+ldm*endP;
	/* ss_diag ss_matA ss_evec 剩下的空间 */
	if (0 == strcmp("bqr", gcg_solver->compP_orth_method)) {
		int length_orth_dbl_ws = sizeC
			+2*(2*sizeV*sizeC-sizeC*sizeC)
			+(sizeV-sizeC)*sizeW
			+9*sizeV+sizeX*sizeP; 
		ops_gcg->DenseMatOrth(mat,nrows,ldm,startP,&endP,
			gcg_solver->compP_orth_zero_tol,
			orth_dbl_ws,length_orth_dbl_ws,int_ws);		
	}
	else {
		LAPACKVEC lapack_vec_P, lapack_vec_ws;
		lapack_vec_P.data   = mat;
		lapack_vec_P.ldd    = ldm;
		lapack_vec_P.ncols  = endP;
		lapack_vec_P.nrows  = nrows;
		
		lapack_vec_ws.data  = orth_dbl_ws;
		lapack_vec_ws.ldd   = ldm;
		lapack_vec_ws.ncols = endP-startP;
		lapack_vec_ws.nrows = nrows;
		if (0 == strcmp("mgs", gcg_solver->compP_orth_method))
			MultiVecOrthSetup_ModifiedGramSchmidt(
				gcg_solver->compP_orth_block_size,
				gcg_solver->compP_orth_max_reorth,
				gcg_solver->compP_orth_zero_tol,
				(void*)&lapack_vec_ws,orth_dbl_ws+ldm*(endP-startP),
				ops_gcg->lapack_ops);
		else if (0 == strcmp("mgs", gcg_solver->compP_orth_method))
			MultiVecOrthSetup_BinaryGramSchmidt(
				gcg_solver->compP_orth_block_size,
				gcg_solver->compP_orth_max_reorth,
				gcg_solver->compP_orth_zero_tol,
				(void*)&lapack_vec_ws,orth_dbl_ws+ldm*(endP-startP),
				ops_gcg->lapack_ops);
		else
			MultiVecOrthSetup_ModifiedGramSchmidt(
				gcg_solver->compP_orth_block_size,
				gcg_solver->compP_orth_max_reorth,
				gcg_solver->compP_orth_zero_tol,
				(void*)&lapack_vec_ws,orth_dbl_ws+ldm*(endP-startP),
				ops_gcg->lapack_ops);
					
		ops_gcg->lapack_ops->MultiVecOrth((void*)&lapack_vec_P,
			startP,&endP,NULL,ops_gcg->lapack_ops);		
	}
	startP += sizeC; endP += sizeC; sizeP = endP-startP;
	
#if DEBUG
	ops_gcg->Printf("startP = %d, endP = %d, sizeP = %d, startW = %d, endW = %d, sizeW = %d, sizeV = %d\n",
			startP,endP,sizeP,startW,endW,sizeW,sizeV);
				
	nrows = sizeV-sizeC; ncols = sizeV-sizeC;
	for (row = 0; row < nrows; ++row) {
		for (col = 0; col < ncols; ++col) {
			ops_gcg->Printf("%6.4e\t",mat[row+ldm*col]);
		}
		ops_gcg->Printf("\n");
	}
	ops_gcg->Printf("matT mat\n");
	ops_gcg->DenseMatQtAP('N','N',sizeV-sizeC,sizeV-sizeC,sizeV-sizeC,sizeV-sizeC,
			1.0,mat   ,sizeV-sizeC,
			    NULL  ,sizeV-sizeC,
			    mat   ,sizeV-sizeC,
			0.0,dbl_ws,sizeV-sizeC,
			dbl_ws+(sizeV-sizeC)*(sizeV-sizeC));
	for (row = 0; row < sizeV-sizeC; ++row) {
		for (col = 0; col < sizeV-sizeC; ++col) {
			ops_gcg->Printf("%6.4e\t",dbl_ws[row+(sizeV-sizeC)*col]);
		}
		ops_gcg->Printf("\n");
	}
#endif
	/* 更新 P */
	start[0] = startN; end[0] = endW ;
	start[1] = 0     ; end[1] = sizeP;
	coef     = ss_evec+(sizeV-sizeC)*(sizeX-sizeC);
	ops_gcg->MultiVecLinearComb(V,mv_ws[0],0,start,end,
			coef,sizeV-sizeC,NULL,0,ops_gcg);
	start[0] = 0     ; end[0] = sizeP;
	start[1] = startP; end[1] = endP ;
	ops_gcg->MultiVecAxpby(1.0,mv_ws[0],0.0,V,start,end,ops_gcg);
	
#if DEBUG
	start[0] = startP; end[0] = endP;
	start[1] = startP; end[1] = endP;
	nrows = end[0]-start[0]; ncols = end[1]-start[1];
	ops_gcg->Printf("PtBP\n");
	ops_gcg->MultiVecQtAP('N','N',V,NULL,V,0,start,end,dbl_ws,nrows,mv_ws[0],ops_gcg);
	for (row = 0; row < nrows; ++row) {
		for (col = 0; col < ncols; ++col) {
			ops_gcg->Printf("%6.4e\t",dbl_ws[row+nrows*col]);
		}
		ops_gcg->Printf("\n");
	}
#endif
#if TIME_GCG
#if USE_MPI
    time_gcg.compP_time += MPI_Wtime();
#else
    time_gcg.compP_time += clock();
#endif
#endif	
	return;	
}
static void ComputeX(void **V, void **ritz_vec)
{
#if TIME_GCG
#if USE_MPI
    time_gcg.compX_time -= MPI_Wtime();
#else
    time_gcg.compX_time -= clock();
#endif
#endif
	int start[2], end[2];
	start[0] = startN; end[0] = endX;
	start[1] = startN; end[1] = endX;
	ops_gcg->MultiVecAxpby(1.0,ritz_vec,0.0,V,start,end,ops_gcg);
#if TIME_GCG
#if USE_MPI
    time_gcg.compX_time += MPI_Wtime();
#else
    time_gcg.compX_time += clock();
#endif
#endif
	return;
}
static void ComputeW(void **V, void *A, void *B, 
	double *ss_eval, void **ritz_vec, int *offset)
{
#if TIME_GCG
#if USE_MPI
    time_gcg.compW_time -= MPI_Wtime();
#else
    time_gcg.compW_time -= clock();
#endif
#endif	
	void **b = ritz_vec;
	int start[2], end[2], block_size, length, inc, idx;
	double *destin = dbl_ws;
	
	/* initialize */
	block_size = 0; startW = endP;	
	for (idx = 0; idx < offset[0]; ++idx) {
		length   = offset[idx*2+2]-offset[idx*2+1];
		/* initialize x */
		start[0] = offset[idx*2+1]  ; end[0] = offset[idx*2+2];
		start[1] = startW+block_size; end[1] = start[1]+length;
		ops_gcg->MultiVecAxpby(1.0,ritz_vec,0.0,V,start,end,ops_gcg);
#if DEBUG
		ops_gcg->Printf("initial W:\n");		
		ops_gcg->MultiVecView(V,start[1],end[1],ops_gcg);	
#endif
		/* set b, b = lambda Bx */
		start[0] = offset[idx*2+1]     ; end[0] = offset[idx*2+2];
		start[1] = offset[1]+block_size; end[1] = start[1]+length;
		ops_gcg->MatDotMultiVec(B,V,b,start,end,ops_gcg);
		ops_gcg->MultiVecLinearComb(NULL,b,0,start,end,
				NULL,0,ss_eval+start[0],1,ops_gcg);
				
		if (gcg_solver->user_defined_multi_linear_solver==0) {			
			inc = 1; 
		    dcopy(&length,ss_eval+start[0],&inc,destin,&inc);
			destin += length;
		}
		
		block_size += length;
#if DEBUG
		ops_gcg->Printf("initial b:\n");		
		ops_gcg->MultiVecView(b,start[1],end[1],ops_gcg);	
#endif	
	}
	endW = startW+block_size;	
	
	/* solve x */
	start[0] = offset[1]; end[0] = start[0]+block_size;
	start[1] = startW   ; end[1] = endW               ;
	/* 是否与 pase 中的线性解法器冲突 */
	if (gcg_solver->user_defined_multi_linear_solver==0) {	        
		MultiLinearSolverSetup_BlockPCG(
			gcg_solver->compW_cg_max_iter,
			gcg_solver->compW_cg_rate,
			gcg_solver->compW_cg_tol,
			gcg_solver->compW_cg_tol_type,
			mv_ws,dbl_ws,int_ws,NULL,ops_gcg);
	}
#if TIME_GCG
#if USE_MPI
    	time_gcg.linsol_time -= MPI_Wtime();
#else
    	time_gcg.linsol_time -= clock();
#endif
#endif
	ops_gcg->MultiLinearSolver(A,b,V,start,end,ops_gcg);
#if TIME_GCG
#if USE_MPI
    	time_gcg.linsol_time += MPI_Wtime();
#else
    	time_gcg.linsol_time += clock();
#endif
#endif


#if DEBUG
	ops_gcg->Printf("W = inv(A) b:\n");		
	ops_gcg->MultiVecView(V,startW,endW,ops_gcg);	
#endif
	/* orth W in V */
	//MultiVecOrthSetup_ModifiedGramSchmidt(2,1e-16,mv_ws[0],dbl_ws,ops_gcg);
	if (0 == strcmp("mgs", gcg_solver->compW_orth_method))
		MultiVecOrthSetup_ModifiedGramSchmidt(
			gcg_solver->compW_orth_block_size,
			gcg_solver->compW_orth_max_reorth,
			gcg_solver->compW_orth_zero_tol,
			mv_ws[0],dbl_ws,ops_gcg);
	else if (0 == strcmp("bgs", gcg_solver->compW_orth_method))
		MultiVecOrthSetup_BinaryGramSchmidt(
			gcg_solver->compW_orth_block_size,
			gcg_solver->compW_orth_max_reorth,
			gcg_solver->compW_orth_zero_tol,
			mv_ws[0],dbl_ws,ops_gcg);
	else
		MultiVecOrthSetup_ModifiedGramSchmidt(
			gcg_solver->compW_orth_block_size,
			gcg_solver->compW_orth_max_reorth,
			gcg_solver->compW_orth_zero_tol,
			mv_ws[0],dbl_ws,ops_gcg);
	
	ops_gcg->MultiVecOrth(V,startW,&endW,B,ops_gcg);
#if DEBUG
	ops_gcg->Printf("Orth W in V, %d, %d\n", startW,endW);		
	ops_gcg->MultiVecView(V,startW,endW,ops_gcg);
	start[0] = startW; end[0] = endW;
	start[1] = startW; end[1] = endW;
	int nrows = end[0]-start[0], ncols = end[1]-start[1], row, col;
	ops_gcg->Printf("WtBW\n");
	ops_gcg->MultiVecQtAP('N','N',V,B,V,0,start,end,dbl_ws,nrows,mv_ws[0],ops_gcg);
	for (row = 0; row < nrows; ++row) {
		for (col = 0; col < ncols; ++col) {
			ops_gcg->Printf("%6.4e\t",dbl_ws[row+nrows*col]);
		}
		ops_gcg->Printf("\n");
	}	
#endif
	//MultiVecOrthSetup_BinaryGramSchmidt(2,1e-16,mv_ws[0],dbl_ws,ops_gcg);
	//ops_gcg->MultiVecOrth(V,startW,&endW,B,ops_gcg);
#if DEBUG
	ops_gcg->Printf("WtBW\n");
	ops_gcg->MultiVecQtAP('N','N',V,B,V,0,start,end,dbl_ws,nrows,mv_ws[0],ops_gcg);
	for (row = 0; row < nrows; ++row) {
		for (col = 0; col < ncols; ++col) {
			ops_gcg->Printf("%6.4e\t",dbl_ws[row+nrows*col]);
		}
		ops_gcg->Printf("\n");
	}
	ops_gcg->Printf("W %d, %d\n",startW,endW);
	ops_gcg->MultiVecView(V,startW,endW,ops_gcg);	
	ops_gcg->Printf("Orth W in V, %d, %d\n", startW,endW);
#endif
	sizeW = endW-startW;
#if 0	
	if (sizeW<block_size) {
		ops_gcg->MultiVecSetRandomValue(V,endW,startW+block_size,ops_gcg);
		endW = startW+block_size;
		ops_gcg->MultiVecOrth(V,startW+sizeW,&endW,B,ops_gcg);
	}
	sizeW = endW-startW;
#endif

#if TIME_GCG
#if USE_MPI
    	time_gcg.compW_time += MPI_Wtime();
#else
    	time_gcg.compW_time += clock();
#endif
#endif	
	return;
}
static void ComputeRayleighRitz(double *ss_matA, double *ss_eval, double *ss_evec, double tol,
		int nevConv, double *ss_diag, void *A, void **V)
{
#if TIME_GCG
#if USE_MPI
    time_gcg.compRR_time -= MPI_Wtime();
#else
    time_gcg.compRR_time -= clock();
#endif
#endif	
	int nrows, ncols, nrowsA, ncolsA, length, incx, incy, idx, start[2], end[2];
	double *source, *destin;
	ops_gcg->Printf("PtAP sizeP = %d\n", sizeP);
	if (sizeP>0) {
		/* 计算 PtAP 部分 */
		nrows  = sizeP      ; ncols  = sizeP      ;
		nrowsA = sizeV-sizeC; ncolsA = sizeV-sizeC;
		/* C = alpha*op(Q)*op(A)*op(P) + beta*C */
		/* dbl_ws: nrows*ncols+nrowA*ncols
		 *       <=(sizeV+sizeP)*sizeP */
		ops_gcg->DenseMatQtAP('L','S',nrowsA,ncolsA,nrows,ncols,
				1.0,ss_evec+(sizeV-sizeC)*(sizeX-sizeC),sizeV-sizeC, /* Q */
				    ss_matA                            ,sizeV-sizeC, /* A */
				    ss_evec+(sizeV-sizeC)*(sizeX-sizeC),sizeV-sizeC, /* P */
				0.0,dbl_ws                             ,nrows      , /* C */
				dbl_ws+nrows*ncols);		
	}	
	
	sizeV  = sizeX+sizeP+sizeW;
	startN = startN + (nevConv-sizeC);
	endN   = endN   + (nevConv-sizeC);
	endN   = (endN<endX)?endN:endX;
	/* test */
	//endN   = (endN<gcg_solver->nevMax+20)?endN:(gcg_solver->nevMax+20);

	sizeN  = endN - startN; 
	sizeC  = nevConv;

	/* 更新 ss_mat ss_evec */
	ss_matA = ss_diag+(sizeV-sizeC);
	ss_evec = ss_matA+(sizeV-sizeC)*(sizeV-sizeC); 
	
	ops_gcg->Printf("WtAW sizeW = %d\n", sizeW);
#if TIME_GCG
#if USE_MPI
    time_gcg.rr_matW_time -= MPI_Wtime();
#else
    time_gcg.rr_matW_time -= clock();
#endif
#endif
	if (sizeW>0) {
		/* 计算 VtAW 部分 */
		start[0] = startN; end[0] = endW;
		start[1] = startW; end[1] = endW;
		destin = ss_matA+(sizeV-sizeC)*(sizeX+sizeP-sizeC);
		/* (endW-startN)*(endW-startW) 个 double 
		 *               (endW-startW) 个 向量 */
		ops_gcg->MultiVecQtAP('S','N',V,A,V,0,start,end,destin,sizeV-sizeC,
				mv_ws[0],ops_gcg);
		/* 对称化 */
		length = sizeX+sizeP-sizeC;
		source = ss_matA+(sizeV-sizeC)*(sizeX+sizeP-sizeC); incx = 1; 
		destin = ss_matA+(sizeX+sizeP-sizeC); incy = sizeV-sizeC;
		for (idx = 0; idx < sizeW; ++idx) {
			dcopy(&length,source,&incx,destin,&incy);
			source += sizeV-sizeC; destin += 1;
		}
	}
#if TIME_GCG
#if USE_MPI
    time_gcg.rr_matW_time += MPI_Wtime();
#else
    time_gcg.rr_matW_time += clock();
#endif
#endif
	
	if (sizeX == sizeV) {
#if DEBUG
		ops_gcg->Printf("V\n");
		ops_gcg->MultiVecView(V,0,sizeX,ops_gcg);
#endif		
		int block_size = gcg_solver->block_size;
		destin     = ss_matA;
		length     = sizeX-sizeC;
		block_size = block_size<length?block_size:length;
		start[0] = sizeC; end[0] = sizeX;
		start[1] = sizeC; end[1] = start[1]+block_size;
		while (length) {			
			ops_gcg->MultiVecQtAP('S','N',V,A,V,0,start,end,
				destin,sizeV-sizeC,mv_ws[0],ops_gcg);
			destin    += (sizeV-sizeC)*block_size;
			length    -= block_size;
			block_size = block_size<length?block_size:length;
			start[1] = end[1]; end[1] = start[1]+block_size;
		}	
	}
	else {
		/* 置零 X P 部分, 忽略 C 部分 */
		length = sizeX+sizeP-sizeC;
		destin = ss_matA;
		for (idx = 0; idx < length; ++idx) {
			memset(destin,0,length*sizeof(double));
			destin += sizeV-sizeC;
		}
		/* 赋值 X 部分的对角线 */
		length = sizeX-sizeC;
		source = ss_eval+sizeC; incx = 1              ; 
		destin = ss_matA      ; incy = (sizeV-sizeC)+1;
		dcopy(&length,source,&incx,destin,&incy);
		/* 更新 PtAP 部分*/
		length = sizeP;
		source = dbl_ws                                           ; incx = 1; 
		destin = ss_matA+(sizeV-sizeC)*(sizeX-sizeC)+(sizeX-sizeC); incy = 1;
		for (idx = 0; idx < length; ++idx) {
			dcopy(&length,source,&incx,destin,&incy);
			source += length; destin += sizeV-sizeC;
		}
	}
	
	/* 记录对角线部分 */
	length = sizeV-sizeC;
	source = ss_matA; incx = (sizeV-sizeC)+1; 
	destin = ss_diag; incy = 1              ;
	dcopy(&length,source,&incx,destin,&incy);
	
#if DEBUG	
	int row, col;
	ops_gcg->Printf("ss_diag:\n");
	for (idx = 0; idx < length; ++idx) ops_gcg->Printf("%f\n",destin[idx]);
#endif	
	/* 计算小规模特征值问题 */
	char   JOBZ, RANGE, UPLO; 
	int    LDA, M, LDZ, INFO, N, LWORK, *IWORK, *IFAIL; 
	double ABSTOL, *AA, *W, *Z, *WORK;
	JOBZ   = 'V'        ; RANGE  = 'A'; UPLO  = 'U'        ;
	LDA    = sizeV-sizeC; ABSTOL = tol; LDZ   = sizeV-sizeC; 
	IWORK  = int_ws; INFO   = 0  ;
	/* 不再计算 C 部分 */
	N      = sizeV-sizeC; M = N;
	IFAIL  = int_ws+5*N; 
	AA     = ss_matA;
	W      = ss_eval+sizeC; 
	Z      = ss_evec; 
	WORK   = Z+LDZ*N;
	/* ss_diag ss_matA ss_evec 剩下的空间 */
	LWORK  = sizeC+2*(2*sizeV*sizeC-sizeC*sizeC)+9*sizeV+sizeX*sizeP; 

#if DEBUG
	ops_gcg->Printf("dsyevx: AA\n"); 		
	for (row = 0; row < N; ++row) {
		for (col = 0; col < N; ++col) {
			ops_gcg->Printf("%6.4e\t",AA[row+col*LDA]);
		}
		ops_gcg->Printf("\n");
	}	
#endif


#if USE_MPI
	/* 当 PAS 调用 GCG 时, 且使用并行怎么办? 
	 * 没关系, PAS 需要保证每个进程都有特征向量 
	 * 同时, 这样的分批计算, 不仅仅是效率的提升
	 * 更重要的是, 保证, 每个进程的特征向量完全一致 */
	int *displs;
	int sendcount, *recvcounts;
	double *recvbuf;
	int IL, IU; int rank, nproc;

	/* 每列多一行, 将特征值拷贝至此, 进行通讯 */
	LDZ  = LDZ+1;
	/* 特征向量不包含 C 的部分 */
	Z    = ss_evec;	
	/* 重置工作空间 */ 
	WORK = Z+LDZ*N; LWORK = LWORK-N;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	/* 分类特征值 */
	destin = ss_eval+sizeC;
	length = N;
	/* 每组至少10个 */
	if (gcg_solver->compRR_min_num <= 0) {
	   gcg_solver->compRR_min_num = N/(nproc+2)>10?N/(nproc+2):10;
	}
	displs = malloc((2*nproc+1)*sizeof(int)); /* 长度需要 2*nproc+1 */
	if (rank == 0) {
	   SplitDoubleArray(destin,length,nproc,
		 gcg_solver->compRR_min_gap,
		 gcg_solver->compRR_min_num,
		 displs,dbl_ws,int_ws);
	}
	MPI_Bcast(displs,nproc+1,MPI_INT,0,MPI_COMM_WORLD);
 	sendcount  = displs[rank+1]-displs[rank];
	recvcounts = displs+nproc+1;
	for (idx = 0; idx < nproc; ++idx) {
		recvcounts[idx] = displs[idx+1]-displs[idx];
	}
	RANGE = 'I';
	/* 1 <= IL <= IU <= N */
	IL = displs[rank]+1; IU = displs[rank+1]  ;
	M  = IU-IL+1;
	/* 不同进程 W Z 不同 */
	W += displs[rank]  ; Z += LDZ*displs[rank];	

#if TIME_GCG
#if USE_MPI
    	time_gcg.dsyevx_time -= MPI_Wtime();
#else
    	time_gcg.dsyevx_time -= clock();
#endif
#endif
	if (sendcount > 0) {
		ops_gcg->Printf("dsyevx: N   = %d, M  = %d, LDA = %d, IL = %d, IU  = %d, LDZ = %d\n", 
		      N, M, LDA, IL, IU, LDZ);
		dsyevx(&JOBZ,&RANGE,&UPLO,&N,AA,&LDA,
				NULL,NULL,&IL,&IU,&ABSTOL,&M,
				W,Z,&LDZ,WORK,&LWORK,IWORK,IFAIL,&INFO);
		assert(M==IU-IL+1);
	}
#if TIME_GCG
#if USE_MPI
    	time_gcg.dsyevx_time += MPI_Wtime();
#else
    	time_gcg.dsyevx_time += clock();
#endif
#endif
	/* 将计算得到的特征值复制到 Z 的最后一行 */
	length  = sendcount;
	source  = W      ; incx    = 1  ;
	destin  = Z+LDZ-1; incy    = LDZ;
	dcopy(&length,source,&incx,destin,&incy);
	recvbuf = ss_evec;
	sendcount *= LDZ;
	for (idx = 0; idx < nproc; ++idx) {
		recvcounts[idx] *= LDZ;
		displs[idx+1]   *= LDZ;
	}
	/* 全聚集特征对, 发送和接收都是连续数据 */

#if DEBUG
	ops_gcg->Printf("before allgaterv sendcount = %d\n", sendcount);
#endif
	MPI_Allgatherv(MPI_IN_PLACE,sendcount,MPI_DOUBLE,
		recvbuf,recvcounts,displs,MPI_DOUBLE,MPI_COMM_WORLD);
#if DEBUG
	ops_gcg->Printf("after  allgaterv sendcount = %d\n", sendcount);
#endif
	free(displs);
	/* 将 Z 的最后一行复制给特征值 */
	length = N;
	source = ss_evec+LDZ-1; incx = LDZ;
	destin = ss_eval+sizeC; incy = 1  ;
	dcopy(&length,source,&incx,destin,&incy);
	/* 移动特征向量 */
#if DEBUG
	ops_gcg->Printf("before memmove length = %d\n", length);
#endif
	length = N; destin = ss_evec; source = ss_evec; 
	for (idx = 0; idx < N; ++idx) {
		/* 保证 source 在被覆盖之前
		 * 将重叠区域的字节拷贝到 destin 中 */
		memmove(destin,source,length*sizeof(double));
		destin += N; source += LDZ;
	}
#if DEBUG
	ops_gcg->Printf("after  memmove length = %d\n", length);
#endif

#else

	ops_gcg->Printf("dsyevx: N = %d, M = %d\n", N, M);
#if TIME_GCG
#if USE_MPI
    time_gcg.dsyevx_time -= MPI_Wtime();
#else
    time_gcg.dsyevx_time -= clock();
#endif
#endif
	/* 保证 ss_evec 是正交归一的 */
	dsyevx(&JOBZ,&RANGE,&UPLO,&N,AA,&LDA,
			NULL,NULL,NULL,NULL,&ABSTOL,&M,
			W,Z,&LDZ,WORK,&LWORK,IWORK,IFAIL,&INFO);
#if TIME_GCG
#if USE_MPI
    time_gcg.dsyevx_time += MPI_Wtime();
#else
    time_gcg.dsyevx_time += clock();
#endif
#endif
	ops_gcg->Printf("dsyevx: N = %d, M = %d\n", N, M);
	assert(M==N);

#endif

	/* 恢复ss_matA对角线部分 */
	length = sizeV-sizeC;
	source = ss_diag; incx = 1              ;
	destin = ss_matA; incy = (sizeV-sizeC)+1;
	dcopy(&length,source,&incx,destin,&incy);
	
#if DEBUG
	ops_gcg->Printf("dsyevx: ss_evec\n");
	for (row = 0; row < N; ++row) {
		for (col = 0; col < N; ++col) {
			ops_gcg->Printf("%6.4e\t",Z[row+col*LDZ]);
		}	
		ops_gcg->Printf("\n");
	}
	ops_gcg->Printf("dsyevx: ss_eval\n");
	for (row = 0; row < M; ++row) ops_gcg->Printf("%6.4e\n",W[row]);
	ops_gcg->Printf("dsyevx: AA\n"); 		
	for (row = 0; row < N; ++row) {
		for (col = 0; col < N; ++col) {
			ops_gcg->Printf("%6.4e\t",AA[row+col*LDA]);
		}
		ops_gcg->Printf("\n");
	}
#endif
#if TIME_GCG
#if USE_MPI
    time_gcg.compRR_time += MPI_Wtime();
#else
    time_gcg.compRR_time += clock();
#endif
#endif
	return;
}
static void GCG(void *A, void *B , double *eval, void **evec,
		int nevGiven, int *nevConv, struct OPS_ *ops)
{	
	/* offsetW[0] 表示有多少个块, 
	 * offsetW[1] <= idx < offsetW[2] 是未收敛的编号 */ 
	int *offsetP, *offsetW, *ptr_tmp;
	gcg_solver = (GCGSolver*)ops->eigen_solver_workspace;
	gcg_solver->nevGiven = nevGiven;	
	ops_gcg = ops;
	int    nevMax, multiMax, block_size, nevInit, nev0, nev;
	int    numIterMax, numIter, numCheck;
	void   **V, **ritz_vec;
	double *ss_matA, *ss_diag, *ss_eval, *ss_evec, *tol;
	int    start[2], end[2], idx; double *coef;

	nevMax     = gcg_solver->nevMax    ; multiMax = gcg_solver->multiMax; 
	block_size = gcg_solver->block_size; 
	/*  gcg_solver->nevInit 表示初始的 sizeX, 可以理解为 最大的 sizeX 
	 *  所谓初始的意思是, 由 Setup 给出的 nevMax+multiMax  
	 *  若 gcg_solver->nevInit <= 当前给出的 nevMax+multiMax, 否则工作空间不够用 */ 
	assert(gcg_solver->nevInit >= nevMax+multiMax); 
	/* 当前 nevInit 其实就等于当前 nevMax+multiMax
	 * 但是 nevInit 可以给出一个值, 使得
	 * block_size <= nevInit <= nevMax+multiMax 
	 * 此时, 初始给出的 sizeX 比最终要计算的 sizeX = nevMax+multiMax 要小
	 * 这样的好处是, dsyevx_ 的规模较小, 但 gcg 整体迭代次数变大, 
	 * 当前测试的例子, 表明其效果不好, 也许当特征值个数真的非常大时会由效果 */
	nevInit    = gcg_solver->nevInit>block_size?gcg_solver->nevInit:block_size;
	nevInit    = gcg_solver->nevInit<(nevMax+multiMax)?nevInit:(nevMax+multiMax);
	numIterMax = gcg_solver->numIterMax; tol = gcg_solver->tol;
	/* 全局变量初始化 */
	sizeC  = 0    ; sizeN = block_size  ; 
	/* sizeX 需要大于 nevGiven */
	sizeX  = nevInit>nevGiven?nevInit:nevGiven;
	sizeP  = 0    ; sizeW = 0           ; sizeV = sizeX+sizeP+sizeW;
	startN = sizeC; endN  = startN+sizeN; endX  = sizeX;
	startP = endX ; endP  = startP+sizeP;
	startW = endP ; endW  = startW+sizeW;
	/* workspace */
	V        = gcg_solver->mv_ws[0]; ritz_vec = evec;
	mv_ws[0] = gcg_solver->mv_ws[1]; mv_ws[1] = gcg_solver->mv_ws[2];
	mv_ws[2] = gcg_solver->mv_ws[3];
	ss_eval  = gcg_solver->dbl_ws; 
	for (idx = 0; idx < (nevMax+multiMax+2*block_size); ++idx) {
	   ss_eval[idx] = 1.0;
	}
	ss_diag  = ss_eval+(nevMax+multiMax)+2*block_size;
	ss_matA  = ss_diag+(sizeV-sizeC);
	ss_evec  = ss_matA+(sizeV-sizeC)*(sizeV-sizeC); 
	int distance = 2*(nevMax+multiMax+2*block_size) /*ss_eval ss_diag */ 
			+(nevMax+multiMax+2*block_size)*(nevMax+multiMax+2*block_size)  /* ss_matA */
			+(nevMax+multiMax+2*block_size)*(nevMax+multiMax+1*block_size); /* ss_evec */ 
	dbl_ws   = gcg_solver->dbl_ws+distance;
	
	offsetP  = gcg_solver->int_ws;
	offsetW  = offsetP + block_size+2; 
	int_ws   = offsetW + block_size+2;

#if TIME_GCG
	time_gcg.checkconv_time = 0.0;
	time_gcg.compP_time     = 0.0; 
	time_gcg.compRR_time    = 0.0; 
	time_gcg.compRV_time    = 0.0;
	time_gcg.compW_time     = 0.0;
	time_gcg.compX_time     = 0.0;
	time_gcg.rr_matW_time   = 0.0;
	time_gcg.dsyevx_time    = 0.0;
	time_gcg.initX_time     = 0.0;
	time_gcg.linsol_time    = 0.0;
#endif
	
	ops_gcg->Printf("initial X\n");
	/* 对 X 赋随机初值且 B 正交归一化 */
	InitializeX(V,ritz_vec,B,nevGiven);	

#if DEBUG
	int row, col;
#endif

	ops_gcg->Printf("ComputeRayleighRitz\n");	
	ComputeRayleighRitz(ss_matA,ss_eval,ss_evec,
		gcg_solver->compRR_tol,0,ss_diag,A,V);	


	for (idx = sizeV; idx < (nevMax+multiMax+2*block_size); ++idx) {
	   ss_eval[idx] = ss_eval[sizeV-1];
	}
	/* 更新 ss_mat ss_evec */
	ss_matA = ss_diag+(sizeV-sizeC);
	ss_evec = ss_matA+(sizeV-sizeC)*(sizeV-sizeC);

	ops_gcg->Printf("ComputeRitzVec\n");	
	ComputeRitzVec(ritz_vec,V,ss_evec);				
	
	*nevConv = (*nevConv)<nevMax?(*nevConv):nevMax;
	nev0 = *nevConv; /* 用户希望收敛的特征对个数 */
	/* 收敛个数达到 nev 后将 P 和 W 部分扩充为 X 部分 */
	nev  = sizeX-multiMax; *nevConv = 0; numIter = 0;
	do {
		ops_gcg->Printf("------------------------------\n");
		ops_gcg->Printf("numIter = %d, sizeC = %d, sizeN = %d, sizeX = %d, sizeP = %d, sizeW = %d, sizeV = %d\n",
				numIter,sizeC,sizeN,sizeX,sizeP,sizeW,sizeV);
		ops_gcg->Printf("CheckConvergence\n");
		if (numIter <= 0) {
		//if (numIter == 0) {
		   //numCheck = nevGiven;
		   numCheck = 0;
		   //numCheck = 20;
		   //numCheck = (nevGiven-10)>(nevGiven/2)?(nevGiven-10):(nevGiven/2);
		} 
		else {
		   //numCheck = (startN+multiMax+sizeN<endX)?(multiMax+sizeN):(endX-startN);
		   numCheck = (startN+sizeN<endX)?(sizeN):(endX-startN);
		   //numCheck = 30; 
		   //numCheck = (numCheck-10)>(numCheck/2)?(numCheck-10):(numCheck/2);
		}
		numCheck = numCheck<gcg_solver->check_conv_max_num?numCheck:gcg_solver->check_conv_max_num;
		*nevConv = CheckConvergence(A,B,ss_eval,ritz_vec,numCheck,tol,offsetW);		
		ops_gcg->Printf("%d\n",*nevConv);		
		if (*nevConv >= nev) {
			if (*nevConv >= nev0) {
				break;
			}
			else {
				nev   += sizeP+sizeW; 
				nev    = nev<nev0?nev:nev0;
				sizeX += sizeP+sizeW; 
				sizeX  = sizeX<(nev0+multiMax)?sizeX:(nev0+multiMax);
				/* 将 P 和 W 部分写入 ritz_vec */
				start[0] = startN; end[0] = endW ;
				start[1] = endX  ; end[1] = sizeX; 
				coef     = ss_evec+(sizeV-sizeC)*(endX-sizeC);
				ops_gcg->MultiVecLinearComb(V,ritz_vec,0, 
						start,end,coef,sizeV-sizeC,NULL,0,ops_gcg);
			
				sizeP  = 0; sizeW = 0; sizeV = sizeX;
				startP = endX ; endP = startP;
				startW = endP ; endW = startW; 
				endX   = sizeX; 
				numIterMax -= numIter; numIter = 0;				
			}
		}
		if (numIter == 0)	{
			sizeP = 0; startP = endX; endP = startP+sizeP;
		}
		else {
			ops_gcg->Printf("ComputeP\n");
			ComputeP(V,ss_evec,offsetP); /* update sizeP startP endP */
		}

		ops_gcg->Printf("ComputeX\n");
		ComputeX(V,ritz_vec);

#if DEBUG		
		ops->MultiVecView(V,0,sizeX,ops);
#endif

		ops_gcg->Printf("ComputeW\n");
		//offsetW = 0;
		ComputeW(V,A,B,ss_eval,ritz_vec,offsetW); /* update sizeW startW endW */
		ptr_tmp = offsetP; offsetP = offsetW; offsetW = ptr_tmp;
		
		ops_gcg->Printf("ComputeRayleighRitz\n");
#if DEBUG	
		ops_gcg->Printf("VtAV\n");
		start[0] = 0; end[0] = sizeX+sizeP+sizeW; start[1] = 0; end[1] = sizeX+sizeP+sizeW;
		ops->MultiVecQtAP('N','N',V,A,V,0,start,end,dbl_ws,sizeX+sizeP+sizeW,mv_ws[0],ops);
		for (row = 0; row < end[0]-start[0]; ++row) {
			for (col = 0; col < end[1]-start[1]; ++col) {
				ops_gcg->Printf("%6.4e\t",dbl_ws[row+col*(end[0]-start[0])]);
			}	
			ops_gcg->Printf("\n");
		}			
		ops_gcg->Printf("VtBV\n");
		ops->MultiVecQtAP('N','N',V,B,V,0,start,end,dbl_ws,sizeX+sizeP+sizeW,mv_ws[0],ops);
		for (row = 0; row < end[0]-start[0]; ++row) {
			for (col = 0; col < end[1]-start[1]; ++col) {
				ops_gcg->Printf("%6.4e\t",dbl_ws[row+col*(end[0]-start[0])]);
			}	
			ops_gcg->Printf("\n");
		}	
#endif
	
		/* 计算完 PtAP 部分后再更新 sizeV */
		ComputeRayleighRitz(ss_matA,ss_eval,ss_evec,
			gcg_solver->compRR_tol,*nevConv,ss_diag,A,V);		

		for (idx = sizeV; idx < (nevMax+multiMax+2*block_size); ++idx) {
		   ss_eval[idx] = ss_eval[sizeV-1];
		}
		ss_matA = ss_diag+(sizeV-sizeC);
		ss_evec = ss_matA+(sizeV-sizeC)*(sizeV-sizeC);
		
		ops_gcg->Printf("ComputeRitzVec\n");
		ComputeRitzVec(ritz_vec,V,ss_evec);
		
		++numIter;
	} while (numIter < numIterMax);
	
	gcg_solver->numIter = numIter+(gcg_solver->numIterMax-numIterMax);	
	/* eval evec 都是 sizeX 长 */
	int inc = 1;
	dcopy(&sizeX,ss_eval,&inc,eval,&inc);
	
#if TIME_GCG
	ops_gcg->Printf("|--GCG----------------------------\n");
	time_gcg.time_total = time_gcg.checkconv_time
		+time_gcg.compP_time
		+time_gcg.compRR_time
		+time_gcg.compRV_time
		+time_gcg.compW_time
		+time_gcg.compX_time
		+time_gcg.initX_time;
	ops_gcg->Printf("|checkconv  compP  compRR (rr_matW  dsyexv)  compRV  compW (linsol)  compX  initX\n");
#if USE_MPI	
	ops_gcg->Printf("|%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
		time_gcg.checkconv_time,		
		time_gcg.compP_time,		
		time_gcg.compRR_time,
		time_gcg.rr_matW_time,
		time_gcg.dsyevx_time,		
		time_gcg.compRV_time,
		time_gcg.compW_time,
		time_gcg.linsol_time,
		time_gcg.compX_time,
		time_gcg.initX_time);	   	
#else	
	ops_gcg->Printf("|%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
		time_gcg.checkconv_time/CLOCKS_PER_SEC,		
		time_gcg.compP_time    /CLOCKS_PER_SEC,
		time_gcg.compRR_time   /CLOCKS_PER_SEC,
		time_gcg.rr_matW_time  /CLOCKS_PER_SEC,
		time_gcg.dsyevx_time   /CLOCKS_PER_SEC,		
		time_gcg.compRV_time   /CLOCKS_PER_SEC,
		time_gcg.compW_time    /CLOCKS_PER_SEC,
		time_gcg.linsol_time   /CLOCKS_PER_SEC,
		time_gcg.compX_time    /CLOCKS_PER_SEC,
		time_gcg.initX_time    /CLOCKS_PER_SEC);	
#endif
	ops_gcg->Printf("|%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\n",
		time_gcg.checkconv_time/time_gcg.time_total*100,
		time_gcg.compP_time    /time_gcg.time_total*100,
		time_gcg.compRR_time   /time_gcg.time_total*100,
		time_gcg.rr_matW_time  /time_gcg.compRR_time*100,
		time_gcg.dsyevx_time   /time_gcg.compRR_time*100,
		time_gcg.compRV_time   /time_gcg.time_total*100,		
		time_gcg.compW_time    /time_gcg.time_total*100,		
		time_gcg.linsol_time   /time_gcg.compW_time*100,
		time_gcg.compX_time    /time_gcg.time_total*100,
		time_gcg.initX_time    /time_gcg.time_total*100);
	ops_gcg->Printf("|--GCG----------------------------\n");
	time_gcg.checkconv_time = 0.0;
	time_gcg.compP_time     = 0.0; 
	time_gcg.compRR_time    = 0.0; 
	time_gcg.compRV_time    = 0.0;
	time_gcg.compW_time     = 0.0;
	time_gcg.compX_time     = 0.0;
	time_gcg.rr_matW_time   = 0.0;
	time_gcg.dsyevx_time    = 0.0;
	time_gcg.initX_time     = 0.0;
	time_gcg.linsol_time    = 0.0;
#endif
	
	return;
}

/* 设定 GCG 的工作空间 */
void EigenSolverSetup_GCG(
	int    nevMax  , int    multiMax  , double gapMin, int block_size, 
	double tol[2]  , int    numIterMax, 
	int user_defined_multi_linear_solver,
	void **mv_ws[4], double *dbl_ws   , int *int_ws,	
	struct OPS_ *ops)
{
	static GCGSolver gcg_solver_static = {
		.nevMax     = 1, .multiMax   = 2   , .gapMin = 0.01, 
		.nevInit    = 3, .nevGiven   = 0   ,
		.block_size = 1, .numIterMax = 4   , .user_defined_multi_linear_solver = 0,
		.mv_ws      = NULL, .dbl_ws  = NULL, .int_ws     = NULL,		
		/* 算法内部参数 */		
		.initX_orth_method     = "mgs",
		.initX_orth_block_size = -1   ,
		.initX_orth_max_reorth = 1    ,
		.initX_orth_zero_tol   = 1e-14,
		.check_conv_max_num    = 15   ,	
		.compP_orth_method     = "mgs", 
		.compP_orth_block_size = -1   ,
		.compP_orth_max_reorth = 1    ,
		.compP_orth_zero_tol   = 1e-14,		
		.compW_orth_method     = "mgs",
		.compW_orth_block_size = -1   ,
		.compW_orth_max_reorth = 1    ,
		.compW_orth_zero_tol   = 1e-14,	
		.compW_cg_max_iter     = 40   ,
		.compW_cg_rate         = 1e-2 , 
		.compW_cg_tol          = 1e-8 ,
		.compW_cg_tol_type     = "abs",		
		.compRR_min_gap        = 0.01 ,
		.compRR_min_num        = -1   ,
		.compRR_tol            = 1e-16,
	};
		
	gcg_solver_static.nevMax     = nevMax;
	gcg_solver_static.multiMax   = multiMax;
	gcg_solver_static.gapMin     = gapMin;	 
	gcg_solver_static.block_size = block_size;
	gcg_solver_static.nevInit    = nevMax+multiMax;
	gcg_solver_static.tol[0]     = tol[0];
	gcg_solver_static.tol[1]     = tol[1];
	gcg_solver_static.numIterMax = numIterMax;
	gcg_solver_static.mv_ws[0]   = mv_ws[0];
	gcg_solver_static.mv_ws[1]   = mv_ws[1];
	gcg_solver_static.mv_ws[2]   = mv_ws[2];
	gcg_solver_static.mv_ws[3]   = mv_ws[3];
	gcg_solver_static.dbl_ws     = dbl_ws;
 	gcg_solver_static.int_ws     = int_ws;
 	
 	gcg_solver_static.compRR_min_gap = gapMin;
 	gcg_solver_static.check_conv_max_num = block_size;
 	gcg_solver_static.user_defined_multi_linear_solver = user_defined_multi_linear_solver;
		
	ops->eigen_solver_workspace = (void *)(&gcg_solver_static);
	ops->EigenSolver            = GCG;
	return;	
}

void EigenSolverCreateWorkspace_GCG(
	int    sizeX, int block_size, void *mat,
	void ***mv_ws, double **dbl_ws, int **int_ws, 
	struct OPS_ *ops)
{
	assert(mv_ws!=NULL);
	int sizeV = sizeX+2*block_size; 
	ops->MultiVecCreateByMat(&mv_ws[0],sizeV,mat,ops);				
	ops->MultiVecSetRandomValue(mv_ws[0],0,sizeV,ops);
	ops->MultiVecCreateByMat(&mv_ws[1],block_size,mat,ops);				
	ops->MultiVecSetRandomValue(mv_ws[1],0,block_size,ops);
	ops->MultiVecCreateByMat(&mv_ws[2],block_size,mat,ops);				
	ops->MultiVecSetRandomValue(mv_ws[2],0,block_size,ops);
	ops->MultiVecCreateByMat(&mv_ws[3],block_size,mat,ops);				
	ops->MultiVecSetRandomValue(mv_ws[3],0,block_size,ops);
	
	int length_dbl_ws = 2*sizeV*sizeV+11*sizeV+sizeX*block_size;
	int length_int_ws = 6*sizeV+2*(block_size+2);
	if (dbl_ws!=NULL) {
		*dbl_ws = malloc(length_dbl_ws*sizeof(double));
		memset(*dbl_ws,0,length_dbl_ws*sizeof(double));
	} 
	if (int_ws!=NULL) {
		*int_ws = malloc(length_int_ws*sizeof(double));	
		memset(*int_ws,0,length_int_ws*sizeof(int));	
	}	
	return;
}
void EigenSolverDestroyWorkspace_GCG(
	int    sizeX, int block_size, void *mat,
	void ***mv_ws, double **dbl_ws, int **int_ws, 
	struct OPS_ *ops)
{
	assert(mv_ws!=NULL);
	int sizeV = sizeX+2*block_size; 
	ops->MultiVecDestroy(&mv_ws[0],sizeV,ops);
	ops->MultiVecDestroy(&mv_ws[1],block_size,ops);
	ops->MultiVecDestroy(&mv_ws[2],block_size,ops);
	ops->MultiVecDestroy(&mv_ws[3],block_size,ops);
	if (dbl_ws!=NULL) {
		free(*dbl_ws); *dbl_ws = NULL;
	}
	if (int_ws!=NULL) {
		free(*int_ws); *int_ws = NULL;
	}
	return;
}


/* 参数设定函数需要在 Setup 之后调用 */
void EigenSolverSetParameters_GCG(
	int check_conv_max_num,
	const char *initX_orth_method, int initX_orth_block_size, int initX_orth_max_reorth, double initX_orth_zero_tol,
	const char *compP_orth_method, int compP_orth_block_size, int compP_orth_max_reorth, double compP_orth_zero_tol,
	const char *compW_orth_method, int compW_orth_block_size, int compW_orth_max_reorth, double compW_orth_zero_tol,
	int compW_cg_max_iter, double compW_cg_rate, double compW_cg_tol, const char *compW_cg_tol_type,
	int compRR_min_num, double compRR_min_gap, double compRR_tol, 
	struct OPS_ *ops)
{
	
	struct GCGSolver_ *gcg_solver = (GCGSolver*)ops->eigen_solver_workspace;
	
	gcg_solver->check_conv_max_num = check_conv_max_num;
	gcg_solver->initX_orth_method  = initX_orth_method;
	if (initX_orth_block_size>0)
		gcg_solver->initX_orth_block_size = initX_orth_block_size;
	if (initX_orth_max_reorth>=0)
		gcg_solver->initX_orth_max_reorth = initX_orth_max_reorth;
	if (initX_orth_zero_tol>0)
		gcg_solver->initX_orth_zero_tol   = initX_orth_zero_tol;
	
	gcg_solver->compP_orth_method = compP_orth_method;
	if (compP_orth_block_size>0)
		gcg_solver->compP_orth_block_size = compP_orth_block_size;
	if (compP_orth_max_reorth>=0)
		gcg_solver->compP_orth_max_reorth = compP_orth_max_reorth;
	if (compP_orth_zero_tol>0)
		gcg_solver->compP_orth_zero_tol   = compP_orth_zero_tol;	
	
	gcg_solver->compW_orth_method = compW_orth_method;
	if (compW_orth_block_size>0)
		gcg_solver->compW_orth_block_size = compW_orth_block_size;
	if (compW_orth_max_reorth>=0)
		gcg_solver->compW_orth_max_reorth = compW_orth_max_reorth;
	if (compW_orth_zero_tol>0)
		gcg_solver->compW_orth_zero_tol   = compW_orth_zero_tol;
	if (compW_cg_max_iter>0)	
		gcg_solver->compW_cg_max_iter = compW_cg_max_iter;
	if (compW_cg_rate>0)
		gcg_solver->compW_cg_rate     = compW_cg_rate;
	if (compW_cg_tol>0)
		gcg_solver->compW_cg_tol      = compW_cg_tol;	
	gcg_solver->compW_cg_tol_type = compW_cg_tol_type;
		
	if (compRR_min_gap>0)
		gcg_solver->compRR_min_gap = compRR_min_gap;
	if (compRR_min_num>0)
		gcg_solver->compRR_min_num = compRR_min_num;
	if (compRR_tol>0)
		gcg_solver->compRR_tol     = compRR_tol;
	
	return;	
}

void EigenSolverSetParametersFromCommandLine_GCG(
	int argc, char* argv[], struct OPS_ *ops)
{
	struct GCGSolver_ *gcg_solver = (GCGSolver*)ops->eigen_solver_workspace;

	ops->GetOptionFromCommandLine("-gcge_nev"        ,'i',
		&gcg_solver->nevMax    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_max_multi"  ,'i',
		&gcg_solver->multiMax  ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_min_gap"    ,'f',
		&gcg_solver->gapMin    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_given_nevec",'i',
		&gcg_solver->nevGiven  ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_block_size" ,'i',
		&gcg_solver->block_size,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_max_niter"  ,'i',
		&gcg_solver->numIterMax,argc,argv, ops);	
	ops->GetOptionFromCommandLine("-gcge_abs_tol"    ,'f',
		&gcg_solver->tol[0]    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_rel_tol"    ,'f',
		&gcg_solver->tol[1]    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_user_defined_multi_lin_sol",'i',
		&gcg_solver->user_defined_multi_linear_solver,argc,argv, ops);
	
	ops->GetOptionFromCommandLine("-gcge_initX_orth_method"    ,'s',
		&gcg_solver->initX_orth_method    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_initX_orth_block_size",'i',
		&gcg_solver->initX_orth_block_size,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_initX_orth_max_reorth",'i',
		&gcg_solver->initX_orth_max_reorth,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_initX_orth_zero_tol"  ,'f',
		&gcg_solver->initX_orth_zero_tol  ,argc,argv, ops);
	
	ops->GetOptionFromCommandLine("-gcge_check_conv_max_num",'i',
		&gcg_solver->check_conv_max_num,argc,argv, ops);
	
	ops->GetOptionFromCommandLine("-gcge_compP_orth_method"    ,'s',
		&gcg_solver->compP_orth_method    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compP_orth_block_size",'i',
		&gcg_solver->compP_orth_block_size,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compP_orth_max_reorth",'i',
		&gcg_solver->compP_orth_max_reorth,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compP_orth_zero_tol"  ,'f',
		&gcg_solver->compP_orth_zero_tol  ,argc,argv, ops);
	
	ops->GetOptionFromCommandLine("-gcge_compW_orth_method"    ,'s',
		&gcg_solver->compW_orth_method    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_orth_block_size",'i',
		&gcg_solver->compW_orth_block_size,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_orth_max_reorth",'i',
		&gcg_solver->compW_orth_max_reorth,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_orth_zero_tol"  ,'f',
		&gcg_solver->compW_orth_zero_tol  ,argc,argv, ops);
	
	ops->GetOptionFromCommandLine("-gcge_compW_cg_max_iter",'i',
		&gcg_solver->compW_cg_max_iter,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_cg_rate"    ,'f',
		&gcg_solver->compW_cg_rate    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_cg_tol"     ,'f',
		&gcg_solver->compW_cg_tol     ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_cg_tol_type",'s',
		&gcg_solver->compW_cg_tol_type,argc,argv, ops);
	
	ops->GetOptionFromCommandLine("-gcge_compRR_min_num",'i',
		&gcg_solver->compRR_min_num,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compRR_min_gap",'i',
		&gcg_solver->compRR_min_gap,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compRR_tol    ",'f',
		&gcg_solver->compRR_tol    ,argc,argv, ops);

	int print_usage = 1;
	ops->GetOptionFromCommandLine("-gcge_print_usage",'i',&print_usage,argc,argv, ops);
	if (print_usage) {
       ops->Printf("\n");
       ops->Printf("Usage: %s [<options>]\n", argv[0]);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_nev          <i>: number of eigenpairs you need               %d (default: 6)\n",gcg_solver->nevMax);
       ops->Printf(" -gcge_max_multi    <i>: maximum of multiplicity of eigenpairs       %d (default: 6)\n",gcg_solver->multiMax);
       ops->Printf(" -gcge_min_gap      <f>: minimum of gap of eigenvalues relatively    %.2e (default: 1e-2)\n",gcg_solver->gapMin);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_max_niter    <i>: maximum of gcg iterations                   %d (default: 100)\n",gcg_solver->numIterMax);
       ops->Printf(" -gcge_block_size   <i>: number of eigenpairs computed in one patch  %d (default: nev/5)\n",gcg_solver->block_size);
       ops->Printf(" -gcge_given_nevec  <i>: number of given initial eigenvectors        %d (default: 0)\n",gcg_solver->nevGiven);
       ops->Printf(" -gcge_abs_tol      <f>: absolute convergence tolerance              %.2e (default: 1e-4)\n",gcg_solver->tol[0]);
       ops->Printf(" -gcge_rel_tol      <f>: relative convergence tolerance              %.2e (default: 1e-4)\n",gcg_solver->tol[1]);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_user_defined_multi_lin_sol  <i>: use user-defined multi linear solver  %d (default: 0[1])\n",gcg_solver->user_defined_multi_linear_solver);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_initX_orth_method  <s>: use which kind of orthogonalization for X  %s (default: mgs[bgs])\n",gcg_solver->initX_orth_method);
       ops->Printf(" -gcge_compP_orth_method  <s>: use which kind of orthogonalization for P  %s (default: bqr[bgs|mgs])\n",gcg_solver->compP_orth_method);
       ops->Printf(" -gcge_compW_orth_method  <s>: use which kind of orthogonalization for W  %s (default: mgs[bgs])\n",gcg_solver->compW_orth_method);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_initX_orth_block_size  <i>: size of vectors orthogonalized in one patch for X  %d (default: -1)\n",gcg_solver->initX_orth_block_size);
       ops->Printf(" -gcge_compP_orth_block_size  <i>: size of vectors orthogonalized in one patch for P  %d (default: -1)\n",gcg_solver->compP_orth_block_size);
       ops->Printf(" -gcge_compW_orth_block_size  <i>: size of vectors orthogonalized in one patch for W  %d (default: -1)\n",gcg_solver->compW_orth_block_size);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_initX_orth_zero_tol  <f>: zero tolerance in orthogonal for X  %.2e (default: 1e-16)\n",gcg_solver->initX_orth_zero_tol);
       ops->Printf(" -gcge_compP_orth_zero_tol  <f>: zero tolerance in orthogonal for P  %.2e (default: 1e-16)\n",gcg_solver->compP_orth_zero_tol);
       ops->Printf(" -gcge_compW_orth_zero_tol  <f>: zero tolerance in orthogonal for W  %.2e (default: 1e-16)\n",gcg_solver->compW_orth_zero_tol);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_initX_orth_max_reorth  <i>: maximum reorthogonal times for X  %d (default: 2)\n",gcg_solver->initX_orth_max_reorth);
       ops->Printf(" -gcge_compP_orth_max_reorth  <i>: maximum reorthogonal times for P  %d (default: 2)\n",gcg_solver->compP_orth_max_reorth);
       ops->Printf(" -gcge_compW_orth_max_reorth  <i>: maximum reorthogonal times for W  %d (default: 2)\n",gcg_solver->compW_orth_max_reorth);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_compW_cg_max_iter  <i>: maximum number of cg iteration       %d (default: 30)\n",gcg_solver->compW_cg_max_iter);
       ops->Printf(" -gcge_compW_cg_rate      <f>: descent rate of residual in cg       %.2e (default: 1e-2)\n",gcg_solver->compW_cg_rate);
       ops->Printf(" -gcge_compW_cg_tol       <f>: convergence tolerance in cg          %.2e (default: 1e-8)\n",gcg_solver->compW_cg_tol);
       ops->Printf(" -gcge_compW_cg_tol_type  <s>: type of convergence tolerance in cg  %s (default: abs[rel|user])\n",gcg_solver->compW_cg_tol_type);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_compRR_min_num  <i>: minimum number for splitting RR eval  %d (default: 10)\n",gcg_solver->compRR_min_num);
       ops->Printf(" -gcge_compRR_min_gap  <f>: minimum gap for splitting RR eval     %.2e (default: 1e-2)\n",gcg_solver->compRR_min_gap);
       ops->Printf(" -gcge_compRR_tol      <f>: convergence tolerance in RR           %.2e (default: 1e-16)\n",gcg_solver->compRR_tol);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_print_orth_zero  <i>: print the zero index in orthogonal      %d (default: 0[1])\n",1);
       ops->Printf(" -gcge_print_split      <i>: print the split information of RR eval  %d (default: 0[1])\n",0);
       ops->Printf(" -gcge_print_conv       <i>: print convergence in each iteration     %d (default: 1[0])\n",1);
       ops->Printf(" -gcge_print_eval       <i>: print the final eigenvalues             %d (default: 1[0])\n",1);
       ops->Printf(" -gcge_print_evec       <i>: print the final eigenvectors            %d (default: 0[1])\n",0);
       ops->Printf(" -gcge_print_time       <i>: print total time of each part           %d (default: 1[0])\n",1);
       ops->Printf(" -gcge_print_usage      <i>: print usage of gcg eigen solver         %d (default: 1[0])\n",1);
       ops->Printf("--------------------------------------------------------------------------------------------------\n");
       //ops->Printf(" -bpcg_print_res        <i>: print residual per five bpcg iteration  (default: 1[0])\n");
    }
	return;	
}
