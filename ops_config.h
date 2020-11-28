/**
 *    @file  ops_config.h
 *   @brief  �����ļ�
 *
 *  �����ļ�
 *
 *  @author  Yu Li, liyu@tjufe.edu.cn
 *
 *       Created:  2020/9/13
 *      Revision:  none
 */
#ifndef  _OPS_CONFIG_H_
#define  _OPS_CONFIG_H_

#define  USE_MPI    0
#define  USE_SLEPC  0
#define  USE_MATLAB 0

#define  PRINT_RANK 0

/* MATLAB �� blas lapack ��ĺ��������� _ */
#if USE_MATLAB
#define FORTRAN_WRAPPER(x) x
#else
#define FORTRAN_WRAPPER(x) x ## _
#endif


#endif  /* -- #ifndef _OPS_CONFIG_H_ -- */
