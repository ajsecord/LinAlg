/*
 * LinAlg 3.0: fixed-sized vector and matrix library for small dimensions with optional LAPACK bindings.
 * Lapack.h
 * Contact: http://www.google.com/search?q=%22Adrian+Secord%22
 * Copyright 2005-2008 Adrian Secord.
 */

#ifndef LINALG_LAPACK_H
#define LINALG_LAPACK_H

extern "C" {
    // -------------------- BLAS routines --------------------
    void daxpy_(const int *, const double *, const double *, const int *, 
                 double *, const int *);
    double dnrm2_(const int *, const double *, const int *);
    double ddot_(const int *, const double *, const int *, const double *, 
              const int *);  

    // -------------------- LAPACK routines --------------------

    int dsyev_(char *jobz, char *uplo, int *n, double *a,
            int *lda, double *w, double *work, int *lwork, 
            int *info);

    int ssyev_(char *jobz, char *uplo, int *n, float *a,
            int *lda, float *w, float *work, int *lwork, 
            int *info);

    int dsyevx_(char *jobz, char *range, char *uplo, int *n, 
             double *a, int *lda, double *vl, double *vu, 
             int * il, int *iu, double *abstol, int *m, double *w, 
             double *z__, int *ldz, double *work, int *lwork, 
             int *iwork, int *ifail, int *info);

    int dgesvd_(char * jobu, char * jobvt, int * m, int * n,
             double * A, int * lda, double * s, double * u,
             int * ldu, double * vt,int * ldvt, double * work,
             int * lwork, int * info);

    int sgesvd_(char * jobu, char * jobvt, int * m, int * n,
             float * A, int * lda, float * s, float * u,
             int * ldu, float * vt,int * ldvt, float * work,
             int * lwork, int * info);

    int dgelss_(int * m, int * n, int * nrhs, double * A, int * lda,
             double * B, int * ldb, double * s,double * rcond,int * rank,
             double * work, int * lwork, int *info);

    void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
    void sgetrf_(int *M, int *N, float *A, int *LDA, int *IPIV, int *INFO);

    int dgetrs_(char *trans, int *n, int *nrhs, 
                double *a, int *lda, int *ipiv, double *b, int * ldb, 
                int *info);

    void dgecon_(char *NORM, int *N, double *A, int *LDA, double *ANORM, 
                double *RCOND, double *WORK, int *IWORK, int *INFO);
    void fgecon_(char *NORM, int *N, float *A, int *LDA, float *ANORM, 
                float *RCOND, float *WORK, int *IWORK, int *INFO);

    void dgetri_(int* N, double* A, int* LDA, int* IPIV, double* WORK,
                 int* LWORK, int* INFO);
    void sgetri_(int* N, float* A, int* LDA, int* IPIV, float* WORK,
                 int* LWORK, int* INFO);

    int ilaenv_(int* IPSPEC, char* NAME, char* OPTS, int* N1, int* N2, 
                int* N3, int* N4);
}

#endif
