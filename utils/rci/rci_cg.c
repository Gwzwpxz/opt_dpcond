/** @file rci\_cg
 *  @brief Mexfile entry function for parallel conjugate gradient method with reverse communicationb interface
 *
 * @author Wenzhi Gao, Shanghai University of Finance and Economics
 * @date Oct, 26th, 2022
 *
 */

#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"

static int rciCGSolve( int nCol, int nRhs, int *colMatBeg, int *colMatIdx, double *colMatElem,
                       double *rhsMat, double *lhsSol, double *dParam, int *iParam,
                       double *tmpArray, double cgTol, double maxIter, int *iterCount ) {
    
    double eone = -1.0;
    MKL_INT ione = 1, method = 1, rci_request;
    
    sparse_matrix_t cscA;
    struct matrix_descr descrA;
    sparse_operation_t transA = SPARSE_OPERATION_NON_TRANSPOSE;
    
    descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;
    
    mkl_sparse_d_create_csc ( &cscA, SPARSE_INDEX_BASE_ZERO, nCol, nCol,
                              colMatBeg, colMatBeg + 1, colMatIdx, colMatElem );
    
    memset(iParam, 0, sizeof(int) * (128 + 2 * nRhs));
    memset(dParam, 0, sizeof(double) * (128 + 2 * nRhs));
    
    dcgmrhs_init(&nCol, lhsSol, &nRhs, rhsMat, &method, &rci_request, iParam, dParam, tmpArray);
    
    if ( rci_request != 0 ) {
        mexErrMsgTxt("Failed to intialize RCI \n");
        goto failure;
    }
    
    iParam[4] = maxIter;
    iParam[8] = 1;
    iParam[9] = 0;
    dParam[0] = cgTol;
    
rci:dcgmrhs(&nCol, lhsSol, &nRhs, rhsMat, &rci_request, iParam, dParam, tmpArray);
    
    if ( rci_request == 0 )
        goto getsln;
    
    if (rci_request == 1) {
        mkl_sparse_d_mv( transA, 1.0, cscA, descrA, tmpArray, 0.0, tmpArray + nCol);
        goto rci;
    }
    goto failure;
    
getsln:dcgmrhs_get(&nCol, lhsSol, &nRhs, rhsMat, &rci_request, iParam, dParam, tmpArray, iterCount);
    mkl_sparse_destroy(cscA);
    MKL_Free_Buffers();
    return 0;
    
failure:
    mkl_sparse_destroy(cscA);
    MKL_Free_Buffers();
    return 1;
}

#define X       plhs[0]
#define A       prhs[0]
#define B       prhs[1]
#define tol     prhs[2]
#define maxIter prhs[3]
/** @brief Matlab entry function
 *  @param[in] nlhs Number of left-hand-side parameters
 *  @param[out] plhs Pointers for left-hand-side parameters
 *  @param[in] nrhs Number of right-hand-side parameters
 *  @param[out] prhs Pointers for left-hand-side parameters
 */
extern void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
    mwSize nCol = 0;
    mwSize nRhs = 0;
    
    mwSize *Ap = NULL;
    mwSize *Ai = NULL;
    
    int *colMatBeg = NULL;
    int *colMatIdx = NULL;
    double *colMatElem = NULL;
    
    double *colObj = NULL;
    double *rowRhs = NULL;
    double *lhsSol = NULL;
    double *rhsMat = NULL;
    
    if ( nrhs != 4 ) {
        mexErrMsgTxt("Invalid number of entries. \n");
    }
    
    /* Get A */
    if ( !mxIsSparse(A) ) {
        mexErrMsgTxt("A must be sparse. \n");
    }
    
    if ( mxGetM(A) != mxGetN(A) || mxGetM(B) != mxGetN(A) ) {
        mexErrMsgTxt("A must be square. B's dimension must agree with A \n");
    }
    
    nRhs = mxGetN(B);
    nCol = mxGetN(A);
    Ap = mxGetJc(A);
    Ai = mxGetIr(A);
    colMatElem = mxGetPr(A);
    
    /* Copy data for MKL RCI */
    colMatBeg = (int *) mxCalloc(nCol + 1, sizeof(int));
    colMatIdx = (int *) mxCalloc(Ap[nCol], sizeof(int));
    
    for ( int i = 0; i < Ap[nCol]; ++i ) {
        colMatIdx[i] = Ai[i];
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        colMatBeg[i + 1] = Ap[i + 1];
    }
    
    if ( !colMatBeg || !colMatIdx || !colMatElem ) {
        mexErrMsgTxt("Failed to extract A. \n");
    }
    
    /* Get B */
    if (!mxIsDouble(B)) {
        mexErrMsgTxt("B must be a double array \n");
    }
    
    rhsMat = mxGetPr(B);
    
    double cgTol = *mxGetPr(tol);
    int maxIteration = (int) *mxGetPr(maxIter);
    
    X = mxCreateDoubleMatrix(nCol, nRhs, mxREAL);
    lhsSol = mxGetPr(X);
    
    /* Working array */
    int *iParam = (int *) mxCalloc(128 + 2 * nRhs, sizeof(int));
    double *dParam = (double *) mxCalloc(128 + 2 * nRhs, sizeof(double));
    double *tmpArray = (double *) mxCalloc(nCol * (3 + nRhs), sizeof(double));
    int *iterCount = (int *) mxCalloc(nRhs, sizeof(int));
    
    int retcode = rciCGSolve(nCol, nRhs, colMatBeg, colMatIdx, colMatElem,
                             rhsMat, lhsSol, dParam, iParam,
                             tmpArray, cgTol, maxIteration, iterCount);
    
    if ( retcode != 0 ) {
        mexErrMsgTxt("CG Failed \n");
    }
    
exit_cleanup:
    
    mxFree(iParam); iParam = NULL;
    mxFree(dParam); dParam = NULL;
    mxFree(tmpArray); tmpArray = NULL;
    mxFree(iterCount); iterCount = NULL;
    mxFree(colMatBeg); colMatBeg = NULL;
    mxFree(colMatIdx); colMatIdx = NULL;
    return;
}
