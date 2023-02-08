### Optimal Diagonal Preconditioning Toolbox

This repository contains implementations and detailed experiment results for the paper [Optimal Diagonal Preconditioning](https://arxiv.org/abs/2209.00809).

**Prerequisites**

The toolbox relies on solving semidefinite programs and needs to invoke optimization solvers through CVX. Users have to install

- CVX  MATLAB toolbox  at http://cvxr.com/cvx/

For testing purposes users need to install

- HDSDP at https://github.com/COPT-Public/HDSDP
- LIBSVM datasets at https://www.csie.ntu.edu.tw/~cjlin/libsvm/
- SuiteSparse matlab interface at https://sparse.tamu.edu/interfaces
- Intel MKL at https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html

and ensure that they are added to MATLAB search path.

**Install**

To install the package, simple clone the repo with

```
git clone https://github.com/Gwzwpxz/opt_dpcond.git
```

and run 

```
optpcd_setup
```

in MATLAB. Upon installation, a toy example will be solved to test the installation.

```
Setting up the optimal diagonal preconditioning repo 
The repo depends on:  

1. CVX   matlab toolbox    at http://cvxr.com/cvx/
2. HDSDP binary (optional) at https://github.com/COPT-Public/HDSDP 

To reproduce the experiments, execute 

    git checkout opt-precond 

and ensure that the following data repos are installed 

3. LIBSVM datasets              at https://www.csie.ntu.edu.tw/~cjlin/libsvm/
4. SuiteSparse matlab interface at https://sparse.tamu.edu/interfaces

Running a toy example.

Solving a two-sided preconditioning problem using bisection 
     kappa    ub - lb 
 5.000e+05  1.000e+06 
 7.500e+05  5.000e+05 
 6.250e+05  2.500e+05 
 5.625e+05  1.250e+05 
 5.313e+05  6.250e+04 
 5.156e+05  3.125e+04 
 5.234e+05  1.562e+04 
Condition number of X  : 2.06e+03 
Condition number of XE : 1.62e+03 
Condition number of DX : 1.28e+03 
Condition number of DXE: 7.23e+02 
Exporting SDP to ./Eprob-R.dat-s 

Installation completes. Check README.md for usage details. 
```

**Usage**

The optimal preconditioning toolbox provides several utilities that allow users to 

- Compute the optimal Left/Right/Two-sided preconditioner
- Export the Left/Right preconditioning SDP to standard SDPA format

The basic usage  can be demonstrated by the following lines of code

```matlab
% Choose preconditioning type
param.ptype = 'R'; 
% Generate preconditioning problem. X is the user data
prob = getoptprob(X, param); 
% Solve the preconditioning problem
sol = optprecond(prob); 
% Get preconditioned matrix
pX = sol.pX;
% Get diagonal preconditioner
E = sol.E; 
% Export problem to SDPA
path = fullfile('.');
pname = 'Eprob';
exportoptprob(prob, path, pname);
```

**Testing**

The experiments from paper `Optimal Diagonal Preconditioning: Theory and Practice` are available at the testing branch, which can be accessed through

```
git checkout opt-precond 
```

in the command line. The tests can be done by modifying and running the testing scripts from `test` directory

- test_precond_libsvm (LIBSVM)
- test_precond_suitesparse (SuiteSparse)
- test_precond_random (Random)

To test preconditioned CG, the users need to install CG from Intel RCI interface using Cmake build system from `utils/rci`

```
mkdir build
cd build
cmake ..
make
```

and obtain the mexfile for CG implementation of multiple RHSs.

**Contact**

Please contact `gwz@163.shufe.edu.cn` for questions on the toolbox.

**Cite as**

> Qu, Z., Gao, W., Hinder, O., Ye, Y., & Zhou, Z. (2022). Optimal Diagonal Preconditioning. *arXiv preprint arXiv:2209.00809*.
