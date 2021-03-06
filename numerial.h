#pragma once

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void pause();

// return with resouce free
#define HANDLE_MUST(varname)                                                   \
  int varname;                                                                 \
  if (0 == 0) {                                                                \
  NUM_OK:                                                                      \
    varname = 1;                                                               \
  } else {                                                                     \
  NUM_ERROR:                                                                   \
    varname = 0;                                                               \
  }

#define MUST_return_ok() goto NUM_OK
#define MUST_return_error()                                                    \
  {                                                                            \
    fprintf(stderr, "%s:%d\t%s()\n", __FILE__, __LINE__, __func__);            \
    goto NUM_ERROR;                                                            \
  }
#define MUST(exp)                                                              \
  if (!(exp)) {                                                                \
    fprintf(stderr, "%s:%d\t%s\n", __FILE__, __LINE__, #exp);                  \
    goto NUM_ERROR;                                                            \
  }

#define newm(type, size) ((type *)malloc(sizeof(type) * (size)))

// i: int, after each{}, i should equal to n
// n: int, MUST greater than 0.
#define each(i, n) for (i = 0; i < (n); i++)

#define range(i, from, to, delta) for (i = (from); i <= (to); i += (delta))

typedef struct {
  // set by caller
  int max_step;
  double tol;

  // set by callee
  int used_step;
} iter_conf_t;

void init_iter_conf(iter_conf_t *conf);

typedef struct {
  int m;
  int n;
  double *data;
} mat_t;

// m*n matrix, use it to read and write.
// i,j are both one-based.
#define mat_v(mat, i, j) ((mat).data[(((i)-1) * ((mat).n) + ((j)-1))])

#define mat_each(mat, i, j) range(i, 1, (mat).m, 1) range(j, 1, (mat).n, 1)

#define init_mat(mat, i, j, exp)                                               \
  mat_each((mat), i, j) { mat_v((mat), i, j) = (exp); }

// MUST call free_mat after using.
// These all will init as 0.
mat_t new_mat(int m, int n);
mat_t new_mat_vec(int n); // vertical vector

mat_t new_mat_clone(mat_t src);

#define vec_v(mat, i) (mat.data[(i)-1])

void free_mat(mat_t *mat);

// format: "12.8" default if format == NULL
void mat_println(const char *format, mat_t mat);

// share same memories.
int mat_is_identical(mat_t A, mat_t B);

// T can be A.
int mat_transpose(mat_t T, mat_t A);

// R = A*B
// R should not be A or B.
// return: 0(error) 1(ok)
// ADVICE: assert(mat_product(R,A,B);
int mat_product(mat_t R, mat_t A, mat_t B);

// R = A+B
// R can be A or B.
int mat_add(mat_t R, mat_t A, mat_t B);

// R = A-B
// R can be A or B.
int mat_sub(mat_t R, mat_t A, mat_t B);

// R = cA
// R can be A.
int mat_scaler(mat_t R, mat_t A, double c);

int mat_assign(mat_t mat, double *data);
int mat_copy(mat_t dst, mat_t src);
int mat_is_same_size(mat_t A, mat_t B);
int mat_is_same_size_3(mat_t A, mat_t B, mat_t C);

int mat_back_substitution(mat_t A, mat_t x, mat_t b);
int mat_back_substitution_L(mat_t x, mat_t L, mat_t b);
int mat_back_substitution_U(mat_t x, mat_t U, mat_t b);

int mat_reduction_qr(mat_t Q, mat_t R, mat_t A);
int mat_reduction_qr_givens(mat_t Q, mat_t R, mat_t A);
int mat_reduction_qr_household(mat_t Q, mat_t R, mat_t A);

// T cannot be A.
int mat_inv_qr(mat_t T, mat_t A,
               int (*reduction_qr)(mat_t Q, mat_t R, mat_t A));

int mat_solve_qr(mat_t x, mat_t A, mat_t b,
                 int (*reduction_qr)(mat_t Q, mat_t R, mat_t A));

int mat_reduction_lu_doolittle(mat_t L, mat_t U, mat_t A);
int mat_reduction_lu_crout(mat_t L, mat_t U, mat_t A);
int mat_det_lu(double *det, mat_t A,
               int (*mat_reduction_lu)(mat_t L, mat_t U, mat_t A));

int mat_reduction_plu_doolittle(mat_t P, mat_t L, mat_t U, mat_t A);

// A = L * transpose(L)
// A must be a positive definite symmetric matrix
int mat_reduction_llt_cholesky(mat_t L, mat_t A);

int mat_inv_L(mat_t T, mat_t L);

// x = Hx + g
// if step <= 0, the iteraion will never stop until error < eps.
// this may cause a dead loop, PLEASE be careful.
int mat_solve_iter_simple(mat_t x, mat_t H, mat_t g, iter_conf_t *conf);
int mat_solve_iter_seidel(mat_t x, mat_t H, mat_t g, iter_conf_t *conf);

int mat_solve_iter_jacobi(mat_t x, mat_t A, mat_t b, iter_conf_t *conf);
int mat_solve_iter_gauss_seidel(mat_t x, mat_t A, mat_t b, iter_conf_t *conf);

int mat_solve_iter_sor(mat_t x, mat_t A, mat_t b, double factor,
                       iter_conf_t *conf);

int mat_solve_iter_steepest_descent(mat_t x, mat_t A, mat_t b,
                                    iter_conf_t *conf);
int mat_solve_iter_conjugate_gradient(mat_t x, mat_t A, mat_t b,
                                      iter_conf_t *conf);

// eigenvalue & eigenvector
int mat_eigen_power_method(double *eigen_v, mat_t eigen_vec, mat_t A,
                           iter_conf_t *conf);

double mat_norm_1(mat_t A);
double mat_norm_F(mat_t A);
double mat_norm_inf(mat_t A);

int mat_cond(double *cond, mat_t A, double (*mat_norm)(mat_t A));

double vec_norm_1(mat_t A);
double vec_norm_2(mat_t A);
double vec_norm_inf(mat_t A);

int vec_inner_product(double *v, mat_t x, mat_t y);

// root
typedef double (*fn_t)(double x);

int root_bisection(double *root, fn_t fn, double a, double b,
                   iter_conf_t *conf);

// You MUST set *x as init value, otherwise, BUG will occur.
int root_iter_fixed_point(double *x, fn_t fn, iter_conf_t *conf);

// fn : f(x)
// fn1: f'(x)
int root_iter_newton(double *x, fn_t fn, fn_t fn1, iter_conf_t *conf);

int root_iter_newton_down(double *x, fn_t f, fn_t f1, double min_lambda,
                          double x_delta, iter_conf_t *conf);

// Pn(x) = a[0]*x^n + a[1]*x^(n-1) + ... + a[n]
// len(P) = n+1
// n: degree of polynomial
double poly_value(double a[], int n, double x);

// MUST free after use its return.
double *exp_taylor(int n);

double vec_dot_product(double a[], double b[], int dim);

// The "itp_" prefix means that interpolation.

// The x and y MUST be vectors and have the same size.
int itp_lagrange(double *y, mat_t X, mat_t Y, double x);
