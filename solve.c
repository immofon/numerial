#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int mat_solve_iter_jacobi(mat_t x, mat_t A, mat_t b, iter_conf_t *conf);
int mat_solve_iter_gauss_seidel(mat_t x, mat_t A, mat_t b, iter_conf_t *conf);

void pause() {
#ifdef _WINDOWS__
  system("PAUSE");
#endif
#ifdef _WIN32
  system("PAUSE");
#endif
}

void init_iter_conf(iter_conf_t *conf) {
  if (conf->max_step <= 1) {
    conf->max_step = 1000;
  }
  if (conf->tol <= 1e-12) {
    conf->tol = 1e-10;
  }
  conf->used_step = 0;
}

// MUST call free_mat after using.
mat_t new_mat(int m, int n) {
  assert(m > 0);
  assert(n > 0);

  int i, j;
  mat_t mat;
  mat.m = m;
  mat.n = n;
  mat.data = newm(double, m *n);
  init_mat(mat, i, j, 0);
  return mat;
}
void free_mat(mat_t *mat) {
  if (mat->data == NULL) {
    return;
  }

  free(mat->data);
  mat->data = NULL;
}

// vertical vector
mat_t new_mat_vec(int n) { return new_mat(n, 1); }

mat_t new_mat_clone(mat_t src) {
  assert(src.data != NULL);

  mat_t R = new_mat(src.m, src.n);
  assert(mat_copy(R, src));
  return R;
}

int mat_is_same_size(mat_t A, mat_t B) { return A.m == B.m && A.n == B.n; }
int mat_is_same_size_3(mat_t A, mat_t B, mat_t C) {
  return mat_is_same_size(A, B) && mat_is_same_size(B, C);
}

// share same memories.
int mat_is_identical(mat_t A, mat_t B) {
  if (!mat_is_same_size(A, B)) {
    return 0;
  }
  if (A.data == NULL || B.data == NULL) {
    return 0;
  }
  if (A.data == B.data) {
    return 1;
  }
  return 0;
}

// R = A*B
// return: 0(error) 1(ok)
// ADVICE: assert(mat_product(R,A,B);
int mat_product(mat_t R, mat_t A, mat_t B) {
  MUST(!mat_is_identical(R, A));
  MUST(!mat_is_identical(R, B));
  MUST(R.m == A.m && R.n == B.n && A.n == B.m);

  int i, j, k;
  double sum;

  range(i, 1, A.m, 1) range(j, 1, B.n, 1) {
    sum = 0;
    range(k, 1, A.n, 1) { sum += mat_v(A, i, k) * mat_v(B, k, j); }
    mat_v(R, i, j) = sum;
  }

  HANDLE_MUST(ret);
  return ret;
}

// R = A+B
// R can be A or B.
int mat_add(mat_t R, mat_t A, mat_t B) {
  int i, j;
  if (!mat_is_same_size_3(R, A, B)) {
    return 0;
  }

  mat_each(R, i, j) { mat_v(R, i, j) = mat_v(A, i, j) + mat_v(B, i, j); }
  return 1;
}

// R = A+B
// R can be A or B.
int mat_sub(mat_t R, mat_t A, mat_t B) {
  if (!mat_is_same_size_3(R, A, B)) {
    return 0;
  }

  int i, j;
  mat_each(R, i, j) { mat_v(R, i, j) = (mat_v(A, i, j)) - (mat_v(B, i, j)); }
  return 1;
}

// R = cA
// R can be A.
int mat_scaler(mat_t R, mat_t A, double c) {
  int i, j;
  if (!mat_is_same_size(R, A)) {
    return 0;
  }

  mat_each(R, i, j) { mat_v(R, i, j) = c * mat_v(A, i, j); }
  return 1;
}

void print_double(const char *format, double v) {
  int i;

  size_t len = strlen(format) + 3 + 1; // "%lf" =>  3
  char *fmt = newm(char, len);
  for (i = 0; i < len; i++) {
    fmt[i] = 0;
  }

  strcat(fmt, "%");
  strcat(fmt, format);
  strcat(fmt, "lf");

  printf(fmt, v);

  free(fmt);
}

// format: "12.8" default if format == NULL
void mat_println(const char *format, mat_t mat) {
  int i, j;

  if (format == NULL) {
    format = "12.8";
  }

  range(i, 1, mat.m, 1) {
    if (i == 1) {
      printf("[");
    } else {
      printf(" ");
    }

    range(j, 1, mat.n, 1) {
      print_double(format, mat_v(mat, i, j));
      printf(" ");
    }

    if (i == mat.m) {
      printf("]");
    }

    printf("\n");
  }
}

int mat_assign(mat_t mat, double *data) {
  assert(mat.data != NULL);
  assert(data != NULL);

  int size = mat.m * mat.n - 1;
  int i;
  range(i, 0, size, 1) { mat.data[i] = data[i]; }
  return 1;
}

int mat_copy(mat_t dst, mat_t src) {
  if (dst.m != dst.m || src.n != src.n) {
    return 0;
  }

  if (dst.data == NULL || src.data == NULL) {
    return 0;
  }

  return mat_assign(dst, src.data);
}

double _mat_inter_product(mat_t A, mat_t B, int i, int j) {
  assert(A.m == B.m);
  assert(1 <= i & i <= A.n);
  assert(1 <= j & j <= B.n);

  int k;
  double sum = 0;
  range(k, 1, A.m, 1) { sum += mat_v(A, k, i) * mat_v(B, k, j); }
  return sum;
}

// T can be A.
int mat_transpose(mat_t T, mat_t A) {
  if (T.m != A.n || T.n != A.m) {
    return 0;
  }

  double tmp;
  int i, j;

  if (mat_is_identical(T, A) && A.m == A.n) {
    range(i, 1, A.n, 1) range(j, i + 1, A.n, 1) {
      tmp = mat_v(A, i, j);
      mat_v(A, i, j) = mat_v(A, j, i);
      mat_v(A, j, i) = tmp;
    }
    return 1;
  }

  mat_each(T, i, j) { mat_v(T, i, j) = mat_v(A, j, i); }
  return 1;
}

int mat_solve_iter_gauss_seidel(mat_t x, mat_t A, mat_t b, iter_conf_t *conf) {
  MUST(A.m == A.n);
  MUST(A.n == b.m);
  MUST(b.n == 1);

  MUST(mat_is_same_size(x, b));

  init_iter_conf(conf);

  int i, j;
  int n = A.m;
  double sum;
  double m;

  range(i, 1, n, 1) { MUST(mat_v(A, i, i) != 0); }

  range(conf->used_step, 1, conf->max_step, 1) {
    m = 0;
    range(i, 1, n, 1) {
      sum = -vec_v(b, i);
      range(j, 1, i - 1, 1) { sum += mat_v(A, i, j) * vec_v(x, j); }
      range(j, i + 1, n, 1) { sum += mat_v(A, i, j) * vec_v(x, j); }
      sum = -(sum / mat_v(A, i, i));

      if (fabs(sum - vec_v(x, i)) > m) {
        m = fabs(sum - vec_v(x, i));
      }
      vec_v(x, i) = sum;
    }

    if (m < conf->tol) {
      MUST_return_ok();
    }
  }
  MUST_return_error();

  HANDLE_MUST(ret);
  return ret;
}

int mat_solve_iter_jacobi(mat_t x, mat_t A, mat_t b, iter_conf_t *conf) {
  mat_t y = new_mat_vec(x.m);

  MUST(A.m == A.n);
  MUST(A.n == b.m);
  MUST(b.n == 1);

  MUST(mat_is_same_size(x, b));

  init_iter_conf(conf);

  int i, j;
  int n = A.m;
  double sum;
  double m;

  range(i, 1, n, 1) { MUST(mat_v(A, i, i) != 0); }

  range(conf->used_step, 1, conf->max_step, 1) {
    m = 0;
    range(i, 1, n, 1) {
      sum = -vec_v(b, i);
      range(j, 1, i - 1, 1) { sum += mat_v(A, i, j) * vec_v(x, j); }
      range(j, i + 1, n, 1) { sum += mat_v(A, i, j) * vec_v(x, j); }
      sum = -(sum / mat_v(A, i, i));

      if (fabs(sum - vec_v(x, i)) > m) {
        m = fabs(sum - vec_v(x, i));
      }
      vec_v(y, i) = sum;
    }
    MUST(mat_copy(x, y));

    if (m < conf->tol) {
      MUST_return_ok();
    }
  }
  MUST_return_error();

  HANDLE_MUST(ret);
  free_mat(&y);
  return ret;
}

int find_index(mat_t X, double x) {
  assert(X.n == 1);
  assert(X.m >= 2);

  int i;

  if (x <= vec_v(X, 1)) {
    return 1;
  }

  range(i, 1, X.m - 1, 1) {
    if (vec_v(X, i) <= x && x <= vec_v(X, i + 1)) {
      return i;
    }
    //
  }

  return X.m - 1;
}

double cubic_spline_intp(mat_t X, mat_t Y, mat_t M, double x) {
  assert(X.n == 1);
  assert(X.m >= 2);
  assert(mat_is_same_size(X, M));

  int i = find_index(X, x);

  double h, t;
  double v = 0;

  h = vec_v(X, i + 1) - vec_v(X, i); // h_{i+1}
  t = x - vec_v(X, i);               // x - x_i
  v += (vec_v(Y, i + 1) - vec_v(M, i + 1) * h * h / 6) * t / h;
  v += t * t * t * vec_v(M, i + 1) / (6 * h);

  t = x - vec_v(X, i + 1); // x - x_{i+1}
  v -= (vec_v(Y, i) - vec_v(M, i) * h * h / 6) * t / h;
  v -= t * t * t * vec_v(M, i) / (6 * h);

  return v;
}

mat_t new_mat_vec_uniform_partition(int n, double a, double b) {
  assert(n >= 2);
  assert(a < b);

  mat_t X = new_mat_vec(n);
  double h = (b - a) / (n - 1);
  int i;
  range(i, 1, n, 1) {
    vec_v(X, i) = a + (i - 1) * h;
    //
  }
  vec_v(X, n) = b;

  return X;
}

mat_t new_mat_function_image(mat_t X, double (*f)(double x)) {
  mat_t Y = new_mat(X.m, X.n);
  int i, j;

  mat_each(X, i, j) { mat_v(Y, i, j) = (*f)(mat_v(X, i, j)); }

  return Y;
}

// You should set boundary condition before solve the equation
mat_t new_mat_cubic_intp_A(mat_t X, mat_t Y) {
  assert(X.m >= 2);
  assert(X.n == 1);
  assert(mat_is_same_size(X, Y));

  int n = X.m;
  mat_t A = new_mat(n, n);
  int i;
  range(i, 1, n - 2, 1) {
    mat_v(A, i + 1, i) =
        (vec_v(X, i + 1) - vec_v(X, i)) / (vec_v(X, i + 2) - vec_v(X, i));
    mat_v(A, i + 1, i + 1) = 2;
    mat_v(A, i + 1, i + 2) = 1 - mat_v(A, i + 1, i);
  }
  return A;
}
mat_t new_mat_vec_cubic_intp_d(mat_t X, mat_t Y) {
  assert(X.m >= 2);
  assert(X.n == 1);
  assert(mat_is_same_size(X, Y));

  int n = X.m;
  mat_t d = new_mat_vec(n);
  int i;

  range(i, 1, n - 2, 1) {
    vec_v(d, i + 1) =
        6 *
        (((vec_v(Y, i + 2) - vec_v(Y, i + 1)) /
          (vec_v(X, i + 2) - vec_v(X, i + 1))) -
         ((vec_v(Y, i + 1) - vec_v(Y, i)) / (vec_v(X, i + 1) - vec_v(X, i)))) /
        (vec_v(X, i + 2) - vec_v(X, i));
  }

  return d;
}

void cubic_intp_set_clamped_boundary(mat_t A, mat_t d, mat_t X, mat_t Y,
                                     double dy1, double dyn) {
  int n = A.m;
  assert(A.m == A.n);
  assert(d.m == n);
  assert(n >= 2);
  assert(d.n == 1);
  assert(mat_is_same_size_3(d, X, Y));

  mat_v(A, 1, 1) = 2;
  mat_v(A, 1, 2) = 1;
  mat_v(A, n, n - 1) = 1;
  mat_v(A, n, n) = 2;

  double h = vec_v(X, 2) - vec_v(X, 1);
  vec_v(d, 1) = 6 * ((vec_v(Y, 2) - vec_v(Y, 1)) / h - dy1) / h;
  h = vec_v(X, n) - vec_v(X, n - 1);
  vec_v(d, n) = 6 * (dyn - (vec_v(Y, n) - vec_v(Y, n - 1)) / h) / h;
}

void cubic_intp_set_natural_boundary(mat_t A, mat_t d, mat_t X, mat_t Y) {
  int n = A.m;
  assert(A.m == A.n);
  assert(d.m == n);
  assert(n >= 2);
  assert(d.n == 1);
  assert(mat_is_same_size_3(d, X, Y));

  mat_v(A, 1, 1) = 2;
  mat_v(A, 1, 2) = 0;
  mat_v(A, n, n - 1) = 0;
  mat_v(A, n, n) = 2;

  vec_v(d, 1) = 0;
  vec_v(d, n) = 0;
}

int test_solve() {
  double A_v[] = {
      2, 0, 0, //
      0, 2, 0, //
      0, 0, 2, //
  };
  mat_t A = new_mat(3, 3);
  assert(mat_assign(A, A_v));

  double b_v[] = {
      1,
      2,
      3,
  };
  mat_t b = new_mat_vec(3);
  assert(mat_assign(b, b_v));

  mat_t x = new_mat_vec(3);

  iter_conf_t conf;
  conf.tol = 1e-6;
  conf.max_step = 100;

  assert(mat_solve_iter_jacobi(x, A, b, &conf));

  mat_println(".8", x);
  return 0;
}

double f(double x) { return x * x * x * x * x; }

int test(int n) {
  int a = -1;
  int b = 1;
  int N = 10000;
  mat_t xx = new_mat_vec_uniform_partition(N, a, b);
  mat_t X = new_mat_vec_uniform_partition(n, a, b);
  mat_t Y = new_mat_function_image(X, &f);
  double dy1 = 0;
  double dyn = 5;

  mat_t A = new_mat_cubic_intp_A(X, Y);
  mat_t d = new_mat_vec_cubic_intp_d(X, Y);
  cubic_intp_set_clamped_boundary(A, d, X, Y, dy1, dyn);
  // cubic_intp_set_natural_boundary(A, d, X, Y);
  mat_t M = new_mat_vec(d.m);

  // mat_println(".6", A);
  // mat_println(".6", d);

  iter_conf_t conf;
  conf.tol = 1e-17;
  conf.max_step = 1000;
  assert(mat_solve_iter_jacobi(M, A, d, &conf));

  double v;
  v = cubic_spline_intp(X, Y, M, 1.0 / 3.0);
  printf("n = %d; err = %.10e\n", n, v - f(1.0 / 3.0));

  FILE *io = fopen("spline.dat", "w");
  int i;
  range(i, 1, N, 1) {
    v = cubic_spline_intp(X, Y, M, vec_v(xx, i));
    fprintf(io, "%e\t%e\n", vec_v(xx, i), v - f(vec_v(xx, i)));
    // fprintf(io, "%e\t%e\n", vec_v(xx, i), v);
  }
  fclose(io);

  return 0;
}

int main() {
  int i;
  range(i, 2, 1000, 1) { test(i); }
}
