#include <time.h>
#include <lua.h>
#include <lauxlib.h>
#include "numerial.h"

#define MAT_T "lua-num.mat"

#define DOUBLE_EQUAL_DELTA 1e-10

static int l_test(lua_State *L) {
	printf("hello luatest\n");
	return 0;
}

static mat_t* l_new_mat(lua_State *L,int m,int n) {
	size_t size = sizeof(mat_t) + m * n * sizeof(double);
	mat_t *A = (mat_t *)lua_newuserdata(L,size);
	A->m = m;
	A->n = n;
	A->data = (double *)(A+1);
	int i;
	range(i,0,m*n-1,1) {
		A->data[i] = 0;
	}
	luaL_getmetatable(L,MAT_T);
	lua_pushvalue(L,-1);
	lua_setfield(L,-2,"__index");
	lua_setmetatable(L,-2);
	return A;
}

static int l_mat(lua_State *L) {
	int m = (int)luaL_checkinteger(L,1);
	int n = (int)luaL_checkinteger(L,2);
	if(m <= 0) {
		luaL_argerror(L,1,"less than 1");
	}
	if(n <= 0) {
		luaL_argerror(L,2,"less than 1");
	}

	mat_t* A = l_new_mat(L,m,n);

	return 1;
}

static int l_mat_get(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	int i = (int)luaL_checkinteger(L,2);
	int j = (int)luaL_checkinteger(L,3);
	luaL_argcheck(L,1<=i && i <= A->m ,2, "out of range");
	luaL_argcheck(L,1<=j && j <= A->n ,3, "out of range");

	lua_pushnumber(L,mat_v(*A,i,j));
	return 1;
}

static int l_mat_set(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	int i = (int)luaL_checkinteger(L,2);
	int j = (int)luaL_checkinteger(L,3);
	double v = (double)luaL_checknumber(L,4);
	luaL_argcheck(L,1<=i && i <= A->m ,2, "out of range");
	luaL_argcheck(L,1<=j && j <= A->n ,3, "out of range");

	mat_v(*A,i,j) = v;
	return 0;
}

//mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
static int l_mat_add(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	mat_t *B = (mat_t *) luaL_checkudata(L,2,MAT_T);
	luaL_argcheck(L,A->m == B->m,2,"not same size");
	luaL_argcheck(L,A->n == B->n,2,"not same size");

	mat_t *R = l_new_mat(L,A->m,A->n);
	if (mat_add(*R,*A,*B)) {
		return 1;
	}

	luaL_error(L,"numerial: mat_add");
	return 0;
}

static int l_mat_sub(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	mat_t *B = (mat_t *) luaL_checkudata(L,2,MAT_T);
	luaL_argcheck(L,A->m == B->m,2,"not same size");
	luaL_argcheck(L,A->n == B->n,2,"not same size");

	mat_t *R = l_new_mat(L,A->m,A->n);
	if (mat_sub(*R,*A,*B)) {
		return 1;
	}

	luaL_error(L,"numerial: mat_sub");
	return 0;
}

static int l_mat_mul(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	mat_t *B = (mat_t *) luaL_checkudata(L,2,MAT_T);
	luaL_argcheck(L,A->n == B->m,2,"can't product");

	mat_t *R = l_new_mat(L,A->m,B->n);
	if (mat_product(*R,*A,*B)) {
		return 1;
	}

	luaL_error(L,"numerial: mat_product");
	return 0;
}

static int l_mat_assign(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_checktype(L,2,LUA_TTABLE);
	lua_len(L,2);
	int size = (int)lua_tointeger(L,-1);
	luaL_argcheck(L,size == A->m*A->n,2,"not same size");
	lua_pop(L,1);

	int i,j;
	double v;
	mat_each(*A,i,j) {
		if (lua_geti(L,2,(i-1)*(A->n) + j) == LUA_TNUMBER) {
			v = (double)lua_tonumber(L,-1);
			lua_pop(L,1);
			mat_v(*A,i,j) = v;
		}else {
			luaL_argerror(L,2,"not number array");
		}
	}

	lua_pushvalue(L,1);
	return 1;
}

static int l_mat_transpose(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	mat_t *T= (mat_t *) l_new_mat(L,A->n,A->m); // T.m == A.n, T.n == A.m

	if(mat_transpose(*T,*A)) {
		return 1;
	}
	luaL_error(L,"mat_transpose");
	return 0;
}

static int l_mat_det(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,A->m==A->n,1,"expect square matrix");

	double det = 0;
	if(mat_det_lu(&det,*A,&mat_reduction_lu_doolittle)) {
		lua_pushnumber(L,det);
		return 1;
	}
	luaL_error(L,"mat_det_lu");
	return 0;
}

static int l_mat_println(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	const char *fmt = luaL_optstring(L,2,"");

	mat_println(fmt,*A);
	return 0;
}

static int l_mat_qr(lua_State *L,int (*reduction_qr)(mat_t Q,mat_t R,mat_t A)) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	mat_t *Q = (mat_t *) l_new_mat(L,A->m,A->m);
	mat_t *R = (mat_t *) l_new_mat(L,A->m,A->n);
	if((reduction_qr)(*Q,*R,*A)) {
		return 2;
	}
	luaL_error(L,"l_mat_qr");
	return 0;
}

static int l_mat_qr_givens(lua_State *L) {
	return l_mat_qr(L,&mat_reduction_qr_givens);
}

static int l_mat_qr_household(lua_State *L) {
	return l_mat_qr(L,&mat_reduction_qr_household);
}

static int l_mat_plu_doolittle(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,A->m==A->n,1,"expect square matrix");

	mat_t *P = (mat_t *) l_new_mat(L,A->m,A->m);
	mat_t *L1 = (mat_t *) l_new_mat(L,A->m,A->m);
	mat_t *U = (mat_t *) l_new_mat(L,A->m,A->m);
	if (mat_reduction_plu_doolittle(*P,*L1,*U,*A)) {
		return 3;
	}
	luaL_error(L,"mat_reduction_plu_doolittle");
	return 0;
}

static int l_mat_llt_cholesky(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,A->m==A->n,1,"expect square matrix");

	mat_t *L1 = (mat_t *) l_new_mat(L,A->m,A->m);
	mat_t *L1t = (mat_t *) l_new_mat(L,A->m,A->m);
	if (mat_reduction_llt_cholesky(*L1,*A)) {
		if (mat_transpose(*L1t,*L1)) {
			return 2;
		}
		luaL_error(L,"mat_transpose");
	}
	luaL_error(L,"mat_reduction_llt_cholesky");
	return 0;
}

static int l_mat_solve_L(lua_State *L) {
	mat_t *L1 = (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,L1->m==L1->n,1,"expect square matrix");
	mat_t *b = (mat_t *) luaL_checkudata(L,2,MAT_T);
	luaL_argcheck(L,L1->n==b->m,2,"worry size");
	luaL_argcheck(L,b->n==1,2,"expect vector");


	mat_t *x = (mat_t *) l_new_mat(L,b->m,b->n);
	if(mat_back_substitution_L(*x,*L1,*b)) {
		return 1;
	}

	luaL_error(L,"mat_back_substitution_L");
	return 0;
}

static int l_mat_solve_U(lua_State *L) {
	mat_t *U= (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,U->m==U->n,1,"expect square matrix");
	mat_t *b = (mat_t *) luaL_checkudata(L,2,MAT_T);
	luaL_argcheck(L,U->n==b->m,2,"worry size");
	luaL_argcheck(L,b->n==1,2,"expect vector");


	mat_t *x = (mat_t *) l_new_mat(L,b->m,b->n);
	if(mat_back_substitution_U(*x,*U,*b)) {
		return 1;
	}

	luaL_error(L,"mat_back_substitution_U");
	return 0;
}

static int l_mat_inv_L(lua_State *L) {
	mat_t *L1= (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,L1->m==L1->n,1,"expect square matrix");

	mat_t *T = (mat_t *) l_new_mat(L,L1->m,L1->n);
	if(mat_inv_L(*T,*L1)) {
		return 1;
	}
	luaL_error(L,"mat_inv_L");
	return 0;
}

// load table at index if exist.
void opt_load_iter_conf(lua_State *L, int index, iter_conf_t *conf) {
	conf->used_step = 0;
	conf->max_step = 0;
	conf->tol = 0;
	if(lua_istable(L,index)) {
		lua_getfield(L,index,"step");
		conf->max_step = (int)lua_tointeger(L,-1);
		lua_pop(L,1);
		lua_getfield(L,index,"tol");
		conf->tol = (double)lua_tonumber(L,-1);
		lua_pop(L,1);
	}
	if (conf->tol == 0) {
		conf->tol = DOUBLE_EQUAL_DELTA/10;
	}
	init_iter_conf(conf);
}

static int l_mat_solve_iter_simple(lua_State *L) {
	mat_t *H= (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,H->m==H->n,1,"expect square matrix");
	mat_t *g= (mat_t *) luaL_checkudata(L,2,MAT_T);
	luaL_argcheck(L,g->n == 1,2,"expect vector");
	luaL_argcheck(L,H->n == g->m,2,"worry size vector");

	iter_conf_t conf;
	opt_load_iter_conf(L,3,&conf);
	mat_t *x = (mat_t *) l_new_mat(L,g->m,g->n);
	if(mat_solve_iter_simple(*x,*H,*g,&conf)) {
		lua_pushinteger(L,conf.used_step);
		return 2;
	}
	luaL_error(L,"mat_solve_iter_simple");
	return 0;
}


static int l_mat_solve_iter_seidel(lua_State *L) {
	mat_t *H= (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,H->m==H->n,1,"expect square matrix");
	mat_t *g= (mat_t *) luaL_checkudata(L,2,MAT_T);
	luaL_argcheck(L,g->n == 1,2,"expect vector");
	luaL_argcheck(L,H->n == g->m,2,"worry size vector");

	iter_conf_t conf;
	opt_load_iter_conf(L,3,&conf);
	mat_t *x = (mat_t *) l_new_mat(L,g->m,g->n);
	if(mat_solve_iter_seidel(*x,*H,*g,&conf)) {
		lua_pushinteger(L,conf.used_step);
		return 2;
	}
	luaL_error(L,"mat_solve_iter_seidel");
	return 0;
}


static int l_mat_solve_iter_jacobi(lua_State *L) {
	mat_t *A= (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,A->m==A->n,1,"expect square matrix");
	mat_t *b= (mat_t *) luaL_checkudata(L,2,MAT_T);
	luaL_argcheck(L,b->n == 1,2,"expect vector");
	luaL_argcheck(L,A->n == b->m,2,"worry size vector");

	iter_conf_t conf;
	opt_load_iter_conf(L,3,&conf);
	mat_t *x = (mat_t *) l_new_mat(L,b->m,b->n);
	if(mat_solve_iter_jacobi(*x,*A,*b,&conf)) {
		lua_pushinteger(L,conf.used_step);
		return 2;
	}
	luaL_error(L,"mat_solve_iter_jacobi");
	return 0;
}

static int l_mat_solve_iter_gauss_seidel(lua_State *L) {
	mat_t *A= (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,A->m==A->n,1,"expect square matrix");
	mat_t *b= (mat_t *) luaL_checkudata(L,2,MAT_T);
	luaL_argcheck(L,b->n == 1,2,"expect vector");
	luaL_argcheck(L,A->n == b->m,2,"worry size vector");

	iter_conf_t conf;
	opt_load_iter_conf(L,3,&conf);
	mat_t *x = (mat_t *) l_new_mat(L,b->m,b->n);
	if(mat_solve_iter_gauss_seidel(*x,*A,*b,&conf)) {
		lua_pushinteger(L,conf.used_step);
		return 2;
	}
	luaL_error(L,"mat_solve_iter_gauss_seidel");
	return 0;
}

static int l_mat_solve_iter_sor(lua_State *L) {
	mat_t *A= (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,A->m==A->n,1,"expect square matrix");
	mat_t *b= (mat_t *) luaL_checkudata(L,2,MAT_T);
	luaL_argcheck(L,b->n == 1,2,"expect vector");
	luaL_argcheck(L,A->n == b->m,2,"worry size vector");
	double factor = (double) luaL_checknumber(L,3);

	iter_conf_t conf;
	opt_load_iter_conf(L,4,&conf);
	mat_t *x = (mat_t *) l_new_mat(L,b->m,b->n);
	if(mat_solve_iter_sor(*x,*A,*b,factor,&conf)) {
		lua_pushinteger(L,conf.used_step);
		return 2;
	}
	luaL_error(L,"mat_solve_iter_sor");
	return 0;
}

static int l_mat_solve_iter_steepest_descent(lua_State *L) {
	mat_t *A= (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,A->m==A->n,1,"expect square matrix");
	mat_t *b= (mat_t *) luaL_checkudata(L,2,MAT_T);
	luaL_argcheck(L,b->n == 1,2,"expect vector");
	luaL_argcheck(L,A->n == b->m,2,"worry size vector");

	iter_conf_t conf;
	opt_load_iter_conf(L,3,&conf);
	mat_t *x = (mat_t *) l_new_mat(L,b->m,b->n);
	if(mat_solve_iter_steepest_descent(*x,*A,*b,&conf)) {
		lua_pushinteger(L,conf.used_step);
		return 2;
	}
	luaL_error(L,"mat_solve_iter_steepest_descent");
	return 0;
}

static int l_mat_solve_iter_conjugate_gradient(lua_State *L) {
	mat_t *A= (mat_t *) luaL_checkudata(L,1,MAT_T);
	luaL_argcheck(L,A->m==A->n,1,"expect square matrix");
	mat_t *b= (mat_t *) luaL_checkudata(L,2,MAT_T);
	luaL_argcheck(L,b->n == 1,2,"expect vector");
	luaL_argcheck(L,A->n == b->m,2,"worry size vector");

	iter_conf_t conf;
	opt_load_iter_conf(L,3,&conf);
	mat_t *x = (mat_t *) l_new_mat(L,b->m,b->n);
	if(mat_solve_iter_conjugate_gradient(*x,*A,*b,&conf)) {
		lua_pushinteger(L,conf.used_step);
		return 2;
	}
	luaL_error(L,"mat_solve_iter_conjugate_gradient");
	return 0;
}

static int l_mat_equal(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	mat_t *B = (mat_t *) luaL_checkudata(L,2,MAT_T);
	if (A->m != B->m || A->n != B->n) {
		lua_pushboolean(L,1);
		return 1;
	}

	double delta;
	int i,j;
	mat_each(*A,i,j) {
		delta = fabs(mat_v(*A,i,j) - mat_v(*B,i,j));
		if ((!(delta == delta)) || delta > DOUBLE_EQUAL_DELTA) {
			lua_pushboolean(L,0);
			return 1;
		}
	}

	lua_pushboolean(L,1);
	return 1;
}

static int l_test_timer(lua_State *L) {
	luaL_argcheck(L,lua_isfunction(L,1),1,"expect function");
	clock_t t;
	t = clock();
	lua_call(L,0,0);
	lua_pushnumber(L,(((double)(clock() - t))/CLOCKS_PER_SEC));
	return 1;
}

static const struct luaL_Reg numlib[] = {
	{"new",l_mat},
	{"get",l_mat_get},
	{"set",l_mat_set},
	{"println",l_mat_println},
	{NULL,NULL},
};

static const struct luaL_Reg mat_metareg[] ={
	{"get",l_mat_get},
	{"set",l_mat_set},
	{"println",l_mat_println},
	{"assign",l_mat_assign},
	{"transpose",l_mat_transpose},
	{"det",l_mat_det},
	{"qr",l_mat_qr_givens},
	{"qr_givens",l_mat_qr_givens},
	{"qr_household",l_mat_qr_household},
	{"plu",l_mat_plu_doolittle},
	{"plu_doolittle",l_mat_plu_doolittle},
	{"llt_cholesky",l_mat_llt_cholesky},
	{"solve_L",l_mat_solve_L},
	{"solve_U",l_mat_solve_U},
	{"inv_L",l_mat_inv_L},
	{"solve_iter_simple",l_mat_solve_iter_simple},
	{"solve_iter_seidel",l_mat_solve_iter_seidel},
	{"solve_iter_jacobi",l_mat_solve_iter_jacobi},
	{"solve_iter_gauss_seidel",l_mat_solve_iter_gauss_seidel},
	{"solve_iter_sor",l_mat_solve_iter_sor},
	{"solve_iter_steepest_descent",l_mat_solve_iter_steepest_descent},
	{"solve_iter_conjugate_gradient",l_mat_solve_iter_conjugate_gradient},
	{"__add",l_mat_add},
	{"__sub",l_mat_sub},
	{"__shl",l_mat_assign},
	{"__eq",l_mat_equal},
	{"__mul",l_mat_mul},
	{NULL,NULL},
};


static const struct luaL_Reg testlib[] ={
	{"timer",l_test_timer},
	{NULL,NULL},
};
int luaopen_num(lua_State *L) {
	luaL_newmetatable(L, MAT_T);
	luaL_setfuncs(L,mat_metareg,0);
	lua_newtable(L);
	luaL_newlib(L,numlib);
	lua_setfield(L,-2,"mat");
	luaL_newlib(L,testlib);
	lua_setfield(L,-2,"test");
	return 1;
}
