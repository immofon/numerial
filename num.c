#include <time.h>
#include <lua.h>
#include <lauxlib.h>
#include "numerial.h"

#define MAT_T "lua-num.mat"

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

#define DOUBLE_EQUAL_DELTA 1e-10
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
	{"qr",l_mat_qr_givens},
	{"qr_givens",l_mat_qr_givens},
	{"qr_household",l_mat_qr_household},
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
