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
	luaL_getmetatable(L,MAT_T);
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
		if (lua_geti(L,2,(i-1)*A->m + j) == LUA_TNUMBER) {
			v = (double)lua_tonumber(L,-1);
			lua_pop(L,1);
			mat_v(*A,i,j) = v;
		}else {
			luaL_argerror(L,2,"not number array");
		}
	}

	return 0;
}

static int l_mat_println(lua_State *L) {
	mat_t *A = (mat_t *) luaL_checkudata(L,1,MAT_T);
	const char *fmt = luaL_optstring(L,2,"");

	mat_println(fmt,*A);
	return 0;
}

static const struct luaL_Reg numlib[] = {
	{"new",l_mat},
	{"println",l_mat_println},
	{"assign",l_mat_assign},
	{NULL,NULL},
};

static const struct luaL_Reg mat_metareg[] ={
	{"__add",l_mat_add},
	{NULL,NULL},
};

int luaopen_num(lua_State *L) {
	luaL_newmetatable(L, MAT_T);
	luaL_setfuncs(L,mat_metareg,0);
	lua_newtable(L);
	luaL_newlib(L,numlib);
	lua_setfield(L,-2,"mat");
	return 1;
}
