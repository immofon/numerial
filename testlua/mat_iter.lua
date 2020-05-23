local mat = (require "num").mat

local H = mat.new(3,3):assign {
	0,1/8,-1/8,
	-1/5,0,1/10,
	1/5,1/5,0,
}

local g= mat.new(3,1):assign {1/8,2/5,-3/5}

local x =  H:solve_iter_simple(g)
assert(x == H*x + g)

x = H:solve_iter_seidel(g, { tol=1e-10 })
assert(x == H*x + g)

local A = mat.new(3,3):assign {
	8,-1,1,
	2,10,-1,
	1,1,-5,
}
local b = mat.new(3,1):assign { 1,4,3 }
x = A:solve_iter_jacobi(b)
assert(A * x == b)

x = A:solve_iter_gauss_seidel(b)
assert(A * x == b)
