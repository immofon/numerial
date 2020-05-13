local mat = (require "num").mat

local H = mat.new(3,3):assign {
	0,1/8,-1/8,
	-1/5,0,1/10,
	1/5,1/5,0,
}

local g= mat.new(3,1):assign {1/8,2/5,-3/5}

local x =  H:solve_iter_simple(g)
assert(x == H*x + g)

x = H:solve_iter_seidel(g)
assert(x == H*x + g)
