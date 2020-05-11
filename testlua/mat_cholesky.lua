local num = require "num"
local mat = num.mat

local A = mat.new(3,3):assign {
	4,2,2,
	2,5,3,
	2,3,6,
}

local L,Lt = A:llt_cholesky()

assert(A == L*Lt)

local B = mat.new(6,6):assign {
	1,1,1,1,1,1,
	1,2,2,2,2,2,
	1,2,3,3,3,3,
	1,2,3,4,4,4,
	1,2,3,4,5,5,
	1,2,3,4,5,6,
}

local X = mat.new(6,1):assign {
	1,8,0,4,3,5,
}

local b = B*X

local L,Lt = B:llt_cholesky()

local y = L:solve_L(b)
local x = Lt:solve_U(y)

assert(x == X)
