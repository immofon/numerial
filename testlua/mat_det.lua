local num = require "num"
local mat = num.mat

local A = mat.new(3,3):assign {
	2,0,3,
	0,-1,2,
	0,0,2,
}

assert(A:det() == -4)
