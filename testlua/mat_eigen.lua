local num = require "num"
local mat = num.mat

local A = mat.new(3,3):assign {
	1,2,3,
	4,3,4,
	2,4,1,
}

local x = mat.new(3,1):assign {
	0.578157,
	1.000000,
	0.717776,
}

local Ax = A*x

Ax:println ""
x:println ""

print(Ax:get(1,1) -  8.183735 * x:get(1,1))

