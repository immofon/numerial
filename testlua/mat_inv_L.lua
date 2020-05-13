local num = require "num"
local mat = num.mat

local A = mat.new(3,3):assign {
	3,0,0,
	1,4,0,
	4,-7,5,
}

local invA = A:inv_L()
local I = mat.new(3,3)
for i=1,3 do
	I:set(i,i,1)
end

assert(invA * A == I)
assert(A * invA == I)
