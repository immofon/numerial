local num = require("num")
local mat = num.mat

local Iexpect = mat.new(4,4):assign {
	1,0,0,0,
	0,1,0,0,
	0,0,1,0,
	0,0,0,1,
}

local I = mat.new(4,4)
for i=1,4 do
	I:set(i,i,1)
end

assert(I == Iexpect)
assert(I:get(1,2) == 0)
assert(I:get(3,3) == 1)
