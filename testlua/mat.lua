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

function mat.new_ident(n)
	local I = mat.new(n,n)
	for i=1,n do
		I:set(i,i,1)
	end
	return I
end

mat.new_ident(10):println ".0"

mat.new(3,3):assign {
	1,2,3,
	1,2,3,
	2,3,4,
} :println ".0"
