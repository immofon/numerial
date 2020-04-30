local num = require("num")
local mat = num.mat

local A = mat.new(3,4)
--[[
A:assign {
1,2,3,4,
2,3,4,5,
3,4,5,6,
}
--]]

A:println ""
A = A+A

A:println ""

local n = 3
local B = mat.new(n,n)

B:assign {
	0,0,1,
	0,0,0,
	1,0,1,
}

B:println ".0"

local Q,R = B:qr_household()
Q:println "4.1"
R:println "4.1"

local C = mat.new(4,4)
C:assign {
	1,2,3,4,
	9,8,7,6,
	4,5,6,7,
	8,4,2,1,
}
Q,R = C:qr_household()
assert(C == Q*R,"qr")
