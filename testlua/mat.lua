local num = require("num")
local mat = num.mat

local A = mat.new(3,4)

A = A+A

mat.println(A,".0")

local n = 3
local B = mat.new(n,n)

mat.assign(B,{
	0,0,1,
	0,0,0,
	1,0,1,
})

mat.println(B,".0")

local Q,R = mat.qr_givens(B)
mat.println(Q,"4.1")
mat.println(R,"4.1")
