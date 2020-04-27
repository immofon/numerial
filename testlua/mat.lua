local num = require("num")
local mat = num.mat

local A = mat.new(3,4)

A = A+A

mat.println(A,".0")

local B = mat.new(2,2)
mat.assign(B,{
	1,2,
	3,4,
})

mat.println(B,".0")
