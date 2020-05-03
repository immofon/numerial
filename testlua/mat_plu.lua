local num = require "num"
local mat = num.mat

local A = mat.new(3,3):assign {
	2,4,-2,
	1,-1,5,
	4,1,-2,
}

local P,L,U = A:plu()

local e1 = mat.new(3,1):assign{1,0,0}
local e2 = mat.new(3,1):assign{0,1,0}
local e3 = mat.new(3,1):assign{0,0,1}


function solve(b)
	local y = L:solve_L(P*b)
	return U:solve_U(y)
end

e1 = solve(e1)
e2 = solve(e2)
e3 = solve(e3)

local B = mat.new(3,3)
for i=1,3 do
	B:set(i,1,e1:get(i,1))
	B:set(i,2,e2:get(i,1))
	B:set(i,3,e3:get(i,1))
end

local I = mat.new(3,3)
for i=1,3 do
	I:set(i,i,1)
end

assert(A*B == I)
assert(B*A == I)

