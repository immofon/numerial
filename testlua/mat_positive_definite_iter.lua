local num = require "num"
local mat = num.mat

local A = mat.new(6,6):assign {
	1.0,0.5,0.4,0.3,0.2,0.1,
	0.5,1.0,0.5,0.4,0.3,0.2,
	0.4,0.5,1.0,0.5,0.4,0.3,
	0.3,0.4,0.5,1.0,0.5,0.4,
	0.2,0.3,0.4,0.5,1.0,0.5,
	0.1,0.2,0.3,0.4,0.5,1.0,
}

local xc = mat.new(6,1):assign {1,8,0,4,3,5}
local b = A*xc

A:println "4.1"

function is_symmetric(A)
	for i=1,6 do
		for j=i+1,6 do
			if (A:get(i,j) ~= A:get(j,i)) then
				return false
			end
		end
	end
	return true
end

print("A is_symmetric",is_symmetric(A))

function solve(A,b)
	P,L,U = A:plu()
	local y = L:solve_L(P*b)
	return U:solve_U(y)
end

local x = solve(A,b)
x:println()

local x,n = A:solve_iter_steepest_descent(b, { tol=1e-6,max_step=2000 })
x:println()
print(n)
