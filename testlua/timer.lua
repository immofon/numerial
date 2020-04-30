local num = require "num"


local sec = num.test.timer(function ()
	local A = num.mat.new(10,10)
	for i = 1,10 do
		for j = 1,10 do
			A:set(i,j,i-j)
		end
	end

	A:println "2.0"
	local Q,R = A:qr()
	Q:println "4.1"
	R:println "4.1"
	num.mat.println(A-Q*R,"")
end)
print(sec)

local mat = num.mat
local a = {1,2,2}

mat.new(3,1):assign {1,2,3}:println ".0"

function mat.vec(t)
	return mat.new(#t,1):assign(t)
end

mat.vec {}:println ".0"


