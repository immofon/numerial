#!/bin/bash

./compile_lua_module.sh
cd testlua

luas="`ls *.lua`"
for file in $luas
do
	lua $file
done
