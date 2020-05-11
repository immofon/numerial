#!/bin/bash

./compile_lua_module.sh
cd testlua
lua mat_plu.lua
lua mat_cholesky.lua
