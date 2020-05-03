#!/bin/bash

./compile_lua_module.sh
cd testlua
lua *.lua
