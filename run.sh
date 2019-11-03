#!/bin/bash

clang -o lattice lattice.c
if [ $? -eq 0 ]; then
	./lattice > run.log
fi
