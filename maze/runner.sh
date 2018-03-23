#!/bin/bash

for f in ./*.c;
  do 
    echo running $f
    clang -emit-llvm -g -c -I ../include/ -o ${f%.c}.bc $f 
       for ((i=1; i<=100; i++))
	do
    		../klee_build/bin/klee --search=smart -max-time=$i ${f%.c}.bc
	done;  
done;
