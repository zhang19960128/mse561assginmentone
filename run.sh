#!/bin/bash
g++ -std=c++11 -o go assignmentone.cpp potential.cpp
for i in `seq 1 100`
do
	./go $i
done
