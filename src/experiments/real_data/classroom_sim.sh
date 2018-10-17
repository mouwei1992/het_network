#! /bin/bash

g++ classroom_sim.cpp -std=c++11 -lopenblas -llapack -pthread -o classroom_sim.out
