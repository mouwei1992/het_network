#! /bin/bash

g++ file_util_test.cpp -std=c++11 -llapack -lopenblas -pthread -o file_test.out
