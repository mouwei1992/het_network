#ifndef DYNAMIC_HET_CPP
#define DYNAMIC_HET_CPP
// this should actually be the head file
// the only file you need to include

// std libraries
#include<cmath>
#include<pthread.h>
#include<cstdlib>
#include<string>
#include<iostream>
#include<regex>

// external dependencies
#include "../dependencies/armadillo"
#include "../dependencies/rapidjson/filereadstream.h"
#include "../dependencies/rapidjson/document.h"

using namespace std;
using namespace rapidjson;
using namespace arma;

// declare all types needed
#include "dynamic_het_types.h"
// declare the class
#include "dynamic_het.h"
// including utility functions
#include "../utils/utils.h"
// including generating positions utility
#include "../sim_generate/sim_generate.h"
// including class methods
#include "dynamic_het_init.cpp"
#include "dynamic_het_utils.cpp"
#include "dynamic_het_mcmc.cpp"


#endif
