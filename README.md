DynamicV6

/date_file

data files reside in this folder, all stored in csv format.

--/simulation
   data simulted from r scripts

--/realdata
   data from real world

/src

program file. with source file written in R or c++11

--/MCMCmodule
c++ source files. the core class is dy_het_network, which could be initialized using adjacency matrix.

--/utils
utility functions for using MCMCmodule, including those that could read csv formatted network raw data.

--/
