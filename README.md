Copyright (C) 2024 Grigorios Loukides and Alessio Conte
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 
Before compiling please install the following libraries:

Eigen: https://gitlab.com/libeigen/eigen/-/releases 
Boost: https://www.boost.org/users/download/
sdsl-lite: https://github.com/simongog/sdsl-lite

Also, edit the LFLFAGS in MAkefile.64-bit.gcc

After this, the program can be compiled by 
make -f Makefile.64-bit.gcc 

We used g++ version (11.4.0) and the program ran on Ubuntu 22.04.1

INFORMATION ABOUT THE INPUT AND OUTPUT
---------------------------------------
Our approach 

  Input parameters (we refer to the parameters using the example in ./compile.txt):

    dataset.txt: This is the input string. It should be a single line of characters.

    z: This is the parameter z (privacy threshold).

    method: zrcbp or zrc or zrcb or zrce: These execute the z-RCBP, z-RC, z-RCB, and z-RCE algorithm from Bernardini et al. JEA 2021, respectively. 
            frontier or tree: AssessET with the frontier or tree-based data structure used in the algorithm from Bernardini et al. 
            frontier_compress or tree_compress: Same as frontier or tree but with chain compression. These are used in the case study.  
            frontier_fixed_d or tree_fixed_d: AssessET with the frontier or tree-based data structure 
            frontier_fixed_d_compress or tree_fixed_d_compress: Same as frontier_fixed_d or tree_fixed_d but with chain compression. These are the AF and AT algorithms.
            best_check: Assessment using the BEST theorem. 
 
    k: This is the parameter k in the paper of Bernardini et al. It should be set to 3.

    output_file: This is the output file for the alternative string that is output by the algorithm. 
    
    d: This is the parameter d. 

Example 
----------------------
./rsds dna.5MB 1000 tree_fixed_d_compress 3 output_file 32

Comments and Questions
----------------------
Grigorios Loukides
grigorios.loukides@kcl.ac.uk

Alessio Conte
alessio.conte@unipi.it
