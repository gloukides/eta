/**
    RSDS: Reverse Safe Data Structure
    Copyright (C) 2019 Grigorios Loukides, Solon P. Pissis

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
**/

typedef int64_t INT;


#include <vector>
#include <unordered_map>
using namespace std;

unsigned int LCParray ( unsigned char * x, INT n, INT * SA, INT * ISA, INT *LCP );
void binarydB ( unsigned char* x, unsigned int z);
int suffix_array_construction(unsigned char* x, INT *&SA, INT *&LCP, INT *&invSA, INT &r_S);
int suffix_array_construction_nomalloc(unsigned char *x, INT *SA, INT *LCP, INT *invSA, INT &r_S);
int suffix_array_construction_fixed_length(unsigned char *x, INT *&SA, INT *&LCP, INT *&invSA, INT &r_S, INT &n);
void binarydB_prefix(unsigned char* x, unsigned int z);
int binarydB_doubling_prefix(unsigned char *x, unsigned int z, unsigned int k);
int binarydB_doubling_prefix_frontier( unsigned char* seq, unsigned int z, unsigned int k);
int binarydB_doubling_prefix_tree( unsigned char* seq, unsigned int z, unsigned int k);
int binarydB_doubling_prefix_no_gab(unsigned char *x, unsigned int z, unsigned int k);

int binarydB_no_prefix(unsigned char *x, unsigned int z, unsigned int k);
int binarydB_no_prefix_gab( unsigned char* seq, unsigned int z, unsigned int k);
int binarydB_no_prefix_exp( unsigned char* seq, unsigned int z, unsigned int k);
int binarydB_doubling_prefix_frontier_compress( unsigned char* seq, unsigned int z, unsigned int k);
int binarydB_doubling_prefix_tree_compress( unsigned char* seq, unsigned int z, unsigned int k);

bool dBgraph_frontier (unsigned int d, INT *SA, INT *LCP, INT *invSA, unsigned int z, INT n);
bool dBgraph_tree (unsigned int d, INT *SA, INT *LCP, INT *invSA, unsigned int z, INT n, bool compress);
bool dBgraph_just_check(  unsigned char * x, unsigned int d, INT *SA, INT *LCP, INT *invSA, double log_z, INT n);
bool dBgraph_frontier_compress(unsigned int d, INT *SA, INT *LCP, INT *invSA, unsigned int z, INT n);

int fixed_d_frontier( unsigned char* seq, unsigned int z, unsigned int k, unsigned int fixed_d);
int fixed_d_frontier_compress( unsigned char* seq, unsigned int z, unsigned int k, unsigned int fixed_d);
int fixed_d_tree( unsigned char* seq, unsigned int z, unsigned int k, unsigned int fixed_d);
int fixed_d_tree_compress( unsigned char* seq, unsigned int z, unsigned int k, unsigned int fixed_d);

void euler_path(std::string &x, int d, string output, INT &line_cnt);
void printCircuit(std::vector<std::vector<int> > adj,  unordered_map<int,string> reverse_B_map, string outputfile, INT &line_cnt);
bool just_check_function(unsigned char* seq, unsigned int z, unsigned int d);
