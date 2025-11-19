/**
    RSDS: Reverse-safe Text Indexing
    Copyright (C) 2020 Grigorios Loukides, Solon P. Pissis, Huiping Chen
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

#ifdef _USE_64
typedef int64_t INT;
#endif

int suffix_array_construction(unsigned char *x, INT *&SA, INT *&LCP, INT *&invSA);
bool dBgraph ( unsigned char * x, unsigned int d, INT *SA, INT *LCP, INT *invSA, double log_z, INT len_x);

void binarydB( unsigned char* x, unsigned int z);
void euler_path(const char *x, int d);



