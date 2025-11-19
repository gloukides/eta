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

#define DEBUG false
#define print_string false  //if it's true it writes a string (eulerian trail) into file

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <sstream>
#include <string>
#include <sys/time.h>
#include "rsds.h"
#include "dbg.h"
#include <chrono>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;
using namespace std::chrono;

int main( int argc, char ** argv )
{

  if(argc!=7)
  {
	//cout<<"Run as ./rsds dataset.txt z zrcbp|zrc|zrcb|zrcbe|frontier k output_file fixed_d\n";
	cout<<"Run as ./rsds dataset.txt z zrcbp|zrc|zrcb|zrce|frontier_fixed_d|frontier|frontier_fixed_d_compress|frontier_compress|tree_fixed_d|tree_fixed_d_compress|tree|tree_compress|best_check k output_file d\n"; 
	exit(-1);
  }  

  std::ifstream is(argv[1]);    

  std::string str(argv[2]);
  
  std::string str3(argv[3]);

  std::string str4(argv[4]);

  
  unsigned int k;             
  std::stringstream(str4)>>k;

  string output(argv[5]);

  std::string str5(argv[6]);

  bool dbl_prefix=false,no_gab=false,no_prefix=false,no_prefix_exp=false;
  bool frontier=false, tree=false, frontier_fixed_d=false,tree_fixed_d=false;
  bool just_check=false;
  bool frontier_compress=false, frontier_fixed_d_compress=false, tree_compress=false, tree_fixed_d_compress=false;      

  if(str3.compare("zrcbp")==0)    //zrcbp  -- z-RCBP in the paper 
  {
        cout << "ZRCBP algorithm"<<endl;
        dbl_prefix=true; 
  }
  else if(str3.compare("frontier")==0)
  {
      cout<<"Frontier algorithm"<<endl;
      frontier=true;
  }
  else if(str3.compare("frontier_compress")==0)
  {
      cout<<"Frontier algorithm with compression"<<endl;
      frontier_compress=true;      
  }
  else if(str3.compare("frontier_fixed_d")==0)
  {
      cout<<"Frontier fixed d algorithm"<<endl;
      frontier_fixed_d=true;

  }
  else if(str3.compare("frontier_fixed_d_compress")==0)
  {
      cout<<"Frontier fixed d algorithm with compression"<<endl;
      frontier_fixed_d_compress=true;      
  }
  else if(str3.compare("tree")==0)
  {
     cout<<"Tree algorithm"<<endl;
     tree=true; 
  }
  else if(str3.compare("tree_compress")==0)
  {
	 cout<<"Tree algorithm with compression"<<endl;
     tree_compress=true;      
  }
  else if(str3.compare("tree_fixed_d")==0)
  {
     cout<<"Tree fixed d algorithm"<<endl;
     tree_fixed_d=true; 
  }
  else if(str3.compare("tree_fixed_d_compress")==0)
  {
     cout<<"Tree fixed d algorithm with compression"<<endl;
     tree_fixed_d_compress=true;      
  }
  else if(str3.compare("zrc")==0)  //zrc   -- z-RC in the paper
  {
        cout << "ZRC algorithm"<<endl;
        dbl_prefix=false;    //no improvement II   
	      no_prefix=true;      
        no_gab=true;         //does not use interval search
  }
  else if(str3.compare("zrce")==0) //zrce  -- z-RCE in the paper
  {
	cout << "ZRCE algorithm"<<endl;
	dbl_prefix=false;
	no_prefix=false;
	no_prefix_exp=true;   //exponential search
        no_gab=true;
  }
  else if(str3.compare("zrcb")==0)  //zrcb  -- z-RCB in the paper 
  {
        cout << "ZRCB algorithm"<<endl;
        dbl_prefix=false;
        no_prefix=true;    //does not use prefix
  }
  else if(str3.compare("best_check")==0)  // just checks using best theorem
  {
	  cout<<"Just check using BEST."<<endl;
	  just_check=true;	  
  }
  else
  {      
	cout<<"Run as ./rsds dataset.txt z zrcbp|zrc|zrcb|zrce|frontier_fixed_d|frontier|frontier_fixed_d_compress|frontier_compress|tree_fixed_d|tree_fixed_d_compress|tree|tree_compress|best_check k output_file d\n"; 
	cout<<"k is the initial prefix size (JEA optimization, default 3)"<<endl;
	 exit(-1);
  }

  unsigned int z;
	
  std::stringstream(str)>>z;
  if(z>=UINT_MAX)
  {
	cout<<"z should be smaller than: "<<UINT_MAX<<endl;
	exit(-1);  
  }
  else if(z<=1)
  {
	cout<<"z should be larger than 1."<<endl;
  }
  cout<<"z is "<<z<<endl;

     
   unsigned int d=INT_MAX;
   if(tree_fixed_d || frontier_fixed_d || frontier_fixed_d_compress || tree_fixed_d_compress || just_check)
   {
       std::stringstream(str5)>>d;
       cout<<"fixed d: "<<d<<endl;
   }
   else
   {   	  
	   cout<<"d plays no role "<<endl;
    }
   

vector<vector<char> > all_input_seqs;

  vector<char> input_seq_vec;
  char c;
  while (is.get(c))     
  {
      if(c=='\n')
      {
		break;
      }
      else
      {
        input_seq_vec.push_back(c);
      }
  }
  is.close();       

cout<<"Text length: "<<input_seq_vec.size()<<"\n"<<endl;

all_input_seqs.push_back(input_seq_vec);

  if(k>=input_seq_vec.size())
  {
        cout<<"k should be smaller than n: "<<input_seq_vec.size()<<endl;
        exit(-1);

  }
 // const char *seq = input_seq_vec.data();

   
  
  srand(time(NULL));   
  clock_t begin = clock();

   if(dbl_prefix)
  {
      INT line_cnt=0;
      for(auto &itt: all_input_seqs)
      {

         string tmpx(itt.begin(),itt.end());
         unsigned char * seq = (unsigned char*)tmpx.c_str();

	  
		      int returned_d;
	 
         high_resolution_clock::time_point t1 = high_resolution_clock::now();

         returned_d= binarydB_doubling_prefix((unsigned char*) seq,z, k);

         high_resolution_clock::time_point t2 = high_resolution_clock::now();

         cout<<"Time for assessment in ms:"<<duration_cast<std::chrono::milliseconds>(t2 - t1).count()<<endl;
         
	        cout<<"Returned: "<<returned_d<<endl;
		      if(returned_d>0)
		      {
			       string string_seq(reinterpret_cast<char*>(seq));

			       euler_path(string_seq, returned_d, output, line_cnt);
			      ++line_cnt;
		      }
       }
  }
  if(frontier)
  {
      INT line_cnt=0;
      for(auto &itt: all_input_seqs)
      {

         string tmpx(itt.begin(),itt.end());
         unsigned char * seq = (unsigned char*)tmpx.c_str();

    
          int returned_d;
          high_resolution_clock::time_point t1 = high_resolution_clock::now();   
          returned_d= binarydB_doubling_prefix_frontier((unsigned char*) seq,z, k);

         high_resolution_clock::time_point t2 = high_resolution_clock::now();

         cout<<"Time for assessment in ms:"<<duration_cast<std::chrono::milliseconds>(t2 - t1).count()<<endl;
         
          cout<<"Returned: "<<returned_d<<endl;
          if(returned_d>0 && print_string)
          {
             string string_seq(reinterpret_cast<char*>(seq));

             euler_path(string_seq, returned_d, output, line_cnt);
            ++line_cnt;
          }
       }   
  }
  if(frontier_compress)
  {
      INT line_cnt=0;
      for(auto &itt: all_input_seqs)
      {

         string tmpx(itt.begin(),itt.end());
         unsigned char * seq = (unsigned char*)tmpx.c_str();

    
          int returned_d;
            
          high_resolution_clock::time_point t1 = high_resolution_clock::now();
         
          returned_d= binarydB_doubling_prefix_frontier_compress((unsigned char*) seq,z, k);
          high_resolution_clock::time_point t2 = high_resolution_clock::now();
          cout<<"Time for assessment in ms:"<<duration_cast<std::chrono::milliseconds>(t2 - t1).count()<<endl;
          
          cout<<"Returned: "<<returned_d<<endl;
          if(returned_d>0 && print_string)
          {
             string string_seq(reinterpret_cast<char*>(seq));

             euler_path(string_seq, returned_d, output, line_cnt);
            ++line_cnt;
          }
       }   
  }      
  if(frontier_fixed_d)
  {
      INT line_cnt=0;
      for(auto &itt: all_input_seqs)
      {

         string tmpx(itt.begin(),itt.end());
         unsigned char * seq = (unsigned char*)tmpx.c_str();

    
          int returned_d;
   
          high_resolution_clock::time_point t1 = high_resolution_clock::now();  
          returned_d= fixed_d_frontier((unsigned char*) seq,z, k, d);
          high_resolution_clock::time_point t2 = high_resolution_clock::now();
          cout<<"Time for assessment in ms:"<<duration_cast<std::chrono::milliseconds>(t2 - t1).count()<<endl;

          cout<<"Returned: "<<returned_d<<endl;
          if(returned_d>0 && print_string)
          {
             string string_seq(reinterpret_cast<char*>(seq));

             euler_path(string_seq, returned_d, output, line_cnt);
            ++line_cnt;
          }
       }   
  }
  if(frontier_fixed_d_compress)
  {
      INT line_cnt=0;
      for(auto &itt: all_input_seqs)
      {

         string tmpx(itt.begin(),itt.end());
         unsigned char * seq = (unsigned char*)tmpx.c_str();

    
          int returned_d;
   
          high_resolution_clock::time_point t1 = high_resolution_clock::now();
          returned_d= fixed_d_frontier_compress((unsigned char*) seq,z, k, d);
          high_resolution_clock::time_point t2 = high_resolution_clock::now();
          cout<<"Time for assessment in ms:"<<duration_cast<std::chrono::milliseconds>(t2 - t1).count()<<endl;

          cout<<"Returned: "<<returned_d<<endl;
          if(returned_d>0 && print_string)
          {
             string string_seq(reinterpret_cast<char*>(seq));

             euler_path(string_seq, returned_d, output, line_cnt);
            ++line_cnt;
          }
       }   
  }
  if(tree)
  {
     INT line_cnt=0;
      for(auto &itt: all_input_seqs)
      {

         string tmpx(itt.begin(),itt.end());
         unsigned char * seq = (unsigned char*)tmpx.c_str();

    
          int returned_d;
         
          high_resolution_clock::time_point t1 = high_resolution_clock::now();
          returned_d= binarydB_doubling_prefix_tree((unsigned char*) seq,z, k);
          high_resolution_clock::time_point t2 = high_resolution_clock::now();
          cout<<"Time for assessment in ms:"<<duration_cast<std::chrono::milliseconds>(t2 - t1).count()<<endl;

          cout<<"Returned: "<<returned_d<<endl;
          if(returned_d>0 && print_string)
          {
             string string_seq(reinterpret_cast<char*>(seq));

             euler_path(string_seq, returned_d, output, line_cnt);
            ++line_cnt;
          }
       } 
  }
  if(tree_compress)
  {
     INT line_cnt=0;
      for(auto &itt: all_input_seqs)
      {

         string tmpx(itt.begin(),itt.end());
         unsigned char * seq = (unsigned char*)tmpx.c_str();

    
          int returned_d;
   
          high_resolution_clock::time_point t1 = high_resolution_clock::now();
          returned_d= binarydB_doubling_prefix_tree_compress((unsigned char*) seq,z, k);
          high_resolution_clock::time_point t2 = high_resolution_clock::now();
          cout<<"Time for assessment in ms:"<<duration_cast<std::chrono::milliseconds>(t2 - t1).count()<<endl;

          cout<<"Returned: "<<returned_d<<endl;
          if(returned_d>0 && print_string)
          {
             string string_seq(reinterpret_cast<char*>(seq));

             euler_path(string_seq, returned_d, output, line_cnt);
            ++line_cnt;
          }
       } 
  }
   if(tree_fixed_d)   
  {
     INT line_cnt=0;
      for(auto &itt: all_input_seqs)
      {

         string tmpx(itt.begin(),itt.end());
         unsigned char * seq = (unsigned char*)tmpx.c_str();

    
          int returned_d;
          high_resolution_clock::time_point t1 = high_resolution_clock::now();

          returned_d= fixed_d_tree((unsigned char*) seq,z, k,d );
          high_resolution_clock::time_point t2 = high_resolution_clock::now();
          cout<<"Time for assessment in ms:"<<duration_cast<std::chrono::milliseconds>(t2 - t1).count()<<endl;

          cout<<"Returned: "<<returned_d<<endl;
          if(returned_d>0 && print_string)
          {
             string string_seq(reinterpret_cast<char*>(seq));

             euler_path(string_seq, returned_d, output, line_cnt);
            ++line_cnt;
          }
       } 
  }
   if(tree_fixed_d_compress)   
  {
     INT line_cnt=0;
      for(auto &itt: all_input_seqs)
      {

         string tmpx(itt.begin(),itt.end());
         unsigned char * seq = (unsigned char*)tmpx.c_str();

    
          int returned_d;
		  cout<<"fixed_d_tree_compress"<<endl;
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
       
          returned_d= fixed_d_tree_compress((unsigned char*) seq,z, k,d );

         high_resolution_clock::time_point t2 = high_resolution_clock::now();
          cout<<"Time for assessment in ms:"<<duration_cast<std::chrono::milliseconds>(t2 - t1).count()<<endl;

          cout<<"Returned: "<<returned_d<<endl;
          if(returned_d>0 && print_string)
          {
             string string_seq(reinterpret_cast<char*>(seq));

             euler_path(string_seq, returned_d, output, line_cnt);
            ++line_cnt;
          }
       } 
  }
  if(no_prefix)
  {
	 INT line_cnt=0;
      for(auto &itt: all_input_seqs)
      {
         string tmpx(itt.begin(),itt.end());
         unsigned char * seq = (unsigned char*)tmpx.c_str();

         int returned_d;
         if(no_gab)
         {
		
	        returned_d= binarydB_no_prefix((unsigned char*) seq,z, k);
		}
	 else
	 {
			
			returned_d = binarydB_no_prefix_gab((unsigned char*) seq,z, k);
         }
        
        cout<<"Returned: "<<returned_d<<endl;
        if(returned_d>0 && print_string)
        {
            string string_seq(reinterpret_cast<char*>(seq));

            euler_path(string_seq, returned_d, output, line_cnt);
            ++line_cnt;
         }

       }
	}
   if(no_prefix_exp)
  {
         INT line_cnt=0;
      for(auto &itt: all_input_seqs)
      {
         string tmpx(itt.begin(),itt.end());
         unsigned char * seq = (unsigned char*)tmpx.c_str();
         int returned_d;

         
         returned_d= binarydB_no_prefix_exp((unsigned char*) seq,z, k);

         cout<<"Returned: "<<returned_d<<endl;
         if(returned_d>0 && print_string)
        {
            string string_seq(reinterpret_cast<char*>(seq));
         
            euler_path(string_seq, returned_d, output, line_cnt);
            ++line_cnt;
         }
       }

  }
  if(just_check)
  {
	  INT line_cnt=0;
      for(auto &itt: all_input_seqs)
      {

         string tmpx(itt.begin(),itt.end());
         unsigned char * seq = (unsigned char*)tmpx.c_str();
                   
	   high_resolution_clock::time_point t1 = high_resolution_clock::now();
       
         bool at_least_z=just_check_function((unsigned char*) seq, z, d);
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
      cout<<"Time for assessment in ms:"<<duration_cast<std::chrono::milliseconds>(t2 - t1).count()<<endl;
       

         (at_least_z==true)?cout<<"TRUE "<<line_cnt<<endl:cout<<"FALSE "<<line_cnt<<endl;        
          
	 }

  }
  clock_t end = clock();

  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  printf("Time needed %lf secs\n", elapsed_secs);


     
	return ( 0 );
}
