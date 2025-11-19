/*
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
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions.hpp>
using namespace boost::math;


#define DEBUG false
#define PRUNE_DET true
#define GAB_OPT_IN_DB false
#define GAB_OPT_IN_NP true

uint64_t timeMs() {
	using namespace std::chrono;
	return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

// #define RANDOM_NODE_EXTRACTION // this does seem to reduce the number of steps


#define NIL -1

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <sys/time.h>
#include <iostream>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <stack>
#include <fstream>
#include <algorithm>
#include "rsds.h"
#include "dbg.h"

#ifdef _USE_64
#include <divsufsort64.h>
#endif

#include <sdsl/bit_vectors.hpp>

using namespace sdsl;
using namespace std;

#include <Eigen/Sparse>
typedef Eigen::Triplet<double> T;
typedef Eigen::SparseMatrix<double> SparseMatrixType;

std::pair<double,double> getMeanVariance(const std::vector<int>& vec) {
    double mean = 0, M2 = 0, variance = 0;

    size_t n = vec.size();
    for(size_t i = 0; i < n; ++i) {
        double delta = vec[i] - mean;
        mean += delta / (i + 1);
        M2 += delta * (vec[i] - mean);
        variance = M2 / (i + 1);
        if (i >= 2) {
            // <-- You can use the running mean and variance here 
        }
    }

    return std::make_pair(mean, variance);
}
struct scc
{
	set<int> nodes;        //vector of an SCC of the graph. Each SCC is represented by a set of nodes (ordered)
	vector<pair<int, int> > edges; //and a vector of edges
	int source;
	int target;
	bool is_path;          //is the SCC comprised of a single node? then this is true
	int lb;             //the lb of the scc, default 1

	/*scc(set<int> nodes_arg, vector<pair<int,int> > edges_arg, int source_arg, int target_arg, bool is_path_arg, int lb_arg)
	{
		for(auto& it : nodes_arg)
			nodes.insert(it); //nodes_arg.begin(),nodes_arg.end());

		for(auto &it : edges_arg)
		{
			edges.push_back(make_pair(it.first, it.second));
		}

		source=source_arg;
		target=target_arg;
		is_path=is_path_arg;
		lb=lb_arg;
	}*/

	scc(set<int>& nodes_arg, vector<pair<int, int> >& edges_arg, int source_arg, int target_arg, bool is_path_arg, int lb_arg)
	{


		for (set<int>::iterator it = nodes_arg.begin(); it != nodes_arg.end(); ++it)
			nodes.insert(*it); //nodes_arg.begin(),nodes_arg.end())

		for (int i = 0; i < edges_arg.size(); ++i)
		{
			edges.push_back(make_pair(edges_arg[i].first, edges_arg[i].second));
		}

		source = source_arg;
		target = target_arg;
		is_path = is_path_arg;
		lb = lb_arg;
	}
};

void print_single_scc(scc &s)
{
	cout << "       ----- scc ---- " << endl;

	cout << "Src: " << s.source << endl;
	cout << "Target: " << s.target << endl;
	cout << "IsPath: " << s.is_path << endl;
	cout << "lb: " << s.lb << endl;
	cout << "Nodes\n";
	set<int> x = s.nodes;

	for (auto & it2 : x)
	{
		cout << it2 << " ";
	}
	cout << endl;
	cout << "Edges=" << s.edges.size() << " :\n";
	//vector<pair<int,int> >y=scc[i].edges
	for (auto &it2 : s.edges)
	{
		cout << it2.first << "->" << it2.second << endl;
	}
	cout << "--------------------\n";
}


enum xtype {SUM, PROD, LEAF, CLOSED}; // SUM=0, PROD=1, LEAF=2, CLOSED=3
// NOTE: to check if root, check parent==NULL;

struct xtnode {
	xtype type;
	unsigned id = 0;
	unsigned lb = 1;

	vector<xtnode*> children;
	xtnode * parent = NULL;

	//STACK allocation : xtnode my_node(param);
	//HEAP allocation:   xtnode * x = new xtnode(param);

	void printme() {
		cout << "ID: " << id << " TYPE: " << type << " LB: " << lb << endl;
	}

	xtnode() {}
	xtnode(xtype mytype, unsigned myid) {
		type = mytype;
		id = myid;
		lb = 1;
		parent = NULL;
		children.clear();
	}
	xtnode(xtype mytype, unsigned myid, unsigned mylb, xtnode * myparent) {
		type = mytype;
		id = myid;
		lb = mylb;
		parent = myparent;
		children.clear();
		// children not initialized since typically nodes are created as LEAF
	}
};

void delete_tree(xtnode * root) {

	for (xtnode * child : root->children) delete_tree(child);
	delete root;
	return;
}



struct scc_vec
{
	xtnode * xn; // node in the XTree (expression-tree)
	unsigned id;
	unordered_set<int> nodes;        //vector of an SCC of the graph. Each SCC is represented by a set of nodes (ordered)
	vector<pair<int, int> > edges; //and a vector of edges
	int source;
	int target;
	bool is_path;          //is the SCC comprised of a single node? then this is true
	unsigned lb;             //the lb of the scc, default 1

	~scc_vec() {
		// nodes.clear();
		// edges.clear();
		xn = NULL;
	}
	scc_vec(const unordered_set<int>& nodes_arg, const vector<pair<int, int> >& edges_arg, int source_arg, int target_arg, bool is_path_arg, int lb_arg)
	{
		for (auto& it : nodes_arg)
			nodes.insert(it); //nodes_arg.begin(),nodes_arg.end());

		for (auto& it : edges_arg)
		{
			edges.push_back(make_pair(it.first, it.second));
		}

		source = source_arg;
		target = target_arg;
		is_path = is_path_arg;
		lb = lb_arg;
	}

	scc_vec* clone() {
		scc_vec* x = new scc_vec(this->nodes, this->edges, this->source, this->target, this->is_path, this->lb);
		x->xn = this->xn;
		x->id = this->id;

		return x;
	}

	scc_vec() {id = 0; xn = NULL;}
};

struct stack_contents
{
	int u;
	int last_i_done;
	//int *disc;
	//int *low;
	//stack<int> *st;
	//bool *stackMember;
	//int scc_id;

	//stack_contents(int u_arg, int disc_arg[], int low_arg[], stack<int> *st_arg, bool stackMember_arg[], int& scc_id_arg)
	stack_contents(int u_arg)//, int disc_arg[], int low_arg[])
	{
		u = u_arg;
		last_i_done = -1;
		//disc=disc_arg;
		//low=low_arg;

		//st=st_arg;
		//stackMember=stackMember_arg;
		//scc_id=scc_id_arg;
	}
};


class Graph_vec
{
public:
	int source;
	int target;
	int V; // No. of vertices
	vector<int> *adj; // A dynamic array of adjacency lists
	vector<int> *adjR; // Same as adj but for reverse edges

	int *node_scc_assoc; //in position i, it has the id of the scc node for node i of the graph. Default value is -1



//	void SCCUtil(int u, int disc[], int low[], stack<int> *st, bool stackMember[],int &);
	void SCCUtil_iter(int u, int disc[], int low[], stack<int> *st, bool stackMember[], int &);
//        void SCCUtil_map(int u, int disc[], int low[], stack<int> *st, bool stackMember[],int &,unordered_map<int,int>& mappings);
	void SCCUtil_map_iter(int u, int disc[], int low[], stack<int> *st, bool stackMember[], int &, unordered_map<int, int>& mappings);

	Graph_vec(int V, int s, int t); // Constructor
	~Graph_vec();
	void addEdge(int v, int w); // function to add an edge to graph
	//void SCC(); // prints strongly connected components
	void SCC_iter(); // prints strongly connected components
	//void SCC_map(unordered_map<int,int>& mappings); // prints strongly connected components
	void SCC_map_iter(unordered_map<int, int>& mappings); // prints strongly connected components
	//void printGraph();
	bool isEulerian();
	bool getV();
	int getSource();
	int getTarget();
	void setSource(int s) {source = s;}
	void setTarget(int t) {target = t;}

	//void removeMultiEdges(Graph& g);
	void copy(Graph_vec& g);
	void quasiChainCompress();
	void quasiMultiChainCompress();

	//void printSccs();
	vector<scc_vec *> sccs;

};
//returns true if the vec (assumed sorted) contains only 1 distinct element
//if length > 1, since it is ASSUMED SORTED
//we just need to see if the first element is the same as the last
bool underlyingDegOne(vector<int> & v){
	if(v.size() == 0) return false;
	if(v.size() == 1 || v[0] == v[v.size()-1]) return true;
	else return false;
}
//returns true if the vec (assumed sorted) contains only 1 distinct element
//if length > 1, since it is ASSUMED SORTED
//we just need to see if the first element is the same as the last
bool underlyingDegOneList(list<int> & v){
	if(v.size() == 0) return false;
	if(v.size() == 1 || v.front() == v.back()) return true;
	else return false;
}

// principle: ->a->b-> => delete b   
void Graph_vec::quasiMultiChainCompress(){
	// removes chains longer than 1 node (preserves node-distinct trails)
	// also compresses multichains (chains with multiplicity larger than 1) longer than 1, same principle
	// node_scc_assoc and sccs will be NOT VALID after this and should be recomputed if necessary	

	unordered_set<int> to_rem; 
	int b;

	//to check multichains easily we need to sort the adjacency lists
	for(int i = 0; i < V; i++){
		sort(adj[i].begin(), adj[i].end());
		sort(adjR[i].begin(), adjR[i].end());
	}


	for(int i = 0; i < V; i++){
		// if(i == this->source || i == this->target) continue; // ACTUALLY THIS IS FINE, we do not compress i

	// ->a->b-> : when looking at a we remove b (both nodes with indeg and outdeg 1). 
	// this eventually finds all except the head of the chain.

	// we also check for multichains 
		if(underlyingDegOne(adj[i]) && underlyingDegOne(adjR[i])){ // a has 1 distinct in-neigh and 1 distinct out-neigh
			b=adj[i][0];
			if(b == this->source || b == this->target) continue; // we cannot remove source or target
			if(underlyingDegOne(adj[b]) && underlyingDegOne(adjR[b])) to_rem.insert(b);
		}
	}

	if(to_rem.size()>0){ // need to remove things

	// remember r has indeg and outdeg 1, just need to re-route edges a->r->c => a->c (multichains, so may be several edges)
		for(int r : to_rem){ 
			int a = adjR[r][0];
			int c = adj[r][0];

			for(int i = 0 ; i < adj[a].size(); i++){
				if(adj[a][i] == r){
					adj[a][i] = c;
				}
			}
			for(int i = 0 ; i < adjR[c].size(); i++){
				if(adjR[c][i] == r){
					adjR[c][i] = a;
				}
			}
		}


		cout << "V="<< V <<", removing " << to_rem.size() << " nodes (s: " << source << ", t:" << target << ") [multichains included]\n";
		//now we removed chains, we have holes and we need to rename to fill
		// we fill hole x by renaming n to x, and decreasing V by 1 (also swap adj lists and free useless)

		if(to_rem.find(source) != to_rem.end() || to_rem.find(target) != to_rem.end() ){
			cout << "ERROR!!! removing source or target\n";
		}

		vector<int> to_rem_v;

		for(int r : to_rem){
			to_rem_v.push_back(r);
		}

		// printvec(to_rem_v);
		sort(to_rem_v.begin(),to_rem_v.end());
		// printvec(to_rem_v);

		// Largest node Usable
		int lnu = V-1; // usable if it is not in to_rem. if the hole we are fixing is larger than this we are done.

		for(int hi = 0; hi < to_rem_v.size() ; hi++){ // fill holes from smallest to largest. each hole once filled is fixed.
			while(to_rem.find(lnu) != to_rem.end() && lnu > 0) lnu--;

			//now use lnu to fill the hole to_rem_v[i] (i.e., rename lnu -> hole)
			int hole = to_rem_v[hi];

			if(hole >= lnu) break;

			if(lnu == source){ source = hole; cout << "SOURCE " << lnu << " -> " << hole << endl;}
			if(lnu == target){ target = hole; cout << "TARGET " << lnu << " -> " << hole << endl;}


			for(int i = 0; i < adj[lnu].size(); i++){ //for every outneigh on of lnu I need to check its inneighs
				int on = adj[lnu][i];
				// we will handle self-loops separately
				if(on != lnu){
					for(int j = 0; j < adjR[on].size(); j++){
						if(adjR[on][j] == lnu){
							adjR[on][j] = hole; // cannot break, there may be multiedges
						} 
					}
				}
				else 
					adj[lnu][i] = hole; // if self-loop, just update adj
			}

			for(int i = 0; i < adjR[lnu].size(); i++){ //for every inneigh in of lnu I need to check its outneighs
				int in = adjR[lnu][i];
				if(in != lnu) { // we will handle self-loops later
					for(int j = 0; j < adj[in].size(); j++){
						if(adj[in][j] == lnu){
							adj[in][j] = hole; // cannot break, there may be multiedges
						} 
					}
				}
				else 
					adjR[lnu][i] = hole; // if self-loop, just update adjR
			}

			//finally swap lists and delete adjs of removed nodes
			// delete adj[r];
			// delete adjR[r];
			adj[hole] = adj[lnu];
			adjR[hole] = adjR[lnu];

			lnu--;
		}

		V -= to_rem.size();
	}

}

void Graph_vec::quasiChainCompress(){
	// removes chains longer than 1 node (preserves node-distinct trails)
	// node_scc_assoc and sccs will be NOT VALID after this and should be recomputed if necessary	

	// int source;
	// int target;
	// int V; // No. of vertices 
	// vector<int> *adj; // A dynamic array of adjacency lists 
	// vector<int> *adjR; // Same as adj but for reverse edges 

	unordered_set<int> to_rem; 
	int b;
	// cout << "V=" << V << " adj->size()=" << adj->size() << endl;
		

	// if(V != adj->size()){
	// 	// cout << "ERROR : V=" << V << " adj->size()=" << adj->size() << endl;
	// 	exit(0);
	// }

	for(int i = 0; i < V; i++){
		// if(i == this->source || i == this->target) continue; // ACTUALLY THIS IS FINE, we do not compress i

	// ->a->b-> : when looking at a we remove b (both nodes with indeg and outdeg 1). 
	// this eventually finds all except the head of the chain.
		if(adj[i].size() == 1 && adjR[i].size() == 1){ // a has indeg and outdeg 1
			b=adj[i][0];
			if(b == this->source || b == this->target) continue; // we cannot remove source or target
			if(adj[b].size() == 1 && adjR[b].size() == 1) to_rem.insert(b);
		}
	}

	if(to_rem.size()>0){ // need to remove things

	// remember r has indeg and outdeg 1, just need to re-route 1 edge a->r->c => a->c (a has outdeg 1, c indeg unlimited)
		for(int r : to_rem){ 
			int a = adjR[r][0];
			int c = adj[r][0];

			adj[a][0] = c;
			for(int i = 0 ; i < adjR[c].size(); i++){
				if(adjR[c][i] == r){
					adjR[c][i] = a;
					break;
				}
			}
		}


		cout << "V="<< V <<", removing " << to_rem.size() << " nodes (s: " << source << ", t:" << target << ")\n";
		//now we removed chains, we have holes and we need to rename to fill
		// we fill hole x by renaming n to x, and decreasing V by 1 (also swap adj lists and free useless)

		if(to_rem.find(source) != to_rem.end() || to_rem.find(target) != to_rem.end() ){
			cout << "ERROR!!! removing source or target\n";
		}

		vector<int> to_rem_v;

		for(int r : to_rem){
			to_rem_v.push_back(r);
		}

		// printvec(to_rem_v);
		sort(to_rem_v.begin(),to_rem_v.end());
		// printvec(to_rem_v);

		// Largest node Usable
		int lnu = V-1; // usable if it is not in to_rem. if the hole we are fixing is larger than this we are done.

		for(int hi = 0; hi < to_rem_v.size() ; hi++){ // fill holes from smallest to largest. each hole once filled is fixed.
			while(to_rem.find(lnu) != to_rem.end() && lnu > 0) lnu--;

			//now use lnu to fill the hole to_rem_v[i] (i.e., rename lnu -> hole)
			int hole = to_rem_v[hi];

			if(hole >= lnu) break;

			if(lnu == source){ source = hole; cout << "SOURCE " << lnu << " -> " << hole << endl;}
			if(lnu == target){ target = hole; cout << "TARGET " << lnu << " -> " << hole << endl;}


			for(int i = 0; i < adj[lnu].size(); i++){ //for every outneigh on of lnu I need to check its inneighs
				int on = adj[lnu][i];
				for(int j = 0; j < adjR[on].size(); j++){
					if(adjR[on][j] == lnu){
						adjR[on][j] = hole; // cannot break, there may be multiedges
					} 
				}
			}

			for(int i = 0; i < adjR[lnu].size(); i++){ //for every inneigh in of lnu I need to check its outneighs
				int in = adjR[lnu][i];
				for(int j = 0; j < adj[in].size(); j++){
					if(adj[in][j] == lnu){
						adj[in][j] = hole; // cannot break, there may be multiedges
					} 
				}
			}

			//finally swap lists and delete adjs of removed nodes
			// delete adj[r];
			// delete adjR[r];
			adj[hole] = adj[lnu];
			adjR[hole] = adjR[lnu];

			lnu--;
		}

		V -= to_rem.size();
	}

}


void Graph_vec::copy(Graph_vec& g)
{
	this->source = g.source;
	this->target = g.target;
	this->V = g.V;

	vector<int> *a(g.adj);
	vector<int> *b(g.adjR);
	this->adj = a;
	this->adjR = b;
}
bool Graph_vec::getV()
{
	return this->V;
}
int Graph_vec::getSource()
{
	return this->source;
}
int Graph_vec::getTarget()
{
	return this->target;
}
Graph_vec::Graph_vec(int V, int s, int t)
{
	this->V = V;
	adj = new vector<int>[V];
	adjR = new vector<int>[V];

	this->source = s;
	this->target = t;

	node_scc_assoc = new int[V];
	for (int i = 0; i < V; ++i)
		node_scc_assoc[i] = -1;
}
Graph_vec::~Graph_vec()
{
	delete[] adj; //greg
	delete[] adjR; //greg
	delete[] node_scc_assoc;
	// sccs.clear();
}
void Graph_vec::addEdge(int v, int w)
{
	adj[v].push_back(w);
	adjR[w].push_back(v);
}
/*
void Graph_vec::SCCUtil_iter(int u2, int disc2[], int low2[], stack<int> *st2,
                             bool stackMember2[], int& scc_id2)
{


	stack<stack_contents > iteration_stack;
	stack<int> *st = st2;

	static int time = 0;
	int u = u2;
	int *disc = disc2;
	int *low = low2;
	bool *stackMember = stackMember2;
	//int scc_id=scc_id2;

	//stack_contents stc(u,disc,low, st,stackMember,scc_id2);
	stack_contents stc(u); //,disc,low);

	disc[u] = low[u] = ++time;
	st->push(u);
	stackMember[u] = true;


	iteration_stack.push(stc);

	while (!iteration_stack.empty())
	{


		stack_contents stkTop = iteration_stack.top();
		u = stkTop.u;
		//disc=stkTop.disc;
		//low=stkTop.low;
		//st=stkTop.st;
		//stackMember=stkTop.stackMember;
		//scc_id2=stkTop.scc_id;
		iteration_stack.pop();

		if (disc[u] == -1)
		{
			// Initialize discovery time and low value
			disc[u] = low[u] = ++time;
			st->push(u);
			stackMember[u] = true;
		}

		if (DEBUG)cout << u << " " << scc_id2 << endl;

		bool new_push = false;
		// Go through all vertices adjacent to this
		// vector<int>::iterator i;
		// for (i = adj[u].begin(); i != adj[u].end(); ++i)

		int i = 0;
		if (stkTop.last_i_done != -1) i = stkTop.last_i_done + 1;

		for ( ; i < adj[u].size(); i++)
		{
			int v = adj[u][i]; // v is current adjacent of 'u'

			// If v is not visited yet, then recur for it
			if (disc[v] == -1)
			{
				//SCCUtil(v, disc, low, st, stackMember, scc_id);

				//stack_contents stc_cur(u,disc,low,st,stackMember,scc_id2);
				stack_contents stc_cur(u); //,disc,low);
				stc_cur.last_i_done = i;
				iteration_stack.push(stc_cur);

				//stack_contents stc_new(v,disc,low,st,stackMember,scc_id2);
				stack_contents stc_new(v); //,disc,low);
				iteration_stack.push(stc_new);

				if (DEBUG) {
					cout << "Added: " << v << endl;
					cout << "Stack size:" << iteration_stack.size() << endl;
				}
				new_push = true;
				break;

				// Check if the subtree rooted with 'v' has a
				// connection to one of the ancestors of 'u'
				// Case 1 (per above discussion on Disc and Low value)
				//low[u] = min(low[u], low[v]);
			}
			// Update low value of 'u' only of 'v' is still in stack
			// (i.e. it's a back edge, not cross edge).
			// Case 2 (per above discussion on Disc and Low value)
			else if (stackMember[v] == true)
				low[u] = min(low[u], disc[v]);


		}

		if (DEBUG) cout << "u: " << u << " " << low[u] << " " << disc[u] << endl;

		if (new_push) continue;


		// head node found, pop the stack and print an SCC
		int w = 0; // To store stack extracted vertices
		if (low[u] == disc[u])
		{

			int sum_scc_lb = 0; //this is the sum in the lb(C_i)=sum(|N+(v)|-1)+1 that we will later increase

			unordered_set<int> nodes_of_scc;
			while (st->top() != u)
			{
				w = (int) st->top();
				//cout <<"Ins:"<< w << " scc_id:"<<scc_id<<endl;

				nodes_of_scc.insert(w);  //add node to SCC

				unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
				if (DEBUG)
					cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << endl;

				if (tmp.size() >= 3)
					sum_scc_lb += tmp.size() - 2;

				//the node w belongs to scc with id scc_id
				node_scc_assoc[w] = scc_id2;

				stackMember[w] = false;
				st->pop();
			}
			w = (int) st->top();
			nodes_of_scc.insert(w);

			unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
			if (DEBUG)
				cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << "  scc_id:" << scc_id2 << endl;

			if (tmp.size() >= 3)
				sum_scc_lb += tmp.size() - 2;

			sum_scc_lb += 1;            //now this becomes the final lb of the component

			//the node w belongs to scc with id scc_id
			node_scc_assoc[w] = scc_id2;


			//empty edge vector just for initialization
			vector<pair<int, int> > edges_of_scc;

			//source=target=-1 for initialization

			bool path_scc = false;
			//if size == 1 it is a path SCC (isPath=true)
			if (nodes_of_scc.size() == 1)
				path_scc = true;

			//the last argument of the CC is the lb of the component, i.e., sum_scc_lb

			scc_vec* scc_obj = new scc_vec(nodes_of_scc, edges_of_scc, -1, -1, path_scc, sum_scc_lb);


			this->sccs.push_back(scc_obj);
			scc_id2++;

			stackMember[w] = false;
			st->pop();
		}

		if (!iteration_stack.empty())
		{
			int v = u;
			u = (iteration_stack.top()).u;
			low[u] = min(low[u], low[v]);
		}

	}
}
*/
void Graph_vec::SCCUtil_iter(int u2, int disc2[], int low2[], stack<int> *st2,
                             bool stackMember2[], int& scc_id2)
{

	stack<stack_contents > iteration_stack;
	stack<int> *st = st2;

	int time = 0;
	int u = u2;
	int *disc = disc2;
	int *low = low2;
	bool *stackMember = stackMember2;
	//int scc_id=scc_id2;

	//stack_contents stc(u,disc,low, st,stackMember,scc_id2);
	stack_contents stc(u); //,disc,low);

	disc[u] = low[u] = ++time;		
	st->push(u);
	
	stackMember[u] = true;
        iteration_stack.push(stc);
	
	//FOR DEBUGGING -- REMOVE AFTERWARDS 
	for(int j=0; j < V; ++j)
	{
		for (int i=0 ; i < adj[j].size(); i++)
		{
			int v = adj[j][i]; // v is current adjacent of 'u'
			if(v>=V)
			{
				cout<<"Node id: "<<v<<" and all nodes in graph are: "<<V<<endl;
				exit(-1);
			}
		}
	}
	
			
	while (!iteration_stack.empty())
	{
		stack_contents stkTop = iteration_stack.top();
		u = stkTop.u;
		
		//if(u==0)cout<<"\n ***************** 0 ************************************* "<<stackMember[u]<<"\n";
		/*if (disc[u] == -1)
		{
			// Initialize discovery time and low value
			disc[u] = low[u] = ++time;
			st->push(u);
			stackMember[u] = true;
		}
		*/
		if (DEBUG)cout << u << " " << scc_id2 << endl;

		// Go through all vertices adjacent to this
		// vector<int>::iterator i;
		// for (i = adj[u].begin(); i != adj[u].end(); ++i)

		int i = 0;
		if (stkTop.last_i_done != -1) i = stkTop.last_i_done + 1;
		
		//cout<<"["<<i<<"] ";

		
		for ( ; i < adj[u].size(); i++)
		{
			int v = adj[u][i]; // v is current adjacent of 'u'

			// If v is not visited yet, then recur for it
			if (disc[v] == -1)
			{	
				disc[v] = low[v] = ++time;
		                st->push(v);
	                        stackMember[v] = true;

				stack_contents stc_new(v); //,disc,low);

				stkTop.last_i_done=i;
				
				iteration_stack.push(stc_new);				
								
				
				break;		
				
			}
			// Update low value of 'u' only of 'v' is still in stack
			// (i.e. it's a back edge, not cross edge).
			// Case 2 (per above discussion on Disc and Low value)
			else if (stackMember[v] == true) //discovered but not in an SCC
			{
			    low[u] = min(low[u], disc[v]);
  			    stkTop.last_i_done=i;
			    
			}
			//else{ //disc[v]!=-1 and stacMember[v]=false i.e., already in an SCC;}
		}
		// cout << "u: " << u << " " << low[u] << " " << disc[u] << endl;

	    if(i==adj[u].size()) // if we have explored all children of u
	   {
		//cout<<"is_size:"<<iteration_stack.size()<<endl;
		//cout<<"Before deletion: "<<(iteration_stack.top()).u<<endl;

		int u2=(iteration_stack.top()).u; //not sure about that 
		if(u2!=u)
		{
			cout<<"u2 ! = u"<<endl;
			exit(-1);
		}
		iteration_stack.pop();
		
		if (!iteration_stack.empty())
		{
			
			int v = (iteration_stack.top()).u;  //top node 
			low[v] = min(low[v], low[u]);
			
		}

		

		//cout<<"cur: "<<u<<" low[u]:"<<low[u]<<" disc[u]:"<<disc[u]<<endl;

		// head node found, pop the stack and print an SCC
		int w = 0; // To store stack extracted vertices
		if (low[u] == disc[u])
		{

			int sum_scc_lb = 0; //this is the sum in the lb(C_i)=sum(|N+(v)|-1)+1 that we will later increase

			//cout<<"st.size():"<<(*st).size()<<endl;

			unordered_set<int> nodes_of_scc;
			/*while (st->top() != u)
			{
				w = (int) st->top();
				//cout <<"Ins:"<< w << " scc_id:"<<scc_id2<<endl;

				nodes_of_scc.insert(w);  //add node to SCC

				unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
				if (DEBUG)
					cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << endl;

				if (tmp.size() >= 3)
					sum_scc_lb += tmp.size() - 2;
				
				//the node w belongs to scc with id scc_id
				node_scc_assoc[w] = scc_id2;
				cout<<"node_scc_assoc["<<w<<"]="<<scc_id2<<endl;

				stackMember[w] = false;				
				st->pop();
				//cout<<"stsz:"<<st->size()<<" component_size:"<<nodes_of_scc.size()<<endl;
			}*/
			while((int)st->top()!=u)
			{
				
				w = (int) st->top();
				st->pop();
				//cout <<"Ins:"<< w << " scc_id:"<<scc_id2<<endl;

				nodes_of_scc.insert(w);  //add node to SCC

				unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
				if (DEBUG)
					cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << endl;

				if (tmp.size() >= 3)
					sum_scc_lb += tmp.size() - 2;
				
				//the node w belongs to scc with id scc_id
				node_scc_assoc[w] = scc_id2;
				

				stackMember[w] = false;				
			
				//cout<<"stsz:"<<st->size()<<" component_size:"<<nodes_of_scc.size()<<endl;
			}
			w = (int) st->top();
			nodes_of_scc.insert(w);

			unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
			if(DEBUG)
				cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << "  scc_id:" << scc_id2 << endl;

			if (tmp.size() >= 3)
				sum_scc_lb += tmp.size() - 2;

			sum_scc_lb += 1;            //now this becomes the final lb of the component

			//the node w belongs to scc with id scc_id
			node_scc_assoc[w] = scc_id2;


			//empty edge vector just for initialization
			vector<pair<int, int> > edges_of_scc;

			//source=target=-1 for initialization

			bool path_scc = false;
			//if size == 1 it is a path SCC (isPath=true)
			if (nodes_of_scc.size() == 1)
				path_scc = true;

			//the last argument of the CC is the lb of the component, i.e., sum_scc_lb

			scc_vec* scc_obj = new scc_vec(nodes_of_scc, edges_of_scc, -1, -1, path_scc, sum_scc_lb);


			this->sccs.push_back(scc_obj);
			scc_id2++;

			stackMember[w] = false;
			st->pop();
		}//end if 

              }//end else 
	   }
	}

void Graph_vec::SCC_iter()
{
	int *disc = new int[V];
	int *low = new int[V];
	bool *stackMember = new bool[V];
	stack<int> *st = new stack<int>();

	int scc_id = 0;

	// Initialize disc and low, and stackMember arrays
	for (int i = 0; i < V; i++)
	{
		disc[i] = NIL;
		low[i] = NIL;
		stackMember[i] = false;
	}

    cout<<"Y: "<<V<<endl;
	// Call the recursive helper function to find strongly
	// connected components in DFS tree with vertex 'i'
	for (int i = 0; i < V; i++)
		if (disc[i] == NIL)
		{
			SCCUtil_iter(i, disc, low, st, stackMember, scc_id);
			if (DEBUG)
				cout << "*<" << i << " scc_id=" << scc_id << endl;
			
		}

	delete[] disc;
	delete[] low;
	delete[] stackMember;
	delete st;

	
	if(DEBUG)
	{	
		for (int i = 0; i < V; ++i)
			cout << "[ " << i << " node_scc_assoc: " << node_scc_assoc[i] << "\n";
	}
	//for each edge
	for (int i = 0; i < V; ++i)
	{
		for (vector<int>::iterator it = adj[i].begin(); it != adj[i].end(); ++it)
		{
			if (DEBUG) cout << "edge: " << i << " " << *it << endl;

			if (node_scc_assoc[i] == node_scc_assoc[*it]) //if both nodes of the edge are in the same SCC
			{
				
				this->sccs[node_scc_assoc[i]]->edges.push_back(make_pair<int, int>((int)i, (int)*it));	//add the edge to the SCC once
				if (DEBUG) cout << "scc:" << node_scc_assoc[i] << " edge:" << i << " " << *it << endl;
			}
			else	//this is a bridge edge
			{
				this->sccs[node_scc_assoc[i]]->target = i;	 //the start node of the edge is the target node of the component of the node
				this->sccs[node_scc_assoc[*it]]->source = *it; //the end node of the ege is the source node of the component of the node

				if (DEBUG)
				{
					cout << "scc1: " << node_scc_assoc[i] << " t=" << this->sccs[node_scc_assoc[i]]->target << " scc2:" << node_scc_assoc[*it] << " s=" << this->sccs[node_scc_assoc[*it]]->source << endl;
				}
			}
		}

	}
	//go to the SCC of the source of the graph and make it the source node in its component
	this->sccs[node_scc_assoc[this->source]]->source = this->source;

	//go to the SCC of the destination node of the graph and make it the target node in its component
	this->sccs[node_scc_assoc[this->target]]->target = this->target;

}

void Graph_vec::SCCUtil_map_iter(int u2, int disc2[], int low2[], stack<int> *st2,
                                 bool stackMember2[], int& scc_id2, unordered_map<int, int>& mappings)
{
	stack<stack_contents > iteration_stack;
	stack<int> *st = st2;

	int time = 0;
	int u = u2;
	int *disc = disc2;
	int *low = low2;
	bool *stackMember = stackMember2;
	//int scc_id=scc_id2;

	//stack_contents stc(u,disc,low, st,stackMember,scc_id2);
	stack_contents stc(u); //,disc,low);

	disc[u] = low[u] = ++time;
	st->push(u);
	
	stackMember[u] = true;
        iteration_stack.push(stc);
	
	while (!iteration_stack.empty())
	{
		stack_contents stkTop = iteration_stack.top();
		u = stkTop.u;
		
		//if(u==0)cout<<"\n ***************** 0 ************************************* "<<stackMember[u]<<"\n";
		/*if (disc[u] == -1)
		{
			// Initialize discovery time and low value
			disc[u] = low[u] = ++time;
			st->push(u);
			stackMember[u] = true;
		}
		*/
		if (DEBUG)cout << u << " " << scc_id2 << endl;

		// Go through all vertices adjacent to this
		// vector<int>::iterator i;
		// for (i = adj[u].begin(); i != adj[u].end(); ++i)

		int i = 0;
		if (stkTop.last_i_done != -1) i = stkTop.last_i_done + 1;
		
		//cout<<"["<<i<<"] ";

		
		for ( ; i < adj[u].size(); i++)
		{
			int v = adj[u][i]; // v is current adjacent of 'u'

			// If v is not visited yet, then recur for it
			if (disc[v] == -1)
			{	
				disc[v] = low[v] = ++time;
		                st->push(v);
	                        stackMember[v] = true;

				stack_contents stc_new(v); //,disc,low);

				stkTop.last_i_done=i;
				
				iteration_stack.push(stc_new);				
				break;		
				
			}
			// Update low value of 'u' only of 'v' is still in stack
			// (i.e. it's a back edge, not cross edge).
			// Case 2 (per above discussion on Disc and Low value)
			else if (stackMember[v] == true) //discovered but not in an SCC
			{
			    low[u] = min(low[u], disc[v]);
  			    stkTop.last_i_done=i;
			    	
			}
			//else{ //disc[v]!=-1 and stacMember[v]=false i.e., already in an SCC;}
		}
		// cout << "u: " << u << " " << low[u] << " " << disc[u] << endl;

	    if(i==adj[u].size()) // if we have explored all children of u
	   {
		//cout<<"is_size:"<<iteration_stack.size()<<endl;
		//cout<<"Before deletion: "<<(iteration_stack.top()).u<<endl;

		int u2=(iteration_stack.top()).u; //not sure about that 
		if(u2!=u)
		{
			cout<<"u2 ! = u"<<endl;
			exit(-1);
		}
		iteration_stack.pop();
		
		if (!iteration_stack.empty())
		{
			
			int v = (iteration_stack.top()).u;  //top node 
			low[v] = min(low[v], low[u]);
			//cout<<"v: "<<v<<" top: "<<u<<endl;
		}

		

		//cout<<"cur: "<<u<<" low[u]:"<<low[u]<<" disc[u]:"<<disc[u]<<endl;

		// head node found, pop the stack and print an SCC
		int w = 0; // To store stack extracted vertices
		if (low[u] == disc[u])
		{

			int sum_scc_lb = 0; //this is the sum in the lb(C_i)=sum(|N+(v)|-1)+1 that we will later increase

			//cout<<"st.size():"<<(*st).size()<<endl;

			unordered_set<int> nodes_of_scc;
			/*while (st->top() != u)
			{
				w = (int) st->top();
				//cout <<"Ins:"<< w << " scc_id:"<<scc_id2<<endl;

				nodes_of_scc.insert(w);  //add node to SCC

				unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
				if (DEBUG)
					cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << endl;

				if (tmp.size() >= 3)
					sum_scc_lb += tmp.size() - 2;
				
				//the node w belongs to scc with id scc_id
				node_scc_assoc[w] = scc_id2;
				cout<<"node_scc_assoc["<<w<<"]="<<scc_id2<<endl;

				stackMember[w] = false;				
				st->pop();
				//cout<<"stsz:"<<st->size()<<" component_size:"<<nodes_of_scc.size()<<endl;
			}*/
			while((int)st->top()!=u)
			{
				
				w = (int) st->top();
				st->pop();
				//cout <<"Ins:"<< w << " scc_id:"<<scc_id2<<endl;

				nodes_of_scc.insert(mappings[w]);  //add node to SCC

				unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
				if (DEBUG)
					cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << endl;

				if (tmp.size() >= 3)
					sum_scc_lb += tmp.size() - 2;
				
				//the node w belongs to scc with id scc_id
				node_scc_assoc[w] = scc_id2;
				//cout<<"node_scc_assoc["<<w<<"]="<<scc_id2<<endl;

				stackMember[w] = false;				
			
				//cout<<"stsz:"<<st->size()<<" component_size:"<<nodes_of_scc.size()<<endl;
			}
			w = (int) st->top();
			nodes_of_scc.insert(mappings[w]);

			unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
			if(DEBUG)
				cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << "  scc_id:" << scc_id2 << endl;

			if (tmp.size() >= 3)
				sum_scc_lb += tmp.size() - 2;

			sum_scc_lb += 1;            //now this becomes the final lb of the component

			//the node w belongs to scc with id scc_id
			node_scc_assoc[w] = scc_id2;


			//empty edge vector just for initialization
			vector<pair<int, int> > edges_of_scc;

			//source=target=-1 for initialization

			bool path_scc = false;
			//if size == 1 it is a path SCC (isPath=true)
			if (nodes_of_scc.size() == 1)
				path_scc = true;

			//the last argument of the CC is the lb of the component, i.e., sum_scc_lb

			scc_vec* scc_obj = new scc_vec(nodes_of_scc, edges_of_scc, -1, -1, path_scc, sum_scc_lb);


			this->sccs.push_back(scc_obj);
			scc_id2++;

			stackMember[w] = false;
			st->pop();
		}//end if 

              }//end else 
	   }


	
    			
}

void Graph_vec::SCC_map_iter(unordered_map<int, int>& mappings)
{
	int *disc = new int[V];
	int *low = new int[V];
	bool *stackMember = new bool[V];
	stack<int> *st = new stack<int>();

	int scc_id = 0;

	// Initialize disc and low, and stackMember arrays
	for (int i = 0; i < V; i++)
	{
		disc[i] = NIL;
		low[i] = NIL;
		stackMember[i] = false;
	}

	unordered_map<int, int> reverse_mappings;
	for (unordered_map<int, int>::iterator it = mappings.begin(); it != mappings.end(); ++it)
	{
		//the value of mappings is the key of reverse_mappings and the key of mappings is the value of reverse_mappings
		reverse_mappings.insert(make_pair<int, int>((int)it->second, (int)it->first));
	}

	// Call the recursive helper function to find strongly
	// connected components in DFS tree with vertex 'i'
	for (int i = 0; i < V; i++)
	{	if (disc[i] == NIL)
		{
			SCCUtil_map_iter(i, disc, low, st, stackMember, scc_id, reverse_mappings);
			if(node_scc_assoc[i]==-1)
			{
				cout<<"danger"<<endl;
			}

			//cout << "*<" << i << " scc_id=" << scc_id << endl;
			//if(scc_id==6143)
			//{
			// cout<< "here 1\n";
			//}
		}
		//cout<<"i: "<<i<<" scc_id'="<<scc_id<<endl;
	}



   	//debug 
   	int minus_one=0;
   	int minus_one_outdegree=0;
	int minus_one_indegree=0;
	int no_edges=0;

	for(int i=0; i<V;++i)
	{
			if(node_scc_assoc[i]==-1)
			{
					//	cout<<"node_scc_assoc["<<i<<"]="<<node_scc_assoc[i]<<endl;
				minus_one++;
				if(adj[i].size()==1)
					minus_one_outdegree++;
				if(adjR[i].size()==1)
					minus_one_indegree++;
			}
			no_edges+=adj[i].size();
	}
	if(minus_one>0)
	{	cout<<"M"<<minus_one<<"|"<<minus_one_indegree<<"|"<<minus_one_outdegree<<"] ";
		cout<<"nedges: "<<no_edges<<"\n";
	}

	//for each edge
	for (int i = 0; i < V; ++i)
	{
		//if (scc_id == 6143)
		//	cout<<"adj[i].size(): "<<adj[i].size()<<" ";

		for (vector<int>::iterator it = adj[i].begin(); it != adj[i].end(); ++it)
		{
			//debug additions
			//if(node_scc_assoc[i]==-1)continue;
			
			//debug additions
		//	if (scc_id == 6143)
		//	{	
		//		cout << "edge: " << i << " " << *it << "->" << reverse_mappings[i] << "," << reverse_mappings[*it] <<endl;
		//	}

			if(node_scc_assoc[i] == node_scc_assoc[*it]) //if both nodes of the edge are in the same SCC
			{
				this->sccs[node_scc_assoc[i]]->edges.push_back(make_pair<int, int>((int)reverse_mappings[i], (int)reverse_mappings[*it]));	//add the edge to the SCC once
				//cout<<"scc:"<<node_scc_assoc[i]<<" edge:"<<i<<" "<<*it<<endl;
			}
			else	//this is a bridge edge
			{
				this->sccs[node_scc_assoc[i]]->target = reverse_mappings[i];	 //the start node of the edge is the target node of the component of the node
				this->sccs[node_scc_assoc[*it]]->source = reverse_mappings[*it]; //the end node of the ege is the source node of the component of the node

				//cout<<"scc1: "<<node_scc_assoc[i]<<" t="<<this->sccs[node_scc_assoc[i]]->target<<" scc2:"<<node_scc_assoc[*it]<<" s="<<this->sccs[node_scc_assoc[*it]]->source<<endl;
			}
		}

	}


	//go to the SCC of the source of the graph and make it the source node in its component
	this->sccs[node_scc_assoc[mappings[this->source]]]->source = this->source;

	//go to the SCC of the destination node of the graph and make it the target node in its component
	this->sccs[node_scc_assoc[mappings[this->target]]]->target = this->target;

	delete[] disc;
	delete[] low;
	delete[] stackMember;
	delete st;

}

bool Graph_vec::isEulerian()
{
	int how_many_Delta_zero = 0;
	bool source_Delta_1 = false, target_Delta_minus_1 = false;

	for (int i = 0; i < this->V; ++i)
	{
		vector<int> l = this->adj[i];	//its size gives outdegree
		vector<int> lR = this->adjR[i]; //its size gives indegree

		//cout<<i<<" "<<lR.size()<<" "<<l.size()<<endl;

		if (l.size() - lR.size() == 0)
			how_many_Delta_zero++;

		if ((i == this->source) && (l.size() - lR.size() == 1))
		{
			source_Delta_1 = true;

		}

		if ((i == this->target) && (l.size() - lR.size() == -1))
		{
			target_Delta_minus_1 = true;

		}

	}
	//cout<<"C:"<<how_many_Delta_zero<<" V="<<V<<" s="<<this->source<<" t="<<this->target<<endl;
	if (how_many_Delta_zero == V)	//Eulerian
		return true;

	if ((how_many_Delta_zero == V - 2) && (source_Delta_1 == true) && (target_Delta_minus_1 == true))
		return true;

	return false;
}

unsigned lb_trivial(scc& s) {
	return 1;
}

unsigned lb_fun(scc& s) {
	return lb_trivial(s);
}

void printMAP(unordered_map<unsigned, scc_vec*>& ID2SCC) {
	for (auto &p : ID2SCC) {
		unsigned id = p.first;
		scc_vec * ss = p.second;

		cout << id << ": [" << ss->id;
		if (ss->xn == NULL)
			cout << " xn==NULL]" << endl;
		else
			cout << ", " << ss->xn->id << ", " << ss->xn->type << "]\n";
	}
}

// call with depth=0;
void printXTREEr(xtnode * root, int depth) {

	for (int i = 0; i < depth - 1 ; i++) cout << "|  ";
	if (depth > 0) cout << "|_ ";

	if (root->type == PROD) cout << "*" << "[" << root->id << "," << root->lb << "]\n";
	if (root->type == SUM) cout << "+" << "[" << root->id << "," << root->lb << "]\n";
	if (root->type == LEAF) cout << "L[" << root->id << "," << root->lb << "]" << "\n";

	for (xtnode * ch : root->children) printXTREEr(ch, depth + 1);
}

void printXTREE(xtnode * root) {
	cout << "======== PRINTING TREE =========\n";

	printXTREEr(root, 0);

	cout << "================================\n";
}

/* ------------------------------------------------------------------------------------------ */

class Graph
{
public:
	int source;
	int target;
	int V; // No. of vertices
	list<int> *adj; // A dynamic array of adjacency lists
	list<int> *adjR; // Same as adj but for reverse edges

	int *node_scc_assoc; //in position i, it has the id of the scc node for node i of the graph. Default value is -1


	// A Recursive DFS based function used by SCC()
	void SCCUtil(int u, int disc[], int low[], stack<int> *st, bool stackMember[], int &);
	void SCCUtil_iter(int u, int disc[], int low[], stack<int> *st, bool stackMember[], int &);
	void SCCUtil_map(int u, int disc[], int low[], stack<int> *st, bool stackMember[], int &, unordered_map<int, int>& mappings);
	void SCCUtil_map_iter(int u, int disc[], int low[], stack<int> *st, bool stackMember[], int &, unordered_map<int, int>& mappings);

	Graph(int V, int s, int t); // Constructor
	~Graph();
	void addEdge(int v, int w); // function to add an edge to graph
	void SCC(); // prints strongly connected components
	void SCC_iter(); // prints strongly connected components
	void SCC_map(unordered_map<int, int>& mappings); // prints strongly connected components
	void SCC_map_iter(unordered_map<int, int>& mappings); // prints strongly connected components
	void printGraph();
	bool isEulerian();
	bool getV();
	int getSource();
	int getTarget();
	void setSource(int s) {source = s;}
	void setTarget(int t) {target = t;}

	void quasiMultiChainCompress();
	
	void removeMultiEdges(Graph& g);
	void copy(Graph& g);
	void printSccs();
	vector<scc > sccs;

};
void Graph::quasiMultiChainCompress()
{
	// removes chains longer than 1 node (preserves node-distinct trails)
	// also compresses multichains (chains with multiplicity larger than 1) longer than 1, same principle
	// node_scc_assoc and sccs will be NOT VALID after this and should be recomputed if necessary	

	unordered_set<int> to_rem; 
	int b;

	//to check multichains easily we need to sort the adjacency lists
	for(int i = 0; i < V; i++){
		adj[i].sort();
		adjR[i].sort();
	}


	for(int i = 0; i < V; i++){
		// if(i == this->source || i == this->target) continue; // ACTUALLY THIS IS FINE, we do not compress i

	// ->a->b-> : when looking at a we remove b (both nodes with indeg and outdeg 1). 
	// this eventually finds all except the head of the chain.

	// we also check for multichains 
		if(underlyingDegOneList(adj[i]) && underlyingDegOneList(adjR[i])){ // a has 1 distinct in-neigh and 1 distinct out-neigh
			b=adj[i].front();
			if(b == this->source || b == this->target) continue; // we cannot remove source or target
			if(underlyingDegOneList(adj[b]) && underlyingDegOneList(adjR[b])) to_rem.insert(b);
		}
	}

	if(to_rem.size()>0){ // need to remove things

	// remember r has indeg and outdeg 1, just need to re-route edges a->r->c => a->c (multichains, so may be several edges)
		for(int r : to_rem)
		{ 
			int a = adjR[r].front();
			int c = adj[r].front();

			//for(int i = 0 ; i < adj[a].size(); i++){
			for(auto &it : adj[a])
			{
				if(it == r){
					it = c;
				}
			}
			for(auto &it : adjR[c])
			{
				if(it == r){
					it = a;
				}
			}
		}


		cout << "V="<< V <<", removing " << to_rem.size() << " nodes (s: " << source << ", t:" << target << ") [multichains included]\n";
		//now we removed chains, we have holes and we need to rename to fill
		// we fill hole x by renaming n to x, and decreasing V by 1 (also swap adj lists and free useless)

		if(to_rem.find(source) != to_rem.end() || to_rem.find(target) != to_rem.end() ){
			cout << "ERROR!!! removing source or target\n";
		}

		vector<int> to_rem_v;

		for(int r : to_rem)
		{
			to_rem_v.push_back(r);
		}

		// printvec(to_rem_v);
		sort(to_rem_v.begin(),to_rem_v.end());
		// printvec(to_rem_v);

		// Largest node Usable
		int lnu = V-1; // usable if it is not in to_rem. if the hole we are fixing is larger than this we are done.

		for(int hi = 0; hi < to_rem_v.size() ; hi++){ // fill holes from smallest to largest. each hole once filled is fixed.
			while(to_rem.find(lnu) != to_rem.end() && lnu > 0) lnu--;

			//now use lnu to fill the hole to_rem_v[i] (i.e., rename lnu -> hole)
			int hole = to_rem_v[hi];

			if(hole >= lnu) break;

			if(lnu == source){ source = hole; cout << "SOURCE " << lnu << " -> " << hole << endl;}
			if(lnu == target){ target = hole; cout << "TARGET " << lnu << " -> " << hole << endl;}


			//for(int i = 0; i < adj[lnu].size(); i++){ //for every outneigh on of lnu I need to check its inneighs
			for(auto &it : adj[lnu]){
				int on = it;
				// we will handle self-loops separately
				if(on != lnu){
					//for(int j = 0; j < adjR[on].size(); j++){
					for(auto &it2 : adjR[on]){
						if(it2 == lnu){
							it2 = hole; // cannot break, there may be multiedges
						} 
					}
				}
				else 
					it = hole; // if self-loop, just update adj
			}

			//for(int i = 0; i < adjR[lnu].size(); i++){ //for every inneigh in of lnu I need to check its outneighs
			for(auto &it : adjR[lnu]){
				int in = it;
				if(in != lnu) { // we will handle self-loops later
					//for(int j = 0; j < adj[in].size(); j++){
					for(auto &it2 : adj[in])
					{
						if(it2 == lnu){
							it2 = hole; // cannot break, there may be multiedges
						} 
					}
				}
				else 
					it = hole; // if self-loop, just update adjR
			}

			//finally swap lists and delete adjs of removed nodes
			// delete adj[r];
			// delete adjR[r];
			adj[hole] = adj[lnu];
			adjR[hole] = adjR[lnu];

			lnu--;
		}

		V -= to_rem.size();
	}

}



void Graph::copy(Graph& g)
{
	this->source = g.source;
	this->target = g.target;
	this->V = g.V;

	list<int> *a(g.adj);
	list<int> *b(g.adjR);
	this->adj = a;
	this->adjR = b;
}
bool Graph::getV()
{
	return this->V;
}
int Graph::getSource()
{
	return this->source;
}
int Graph::getTarget()
{
	return this->target;
}
Graph::Graph(int V, int s, int t)
{
	this->V = V;
	adj = new list<int>[V];
	adjR = new list<int>[V];

	this->source = s;
	this->target = t;

	node_scc_assoc = new int[V];
	for (int i = 0; i < V; ++i)
		node_scc_assoc[i] = -1;
}
Graph::~Graph()
{
	delete[] node_scc_assoc;
	delete[] adj;
	delete[] adjR;
}
void Graph::addEdge(int v, int w)
{
	adj[v].push_back(w);
	adjR[w].push_back(v);
}

// A recursive function that finds and prints strongly connected
// components using DFS traversal
// u --> The vertex to be visited next
// disc[] --> Stores discovery times of visited vertices
// low[] -- >> earliest visited vertex (the vertex with minimum
//			 discovery time) that can be reached from subtree
//			 rooted with current vertex
// *st -- >> To store all the connected ancestors (could be part
//		 of SCC)
// stackMember[] --> bit/index array for faster check whether
//				 a node is in stack
void Graph::SCCUtil(int u, int disc[], int low[], stack<int> *st,
                    bool stackMember[], int& scc_id)
{
	// A static variable is used for simplicity, we can avoid use
	// of static variable by passing a pointer.
	static int time = 0;

	// Initialize discovery time and low value
	disc[u] = low[u] = ++time;
	st->push(u);
	stackMember[u] = true;

	// Go through all vertices adjacent to this
	list<int>::iterator i;
	for (i = adj[u].begin(); i != adj[u].end(); ++i)
	{
		int v = *i; // v is current adjacent of 'u'

		// If v is not visited yet, then recur for it
		if (disc[v] == -1)
		{
			SCCUtil(v, disc, low, st, stackMember, scc_id);

			// Check if the subtree rooted with 'v' has a
			// connection to one of the ancestors of 'u'
			// Case 1 (per above discussion on Disc and Low value)
			low[u] = min(low[u], low[v]);
		}

		// Update low value of 'u' only of 'v' is still in stack
		// (i.e. it's a back edge, not cross edge).
		// Case 2 (per above discussion on Disc and Low value)
		else if (stackMember[v] == true)
			low[u] = min(low[u], disc[v]);
	}

	// head node found, pop the stack and print an SCC
	int w = 0; // To store stack extracted vertices
	if (low[u] == disc[u])
	{
		int sum_scc_lb = 0; //this is the sum in the lb(C_i)=sum(|N+(v)|-1)+1 that we will later increase

		set<int> nodes_of_scc;
		while (st->top() != u)
		{
			w = (int) st->top();
			//cout <<"Ins:"<< w << " scc_id:"<<scc_id<<endl;

			nodes_of_scc.insert(w);  //add node to SCC

			unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
			if (DEBUG)
				cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << endl;

			if (tmp.size() >= 3)
				sum_scc_lb += tmp.size() - 2;

			//the node w belongs to scc with id scc_id
			node_scc_assoc[w] = scc_id;

			stackMember[w] = false;
			st->pop();
		}
		w = (int) st->top();
		nodes_of_scc.insert(w);

		unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
		if (DEBUG)
			cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << "  scc_id:" << scc_id << endl;

		if (tmp.size() >= 3)
			sum_scc_lb += tmp.size() - 2;

		sum_scc_lb += 1;            //now this becomes the final lb of the component

		//the node w belongs to scc with id scc_id
		node_scc_assoc[w] = scc_id;


		//empty edge vector just for initialization
		vector<pair<int, int> > edges_of_scc;

		//source=target=-1 for initialization

		bool path_scc = false;
		//if size == 1 it is a path SCC (isPath=true)
		if (nodes_of_scc.size() == 1)
			path_scc = true;

		//the last argument of the CC is the lb of the component, i.e., sum_scc_lb

		scc scc_obj(nodes_of_scc, edges_of_scc, -1, -1, path_scc, sum_scc_lb);


		this->sccs.push_back(scc_obj);
		scc_id++;

		stackMember[w] = false;
		st->pop();
	}
	//for each SCC record its edges
	/*for(int i=0;i<this->V ; ++i)
	{
		cout<<endl;
	}*/
	//cout<<"!!\n";
}

void Graph::SCCUtil_iter(int u2, int disc2[], int low2[], stack<int> *st2,
                         bool stackMember2[], int& scc_id2)
{


	stack<stack_contents > iteration_stack;
	stack<int> *st = st2;

	static int time = 0;
	int u = u2;
	int *disc = disc2;
	int *low = low2;
	bool *stackMember = stackMember2;
	//int scc_id=scc_id2;

	//stack_contents stc(u,disc,low, st,stackMember,scc_id2);
	stack_contents stc(u); //,disc,low);

	disc[u] = low[u] = ++time;
	st->push(u);
	stackMember[u] = true;


	iteration_stack.push(stc);

	while (!iteration_stack.empty())
	{


		stack_contents stkTop = iteration_stack.top();
		u = stkTop.u;

		int i=0;		
		if (stkTop.last_i_done != -1) i = stkTop.last_i_done + 1;

		list<int>::iterator it=adj[u].begin();

		advance(it,i);  // goes to where it stopped (needed for dfs)
 
		// Go through all vertices adjacent to this
		for (; it != adj[u].end(); ++it)
		{
			int v = *it; // v is current adjacent of 'u'

			// If v is not visited yet, then recur for it
			if (disc[v] == -1)
			{
				disc[v]=low[v]=++time;
				st->push(v);
				stackMember[v]=true;
			
				stack_contents stc_new(v); 
				
				stkTop.last_i_done=i;

				iteration_stack.push(stc_new);

				break;

			}
			// Update low value of 'u' only of 'v' is still in stack
			// (i.e. it's a back edge, not cross edge).
			// Case 2 (per above discussion on Disc and Low value)
			else if (stackMember[v] == true)
			{
				low[u] = min(low[u], disc[v]);
				stkTop.last_i_done=i;
			}

			i++;
		}


	if(i==adj[u].size())
       {
		int u2=(iteration_stack.top()).u; //not sure about that 
		if(u2!=u)
		{
			cout<<"u2 ! = u"<<endl;
			exit(-1);
		}
		iteration_stack.pop();
		
		if (!iteration_stack.empty())
		{
			
			int v = (iteration_stack.top()).u;  //top node 
			low[v] = min(low[v], low[u]);			
		}

		// head node found, pop the stack and print an SCC
		int w = 0; // To store stack extracted vertices

		if (low[u] == disc[u])
		{

			int sum_scc_lb = 0; //this is the sum in the lb(C_i)=sum(|N+(v)|-1)+1 that we will later increase

			set<int> nodes_of_scc;
			while (st->top() != u)
			{
				w = (int) st->top();
				st->pop();
				

				nodes_of_scc.insert(w);  //add node to SCC

				unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
				if (DEBUG)
					cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << endl;

				if (tmp.size() >= 3)
					sum_scc_lb += tmp.size() - 2;

				//the node w belongs to scc with id scc_id
				node_scc_assoc[w] = scc_id2;

				stackMember[w] = false;
				
			}
			w = (int) st->top();
			nodes_of_scc.insert(w);

			unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
			if (DEBUG)
				cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << "  scc_id:" << scc_id2 << endl;

			if (tmp.size() >= 3)
				sum_scc_lb += tmp.size() - 2;

			sum_scc_lb += 1;            //now this becomes the final lb of the component

			//the node w belongs to scc with id scc_id
			node_scc_assoc[w] = scc_id2;


			//empty edge vector just for initialization
			vector<pair<int, int> > edges_of_scc;

			//source=target=-1 for initialization

			bool path_scc = false;
			//if size == 1 it is a path SCC (isPath=true)
			if (nodes_of_scc.size() == 1)
				path_scc = true;

			//the last argument of the CC is the lb of the component, i.e., sum_scc_lb

			scc scc_obj(nodes_of_scc, edges_of_scc, -1, -1, path_scc, sum_scc_lb);


			this->sccs.push_back(scc_obj);
			scc_id2++;

			stackMember[w] = false;
			st->pop();
		}//end if(low[u]=disc[u])

		
          }//end if(i==adj[i].size())

	}
}
void Graph::SCC_iter()
{
	int *disc = new int[V];
	int *low = new int[V];
	bool *stackMember = new bool[V];
	stack<int> *st = new stack<int>();

	int scc_id = 0;

	// Initialize disc and low, and stackMember arrays
	for (int i = 0; i < V; i++)
	{
		disc[i] = NIL;
		low[i] = NIL;
		stackMember[i] = false;
	}

	// Call the recursive helper function to find strongly
	// connected components in DFS tree with vertex 'i'
	for (int i = 0; i < V; i++)
		if (disc[i] == NIL)
		{
			SCCUtil_iter(i, disc, low, st, stackMember, scc_id);
			if (DEBUG)
				cout << "*<" << i << " scc_id=" << scc_id << endl;
			//scc_id++;
		}

	delete[] disc;
	delete[] low;
	delete[] stackMember;
	delete st;

	if (DEBUG)
	{	for (int i = 0; i < V; ++i)
			cout << "[ " << i << " node_scc_assoc: " << node_scc_assoc[i] << "\n";
	}
	//for each edge
	for (int i = 0; i < V; ++i)
	{
		for (list<int>::iterator it = adj[i].begin(); it != adj[i].end(); ++it)
		{
			if (DEBUG) cout << "edge: " << i << " " << *it << endl;

			if (node_scc_assoc[i] == node_scc_assoc[*it]) //if both nodes of the edge are in the same SCC
			{
				this->sccs[node_scc_assoc[i]].edges.push_back(make_pair<int, int>((int)i, (int)*it));	//add the edge to the SCC once
				if (DEBUG) cout << "scc:" << node_scc_assoc[i] << " edge:" << i << " " << *it << endl;
			}
			else	//this is a bridge edge
			{
				this->sccs[node_scc_assoc[i]].target = i;	 //the start node of the edge is the target node of the component of the node
				this->sccs[node_scc_assoc[*it]].source = *it; //the end node of the ege is the source node of the component of the node

				if (DEBUG)
				{
					cout << "scc1: " << node_scc_assoc[i] << " t=" << this->sccs[node_scc_assoc[i]].target << " scc2:" << node_scc_assoc[*it] << " s=" << this->sccs[node_scc_assoc[*it]].source << endl;
				}
			}
		}

	}
	//go to the SCC of the source of the graph and make it the source node in its component
	this->sccs[node_scc_assoc[this->source]].source = this->source;

	//go to the SCC of the destination node of the graph and make it the target node in its component
	this->sccs[node_scc_assoc[this->target]].target = this->target;

}

// The function to do DFS traversal. It uses SCCUtil()
void Graph::SCC()
{
	int *disc = new int[V];
	int *low = new int[V];
	bool *stackMember = new bool[V];
	stack<int> *st = new stack<int>();

	int scc_id = 0;

	// Initialize disc and low, and stackMember arrays
	for (int i = 0; i < V; i++)
	{
		disc[i] = NIL;
		low[i] = NIL;
		stackMember[i] = false;
	}

	// Call the recursive helper function to find strongly
	// connected components in DFS tree with vertex 'i'
	for (int i = 0; i < V; i++)
		if (disc[i] == NIL)
		{
			SCCUtil(i, disc, low, st, stackMember, scc_id);
			if (DEBUG)
				cout << "*<" << i << " scc_id=" << scc_id << endl;
			//scc_id++;
		}
	if (DEBUG)
	{	for (int i = 0; i < V; ++i)
			cout << "[ " << i << " node_scc_assoc: " << node_scc_assoc[i] << "\n";
	}
	//for each edge
	for (int i = 0; i < V; ++i)
	{
		for (list<int>::iterator it = adj[i].begin(); it != adj[i].end(); ++it)
		{
			if (DEBUG) cout << "edge: " << i << " " << *it << endl;

			if (node_scc_assoc[i] == node_scc_assoc[*it]) //if both nodes of the edge are in the same SCC
			{
				this->sccs[node_scc_assoc[i]].edges.push_back(make_pair<int, int>((int)i, (int)*it));	//add the edge to the SCC once
				if (DEBUG) cout << "scc:" << node_scc_assoc[i] << " edge:" << i << " " << *it << endl;
			}
			else	//this is a bridge edge
			{
				this->sccs[node_scc_assoc[i]].target = i;	 //the start node of the edge is the target node of the component of the node
				this->sccs[node_scc_assoc[*it]].source = *it; //the end node of the ege is the source node of the component of the node

				if (DEBUG)
				{
					cout << "scc1: " << node_scc_assoc[i] << " t=" << this->sccs[node_scc_assoc[i]].target << " scc2:" << node_scc_assoc[*it] << " s=" << this->sccs[node_scc_assoc[*it]].source << endl;
				}
			}
		}

	}
	//go to the SCC of the source of the graph and make it the source node in its component
	this->sccs[node_scc_assoc[this->source]].source = this->source;

	//go to the SCC of the destination node of the graph and make it the target node in its component
	this->sccs[node_scc_assoc[this->target]].target = this->target;


	delete[] disc;
	delete[] low;
	delete[] stackMember;
	delete st;

}
void Graph::SCCUtil_map(int u, int disc[], int low[], stack<int> *st,
                        bool stackMember[], int& scc_id, unordered_map<int, int>& mappings)
{
	// A static variable is used for simplicity, we can avoid use
	// of static variable by passing a pointer.
	static int time = 0;

	// Initialize discovery time and low value
	disc[u] = low[u] = ++time;
	st->push(u);
	stackMember[u] = true;

	// Go through all vertices adjacent to this
	list<int>::iterator i;
	for (i = adj[u].begin(); i != adj[u].end(); ++i)
	{
		int v = *i; // v is current adjacent of 'u'

		// If v is not visited yet, then recur for it
		if (disc[v] == -1)
		{
			SCCUtil_map(v, disc, low, st, stackMember, scc_id, mappings);

			// Check if the subtree rooted with 'v' has a
			// connection to one of the ancestors of 'u'
			// Case 1 (per above discussion on Disc and Low value)
			low[u] = min(low[u], low[v]);
		}

		// Update low value of 'u' only of 'v' is still in stack
		// (i.e. it's a back edge, not cross edge).
		// Case 2 (per above discussion on Disc and Low value)
		else if (stackMember[v] == true)
			low[u] = min(low[u], disc[v]);
	}

	// head node found, pop the stack and print an SCC
	int w = 0; // To store stack extracted vertices
	if (low[u] == disc[u])
	{
		int sum_scc_lb = 0; //this is the sum in the lb(C_i)=sum(|N+(v)|-1)+1 that we will later increase

		set<int> nodes_of_scc;
		while (st->top() != u)
		{
			w = (int) st->top();
			//cout <<"Ins:"<< w << "->"<<mappings[w]<<" scc_id:"<<scc_id<<endl;

			nodes_of_scc.insert(mappings[w]);  //add node to SCC

			unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
			//cout<<"w="<<w<<"->"<<mappings[w]<<" # unique out-neighbors: "<<tmp.size()<<endl;

			if (tmp.size() >= 3)
				sum_scc_lb += tmp.size() - 2;

			//the node w belongs to scc with id scc_id
			node_scc_assoc[w] = scc_id;

			stackMember[w] = false;
			st->pop();
		}
		w = (int) st->top();
		//cout <<"$"<< w <<"->"<<mappings[w]<< "\n";
		nodes_of_scc.insert(mappings[w]);

		unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
		//cout<<"w="<<w<<"->"<<mappings[w]<<" # unique out-neighbors: "<<tmp.size()<<"  scc_id:"<<scc_id<<endl;
		if (tmp.size() >= 3)
			sum_scc_lb += tmp.size() - 2;

		sum_scc_lb += 1;            //now this becomes the final lb of the component

		//the node w belongs to scc with id scc_id
		node_scc_assoc[w] = scc_id;


		//empty edge vector just for initialization
		vector<pair<int, int> > edges_of_scc;

		//source=target=-1 for initialization

		bool path_scc = false;
		//if size == 1 it is a path SCC (isPath=true)
		if (nodes_of_scc.size() == 1)
			path_scc = true;

		//the last argument of the CC is the lb of the component, i.e., sum_scc_lb

		scc scc_obj(nodes_of_scc, edges_of_scc, -1, -1, path_scc, sum_scc_lb);


		this->sccs.push_back(scc_obj);
		scc_id++;

		stackMember[w] = false;
		st->pop();
	}
	//for each SCC record its edges
	/*for(int i=0;i<this->V ; ++i)
	{
		cout<<endl;
	}*/
	//cout<<"!!\n";
}

void Graph::SCCUtil_map_iter(int u2, int disc2[], int low2[], stack<int> *st2,
                             bool stackMember2[], int& scc_id2, unordered_map<int, int>& mappings)
{

stack<stack_contents > iteration_stack;
	stack<int> *st = st2;

	static int time = 0;
	int u = u2;
	int *disc = disc2;
	int *low = low2;
	bool *stackMember = stackMember2;
	//int scc_id=scc_id2;

	//stack_contents stc(u,disc,low, st,stackMember,scc_id2);
	stack_contents stc(u); //,disc,low);

	disc[u] = low[u] = ++time;
	st->push(u);
	stackMember[u] = true;


	iteration_stack.push(stc);

	while (!iteration_stack.empty())
	{


		stack_contents stkTop = iteration_stack.top();
		u = stkTop.u;

		int i=0;		
		if (stkTop.last_i_done != -1) i = stkTop.last_i_done + 1;

		list<int>::iterator it=adj[u].begin();

		advance(it,i);  // goes to where it stopped (needed for dfs)
 
		// Go through all vertices adjacent to this
		for (; it != adj[u].end(); ++it)
		{
			int v = *it; // v is current adjacent of 'u'

			// If v is not visited yet, then recur for it
			if (disc[v] == -1)
			{
				disc[v]=low[v]=++time;
				st->push(v);
				stackMember[v]=true;
			
				stack_contents stc_new(v); 
				
				stkTop.last_i_done=i;

				iteration_stack.push(stc_new);

				break;

			}
			// Update low value of 'u' only of 'v' is still in stack
			// (i.e. it's a back edge, not cross edge).
			// Case 2 (per above discussion on Disc and Low value)
			else if (stackMember[v] == true)
			{
				low[u] = min(low[u], disc[v]);
				stkTop.last_i_done=i;
			}

			i++;
		}


	if(i==adj[u].size())
       {
		int u2=(iteration_stack.top()).u; //not sure about that 
		if(u2!=u)
		{
			cout<<"u2 ! = u"<<endl;
			exit(-1);
		}
		iteration_stack.pop();
		
		if (!iteration_stack.empty())
		{
			
			int v = (iteration_stack.top()).u;  //top node 
			low[v] = min(low[v], low[u]);			
		}

		// head node found, pop the stack and print an SCC
		int w = 0; // To store stack extracted vertices

		if (low[u] == disc[u])
		{

			int sum_scc_lb = 0; //this is the sum in the lb(C_i)=sum(|N+(v)|-1)+1 that we will later increase

			set<int> nodes_of_scc;
			while (st->top() != u)
			{
				w = (int) st->top();
				st->pop();
				

				nodes_of_scc.insert(mappings[w]);  //add node to SCC

				unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
				if (DEBUG)
					cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << endl;

				if (tmp.size() >= 3)
					sum_scc_lb += tmp.size() - 2;

				//the node w belongs to scc with id scc_id
				node_scc_assoc[w] = scc_id2;

				stackMember[w] = false;
				
			}
			w = (int) st->top();
			nodes_of_scc.insert(mappings[w]);

			unordered_set<int> tmp(adj[w].begin(), adj[w].end()); //contains unique out-neighbors of the node
			if (DEBUG)
				cout << "w=" << w << " # unique out-neighbors: " << tmp.size() << "  scc_id:" << scc_id2 << endl;

			if (tmp.size() >= 3)
				sum_scc_lb += tmp.size() - 2;

			sum_scc_lb += 1;            //now this becomes the final lb of the component

			//the node w belongs to scc with id scc_id
			node_scc_assoc[w] = scc_id2;


			//empty edge vector just for initialization
			vector<pair<int, int> > edges_of_scc;

			//source=target=-1 for initialization

			bool path_scc = false;
			//if size == 1 it is a path SCC (isPath=true)
			if (nodes_of_scc.size() == 1)
				path_scc = true;

			//the last argument of the CC is the lb of the component, i.e., sum_scc_lb

			scc scc_obj(nodes_of_scc, edges_of_scc, -1, -1, path_scc, sum_scc_lb);


			this->sccs.push_back(scc_obj);
			scc_id2++;

			stackMember[w] = false;
			st->pop();
		}//end if(low[u]=disc[u])

		
          }//end if(i==adj[i].size())

	}

}
// The function to do DFS traversal. It uses SCCUtil()
void Graph::SCC_map(unordered_map<int, int>& mappings)
{
	int *disc = new int[V];
	int *low = new int[V];
	bool *stackMember = new bool[V];
	stack<int> *st = new stack<int>();

	int scc_id = 0;

	// Initialize disc and low, and stackMember arrays
	for (int i = 0; i < V; i++)
	{
		disc[i] = NIL;
		low[i] = NIL;
		stackMember[i] = false;
	}

	unordered_map<int, int> reverse_mappings;
	for (unordered_map<int, int>::iterator it = mappings.begin(); it != mappings.end(); ++it)
	{
		//the value of mappings is the key of reverse_mappings and the key of mappings is the value of reverse_mappings
		reverse_mappings.insert(make_pair<int, int>((int)it->second, (int)it->first));
	}

	// Call the recursive helper function to find strongly
	// connected components in DFS tree with vertex 'i'
	for (int i = 0; i < V; i++)
		if (disc[i] == NIL)
		{
			SCCUtil_map(i, disc, low, st, stackMember, scc_id, reverse_mappings);
			//cout<<"*<"<<i<<" scc_id="<<scc_id<<endl;

		}

	//for(int i=0;i<V;++i)
	//	cout<<"[ "<<i<<"->"<<reverse_mappings[i]<<" node_scc_assoc: "<<node_scc_assoc[i]<<"\n";

	//for each edge
	for (int i = 0; i < V; ++i)
	{
		for (list<int>::iterator it = adj[i].begin(); it != adj[i].end(); ++it)
		{
			//cout<<"edge: "<<i<<" "<<*it<<"->"<<reverse_mappings[i]<<","<<reverse_mappings[*it]<<endl;

			if (node_scc_assoc[i] == node_scc_assoc[*it]) //if both nodes of the edge are in the same SCC
			{
				this->sccs[node_scc_assoc[i]].edges.push_back(make_pair<int, int>((int)reverse_mappings[i], (int)reverse_mappings[*it]));	//add the edge to the SCC once
				//cout<<"scc:"<<node_scc_assoc[i]<<" edge:"<<i<<" "<<*it<<endl;
			}
			else	//this is a bridge edge
			{
				this->sccs[node_scc_assoc[i]].target = reverse_mappings[i];	 //the start node of the edge is the target node of the component of the node
				this->sccs[node_scc_assoc[*it]].source = reverse_mappings[*it]; //the end node of the ege is the source node of the component of the node

				//cout<<"scc1: "<<node_scc_assoc[i]<<" t="<<this->sccs[node_scc_assoc[i]].target<<" scc2:"<<node_scc_assoc[*it]<<" s="<<this->sccs[node_scc_assoc[*it]].source<<endl;
			}
		}

	}
	//go to the SCC of the source of the graph and make it the source node in its component
	this->sccs[node_scc_assoc[mappings[this->source]]].source = this->source;

	//go to the SCC of the destination node of the graph and make it the target node in its component
	this->sccs[node_scc_assoc[mappings[this->target]]].target = this->target;

	delete[] disc;
	delete[] low;
	delete[] stackMember;
	delete st;

}
void Graph::SCC_map_iter(unordered_map<int, int>& mappings)
{
	int *disc = new int[V];
	int *low = new int[V];
	bool *stackMember = new bool[V];
	stack<int> *st = new stack<int>();

	int scc_id = 0;

	// Initialize disc and low, and stackMember arrays
	for (int i = 0; i < V; i++)
	{
		disc[i] = NIL;
		low[i] = NIL;
		stackMember[i] = false;
	}

	unordered_map<int, int> reverse_mappings;
	for (unordered_map<int, int>::iterator it = mappings.begin(); it != mappings.end(); ++it)
	{
		//the value of mappings is the key of reverse_mappings and the key of mappings is the value of reverse_mappings
		reverse_mappings.insert(make_pair<int, int>((int)it->second, (int)it->first));
	}

	// Call the recursive helper function to find strongly
	// connected components in DFS tree with vertex 'i'
	for (int i = 0; i < V; i++)
		if (disc[i] == NIL)
		{
			SCCUtil_map_iter(i, disc, low, st, stackMember, scc_id, reverse_mappings);
			//cout<<"*<"<<i<<" scc_id="<<scc_id<<endl;

		}

	//for(int i=0;i<V;++i)
	//	cout<<"[ "<<i<<"->"<<reverse_mappings[i]<<" node_scc_assoc: "<<node_scc_assoc[i]<<"\n";

	//for each edge
	for (int i = 0; i < V; ++i)
	{
		for (list<int>::iterator it = adj[i].begin(); it != adj[i].end(); ++it)
		{
			//cout<<"edge: "<<i<<" "<<*it<<"->"<<reverse_mappings[i]<<","<<reverse_mappings[*it]<<endl;

			if (node_scc_assoc[i] == node_scc_assoc[*it]) //if both nodes of the edge are in the same SCC
			{
				this->sccs[node_scc_assoc[i]].edges.push_back(make_pair<int, int>((int)reverse_mappings[i], (int)reverse_mappings[*it]));	//add the edge to the SCC once
				//cout<<"scc:"<<node_scc_assoc[i]<<" edge:"<<i<<" "<<*it<<endl;
			}
			else	//this is a bridge edge
			{
				this->sccs[node_scc_assoc[i]].target = reverse_mappings[i];	 //the start node of the edge is the target node of the component of the node
				this->sccs[node_scc_assoc[*it]].source = reverse_mappings[*it]; //the end node of the ege is the source node of the component of the node

				//cout<<"scc1: "<<node_scc_assoc[i]<<" t="<<this->sccs[node_scc_assoc[i]].target<<" scc2:"<<node_scc_assoc[*it]<<" s="<<this->sccs[node_scc_assoc[*it]].source<<endl;
			}
		}

	}
	//go to the SCC of the source of the graph and make it the source node in its component
	this->sccs[node_scc_assoc[mappings[this->source]]].source = this->source;

	//go to the SCC of the destination node of the graph and make it the target node in its component
	this->sccs[node_scc_assoc[mappings[this->target]]].target = this->target;

	delete[] disc;
	delete[] low;
	delete[] stackMember;
	delete st;

}

void Graph::printSccs()
{
	cout << "\n Graph Print Sccs:\n";

	for (int i = 0; i < sccs.size(); ++i)
	{
		cout << "--------------------\n";
		cout << "SCC_id: " << i << endl;
		cout << "Src: " << sccs[i].source << endl;
		cout << "Target: " << sccs[i].target << endl;
		cout << "IsPath: " << sccs[i].is_path << endl;
		cout << "lb: " << sccs[i].lb << endl;
		cout << "Nodes\n";
		set<int> x = sccs[i].nodes;

		for (auto & it2 : x)
		{
			cout << it2 << " ";
		}
		cout << endl;
		cout << "Edges=" << sccs[i].edges.size() << " :\n";
		//vector<pair<int,int> >y=scc[i].edges
		for (auto &it2 : sccs[i].edges)
		{
			cout << it2.first << "->" << it2.second << endl;
		}
		cout << "--------------------\n";
	}
}

void Graph::printGraph()
{
	cout << "\n ---------- printGraph ----------- \n";

	for (int i = 0; i < this->V; ++i)
	{
		list<int> l = this->adj[i];
		cout << i << "|";
		for (auto &it : l)
			cout << it << ", ";
		cout << endl;
	}
	/* cout<<endl<<"Reverse: "<<endl;
	for(int i=0;i<this->V;++i)
	{
		list<int> l=this->adjR[i];
		cout<<i<<"|";
		for(auto &it:l)
			cout<<it<<", ";
		cout<<endl;
	}*/
	cout << "\n ---------------------------\n";
}
bool Graph::isEulerian()
{
	int how_many_Delta_zero = 0;
	bool source_Delta_1 = false, target_Delta_minus_1 = false;

	for (int i = 0; i < this->V; ++i)
	{
		list<int> l = this->adj[i];	//its size gives outdegree
		list<int> lR = this->adjR[i]; //its size gives indegree

		//cout<<i<<" "<<lR.size()<<" "<<l.size()<<endl;

		if (l.size() - lR.size() == 0)
			how_many_Delta_zero++;

		if ((i == this->source) && (l.size() - lR.size() == 1))
		{
			source_Delta_1 = true;

		}

		if ((i == this->target) && (l.size() - lR.size() == -1))
		{
			target_Delta_minus_1 = true;

		}

	}
	//cout<<"C:"<<how_many_Delta_zero<<" V="<<V<<" s="<<this->source<<" t="<<this->target<<endl;
	if (how_many_Delta_zero == V)	//Eulerian
		return true;

	if ((how_many_Delta_zero == V - 2) && (source_Delta_1 == true) && (target_Delta_minus_1 == true))
		return true;

	return false;
}

void Graph::removeMultiEdges(Graph& g2)
{
	cout << "\n Before\n";
	g2.printGraph();
	for (int i = 0; i < this->V; ++i)
	{
		g2.adj[i].sort();
		g2.adj[i].unique();
		g2.adjR[i].sort();
		g2.adjR[i].unique();
	}
	cout << "\n After\n";
	g2.printGraph();
}

void free_arrays(INT *SA, INT *LCP, INT *invSA)
{
	free(invSA);
	free(SA);
	free(LCP);
}
bool dBgraph_just_check(  unsigned char * x, unsigned int d, INT *SA, INT *LCP, INT *invSA, double log_z, INT n)
	{

	
			if(d==0)
			return true;

		INT cluster_id=(INT) 0;
		INT *C=new INT[n];

		for (INT i=0; i<n ; i++)
		{
			if(LCP[i]>=d-1)
			{
				C[i]=(INT)cluster_id;			
				
			}
			else
			{
				cluster_id++;
				C[i]=cluster_id;
				
			}	
		
		}


		multimap<INT, INT> map_for_outdegree;


		multimap<INT, INT> map_for_indegree;


		unordered_map<INT, INT> node_id_consec; 
		INT node_id_consec_ind=0;

			unordered_set<INT> distinct_nodes;
	 
		for(INT i=0; i<=n-d ; ++i)
		{
			INT cluster_id_pref=C[invSA[i]]; 
			INT cluster_id_suff=C[invSA[i+1]]; 

			distinct_nodes.insert(cluster_id_pref);
					distinct_nodes.insert(cluster_id_suff);


			map_for_outdegree.insert(std::pair<INT,INT>(cluster_id_pref,cluster_id_suff));
			map_for_indegree.insert(std::pair<INT,INT>(cluster_id_suff,cluster_id_pref));
			
					unordered_map<INT, INT>::const_iterator nic_find=node_id_consec.find(cluster_id_pref);

			if(nic_find==node_id_consec.end())
			{	node_id_consec[cluster_id_pref]=node_id_consec_ind;
						node_id_consec_ind++; 
				}

			unordered_map<INT, INT>::const_iterator nic_find2=node_id_consec.find(cluster_id_suff);

					if(nic_find2==node_id_consec.end())
					{       node_id_consec[cluster_id_suff]=node_id_consec_ind;
							node_id_consec_ind++;
					}

		}


auto algo_start=chrono::steady_clock::now();

		 INT matrix_dim=distinct_nodes.size();

		 std::vector<T> tripletList;
	  
		INT node_t=C[invSA[n-d+1]];
	 
		std::pair <std::multimap<INT,INT>::iterator, std::multimap<INT,INT>::iterator> ret_outdegree;
		std::unordered_map<INT,INT> node_outdegree;
	   
		unordered_set<INT> used_keys; 

		unordered_map<INT,INT> a_uu;

	   
		
		double denominator_log_sum=0.0; 
		double nominator_log_sum=0.0;

	  
		for(multimap<INT,INT>::const_iterator it=map_for_outdegree.begin();it!=map_for_outdegree.end();++it)
		{
			ret_outdegree = map_for_outdegree.equal_range(it->first);
		   
			   
			std::unordered_map<INT,INT>::const_iterator it3=node_outdegree.find(it->first);
				if(it3==node_outdegree.end())
			{
			node_outdegree.insert(std::pair<INT,INT>(it->first,1));
			}    	    
			else
			node_outdegree[it->first]++;

			  
			   unordered_set<INT>::const_iterator used_keys_it=used_keys.find(it->first);
			   
			   if(used_keys_it==used_keys.end())
			   {
		  

					used_keys.insert(it->first);


				  unordered_map<INT,INT> edge_multiplicity;

			 for (std::multimap<INT,INT>::iterator it2=ret_outdegree.first; it2!=ret_outdegree.second; ++it2)
			 { 		
				unordered_map<INT,INT>::const_iterator it3=edge_multiplicity.find(it2->second);
				if(it3==edge_multiplicity.end())
					edge_multiplicity.insert(std::pair<INT,INT>(it2->second,1));
				else
					edge_multiplicity[it2->second]++;
			 }
			 

				 
			 for(unordered_map<INT,INT>::const_iterator it4=edge_multiplicity.begin();it4!=edge_multiplicity.end();++it4)
			{
					
				  for(INT x=1;x<=it4->second;++x)
					{    
					  denominator_log_sum+=log(x);
					}
	  
				
					if(it->first != it4->first)
					{
						tripletList.push_back(T(node_id_consec[it->first],node_id_consec[it4->first],-1*it4->second));
				   }
			   else  
			   {
				a_uu[it->first]=it4->second;
			   }
				   
								 
			}
			 }
			
		}

		for(unordered_map<INT,INT>::const_iterator itx=node_id_consec.begin(); itx!=node_id_consec.end();++itx)
		{
		

			if(itx->first == node_t)
			{
				tripletList.push_back(T(itx->second,itx->second,node_outdegree[node_t]+1-a_uu[node_t]));
			}
			else
			{
				tripletList.push_back(T(itx->second,itx->second,node_outdegree[itx->first]-a_uu[itx->first]));
			}
		}

	   
		for(unordered_map<INT,INT>::const_iterator it=node_outdegree.begin();it!=node_outdegree.end();++it)
		{
			 
			  if(it->first!=node_t)
			 {
				 for(INT x=1;x<=it->second-1;++x)
					nominator_log_sum+=log(x);
			 }
			  else
			  {
				   for(INT x=1;x<=it->second;++x)
					 nominator_log_sum+=log(x);
			  }
		}

		std::pair <std::multimap<INT,INT>::iterator, std::multimap<INT,INT>::iterator> ret_indegree;
		std::unordered_map<INT,INT> node_indegree;

		for(multimap<INT,INT>::const_iterator it=map_for_indegree.begin();it!=map_for_indegree.end();++it)
		{

			ret_indegree = map_for_indegree.equal_range(it->first);
			  
			std::unordered_map<INT,INT>::const_iterator it3=node_indegree.find(it->first);
			if(it3==node_indegree.end())
			 {
				node_indegree.insert(std::pair<INT,INT>(it->first,1));
			 }
			else
				node_indegree[it->first]++;

		}


		

	 
		SparseMatrixType m((INT)matrix_dim,(INT)matrix_dim);
		m.setFromTriplets(tripletList.begin(), tripletList.end());

	delete []C; 
	free(invSA);
	free(SA);
	free(LCP);

	if(PRUNE_DET)
	{
	   if((nominator_log_sum - denominator_log_sum) >= log_z)
		{
			
		auto algo_end_no_det = chrono::steady_clock::now();
		cout<<"log(numerator): "<<nominator_log_sum<<" log(denominator): "<<denominator_log_sum<<" log(num)-log(denom):"<<nominator_log_sum-denominator_log_sum<<" log(z)"<<log_z<<endl;
		cout<<"exp(numerator): "<<exp(nominator_log_sum)<<" exp(denominator): "<<exp(denominator_log_sum)<<" ratio: "<<exp(nominator_log_sum)/exp(denominator_log_sum)<<endl;
		cout<<"Time for algorithm (no determinant) in ms: "<<chrono::duration_cast<chrono::milliseconds>(algo_end_no_det - algo_start).count()<<endl;

			return true;
		}
	}

		Eigen::SparseLU<Eigen::SparseMatrix<double>  > solver;
		solver.analyzePattern(m);
		solver.factorize(m);

		double logdet=solver.logAbsDeterminant();

	
	   
		double log_N_gt=logdet+nominator_log_sum-denominator_log_sum;


auto algo_end = chrono::steady_clock::now();
cout<<"Time for algorithm in ms: "<<chrono::duration_cast<chrono::milliseconds>(algo_end - algo_start).count()<<endl;
			
		cout<<"log(det):"<<logdet<<" log(numerator): "<<nominator_log_sum<<" log(denominator): "<<denominator_log_sum<<" log(det)+log(numerator)-log(denominator): "<<log_N_gt<<" log(z):"<<log_z<<endl;
		cout<<exp(logdet)*exp(nominator_log_sum)/exp(denominator_log_sum)<<endl;
		if(log_N_gt >= log_z)
		{
			cout<<"#ET: "<<exp(log_N_gt)<<endl;
			return true;
		}
		cout<<"#ET nodet: "<<exp(log_N_gt)<<endl;	
		return false;
}

bool just_check_function(unsigned char* seq, unsigned int z, unsigned int d)
{
        //clock_t begin=clock();
	auto begin= chrono::steady_clock::now();

	INT *SA;
        INT *LCP;
        INT *invSA;
        INT n=strlen ( ( char * ) seq );
	

        string tmpx((const char*)seq);

	INT r_S=-1;
        if(DEBUG)
                cout<<"Suffix array construction\n";

		SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
			invSA = ( INT * ) calloc( n , sizeof( INT ) );
			LCP = ( INT * ) calloc  ( n, sizeof( INT ) );			
			
			if( ( SA == NULL) )
			{
					return ( 0 );
			}

			if( ( invSA == NULL) )
			{
					fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
					return ( 0 );
			}
			
			if( ( LCP == NULL) )
			{
					fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
					return ( 0 );
			}
			
        suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);

        //clock_t end = clock();
        //double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	auto end=chrono::steady_clock::now();
	cout<<"Time :"<<chrono::duration_cast<chrono::milliseconds>(end - begin	).count()<<" ms for suffix array construction."<<endl;

        cout<<"r_S="<<r_S<<endl;
	
        double log_z=log(z);
	
	
       
	bool ret=dBgraph_just_check(seq, d, SA, LCP, invSA, log_z,n );
	
	return ret;
}

int binarydB_doubling_prefix( unsigned char* seq, unsigned int z, unsigned int k)
{
	INT *SA;
	INT *LCP;
	INT *invSA;
	INT n = strlen ( ( char * ) seq );

	string tmpx((const char*)seq);

	INT r_S = -1;

	cout << "Initial suffix array construction: " << endl;


	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}

	suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);

//			suffix_array_construction(seq, SA, LCP, invSA, r_S);


	cout << "r(S): " << r_S - 1 << endl;

	double log_z = log(z);

	INT l = 0;
	INT r;

	if (GAB_OPT_IN_DB)
		r = r_S;
	else
		r = n;


	srand(time(NULL));

	if (l < 0)
		l = 0;



	double time_for_doubling = 0.0;

	int cnt = 0;
	//int no_iter=0;

	int prefix_sz = 2 * k < n ? 2 * k : n;

	int largest_prefix_iter = ceil(log((double)n / (double)k) / log(2));


	bool largest_prefix_found = false;
	bool old_prefix = false;
	unsigned char *prefix = NULL;
	INT r_S_pref = -1;

	bool unique_prefix = false;
	int no_of_unique = 0, no_of_non_unique = 0;


	INT d = floor((double)(l + r) / 2.0);

	while (l < r)
	{
		if (!old_prefix)
		{

			prefix = (unsigned char *)strndup((const char *)seq, prefix_sz);

			r_S_pref = -1;

			free_arrays(SA, LCP, invSA);

			SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
			invSA = ( INT * ) calloc( n , sizeof( INT ) );
			LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

			suffix_array_construction_nomalloc(prefix, SA, LCP, invSA, r_S_pref);
			//suffix_array_construction(prefix, SA, LCP, invSA, r_S_pref);

			if (r_S_pref < d)
			{
				unique_prefix = true;
				no_of_unique++;
			}
			else
			{
				unique_prefix = false;
				no_of_non_unique++;
			}
			old_prefix = false;
		}


		d = floor((double)(l + r) / 2.0);
		//cout<<"d: "<<d<<endl;
		if (!largest_prefix_found)
		{
			while (cnt < largest_prefix_iter)
			{
				if (!unique_prefix && dBgraph (prefix, d, SA, LCP, invSA, log_z, prefix_sz ))
				{
					l = d + 1;

					old_prefix = true;

					break;
				}
				else
				{

					if (cnt == largest_prefix_iter)
					{

						old_prefix = true;
						r = d;
						break;
					}
					else
					{
						prefix_sz = prefix_sz * 2 < n ? prefix_sz * 2 : n;
						cnt++;

						if (cnt == largest_prefix_iter)
							largest_prefix_found = true;

						prefix = (unsigned char *)strndup((const char *)seq, prefix_sz);

						free ( invSA );
						free ( SA );
						free ( LCP );

						r_S_pref = -1;
						clock_t begin_sac = clock();

						SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
						invSA = ( INT * ) calloc( n , sizeof( INT ) );
						LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

						suffix_array_construction_nomalloc(prefix, SA, LCP, invSA, r_S_pref);
						//suffix_array_construction(prefix, SA, LCP, invSA, r_S_pref);
						clock_t end_sac = clock();

						time_for_doubling += double(end_sac - begin_sac) / CLOCKS_PER_SEC;

						if (r_S_pref < d)
						{
							unique_prefix = true;
							no_of_unique++;
						}
						else
						{
							unique_prefix = false;
							no_of_non_unique++;

						}

						if (prefix_sz == n)
						{
							//old_prefix==true;
							break;  //new addition
						}
					}
				}
			}

		}
		else
		{

			if (r_S_pref >= d && dBgraph (prefix, d, SA, LCP, invSA, log_z, prefix_sz ))
			{
				l = d + 1;

			}
			else
				r = d;
		}
	}

	if (l > 0)
	{
		cout << "Prefix size: " << prefix_sz << " Answer d = " << l - 1 << endl;
		d = l - 1;
	}

	if (l == 0)
		cout << "FAIL\n";

	cout << "Unique_prefixes: " << no_of_unique << " non_unique_prefixes: " << no_of_non_unique << " Ratio unique/all: " << (double)no_of_unique / ((double)(no_of_unique + no_of_non_unique)) << endl;
	free(prefix);
	free ( invSA );
	free ( SA );
	free ( LCP );


	return d;

}

int fixed_d_frontier( unsigned char* seq, unsigned int z, unsigned int k, unsigned int fixed_d)
{

	INT *SA;
	INT *LCP;
	INT *invSA;
	INT n = strlen ( ( char * ) seq );

	string tmpx((const char*)seq);

	INT r_S = -1;

	cout << "Initial suffix array construction: " << endl;


	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}

	suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);

	double log_z = log(z);


	bool good_d = dBgraph_frontier (fixed_d, SA, LCP, invSA, z, n );
	if (good_d == true)
	{
		return fixed_d;
	}

	return 0;
}



int fixed_d_frontier_compress( unsigned char* seq, unsigned int z, unsigned int k, unsigned int fixed_d)
{

	INT *SA;
	INT *LCP;
	INT *invSA;
	INT n = strlen ( ( char * ) seq );

	string tmpx((const char*)seq);

	INT r_S = -1;

	cout << "Initial suffix array construction: " << endl;


	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}

	suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);

	double log_z = log(z);


	bool good_d = dBgraph_frontier_compress (fixed_d, SA, LCP, invSA, z, n );
	if (good_d == true)
	{
		return fixed_d;
	}

	return 0;
}


	
int binarydB_doubling_prefix_frontier( unsigned char* seq, unsigned int z, unsigned int k)
{
	INT *SA;
	INT *LCP;
	INT *invSA;
	INT n = strlen ( ( char * ) seq );

	string tmpx((const char*)seq);

	INT r_S = -1;

	cout << "Initial suffix array construction: " << endl;


	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}

	suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);

//			suffix_array_construction(seq, SA, LCP, invSA, r_S);


	cout << "r(S): " << r_S - 1 << endl;

	double log_z = log(z);

	INT l = 0;
	INT r;

	if (GAB_OPT_IN_DB)
		r = r_S;
	else
		r = n;


	srand(time(NULL));

	if (l < 0)
		l = 0;



	double time_for_doubling = 0.0;

	int cnt = 0;
	//int no_iter=0;

	int prefix_sz = 2 * k < n ? 2 * k : n;

	int largest_prefix_iter = ceil(log((double)n / (double)k) / log(2));


	bool largest_prefix_found = false;
	bool old_prefix = false;
	unsigned char *prefix = NULL;
	INT r_S_pref = -1;

	bool unique_prefix = false;
	int no_of_unique = 0, no_of_non_unique = 0;


	INT d = floor((double)(l + r) / 2.0);

	while (l < r)
	{
		if (!old_prefix)
		{

			prefix = (unsigned char *)strndup((const char *)seq, prefix_sz);

			r_S_pref = -1;

			free_arrays(SA, LCP, invSA);

			SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
			invSA = ( INT * ) calloc( n , sizeof( INT ) );
			LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

			suffix_array_construction_nomalloc(prefix, SA, LCP, invSA, r_S_pref);
			//suffix_array_construction(prefix, SA, LCP, invSA, r_S_pref);

			if (r_S_pref < d)
			{
				unique_prefix = true;
				no_of_unique++;
			}
			else
			{
				unique_prefix = false;
				no_of_non_unique++;
			}
			old_prefix = false;
		}


		d = floor((double)(l + r) / 2.0);

		if (!largest_prefix_found)
		{
			while (cnt < largest_prefix_iter)
			{
				if (!unique_prefix && dBgraph_frontier (d, SA, LCP, invSA, z, prefix_sz ))
				{
					l = d + 1;

					old_prefix = true;

					break;
				}
				else
				{

					if (cnt == largest_prefix_iter)
					{

						old_prefix = true;
						r = d;
						break;
					}
					else
					{
						prefix_sz = prefix_sz * 2 < n ? prefix_sz * 2 : n;
						cnt++;

						if (cnt == largest_prefix_iter)
							largest_prefix_found = true;

						prefix = (unsigned char *)strndup((const char *)seq, prefix_sz);

						free ( invSA );
						free ( SA );
						free ( LCP );

						r_S_pref = -1;
						clock_t begin_sac = clock();

						SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
						invSA = ( INT * ) calloc( n , sizeof( INT ) );
						LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

						suffix_array_construction_nomalloc(prefix, SA, LCP, invSA, r_S_pref);
						//suffix_array_construction(prefix, SA, LCP, invSA, r_S_pref);
						clock_t end_sac = clock();

						time_for_doubling += double(end_sac - begin_sac) / CLOCKS_PER_SEC;

						if (r_S_pref < d)
						{
							unique_prefix = true;
							no_of_unique++;
						}
						else
						{
							unique_prefix = false;
							no_of_non_unique++;

						}

						if (prefix_sz == n)
						{
							//old_prefix==true;
							break;  //new addition
						}
					}
				}
			}

		}
		else
		{

			if (r_S_pref >= d && dBgraph_frontier ( d, SA, LCP, invSA, z, prefix_sz ))
			{
				l = d + 1;

			}
			else
				r = d;
		}

	}

	if (l > 0)
	{
		cout << "Prefix size: " << prefix_sz << " Answer d = " << l - 1 << endl;
		d = l - 1;
	}

	if (l == 0)
		cout << "FAIL\n";

	cout << "Unique_prefixes: " << no_of_unique << " non_unique_prefixes: " << no_of_non_unique << " Ratio unique/all: " << (double)no_of_unique / ((double)(no_of_unique + no_of_non_unique)) << endl;
	free(prefix);
	free ( invSA );
	free ( SA );
	free ( LCP );


	return d;

}
int binarydB_doubling_prefix_frontier_compress( unsigned char* seq, unsigned int z, unsigned int k)
{
	INT *SA;
	INT *LCP;
	INT *invSA;
	INT n = strlen ( ( char * ) seq );

	string tmpx((const char*)seq);

	INT r_S = -1;

	cout << "Initial suffix array construction: " << endl;


	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}

	suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);

//			suffix_array_construction(seq, SA, LCP, invSA, r_S);


	cout << "r(S): " << r_S - 1 << endl;

	double log_z = log(z);

	INT l = 0;
	INT r;

	if (GAB_OPT_IN_DB)
		r = r_S;
	else
		r = n;


	srand(time(NULL));

	if (l < 0)
		l = 0;



	double time_for_doubling = 0.0;

	int cnt = 0;
	//int no_iter=0;

	int prefix_sz = 2 * k < n ? 2 * k : n;

	int largest_prefix_iter = ceil(log((double)n / (double)k) / log(2));


	bool largest_prefix_found = false;
	bool old_prefix = false;
	unsigned char *prefix = NULL;
	INT r_S_pref = -1;

	bool unique_prefix = false;
	int no_of_unique = 0, no_of_non_unique = 0;


	INT d = floor((double)(l + r) / 2.0);

	while (l < r)
	{
		if (!old_prefix)
		{

			prefix = (unsigned char *)strndup((const char *)seq, prefix_sz);

			r_S_pref = -1;

			free_arrays(SA, LCP, invSA);

			SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
			invSA = ( INT * ) calloc( n , sizeof( INT ) );
			LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

			suffix_array_construction_nomalloc(prefix, SA, LCP, invSA, r_S_pref);
			//suffix_array_construction(prefix, SA, LCP, invSA, r_S_pref);

			if (r_S_pref < d)
			{
				unique_prefix = true;
				no_of_unique++;
			}
			else
			{
				unique_prefix = false;
				no_of_non_unique++;
			}
			old_prefix = false;
		}


		d = floor((double)(l + r) / 2.0);

		if (!largest_prefix_found)
		{
			while (cnt < largest_prefix_iter)
			{
				if (!unique_prefix && dBgraph_frontier_compress (d, SA, LCP, invSA, z, prefix_sz ))
				{
					l = d + 1;

					old_prefix = true;

					break;
				}
				else
				{

					if (cnt == largest_prefix_iter)
					{

						old_prefix = true;
						r = d;
						break;
					}
					else
					{
						prefix_sz = prefix_sz * 2 < n ? prefix_sz * 2 : n;
						cnt++;

						if (cnt == largest_prefix_iter)
							largest_prefix_found = true;

						prefix = (unsigned char *)strndup((const char *)seq, prefix_sz);

						free ( invSA );
						free ( SA );
						free ( LCP );

						r_S_pref = -1;
						clock_t begin_sac = clock();

						SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
						invSA = ( INT * ) calloc( n , sizeof( INT ) );
						LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

						suffix_array_construction_nomalloc(prefix, SA, LCP, invSA, r_S_pref);
						//suffix_array_construction(prefix, SA, LCP, invSA, r_S_pref);
						clock_t end_sac = clock();

						time_for_doubling += double(end_sac - begin_sac) / CLOCKS_PER_SEC;

						if (r_S_pref < d)
						{
							unique_prefix = true;
							no_of_unique++;
						}
						else
						{
							unique_prefix = false;
							no_of_non_unique++;

						}

						if (prefix_sz == n)
						{
							//old_prefix==true;
							break;  //new addition
						}
					}
				}
			}

		}
		else
		{

			if (r_S_pref >= d && dBgraph_frontier_compress ( d, SA, LCP, invSA, z, prefix_sz ))
			{
				l = d + 1;

			}
			else
				r = d;
		}

	}

	if (l > 0)
	{
		cout << "Prefix size: " << prefix_sz << " Answer d = " << l - 1 << endl;
		d = l - 1;
	}

	if (l == 0)
		cout << "FAIL\n";

	cout << "Unique_prefixes: " << no_of_unique << " non_unique_prefixes: " << no_of_non_unique << " Ratio unique/all: " << (double)no_of_unique / ((double)(no_of_unique + no_of_non_unique)) << endl;
	free(prefix);
	free ( invSA );
	free ( SA );
	free ( LCP );


	return d;
}
int fixed_d_tree( unsigned char* seq, unsigned int z, unsigned int k, unsigned int fixed_d)
{

	INT *SA;
	INT *LCP;
	INT *invSA;
	INT n = strlen ( ( char * ) seq );

	string tmpx((const char*)seq);

	INT r_S = -1;

	cout << "Initial suffix array construction: " << endl;


	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}

	suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);

	double log_z = log(z);


	bool good_d = dBgraph_tree (fixed_d, SA, LCP, invSA, z, n, false );
	if (good_d == true)
	{
		return fixed_d;
	}

	return 0;
}
int fixed_d_tree_compress( unsigned char* seq, unsigned int z, unsigned int k, unsigned int fixed_d)
{

	INT *SA;
	INT *LCP;
	INT *invSA;
	INT n = strlen ( ( char * ) seq );

	string tmpx((const char*)seq);

	INT r_S = -1;

	cout << "Initial suffix array construction: " << endl;


	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}

	suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);

	double log_z = log(z);


	bool good_d = dBgraph_tree (fixed_d, SA, LCP, invSA, z, n, true );
	if (good_d == true)
	{
		return fixed_d;
	}

	return 0;
}
int binarydB_doubling_prefix_tree( unsigned char* seq, unsigned int z, unsigned int k)
{
	INT *SA;
	INT *LCP;
	INT *invSA;
	INT n = strlen ( ( char * ) seq );

	string tmpx((const char*)seq);

	INT r_S = -1;

	cout << "Initial suffix array construction: " << endl;


	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}

	suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);

//			suffix_array_construction(seq, SA, LCP, invSA, r_S);


	cout << "r(S): " << r_S - 1 << endl;

	double log_z = log(z);

	INT l = 0;
	INT r;

	if (GAB_OPT_IN_DB)
		r = r_S;
	else
		r = n;


	srand(time(NULL));

	if (l < 0)
		l = 0;



	double time_for_doubling = 0.0;

	int cnt = 0;
	//int no_iter=0;

	int prefix_sz = 2 * k < n ? 2 * k : n;

	int largest_prefix_iter = ceil(log((double)n / (double)k) / log(2));


	bool largest_prefix_found = false;
	bool old_prefix = false;
	unsigned char *prefix = NULL;
	INT r_S_pref = -1;

	bool unique_prefix = false;
	int no_of_unique = 0, no_of_non_unique = 0;


	INT d = floor((double)(l + r) / 2.0);

	while (l < r)
	{
		if (!old_prefix)
		{

			prefix = (unsigned char *)strndup((const char *)seq, prefix_sz);

			r_S_pref = -1;

			free_arrays(SA, LCP, invSA);

			SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
			invSA = ( INT * ) calloc( n , sizeof( INT ) );
			LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

			suffix_array_construction_nomalloc(prefix, SA, LCP, invSA, r_S_pref);
			//suffix_array_construction(prefix, SA, LCP, invSA, r_S_pref);

			if (r_S_pref < d)
			{
				unique_prefix = true;
				no_of_unique++;
			}
			else
			{
				unique_prefix = false;
				no_of_non_unique++;
			}
			old_prefix = false;
		}


		d = floor((double)(l + r) / 2.0);

		if (!largest_prefix_found)
		{
			while (cnt < largest_prefix_iter)
			{
				if (!unique_prefix && dBgraph_tree (d, SA, LCP, invSA, z, prefix_sz,false ))
				{
					l = d + 1;

					old_prefix = true;

					break;
				}
				else
				{

					if (cnt == largest_prefix_iter)
					{

						old_prefix = true;
						r = d;
						break;
					}
					else
					{
						prefix_sz = prefix_sz * 2 < n ? prefix_sz * 2 : n;
						cnt++;

						if (cnt == largest_prefix_iter)
							largest_prefix_found = true;

						prefix = (unsigned char *)strndup((const char *)seq, prefix_sz);

						free ( invSA );
						free ( SA );
						free ( LCP );

						r_S_pref = -1;
						clock_t begin_sac = clock();

						SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
						invSA = ( INT * ) calloc( n , sizeof( INT ) );
						LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

						suffix_array_construction_nomalloc(prefix, SA, LCP, invSA, r_S_pref);
						//suffix_array_construction(prefix, SA, LCP, invSA, r_S_pref);
						clock_t end_sac = clock();

						time_for_doubling += double(end_sac - begin_sac) / CLOCKS_PER_SEC;

						if (r_S_pref < d)
						{
							unique_prefix = true;
							no_of_unique++;
						}
						else
						{
							unique_prefix = false;
							no_of_non_unique++;

						}

						if (prefix_sz == n)
						{
							//old_prefix==true;
							break;  //new addition
						}
					}
				}
			}

		}
		else
		{

			if (r_S_pref >= d && dBgraph_tree ( d, SA, LCP, invSA, z, prefix_sz,false ))
			{
				l = d + 1;

			}
			else
				r = d;
		}

	}

	if (l > 0)
	{
		cout << "Prefix size: " << prefix_sz << " Answer d = " << l - 1 << endl;
		d = l - 1;
	}

	if (l == 0)
		cout << "FAIL\n";

	cout << "Unique_prefixes: " << no_of_unique << " non_unique_prefixes: " << no_of_non_unique << " Ratio unique/all: " << (double)no_of_unique / ((double)(no_of_unique + no_of_non_unique)) << endl;
	free(prefix);
	free ( invSA );
	free ( SA );
	free ( LCP );


	return d;

}
int binarydB_doubling_prefix_tree_compress( unsigned char* seq, unsigned int z, unsigned int k)
{
	INT *SA;
	INT *LCP;
	INT *invSA;
	INT n = strlen ( ( char * ) seq );

	string tmpx((const char*)seq);

	INT r_S = -1;

	cout << "Initial suffix array construction: " << endl;


	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}

	suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);

//			suffix_array_construction(seq, SA, LCP, invSA, r_S);


	cout << "r(S): " << r_S - 1 << endl;

	double log_z = log(z);

	INT l = 0;
	INT r;

	if (GAB_OPT_IN_DB)
		r = r_S;
	else
		r = n;


	srand(time(NULL));

	if (l < 0)
		l = 0;



	double time_for_doubling = 0.0;

	int cnt = 0;
	//int no_iter=0;

	int prefix_sz = 2 * k < n ? 2 * k : n;

	int largest_prefix_iter = ceil(log((double)n / (double)k) / log(2));


	bool largest_prefix_found = false;
	bool old_prefix = false;
	unsigned char *prefix = NULL;
	INT r_S_pref = -1;

	bool unique_prefix = false;
	int no_of_unique = 0, no_of_non_unique = 0;


	INT d = floor((double)(l + r) / 2.0);

	while (l < r)
	{
		if (!old_prefix)
		{

			prefix = (unsigned char *)strndup((const char *)seq, prefix_sz);

			r_S_pref = -1;

			free_arrays(SA, LCP, invSA);

			SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
			invSA = ( INT * ) calloc( n , sizeof( INT ) );
			LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

			suffix_array_construction_nomalloc(prefix, SA, LCP, invSA, r_S_pref);
			//suffix_array_construction(prefix, SA, LCP, invSA, r_S_pref);

			if (r_S_pref < d)
			{
				unique_prefix = true;
				no_of_unique++;
			}
			else
			{
				unique_prefix = false;
				no_of_non_unique++;
			}
			old_prefix = false;
		}


		d = floor((double)(l + r) / 2.0);

		if (!largest_prefix_found)
		{
			while (cnt < largest_prefix_iter)
			{
				if (!unique_prefix && dBgraph_tree (d, SA, LCP, invSA, z, prefix_sz,true))
				{
					l = d + 1;

					old_prefix = true;

					break;
				}
				else
				{

					if (cnt == largest_prefix_iter)
					{

						old_prefix = true;
						r = d;
						break;
					}
					else
					{
						prefix_sz = prefix_sz * 2 < n ? prefix_sz * 2 : n;
						cnt++;

						if (cnt == largest_prefix_iter)
							largest_prefix_found = true;

						prefix = (unsigned char *)strndup((const char *)seq, prefix_sz);

						free ( invSA );
						free ( SA );
						free ( LCP );

						r_S_pref = -1;
						clock_t begin_sac = clock();

						SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
						invSA = ( INT * ) calloc( n , sizeof( INT ) );
						LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

						suffix_array_construction_nomalloc(prefix, SA, LCP, invSA, r_S_pref);
						//suffix_array_construction(prefix, SA, LCP, invSA, r_S_pref);
						clock_t end_sac = clock();

						time_for_doubling += double(end_sac - begin_sac) / CLOCKS_PER_SEC;

						if (r_S_pref < d)
						{
							unique_prefix = true;
							no_of_unique++;
						}
						else
						{
							unique_prefix = false;
							no_of_non_unique++;

						}

						if (prefix_sz == n)
						{
							//old_prefix==true;
							break;  //new addition
						}
					}
				}
			}

		}
		else
		{

			if (r_S_pref >= d && dBgraph_tree ( d, SA, LCP, invSA, z, prefix_sz, true ))
			{
				l = d + 1;

			}
			else
				r = d;
		}

	}

	if (l > 0)
	{
		cout << "Prefix size: " << prefix_sz << " Answer d = " << l - 1 << endl;
		d = l - 1;
	}

	if (l == 0)
		cout << "FAIL\n";

	cout << "Unique_prefixes: " << no_of_unique << " non_unique_prefixes: " << no_of_non_unique << " Ratio unique/all: " << (double)no_of_unique / ((double)(no_of_unique + no_of_non_unique)) << endl;
	free(prefix);
	free ( invSA );
	free ( SA );
	free ( LCP );


	return d;

}
int binarydB_no_prefix( unsigned char* seq, unsigned int z, unsigned int k)
{
	INT *SA;
	INT *LCP;
	INT *invSA;
	INT n = strlen ( ( char * ) seq );
	string tmpx((const char*)seq);
	INT r_S = -1;

	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}


	suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);
	double log_z = log(z);

	INT l = 0;
	INT r = n;

	srand(time(NULL));

	INT d = 0;

	while (l < r)
	{
		d = floor((double)(l + r) / 2.0);
		if (dBgraph (seq, d, SA, LCP, invSA, log_z, n ))
		{
			l = d + 1;
		}
		else
		{
			r = d;
		}
	}
	if (l > 0) {
		cout << "Answer d = " << l - 1 << endl;
		d = l - 1;
	}
	if (l == 0)  cout << "FAIL\n";
	free ( invSA );
	free ( SA );
	free ( LCP );
	return d;
}

int binarydB_no_prefix_gab( unsigned char* seq, unsigned int z, unsigned int k)
{
	INT *SA;
	INT *LCP;
	INT *invSA;
	INT n = strlen ( ( char * ) seq );
	string tmpx((const char*)seq);
	INT r_S = -1;


	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}

	suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);


	double log_z = log(z);

	INT l = 0;
	INT r = 0;
	if (GAB_OPT_IN_NP)
	{
		r = r_S;
	}
	else
		r = n;


	srand(time(NULL));
	INT d = 0;

	while (l < r)
	{
		d = floor((double)(l + r) / 2.0);
		if (dBgraph (seq, d, SA, LCP, invSA, log_z, n ))
		{
			l = d + 1;
		}
		else
		{
			r = d;
		}
	}
	if (l > 0) {
		cout << "Answer d = " << l - 1 << endl;
		d = l - 1;
	}
	if (l == 0)  cout << "FAIL\n";
	free ( invSA );
	free ( SA );
	free ( LCP );
	return d;
}

int binarydB_no_prefix_exp( unsigned char* seq, unsigned int z, unsigned int k)
{
	INT *SA;
	INT *LCP;
	INT *invSA;
	INT n = strlen ( ( char * ) seq );
	string tmpx((const char*)seq);
	INT r_S = -1;

	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );

	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}

	suffix_array_construction_nomalloc(seq, SA, LCP, invSA, r_S);

	double log_z = log(z);

	INT l = 0;
	INT r;
	r = n;
	srand(time(NULL));
	INT i = 2;
	while (i < n && dBgraph (seq, i, SA, LCP, invSA, log_z, n )) {
		i = i * 2;
	}
	l = i / 2;
	r = (i < r) ? i : r;
	//int no_iter=0;
	INT d;

	while (l < r)
	{

		d = floor((double)(l + r) / 2.0);

		if (dBgraph (seq, d, SA, LCP, invSA, log_z, n ))
		{
			l = d + 1;
		}
		else  //a_d < z
		{
			r = d;
		}
	}
	if (l > 0)
	{
		cout << "Answer d = " << l - 1 << endl;
		d = l - 1;
	}
	if (l == 0) cout << "FAIL\n";

	free ( invSA );
	free ( SA );
	free ( LCP );
	return d;
}
void binarydB( unsigned char* seq, unsigned int z)
{
	INT *SA;
	INT *LCP;
	INT *invSA;

	INT n = strlen ( ( char * ) seq );

	INT r_S = -1;

	suffix_array_construction(seq, SA, LCP, invSA, r_S);

	double log_z = log(z);

	INT l = 0;
	INT r;

	if (GAB_OPT_IN_DB)
		r = r_S;
	else
		r = n;


	srand(time(NULL));


	if (l < 0)
		l = 0;



	while (l < r)
	{
		INT d = floor((double)(l + r) / 2.0);

		if (d == 0)
		{
			l = 0;
			break;
		}
		if (dBgraph (seq, d, SA, LCP, invSA, log_z, n ))
			l = d + 1;
		else
			r = d;
	}
	if (l > 0)
	{
		cout << "Answer d = " << l - 1 << endl;
	}
	if (l == 0)
		cout << "FAIL\n";



	free ( invSA );
	free ( SA );
	free ( LCP );

}

double factorial_compute(unsigned int x)
{
	try {
		return factorial<double>(x);
	}
	catch (exception& e)
	{
		cout << e.what() << endl;
	}

	return 0;
}


unsigned int LCParray ( unsigned char * x, INT n, INT * SA, INT * ISA, INT * LCP )
{
	INT i = 0, j = 0;

	LCP[0] = 0;
	for ( i = 0; i < n; i++ )
		if ( ISA[i] != 0 )
		{
			if ( i == 0) j = 0;
			else j = (LCP[ISA[i - 1]] >= 2) ? LCP[ISA[i - 1]] - 1 : 0;
			while ( x[i + j] == x[SA[ISA[i] - 1] + j] )
				j++;
			LCP[ISA[i]] = j;
		}

	return ( 1 );
}

int suffix_array_construction(unsigned char *x, INT *&SA, INT *&LCP, INT *&invSA, INT &r_S)
{

	INT n = strlen ( ( char * ) x );

	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	if ( ( SA == NULL) )
	{
		return ( 0 );
	}

#ifdef _USE_64
	if ( divsufsort64( x, SA,  n ) != 0 )
	{
		exit( EXIT_FAILURE );
	}
#endif

	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	if ( ( invSA == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
		return ( 0 );
	}

	for ( INT i = 0; i < n; i ++ )
		invSA [SA[i]] = i;


	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
	if ( ( LCP == NULL) )
	{
		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
		return ( 0 );
	}

	if ( LCParray( x, n, SA, invSA, LCP ) != 1 )
	{
		fprintf(stderr, " Error: LCP computation failed.\n" );
		exit( EXIT_FAILURE );
	}


	if (GAB_OPT_IN_DB || GAB_OPT_IN_NP)
	{	INT maxLCP = -1;

		for (int i = 0; i < n; ++i)
		{
			if (LCP[i] > maxLCP)
				maxLCP = LCP[i];
		}

		r_S = maxLCP + 1;
	}
	return 0;
}
int suffix_array_construction_nomalloc(unsigned char *x, INT *SA, INT *LCP, INT *invSA, INT &r_S)
{

	INT n = strlen ( ( char * ) x );



#ifdef _USE_64
	if ( divsufsort64( x, SA,  n ) != 0 )
	{
		exit( EXIT_FAILURE );
	}
#endif




	for ( INT i = 0; i < n; i ++ )
		invSA [SA[i]] = i;




	if ( LCParray( x, n, SA, invSA, LCP ) != 1 )
	{
		fprintf(stderr, " Error: LCP computation failed.\n" );
		exit( EXIT_FAILURE );
	}


	if (GAB_OPT_IN_DB || GAB_OPT_IN_NP)
	{	INT maxLCP = -1;

		for (int i = 0; i < n; ++i)
		{
			if (LCP[i] > maxLCP)
				maxLCP = LCP[i];
		}

		r_S = maxLCP + 1;
	}
	return 0;
}

bool dBgraph (  unsigned char * x, unsigned int d, INT *SA, INT *LCP, INT *invSA, double log_z, INT n)
{

	if (d == 0)
		return true;

	INT cluster_id = (INT) 0;
	INT *C = new INT[n];

	for (INT i = 0; i < n ; i++)
	{
		if (LCP[i] >= d - 1)
		{
			C[i] = (INT)cluster_id;

		}
		else
		{
			cluster_id++;
			C[i] = cluster_id;

		}

	}


	multimap<int, int> map_for_outdegree;


	multimap<int, int> map_for_indegree;


	unordered_map<int, int> node_id_consec;
	int node_id_consec_ind = 0;

	unordered_set<INT> distinct_nodes;

	for (INT i = 0; i <= n - d ; ++i)
	{
		INT cluster_id_pref = C[invSA[i]];
		INT cluster_id_suff = C[invSA[i + 1]];

		distinct_nodes.insert(cluster_id_pref);
		distinct_nodes.insert(cluster_id_suff);


		map_for_outdegree.insert(std::pair<int, int>(cluster_id_pref, cluster_id_suff));
		map_for_indegree.insert(std::pair<int, int>(cluster_id_suff, cluster_id_pref));

		unordered_map<int, int>::const_iterator nic_find = node_id_consec.find(cluster_id_pref);

		if (nic_find == node_id_consec.end())
		{	node_id_consec[cluster_id_pref] = node_id_consec_ind;
			node_id_consec_ind++;
		}

		unordered_map<int, int>::const_iterator nic_find2 = node_id_consec.find(cluster_id_suff);

		if (nic_find2 == node_id_consec.end())
		{	node_id_consec[cluster_id_suff] = node_id_consec_ind;
			node_id_consec_ind++;
		}

	}



	INT matrix_dim = distinct_nodes.size();

	std::vector<T> tripletList;

	int node_t = C[invSA[n - d + 1]];

	std::pair <std::multimap<int, int>::iterator, std::multimap<int, int>::iterator> ret_outdegree;
	std::unordered_map<int, int> node_outdegree;

	unordered_set<int> used_keys;

	unordered_map<int, int> a_uu;



	double denominator_log_sum = 0.0;
	double nominator_log_sum = 0.0;


	for (multimap<int, int>::const_iterator it = map_for_outdegree.begin(); it != map_for_outdegree.end(); ++it)
	{
		ret_outdegree = map_for_outdegree.equal_range(it->first);


		std::unordered_map<int, int>::const_iterator it3 = node_outdegree.find(it->first);
		if (it3 == node_outdegree.end())
		{
			node_outdegree.insert(std::pair<int, int>(it->first, 1));
		}
		else
			node_outdegree[it->first]++;


		unordered_set<int>::const_iterator used_keys_it = used_keys.find(it->first);

		if (used_keys_it == used_keys.end())
		{


			used_keys.insert(it->first);


			unordered_map<int, int> edge_multiplicity;

			for (std::multimap<int, int>::iterator it2 = ret_outdegree.first; it2 != ret_outdegree.second; ++it2)
			{
				unordered_map<int, int>::const_iterator it3 = edge_multiplicity.find(it2->second);
				if (it3 == edge_multiplicity.end())
					edge_multiplicity.insert(std::pair<int, int>(it2->second, 1));
				else
					edge_multiplicity[it2->second]++;
			}



			for (unordered_map<int, int>::const_iterator it4 = edge_multiplicity.begin(); it4 != edge_multiplicity.end(); ++it4)
			{

				for (int x = 1; x <= it4->second; ++x)
				{
					denominator_log_sum += log(x);
				}


				if (it->first != it4->first)
				{
					tripletList.push_back(T(node_id_consec[it->first], node_id_consec[it4->first], -1 * it4->second));
				}
				else
				{
					a_uu[it->first] = it4->second;
				}


			}
		}

	}

	for (unordered_map<int, int>::const_iterator itx = node_id_consec.begin(); itx != node_id_consec.end(); ++itx)
	{


		if (itx->first == node_t)
		{
			tripletList.push_back(T(itx->second, itx->second, node_outdegree[node_t] + 1 - a_uu[node_t]));
		}
		else
		{
			tripletList.push_back(T(itx->second, itx->second, node_outdegree[itx->first] - a_uu[itx->first]));
		}
	}


	for (unordered_map<int, int>::const_iterator it = node_outdegree.begin(); it != node_outdegree.end(); ++it)
	{

		if (it->first != node_t)
		{
			for (int x = 1; x <= it->second - 1; ++x)
				nominator_log_sum += log(x);
		}
		else
		{
			for (int x = 1; x <= it->second; ++x)
				nominator_log_sum += log(x);
		}

	}

	std::pair <std::multimap<int, int>::iterator, std::multimap<int, int>::iterator> ret_indegree;
	std::unordered_map<int, int> node_indegree;

	for (multimap<int, int>::const_iterator it = map_for_indegree.begin(); it != map_for_indegree.end(); ++it)
	{

		ret_indegree = map_for_indegree.equal_range(it->first);

		std::unordered_map<int, int>::const_iterator it3 = node_indegree.find(it->first);
		if (it3 == node_indegree.end())
		{
			node_indegree.insert(std::pair<int, int>(it->first, 1));
		}
		else
			node_indegree[it->first]++;

	}





	SparseMatrixType m((INT)matrix_dim, (INT)matrix_dim);
	m.setFromTriplets(tripletList.begin(), tripletList.end());

	delete []C;
	cout << "d=" << d << endl;

	if (PRUNE_DET)
	{
		if ((nominator_log_sum - denominator_log_sum) >= log_z)
		{
			cout << "ratio num: " << exp(nominator_log_sum) << " ratio denom: " << exp(denominator_log_sum) << endl;
			cout << "ratio: " << exp(nominator_log_sum) / exp(denominator_log_sum) << endl;
			cout << "prune det YES" << endl;
			return true;
		}
	}

	Eigen::SparseLU<Eigen::SparseMatrix<double>  > solver;
	solver.analyzePattern(m);
	solver.factorize(m);

	double logdet = solver.logAbsDeterminant();



	double log_N_gt = logdet + nominator_log_sum - denominator_log_sum;

	if (log_N_gt >= log_z)
	{
		cout << "det YES" << endl;
		return true;
	}
	cout << "ratio: " << exp(nominator_log_sum) << " " << exp(denominator_log_sum) << " " << exp(log_z) << endl;
	cout << "det NO " << exp(log_N_gt) << " " << exp(log_z) << endl;
	return false;
}
bool dBgraph_frontier (unsigned int d, INT *SA, INT *LCP, INT *invSA, unsigned int z, INT n)
{

	if (d == 0)
		return true;

	INT cluster_id = (INT) 0;
	INT *C = new INT[n];

	for (INT i = 0; i < n ; i++)
	{
		if (LCP[i] >= d - 1)
		{
			C[i] = (INT)cluster_id;

		}
		else
		{
			cluster_id++;
			C[i] = cluster_id;

		}

	}

	//maps each cluster_id (i.e., node) to a consecutive integer starting from 0
	//this is needed because the matrix L needs to have co-ordinates that are consecutive integers starting from 0
	unordered_map<int, int> node_id_consec;
	int node_id_consec_ind = 0;

	vector<pair<INT, INT> > edges;

	for (INT i = 0; i <= n - d ; ++i)
	{
		//INT cluster_id_pref=C[invSA[i]]; //this is cluster id of prefix
		//INT cluster_id_suff=C[invSA[i+1]]; // this is cluster if of suffix
		edges.push_back(make_pair<INT, INT>((INT) C[invSA[i]], (INT) C[invSA[i + 1]]));

		//cout<<C[invSA[i]]<<" "<<C[invSA[i+1]]<<endl;
	}

	vector<pair<INT, INT> > edges_from_file;

	int line_cnt = 0;
	int source_from_file = -1;
	int target_from_file = -1;

	//stores for each node-id in input_file its corresponding node-id in the graph here (graph nodes need to be consecutive and start from 0


	unordered_map<int, int> consec_nodes;
	int start_cnt = 0;

	for (auto & it_edges : edges)
	{
		int start_node = it_edges.first;
		int end_node = it_edges.second;

		unordered_map<int, int>::const_iterator it = consec_nodes.find(start_node);

		int start_node_consec = -1;

		if (it == consec_nodes.end())
		{
			consec_nodes[start_node] = start_cnt;
			start_node_consec = start_cnt;
			start_cnt++;
		}
		else
			start_node_consec = consec_nodes[start_node];

		it = consec_nodes.find(end_node);
		int end_node_consec = -1;

		if (it == consec_nodes.end())
		{
			consec_nodes[end_node] = start_cnt;
			end_node_consec = start_cnt;
			start_cnt++;
		}
		else
			end_node_consec = consec_nodes[end_node];

		edges_from_file.push_back(make_pair<int, int>((int)start_node_consec, (int)end_node_consec));

	}

	Graph g1(start_cnt, edges_from_file.front().first, edges_from_file.back().second);

	if (DEBUG)
	{
		cout << "source: " << edges_from_file.front().first << endl;
		cout << "target: " << edges_from_file.back().second << endl;

		cout << "Read edges\n";
	}

	map<pair<int, int>, int > edge_multiplicities;
	unordered_map<int, int> node_outdegree;

	for (auto &it : edges_from_file)
	{
		if (DEBUG)cout << it.first << " " << it.second << endl;
		g1.addEdge(it.first, it.second);

		auto it_find = edge_multiplicities.find(make_pair(it.first, it.second));
		if (it_find == edge_multiplicities.end())
			edge_multiplicities[make_pair<int, int>(it.first, it.second)] = 1;
		else
			edge_multiplicities[make_pair<int, int>(it.first, it.second)]++;

		auto it_find2 = node_outdegree.find(it.first);
		if (it_find2 == node_outdegree.end())
			node_outdegree[it.first] = 1;
		else
			node_outdegree[it.first]++;

	}
	double numerator_log_sum = 0.0;
	double denominator_log_sum = 0.0;

	//cout<<"Node outdegrees:\n";
	for (auto &it : node_outdegree)
	{
		//cout<<it.first<<" "<<it.second<<endl;
		if (it.first != edges_from_file.back().second)
		{
			for (int x = 1; x <= it.second - 1; ++x)
				numerator_log_sum += log(x);
		}
		else
		{
			for (int x = 1; x <= it.second; ++x)
				numerator_log_sum += log(x);
		}

	}

	//cout<<"Edge multi:\n";
	for (auto &it : edge_multiplicities)
	{
		//cout<<it.first.first<<","<<it.first.second<<" "<<it.second<<endl;
		for (int x = 1; x <= it.second; ++x)
		{
			denominator_log_sum += log(x);
		}
	}
	if (numerator_log_sum - denominator_log_sum >= log(z))
	{
		cout << "based on ratio d: " << d << " YES" << endl;
		return true;
	}
//  cout<<"Numerator: "<<exp(numerator_log_sum)<<" Denominator: "<<exp(denominator_log_sum)<<endl;




	cout << "\n Current Graph. Number of nodes: " << g1.V << " number of edges: " << edges_from_file.size() << endl;

	if (!g1.isEulerian())
	{
		cout << "Not eulerian: NO\n";
		exit(-1);
	}
	//g1.printGraph();

	g1.SCC_iter();

	//creates stack and adds sccs
	stack<vector<scc > > ST;


	//remove trivial sccs
	for (vector<scc>::iterator it = g1.sccs.begin(); it != g1.sccs.end();)
	{
		//scc x = *it;
		//if (x.nodes.size() == 1)
		if (it->nodes.size() == 1)
		{

			it = g1.sccs.erase(it);

		}
		else
			++it;
	}

	if (g1.sccs.empty())
	{
		cout << "All SCCs of input graph are trivial.\n";
		cout << "One distinct Eulerian tour.\n";
		if (z > 1)
			cout << "z>1, so NO.\n";
		else if (z == 1)
			cout << "z=1, so YES.\n";

		return 0;
	}

	ST.push(g1.sccs);

	int bound = 1;
	for (auto &it : g1.sccs)
	{
		bound *= it.lb;
		//cout<<it.lb<<", ";
	}

	//if(DEBUG)
	//if(d==4)
	cout << "\n bound= " << bound << " " << z << endl;

	//fixed seed fro reproducibility
	/*unsigned int seed1=1234567;
	    std::mt19937 gen (seed1);*/

	//replace with the folloiwing if you need different seeds
	std::random_device rd;
	std::mt19937 gen(rd());

	int parameter_t = 0;      //how many times times were popped from ST;
	int parameter_for = 0;    //how many times the for loop in step 11 is executed;
	int unary_iterations=0;

	bool ST_was_empty = false;



	while (bound < z)
	{

		if (ST.empty())
		{
			//cout<<"NO"<<endl;
			ST_was_empty = true;
			break;
		}

		//GL for debuginning
		if (DEBUG)
		{

			//GL: this is for debugging
			stack<vector<scc > >  tmp_ST(ST);
			cout << " ---- Contents of ST before pop ------" << endl;
			while (!tmp_ST.empty())
			{
				cout << " ------- Vector of sccs -----" << endl;
				vector<scc> c = tmp_ST.top();
				for (auto &it : c)
					print_single_scc(it);
				tmp_ST.pop();
			}
			cout << " --------- End of tmp_ST ------- " << endl;
		}

		vector<scc> f_copy = ST.top();  //get the top element (vector of sccs) and remove it from the list
		vector<scc> f;

		for (int xxx = 0; xxx < f_copy.size(); ++xxx)
		{
			f.push_back(f_copy[xxx]);
		}

		//vector<scc> f(f_copy);  //this is because f_copy is just a reference and is lost after pop
		ST.pop();
		parameter_t++;


		if (DEBUG)
		{
			cout << "  --- what I popped --- \n";
			for (auto &it : f)
			{
				print_single_scc(it);
			}
			cout << " ---  end of what I popped --- \n";

			stack<vector<scc > >  tmp_ST(ST);
			cout << " ---- Contents of ST after pop ------" << endl;
			while (!tmp_ST.empty())
			{
				cout << " ------- Vector of sccs -----" << endl;
				vector<scc> c = tmp_ST.top();
				for (auto &it : c)
					print_single_scc(it);
				tmp_ST.pop();
			}
			cout << " --------- End of tmp_ST ------- " << endl;
		}

		//get a random scc from f  with index 0 to f.size()-1

		std::uniform_int_distribution<int> dist(0, f.size() - 1);
		int scc_rand_id = dist(gen);
		if (DEBUG)
			cout << "SCC rand id: " << scc_rand_id << endl;


		int f_bound = 1;
		for (auto &it : f) //f is a vector of scc structs
		{
			//scc x = it;
			//if (DEBUG)
			//cout << "# " << it.nodes.size() << " " << it.lb << endl;
			//f_bound *= x.lb;
			f_bound *= it.lb;
		}
		bound = bound - f_bound;

		/*if(d==4)
		{	cout<<f.size()<<" fb:"<<f_bound<<" b:"<<bound<<endl;

		}*/
		scc x = f[scc_rand_id]; //my current C_i from the vector of sccs f

		if (DEBUG)
			cout << "\n I chose " << scc_rand_id << endl;

		//19-10-2022 addition
		f.erase(f.begin() + scc_rand_id); //remove the C_i you chose randomly above

		//for each distinct out-neighbor of the source node
		unordered_set<int> seen_outneighbors;
		//for(list<int>::const_iterator it=g1.adj[x.source].begin(); it!=g1.adj[x.source].end();++it)

		if(x.edges.size()<=1)
			unary_iterations++;

		for (vector<pair<int, int> >::const_iterator it = x.edges.begin(); it != x.edges.end(); ++it)
		{
			//19-10-22
			vector<scc> f_neighbor(f);


			parameter_for++;

			//if the edge does not start from x.source continue, so that you get only outneighbors of x.source (non-distinct)
			if (it->first != x.source) continue;

			int outneighbor = it->second;


			if (DEBUG)
			{
				cout << "out: " << outneighbor << endl;
				cout << "x.source=" << x.source << " outneighbor: " << outneighbor << endl;
			}
			//if out-neighbor of source in C_i
			//if(g1.node_scc_assoc[*it]==g1.node_scc_assoc[x.source])
			//{
			//find the out-neighbor in seen
			unordered_set<int>::const_iterator it2 = seen_outneighbors.find(outneighbor);

			//if you have not seen the outneighbor before
			if (it2 == seen_outneighbors.end())
			{
				//compute SCC in the graph of the edges of the component except (x.source,*it)
				//this graph has the same number of nodes as the currect C_i
				//set source and target to -1
				if (DEBUG)
				{
					cout << "C\(s,u) node-set size:" << x.nodes.size() << endl;
					cout << "I will delete: " << x.source << "->" << outneighbor << endl;

				}
				Graph cur_C_except_edge(x.nodes.size(), -1, -1);

				unordered_map<int, int> consec;
				int consec_id = 0;

				int multiplicity_of_removed = 0;
				for (auto & it3 : x.edges)
				{
					//cout<<"edge:"<<it3.first<<" "<<it3.second<<endl;
					if (DEBUG)
					{	for (auto &itx : consec)
						{
							cout << "<" << itx.first << " " << itx.second << ">" << endl;
						}
					}
					if (it3.first != x.source || it3.second != outneighbor)			//add everything except edge (x.source,*it)
					{
						unordered_map<int, int>::const_iterator it_cons = consec.find(it3.first);

						unordered_map<int, int>::const_iterator it_cons2 = consec.find(it3.second);

						if (it_cons == consec.end() && it_cons2 == consec.end())	//I didn't find it3.first nor it3.second, I need to put consec_id
						{
							int start_node = 1;
							int end_node = -1;

							consec[it3.first] = consec_id;
							start_node = consec_id;
							consec_id++;

							it_cons2 = consec.find(it3.second); //because the addition of it3.first updates consec and it3.first may be it3.second

							if (it_cons2 == consec.end()) //
							{
								consec[it3.second] = consec_id;
								end_node = consec_id;
								consec_id++;
							}
							else                        //it3.first == it3.second
								end_node = start_node;


							cur_C_except_edge.addEdge(start_node, end_node);

							//cout<<"1had:"<<it3.first<<" "<<it3.second<<" added: "<<start_node<<","<<end_node<<endl;

						}
						else if (it_cons != consec.end() && it_cons2 == consec.end())
						{
							consec[it3.second] = consec_id;
							cur_C_except_edge.addEdge(consec[it3.first], consec_id);
							//cout<<"2had:"<<it3.first<<" "<<it3.second<<"added: "<<consec[it3.first]<<","<<consec_id<<endl;
							consec_id++;
						}
						else if (it_cons == consec.end() && it_cons2 != consec.end())
						{
							consec[it3.first] = consec_id;
							cur_C_except_edge.addEdge(consec_id, consec[it3.second]);
							//cout<<"3had:"<<it3.first<<" "<<it3.second<<"added: "<<consec_id<<","<<consec[it3.second]<<endl;
							consec_id++;
						}
						else if (it_cons != consec.end() && it_cons2 != consec.end())
						{
							cur_C_except_edge.addEdge(consec[it3.first], consec[it3.second]);
							//cout<<"4had:"<<it3.first<<" "<<it3.second<<"added: "<<consec[it3.first]<<","<<consec[it3.second]<<endl;
						}
					}
					if (it3.first == x.source && it3.second == outneighbor)
					{
						multiplicity_of_removed++;
						//cout<<x.source<<"->"<<outneighbor<<" has multiplicity: "<<multiplicity_of_removed<<endl;
					}
				}

				if (multiplicity_of_removed > 1)
				{
					unordered_map<int, int>::const_iterator it_cons = consec.find(x.source);
					unordered_map<int, int>::const_iterator it_cons2 = consec.find(outneighbor);
					int mapped_1 = -1;
					int mapped_2 = -1;

					if (it_cons == consec.end())
					{	mapped_1 = consec_id;
						consec[x.source] = consec_id; //new addition
						consec_id++;
					}
					else mapped_1 = consec[x.source];

					if (it_cons2 == consec.end())
					{
						mapped_2 = consec_id;
						consec[outneighbor] = consec_id; //new addition
						consec_id++;
					}
					else mapped_2 = consec[outneighbor];

					//cout<<mapped_1<<" ! "<<mapped_2<<endl;

					for (int i = 1; i < multiplicity_of_removed; ++i)
					{

						cur_C_except_edge.addEdge(mapped_1, mapped_2);
						//cout<<"had:"<<x.source<<","<<outneighbor<<" added:"<<mapped_1<<","<<mapped_2<<endl;
					}
				}

				cur_C_except_edge.source = outneighbor; //the node *it is the new source of graph C\(s,u)
				cur_C_except_edge.target = x.target; //the target of C_i is the target of C\(s,u)

				if (DEBUG)
				{
					cout << "C_i\(s,u) source :" << x.source << " outneighbor: " << outneighbor << endl;

					cout << "C_i\(s,u) graph:";
					cur_C_except_edge.printGraph();
				}



				/*if (!f.empty())
				{
				cout<<" ------ in !f.empty() BEFORE ---- \n";
				for(auto &it: f)
				{
						print_single_scc(it);

				}
				 cout<<" ------ end in !f.empty() ---- \n";
				cout<<"scc_rand_id: "<<scc_rand_id<<endl;


				//-->					f.erase(f.begin() + scc_rand_id); //remove the C_i you chose randomly above


				if(d==4 && bound>0)
					cout << "Removed: " << scc_rand_id << endl;

				cout<<" ------ in !f.empty() AFTER ---- \n";
				for(auto &it: f)
				{
						print_single_scc(it);

				}
				 cout<<" ------ end in !f.empty() ---- \n";
				cout<<"scc_rand_id: "<<scc_rand_id<<endl;

				}

				//these SCCs will be annotated with source and target nodes belonging only to the nodes of this graph -- maybe wrong
				for(auto &itz : consec)
				{
				cout<<"Map: "<<itz.first<<"->"<<itz.second<<endl;
				}
				*/

				cur_C_except_edge.SCC_map_iter(consec);




				//	cur_C_except_edge.printSccs();

				// add these new components that have size >1 to f with f.push_back();
				for (auto &itx : cur_C_except_edge.sccs)
				{
					//scc x2 = itx;
					//if (x2.nodes.size() > 1)
					if (itx.nodes.size() > 1)
						f_neighbor.push_back(itx);
				}



				if (!f_neighbor.empty())
				{
					ST.push(f_neighbor);
					f_bound = 1;
					for (auto &ity : f_neighbor) //f is a vector of scc structs
					{
						f_bound *= ity.lb;
					}

					bound += f_bound;
				}
				else
				{
					bound++;
				}

				seen_outneighbors.insert(outneighbor);   //add the out-neighbor to seen
				//GL for debugging
				/*if(d==4)
				{
					cout<<"Seen out-neighbors\n";
					for(auto& itg : seen_outneighbors)
						cout<<"["<<itg<<"]";
					cout<<endl;
				}*/

			}// end of  if(it2==seen_outneighbors.end())

			//}  // end of if(g1.node_scc_assoc[*it]==g1.node_scc_assoc[x.source])
		} //end of for(vector<pair<int,int> >::const_iterator it=x.edges.begin(); it!=x.edges.end(); ++it)

	}
	cout << "parameter t: " << parameter_t << endl;
	cout << "parameter_for: " << parameter_for << endl;
	cout <<"unary iterations: "<<unary_iterations<<endl;
	cout << "bound: " << bound << endl;
	
	delete []C;
	if (ST_was_empty)
	{
		cout << "d: " << d << " NO\n";
		return false;
	}
	else
	{	cout << "d: " << d << " YES \n\n";
		return true;
	}


	return false;
}

bool dBgraph_frontier_compress(unsigned int d, INT *SA, INT *LCP, INT *invSA, unsigned int z, INT n)
{

	if (d == 0)
		return true;

	INT cluster_id = (INT) 0;
	INT *C = new INT[n];

	for (INT i = 0; i < n ; i++)
	{
		if (LCP[i] >= d - 1)
		{
			C[i] = (INT)cluster_id;

		}
		else
		{
			cluster_id++;
			C[i] = cluster_id;

		}

	}

	//maps each cluster_id (i.e., node) to a consecutive integer starting from 0
	//this is needed because the matrix L needs to have co-ordinates that are consecutive integers starting from 0
	unordered_map<int, int> node_id_consec;
	int node_id_consec_ind = 0;

	vector<pair<INT, INT> > edges;

	for (INT i = 0; i <= n - d ; ++i)
	{
		//INT cluster_id_pref=C[invSA[i]]; //this is cluster id of prefix
		//INT cluster_id_suff=C[invSA[i+1]]; // this is cluster if of suffix
		edges.push_back(make_pair<INT, INT>((INT) C[invSA[i]], (INT) C[invSA[i + 1]]));

		//cout<<C[invSA[i]]<<" "<<C[invSA[i+1]]<<endl;
	}

	vector<pair<INT, INT> > edges_from_file;

	int line_cnt = 0;
	int source_from_file = -1;
	int target_from_file = -1;

	//stores for each node-id in input_file its corresponding node-id in the graph here (graph nodes need to be consecutive and start from 0


	unordered_map<int, int> consec_nodes;
	int start_cnt = 0;

	for (auto & it_edges : edges)
	{
		int start_node = it_edges.first;
		int end_node = it_edges.second;

		unordered_map<int, int>::const_iterator it = consec_nodes.find(start_node);

		int start_node_consec = -1;

		if (it == consec_nodes.end())
		{
			consec_nodes[start_node] = start_cnt;
			start_node_consec = start_cnt;
			start_cnt++;
		}
		else
			start_node_consec = consec_nodes[start_node];

		it = consec_nodes.find(end_node);
		int end_node_consec = -1;

		if (it == consec_nodes.end())
		{
			consec_nodes[end_node] = start_cnt;
			end_node_consec = start_cnt;
			start_cnt++;
		}
		else
			end_node_consec = consec_nodes[end_node];

		edges_from_file.push_back(make_pair<int, int>((int)start_node_consec, (int)end_node_consec));

	}

	Graph g1(start_cnt, edges_from_file.front().first, edges_from_file.back().second);

	if (DEBUG)
	{
		cout << "source: " << edges_from_file.front().first << endl;
		cout << "target: " << edges_from_file.back().second << endl;

		cout << "Read edges\n";
	}

	map<pair<int, int>, int > edge_multiplicities;
	unordered_map<int, int> node_outdegree;

	for (auto &it : edges_from_file)
	{
		if (DEBUG)cout << it.first << " " << it.second << endl;
		g1.addEdge(it.first, it.second);

		auto it_find = edge_multiplicities.find(make_pair(it.first, it.second));
		if (it_find == edge_multiplicities.end())
			edge_multiplicities[make_pair<int, int>(it.first, it.second)] = 1;
		else
			edge_multiplicities[make_pair<int, int>(it.first, it.second)]++;

		auto it_find2 = node_outdegree.find(it.first);
		if (it_find2 == node_outdegree.end())
			node_outdegree[it.first] = 1;
		else
			node_outdegree[it.first]++;

	}
	double numerator_log_sum = 0.0;
	double denominator_log_sum = 0.0;

	//cout<<"Node outdegrees:\n";
	for (auto &it : node_outdegree)
	{
		//cout<<it.first<<" "<<it.second<<endl;
		if (it.first != edges_from_file.back().second)
		{
			for (int x = 1; x <= it.second - 1; ++x)
				numerator_log_sum += log(x);
		}
		else
		{
			for (int x = 1; x <= it.second; ++x)
				numerator_log_sum += log(x);
		}

	}

	//cout<<"Edge multi:\n";
	for (auto &it : edge_multiplicities)
	{
		//cout<<it.first.first<<","<<it.first.second<<" "<<it.second<<endl;
		for (int x = 1; x <= it.second; ++x)
		{
			denominator_log_sum += log(x);
		}
	}
	if (numerator_log_sum - denominator_log_sum >= log(z))
	{
		cout << "based on ratio d: " << d << " YES" << endl;
		return true;
	}
//  cout<<"Numerator: "<<exp(numerator_log_sum)<<" Denominator: "<<exp(denominator_log_sum)<<endl;

  cout << "\n Current Graph. Number of nodes: " << g1.V << " number of edges: " << edges_from_file.size() << endl;
	g1.quasiMultiChainCompress();
	
       int  numedges = 0;
                for (int i=0; i < g1.V; i++)
                {
                        numedges += g1.adj[i].size();
                }       


	cout << "\n After quasiMultiChainCompression: Current Graph. Number of nodes: " << g1.V << " number of edges: " << numedges << endl;

	if (!g1.isEulerian())
	{
		cout << "Not eulerian: NO\n";
		exit(-1);
	}
	//g1.printGraph();

	g1.SCC_iter();

	//creates stack and adds sccs
	stack<vector<scc > > ST;


	//remove trivial sccs
	for (vector<scc>::iterator it = g1.sccs.begin(); it != g1.sccs.end();)
	{
		//scc x = *it;
		//if (x.nodes.size() == 1)
		if (it->nodes.size() == 1)
		{

			it = g1.sccs.erase(it);

		}
		else
			++it;
	}

	if (g1.sccs.empty())
	{
		cout << "All SCCs of input graph are trivial.\n";
		cout << "One distinct Eulerian tour.\n";
		if (z > 1)
			cout << "z>1, so NO.\n";
		else if (z == 1)
			cout << "z=1, so YES.\n";

		return 0;
	}

	ST.push(g1.sccs);

	int bound = 1;
	for (auto &it : g1.sccs)
	{
		bound *= it.lb;
		//cout<<it.lb<<", ";
	}

	//if(DEBUG)
	//if(d==4)
	cout << "\n bound= " << bound << " " << z << endl;
        double min_bound=bound;

	//fixed seed fro reproducibility
	/*unsigned int seed1=1234567;
	    std::mt19937 gen (seed1);*/

	//replace with the folloiwing if you need different seeds
	std::random_device rd;
	std::mt19937 gen(rd());

	int parameter_t = 0;      //how many times times were popped from ST;
	int parameter_for = 0;    //how many times the for loop in step 11 is executed;
	int unary_iterations=0;

	bool ST_was_empty = false;


        vector<int> bound_vec;
	bound_vec.push_back(bound);

	while (bound < z)
	{

		if (ST.empty())
		{
			//cout<<"NO"<<endl;
			ST_was_empty = true;
			break;
		}

		//GL for debuginning
		if (DEBUG)
		{

			//GL: this is for debugging
			stack<vector<scc > >  tmp_ST(ST);
			cout << " ---- Contents of ST before pop ------" << endl;
			while (!tmp_ST.empty())
			{
				cout << " ------- Vector of sccs -----" << endl;
				vector<scc> c = tmp_ST.top();
				for (auto &it : c)
					print_single_scc(it);
				tmp_ST.pop();
			}
			cout << " --------- End of tmp_ST ------- " << endl;
		}

		vector<scc> f_copy = ST.top();  //get the top element (vector of sccs) and remove it from the list
		vector<scc> f;

		for (int xxx = 0; xxx < f_copy.size(); ++xxx)
		{
			f.push_back(f_copy[xxx]);
		}

		//vector<scc> f(f_copy);  //this is because f_copy is just a reference and is lost after pop
		ST.pop();
		parameter_t++;

		if(parameter_t % 1000==0)
			cout<<parameter_t<<endl;

		if (DEBUG)
		{
			cout << "  --- what I popped --- \n";
			for (auto &it : f)
			{
				print_single_scc(it);
			}
			cout << " ---  end of what I popped --- \n";

			stack<vector<scc > >  tmp_ST(ST);
			cout << " ---- Contents of ST after pop ------" << endl;
			while (!tmp_ST.empty())
			{
				cout << " ------- Vector of sccs -----" << endl;
				vector<scc> c = tmp_ST.top();
				for (auto &it : c)
					print_single_scc(it);
				tmp_ST.pop();
			}
			cout << " --------- End of tmp_ST ------- " << endl;
		}

		//get a random scc from f  with index 0 to f.size()-1

		std::uniform_int_distribution<int> dist(0, f.size() - 1);
		int scc_rand_id = dist(gen);
		if (DEBUG)
			cout << "SCC rand id: " << scc_rand_id << endl;


		int f_bound = 1;
		for (auto &it : f) //f is a vector of scc structs
		{
			//scc x = it;
			//if (DEBUG)
			//cout << "# " << it.nodes.size() << " " << it.lb << endl;
			//f_bound *= x.lb;
			f_bound *= it.lb;
		}
		bound = bound - f_bound;

		/*if(d==4)
		{	cout<<f.size()<<" fb:"<<f_bound<<" b:"<<bound<<endl;

		}*/
		scc x = f[scc_rand_id]; //my current C_i from the vector of sccs f

		if (DEBUG)
			cout << "\n I chose " << scc_rand_id << endl;

		//19-10-2022 addition
		f.erase(f.begin() + scc_rand_id); //remove the C_i you chose randomly above

		//for each distinct out-neighbor of the source node
		unordered_set<int> seen_outneighbors;
		//for(list<int>::const_iterator it=g1.adj[x.source].begin(); it!=g1.adj[x.source].end();++it)

		if(x.edges.size()<=1)
			unary_iterations++;
		 
		for (vector<pair<int, int> >::const_iterator it = x.edges.begin(); it != x.edges.end(); ++it)
		{
			//19-10-22
			vector<scc> f_neighbor(f);


			parameter_for++;

			//if the edge does not start from x.source continue, so that you get only outneighbors of x.source (non-distinct)
			if (it->first != x.source) continue;

			int outneighbor = it->second;


			if (DEBUG)
			{
				cout << "out: " << outneighbor << endl;
				cout << "x.source=" << x.source << " outneighbor: " << outneighbor << endl;
			}
			//if out-neighbor of source in C_i
			//if(g1.node_scc_assoc[*it]==g1.node_scc_assoc[x.source])
			//{
			//find the out-neighbor in seen
			unordered_set<int>::const_iterator it2 = seen_outneighbors.find(outneighbor);

			//if you have not seen the outneighbor before
			if (it2 == seen_outneighbors.end())
			{
				//compute SCC in the graph of the edges of the component except (x.source,*it)
				//this graph has the same number of nodes as the currect C_i
				//set source and target to -1
				if (DEBUG)
				{
					cout << "C\(s,u) node-set size:" << x.nodes.size() << endl;
					cout << "I will delete: " << x.source << "->" << outneighbor << endl;

				}
				Graph cur_C_except_edge(x.nodes.size(), -1, -1);

				unordered_map<int, int> consec;
				int consec_id = 0;

				int multiplicity_of_removed = 0;
				for (auto & it3 : x.edges)
				{
					//cout<<"edge:"<<it3.first<<" "<<it3.second<<endl;
					if (DEBUG)
					{	for (auto &itx : consec)
						{
							cout << "<" << itx.first << " " << itx.second << ">" << endl;
						}
					}
					if (it3.first != x.source || it3.second != outneighbor)			//add everything except edge (x.source,*it)
					{
						unordered_map<int, int>::const_iterator it_cons = consec.find(it3.first);

						unordered_map<int, int>::const_iterator it_cons2 = consec.find(it3.second);

						if (it_cons == consec.end() && it_cons2 == consec.end())	//I didn't find it3.first nor it3.second, I need to put consec_id
						{
							int start_node = 1;
							int end_node = -1;

							consec[it3.first] = consec_id;
							start_node = consec_id;
							consec_id++;

							it_cons2 = consec.find(it3.second); //because the addition of it3.first updates consec and it3.first may be it3.second

							if (it_cons2 == consec.end()) //
							{
								consec[it3.second] = consec_id;
								end_node = consec_id;
								consec_id++;
							}
							else                        //it3.first == it3.second
								end_node = start_node;


							cur_C_except_edge.addEdge(start_node, end_node);

							//cout<<"1had:"<<it3.first<<" "<<it3.second<<" added: "<<start_node<<","<<end_node<<endl;

						}
						else if (it_cons != consec.end() && it_cons2 == consec.end())
						{
							consec[it3.second] = consec_id;
							cur_C_except_edge.addEdge(consec[it3.first], consec_id);
							//cout<<"2had:"<<it3.first<<" "<<it3.second<<"added: "<<consec[it3.first]<<","<<consec_id<<endl;
							consec_id++;
						}
						else if (it_cons == consec.end() && it_cons2 != consec.end())
						{
							consec[it3.first] = consec_id;
							cur_C_except_edge.addEdge(consec_id, consec[it3.second]);
							//cout<<"3had:"<<it3.first<<" "<<it3.second<<"added: "<<consec_id<<","<<consec[it3.second]<<endl;
							consec_id++;
						}
						else if (it_cons != consec.end() && it_cons2 != consec.end())
						{
							cur_C_except_edge.addEdge(consec[it3.first], consec[it3.second]);
							//cout<<"4had:"<<it3.first<<" "<<it3.second<<"added: "<<consec[it3.first]<<","<<consec[it3.second]<<endl;
						}
					}
					if (it3.first == x.source && it3.second == outneighbor)
					{
						multiplicity_of_removed++;
						//cout<<x.source<<"->"<<outneighbor<<" has multiplicity: "<<multiplicity_of_removed<<endl;
					}
				}

				if (multiplicity_of_removed > 1)
				{
					unordered_map<int, int>::const_iterator it_cons = consec.find(x.source);
					unordered_map<int, int>::const_iterator it_cons2 = consec.find(outneighbor);
					int mapped_1 = -1;
					int mapped_2 = -1;

					if (it_cons == consec.end())
					{	mapped_1 = consec_id;
						consec[x.source] = consec_id; //new addition
						consec_id++;
					}
					else mapped_1 = consec[x.source];

					if (it_cons2 == consec.end())
					{
						mapped_2 = consec_id;
						consec[outneighbor] = consec_id; //new addition
						consec_id++;
					}
					else mapped_2 = consec[outneighbor];

					//cout<<mapped_1<<" ! "<<mapped_2<<endl;

					for (int i = 1; i < multiplicity_of_removed; ++i)
					{

						cur_C_except_edge.addEdge(mapped_1, mapped_2);
						//cout<<"had:"<<x.source<<","<<outneighbor<<" added:"<<mapped_1<<","<<mapped_2<<endl;
					}
				}

				cur_C_except_edge.source = outneighbor; //the node *it is the new source of graph C\(s,u)
				cur_C_except_edge.target = x.target; //the target of C_i is the target of C\(s,u)

				if (DEBUG)
				{
					cout << "C_i\(s,u) source :" << x.source << " outneighbor: " << outneighbor << endl;

					cout << "C_i\(s,u) graph:";
					cur_C_except_edge.printGraph();
				}



				/*if (!f.empty())
				{
				cout<<" ------ in !f.empty() BEFORE ---- \n";
				for(auto &it: f)
				{
						print_single_scc(it);

				}
				 cout<<" ------ end in !f.empty() ---- \n";
				cout<<"scc_rand_id: "<<scc_rand_id<<endl;


				//-->					f.erase(f.begin() + scc_rand_id); //remove the C_i you chose randomly above


				if(d==4 && bound>0)
					cout << "Removed: " << scc_rand_id << endl;

				cout<<" ------ in !f.empty() AFTER ---- \n";
				for(auto &it: f)
				{
						print_single_scc(it);

				}
				 cout<<" ------ end in !f.empty() ---- \n";
				cout<<"scc_rand_id: "<<scc_rand_id<<endl;

				}

				//these SCCs will be annotated with source and target nodes belonging only to the nodes of this graph -- maybe wrong
				for(auto &itz : consec)
				{
				cout<<"Map: "<<itz.first<<"->"<<itz.second<<endl;
				}
				*/

				cur_C_except_edge.SCC_map_iter(consec);




				//	cur_C_except_edge.printSccs();

				// add these new components that have size >1 to f with f.push_back();
				for (auto &itx : cur_C_except_edge.sccs)
				{
					//scc x2 = itx;
					//if (x2.nodes.size() > 1)
					if (itx.nodes.size() > 1)
						f_neighbor.push_back(itx);
				}



				if (!f_neighbor.empty())
				{
					ST.push(f_neighbor);
					f_bound = 1;
					for (auto &ity : f_neighbor) //f is a vector of scc structs
					{
						f_bound *= ity.lb;
					}

					bound += f_bound;
				}
				else
				{
					bound++;
				}

				seen_outneighbors.insert(outneighbor);   //add the out-neighbor to seen
				//GL for debugging
				/*if(d==4)
				{
					cout<<"Seen out-neighbors\n";
					for(auto& itg : seen_outneighbors)
						cout<<"["<<itg<<"]";
					cout<<endl;
				}*/

			}// end of  if(it2==seen_outneighbors.end())

			//}  // end of if(g1.node_scc_assoc[*it]==g1.node_scc_assoc[x.source])
		} //end of for(vector<pair<int,int> >::const_iterator it=x.edges.begin(); it!=x.edges.end(); ++it)

		//cout << "it="<<parameter_t<<" bound="<<bound<<" "<< z << endl;
		bound_vec.push_back(bound);

	}
        double mean=0.0;
        auto p=getMeanVariance(bound_vec);
        cout<<"bound mean:"<<p.first<<" stdev: "<<sqrt(p.second)<<" min:"<<min_bound<<" max:"<<bound<<endl;
             
	cout << "parameter t: " << parameter_t << endl;
	//cout << "parameter_for: " << parameter_for << endl;
	cout << "unary iterations: "<<unary_iterations<<endl;
	cout << "bound: " << bound << endl;
	delete []C;
	if (ST_was_empty)
	{
		cout << "d: " << d << " NO\n";
		return false;
	}
	else
	{	cout << "d: " << d << " YES \n\n";
		return true;
	}


	return false;
}

bool dBgraph_tree (unsigned int d, INT *SA, INT *LCP, INT *invSA, unsigned int z, INT n, bool compress)
{    
    
	if (d == 0)
		return true;

	INT cluster_id = (INT) 0;
	INT *C = new INT[n];

	for (INT i = 0; i < n ; i++)
	{
		if (LCP[i] >= d - 1)
		{
			C[i] = (INT)cluster_id;

		}
		else
		{
			cluster_id++;
			C[i] = cluster_id;

		}

	}

	//maps each cluster_id (i.e., node) to a consecutive integer starting from 0
	//this is needed because the matrix L needs to have co-ordinates that are consecutive integers starting from 0
	unordered_map<int, int> node_id_consec;
	int node_id_consec_ind = 0;

	vector<pair<INT, INT> > edges;

	for (INT i = 0; i <= n - d ; ++i)
	{
		//INT cluster_id_pref=C[invSA[i]]; //this is cluster id of prefix
		//INT cluster_id_suff=C[invSA[i+1]]; // this is cluster if of suffix
		edges.push_back(make_pair<INT, INT>((INT) C[invSA[i]], (INT) C[invSA[i + 1]]));

		//cout<<C[invSA[i]]<<" "<<C[invSA[i+1]]<<endl;
	}

	vector<pair<INT, INT> > edges_from_file;

	int line_cnt = 0;
	int source_from_file = -1;
	int target_from_file = -1;

	//stores for each node-id in input_file its corresponding node-id in the graph here (graph nodes need to be consecutive and start from 0


	unordered_map<int, int> consec_nodes;
	int start_cnt = 0;

	for (auto & it_edges : edges)
	{
		int start_node = it_edges.first;
		int end_node = it_edges.second;

		unordered_map<int, int>::const_iterator it = consec_nodes.find(start_node);

		int start_node_consec = -1;

		if (it == consec_nodes.end())
		{
			consec_nodes[start_node] = start_cnt;
			start_node_consec = start_cnt;
			start_cnt++;
		}
		else
			start_node_consec = consec_nodes[start_node];

		it = consec_nodes.find(end_node);
		int end_node_consec = -1;

		if (it == consec_nodes.end())
		{
			consec_nodes[end_node] = start_cnt;
			end_node_consec = start_cnt;
			start_cnt++;
		}
		else
			end_node_consec = consec_nodes[end_node];

		edges_from_file.push_back(make_pair<int, int>((int)start_node_consec, (int)end_node_consec));

	}

	Graph_vec g1(start_cnt, edges_from_file.front().first, edges_from_file.back().second);
  
	if (DEBUG)
	{
		cout << "source: " << edges_from_file.front().first << endl;
		cout << "target: " << edges_from_file.back().second << endl;

		cout << "Read edges\n";
	}

	map<pair<int, int>, int > edge_multiplicities;
	unordered_map<int, int> node_outdegree;

	for (auto &it : edges_from_file)
	{
		if (DEBUG)cout << it.first << " " << it.second << endl;
		g1.addEdge(it.first, it.second);

		auto it_find = edge_multiplicities.find(make_pair(it.first, it.second));
		if (it_find == edge_multiplicities.end())
			edge_multiplicities[make_pair<int, int>(it.first, it.second)] = 1;
		else
			edge_multiplicities[make_pair<int, int>(it.first, it.second)]++;

		auto it_find2 = node_outdegree.find(it.first);
		if (it_find2 == node_outdegree.end())
			node_outdegree[it.first] = 1;
		else
			node_outdegree[it.first]++;

	}
	double numerator_log_sum = 0.0;
	double denominator_log_sum = 0.0;

	//cout<<"Node outdegrees:\n";
	for (auto &it : node_outdegree)
	{
		//cout<<it.first<<" "<<it.second<<endl;
		if (it.first != edges_from_file.back().second)
		{
			for (int x = 1; x <= it.second - 1; ++x)
				numerator_log_sum += log(x);
		}
		else
		{
			for (int x = 1; x <= it.second; ++x)
				numerator_log_sum += log(x);
		}

	}

	//cout<<"Edge multi:\n";
	for (auto &it : edge_multiplicities)
	{
		//cout<<it.first.first<<","<<it.first.second<<" "<<it.second<<endl;
		for (int x = 1; x <= it.second; ++x)
		{
			denominator_log_sum += log(x);
		}
	}
       cout<<"log(numerator): "<<numerator_log_sum<<" log(denominator): "<<denominator_log_sum<<" log(numer)-log(denom):"<<numerator_log_sum-denominator_log_sum<<" log(z):"<<log(z)<<endl;

	if (numerator_log_sum - denominator_log_sum >= log(z))
	{
		cout << "based on ratio d: " << d << " YES" << endl;
		return true;
	}


  if(compress)
  {
	    cout << "\n Current Graph. Number of nodes: " << g1.V << " number of edges: " << edges_from_file.size() << endl;
	  		g1.quasiMultiChainCompress();

		int  numedges = 0;
		for (int i=0; i < g1.V; i++)
		{
			numedges += g1.adj[i].size();
		}	

	  	cout << "\n After quasiMultiChainCompress: Current Graph. Number of nodes: " << g1.V << " number of edges: " << numedges << endl;
  }
  else
     cout << "\n Current Graph. Number of nodes: " << g1.V << " number of edges: " << edges_from_file.size() << endl;

	if (!g1.isEulerian())
	{
		cout << "Not eulerian: NO\n";
		exit(-1);
	}
	//g1.printGraph();


	uint64_t t0 = timeMs(); // START OF ALGORITHM


	g1.SCC_iter();

	cout<<"# of SCCs: "<<g1.sccs.size()<<endl;
	
	//creates stack and adds sccs
	stack<vector<scc_vec > > ST;

	int total_size_of_sccs=0;
	//remove trivial sccs
	for (int iscc = 0; iscc < g1.sccs.size(); iscc++)
	{
		scc_vec *xx = g1.sccs[iscc];
		if (xx->nodes.size() == 1)
		{
			// delete &(*it); // TODO -- how do we delete them?
			delete g1.sccs[iscc];
			g1.sccs.erase(g1.sccs.begin() + iscc);
			iscc--;
			
		}
		else
			total_size_of_sccs=+xx->nodes.size();
	}
	cout<<"Total size of sccs: "<<total_size_of_sccs<<endl;
	if (g1.sccs.empty())
	{
		cout << "All SCCs of input graph are trivial.\n";
		cout << "One distinct Eulerian tour.\n";
		if (z > 1)
			cout << "z>1, so NO.\n";
		else if (z == 1)
			cout << "z=1, so YES.\n";

		return 0;
	}


	// ========================================================
	// 					TREE / ID2SCC
	// ========================================================
	// ID counter for SCC, 0 used for internal nodes if necessary (although internal nodes can keep the id of the component they represented originally)
	// curr is assumed unused (incremented after each use)
	unsigned currSccID = 1;

	// SUM/PROD tree : internal nodes have SUM or PROD xtype, and leaves LEAF type,
	// and they contain the ID of an SCC in the map
	// futhermore each node has the LB of its subtree (which is kept updated)
	xtnode * TREE = new xtnode(PROD, 777);
	TREE->lb = 1;

	//vector of sccs
	vector<scc_vec *> SCCS;

	xtnode * xtemp = NULL;

	for (auto &c : g1.sccs) {
		c->id = currSccID;
		xtemp = new xtnode(LEAF, currSccID, c->lb, TREE); // TODO : compute the lb with lb_fun ?
		currSccID++; //curr must be left to an unused value

		c->xn = xtemp; //link scc to its xtnode

		TREE->children.push_back(xtemp); // put in the tree

		TREE->lb = TREE->lb * c->lb;

		SCCS.push_back(c);

		if (DEBUG) cout << "SCC " << c->id << " has LB : " << c->lb << endl;
	}

	//if (DEBUG) 
		cout << "TREE has LB : " << TREE->lb << endl;

	// ========================================================
//fixed seed fro reproducibility
	/*unsigned int seed1=1234567;
	      std::mt19937 gen (seed1);*/

	//replace with the folloiwing if you need different seeds
	std::random_device rd;
	std::mt19937 gen(rd());

	bool MAP_is_empty = false;
	unsigned while_count = 0;
	// ============= NEW STATS CODE ===============  counter
	unsigned unaryIterations = 0;
	// ============= END NEW STATS CODE ===============  

    int min_bound=TREE->lb;
    vector<int> bound_vec;
    bound_vec.push_back(TREE->lb);

	while (TREE->lb < z)
	{
		if (DEBUG) cout << "============= ITER " << while_count << endl;
		while_count++;

		if (SCCS.empty())
		{
			//cout<<"NO"<<endl;
			MAP_is_empty = true;
			break;
		}

		// vector<scc> f_copy=ST.top();    //get the top element (vector of sccs) and remove it from the vector
		// vector<scc> f(f_copy);  //this is because f_copy is just a reference and is lost after pop
		// ST.pop();

		//get a random scc from f  with index 0 to f.size()-1



// SELECT SCC FROM TREE/MAP (TODO finally, we want this to be random or subject to strategies)
		// std::uniform_int_distribution<int> dist(0, f.size()-1);
		// int scc_rand_id=dist(gen);

// useful for random extraction:
		// std::uniform_int_distribution<int> dist(0, f.size()-1);
		// int scc_rand_id=dist(gen);
		// if(DEBUG)
		// cout<<"SCC rand id: "<<scc_rand_id<<endl;


		// scc  x=f[scc_rand_id]; //my current C_i from the vector of sccs f

		scc_vec * x = NULL;
		xtnode * node_of_x = NULL;
		int id_to_remove = -1;

		// SMALLEST NODE EXTRACTION
		int minval = INT32_MAX;
		int pos = 0;
		for (int i = 0; i < SCCS.size(); i++) {
			if (SCCS[i]->nodes.size() < minval) {
				id_to_remove = SCCS[i]->id;
				x = SCCS[i];
				node_of_x = x->xn;
				pos = i;
			}
		}
		// id_to_remove = p->first;
		// x = p->second;
		// node_of_x = x->xn;

		if (DEBUG || node_of_x->type != LEAF) {
			cout << "[FATAL] x (" << x->id << " / " << node_of_x->id << ") is not a leaf node (type " << node_of_x->type << ")\n";

			printXTREE(TREE);

		}

		// must remove SCC from MAP and changed in TREE
		SCCS.erase(SCCS.begin() + pos);

		if (DEBUG) cout << "DELETED FROM VEC\n";

		// xtnode of x needs to be turned from leaf to SUM (may be compressed if it has only 1 child)
		// TODO free scc x?
		node_of_x->type = SUM;

		unsigned old_x_lb = node_of_x->lb; // for checking, the result should never decrease
		node_of_x->lb = 0;

		if (DEBUG)
			cout << "\n I chose " << x->id << "(source:" << x->source << " target:" << x->target << " lb:" << x->lb << ")" << endl; // sou ta lb ispath


		//for each distinct out-neighbor of the source node
		unordered_set<int> seen_outneighbors;
		//for(vector<int>::const_iterator it=g1.adj[x->source].begin(); it!=g1.adj[x->source].end();++it)
		if (DEBUG) {
			for ( auto& E : x->edges) {
				cout << "EDGE " << E.first << " " << E.second << endl;
			}
		}

		// ================================
		// We are exploring all outedges of the source
		// ================================
		unsigned nonTrivialChildrenOfX = 0;

		for (const auto& it : x->edges)
		{
			if (DEBUG) cout << "PAIR : " << it.first << " - " << it.second;
			//if the edge does not start from x->source continue, so that you get only outneighbors of x->source (non-distinct)
			if (it.first != x->source) {
				if (DEBUG) cout << " -not from source-\n";
				continue;
			} else {
				if (DEBUG) cout << " -----------> PROCESSING\n";
			}

			int outneighbor = it.second; // THIS is one possible step of the tour, the for should try all

			if (DEBUG)
				cout << "x->source=" << x->source << " outneighbor: " << outneighbor << endl;

			//if out-neighbor of source in C_i
			//if(g1.node_scc_assoc[*it]==g1.node_scc_assoc[x->source])
			//{
			//find the out-neighbor in seen
			unordered_set<int>::iterator it2 = seen_outneighbors.find(outneighbor);

			//if you have not seen the outneighbor before
			if (it2 == seen_outneighbors.end())
			{
				//cout<<"seen_outneighbrs.size():"<<seen_outneighbors.size()<<" I didn't find:"<<outneighbor<<" "<<endl;

				//compute SCC in the graph of the edges of the component except (x->source,*it)
				//this graph has the same number of nodes as the currect C_i
				//set source and target to -1
				if (DEBUG)
				{
					cout << "C\\(s,u) node-set size:" << x->nodes.size() << endl;
					cout << "I will delete: " << x->source << "->" << outneighbor << endl;
				}
				Graph_vec cur_C_except_edge(x->nodes.size(), -1, -1);

				unordered_map<int, int> consec;
				int consec_id = 0;

				int multiplicity_of_removed = 0;
				for (const auto & it3 : x->edges)
				{
					if (it3.first != x->source || it3.second != outneighbor)			//add everything except edge (x->source,*it)
					{
						unordered_map<int, int>::const_iterator it_cons = consec.find(it3.first);

						unordered_map<int, int>::const_iterator it_cons2 = consec.find(it3.second);

						if (it_cons == consec.end() && it_cons2 == consec.end())	//I didn't find it3.first nor it3.second, I need to put consec_id
						{
							int start_node = 1;
							int end_node = -1;

							consec[it3.first] = consec_id;
							start_node = consec_id;
							consec_id++;

							it_cons2 = consec.find(it3.second); //because the addition of it3.first updates consec and it3.first may be it3.second

							if (it_cons2 == consec.end()) //
							{
								consec[it3.second] = consec_id;
								end_node = consec_id;
								consec_id++;
							}
							else                        //it3.first == it3.second
								end_node = start_node;


							cur_C_except_edge.addEdge(start_node, end_node);

							if (DEBUG) cout << "1had:" << it3.first << " " << it3.second << " added: " << start_node << "," << end_node << endl;
						}
						else if (it_cons != consec.end() && it_cons2 == consec.end())
						{
							consec[it3.second] = consec_id;
							cur_C_except_edge.addEdge(consec[it3.first], consec_id);
							//cout<<"had:"<<it3.first<<" "<<it3.second<<"added: "<<consec[it3.first]<<","<<consec_id<<endl;
							consec_id++;
						}
						else if (it_cons == consec.end() && it_cons2 != consec.end())
						{
							consec[it3.first] = consec_id;
							cur_C_except_edge.addEdge(consec_id, consec[it3.second]);
							//cout<<"had:"<<it3.first<<" "<<it3.second<<"added: "<<consec_id<<","<<consec[it3.second]<<endl;
							consec_id++;
						}
						else
						{
							cur_C_except_edge.addEdge(consec[it3.first], consec[it3.second]);
							//cout<<"had:"<<it3.first<<" "<<it3.second<<"added: "<<consec[it3.first]<<","<<consec[it3.second]<<endl;
						}
					}
					if (it3.first == x->source && it3.second == outneighbor)
					{
						multiplicity_of_removed++;
						//cout<<x->source<<"->"<<outneighbor<<" has multiplicity: "<<multiplicity_of_removed<<endl;
					}
				}

				if (multiplicity_of_removed > 1)
				{
					unordered_map<int, int>::const_iterator it_cons = consec.find(x->source);
					unordered_map<int, int>::const_iterator it_cons2 = consec.find(outneighbor);
					int mapped_1 = -1;
					int mapped_2 = -1;

					if (it_cons == consec.end())
					{	mapped_1 = consec_id;
						consec_id++;
					}
					else mapped_1 = consec[x->source];

					if (it_cons2 == consec.end())
					{
						mapped_2 = consec_id;
						consec_id++;
					}
					else mapped_2 = consec[outneighbor];

					for (int i = 1; i < multiplicity_of_removed; ++i)
					{

						cur_C_except_edge.addEdge(mapped_1, mapped_2);
						//cout<<"had:"<<x->source<<","<<outneighbor<<" added:"<<mapped_1<<","<<mapped_2<<endl;
					}
				}

				cur_C_except_edge.source = outneighbor; //the node *it is the new source of graph C\(s,u)
				cur_C_except_edge.target = x->target; //the target of C_i is the target of C\(s,u)
				//cout << "C_i\(s,u) source :" << x->source << " outneighbor: " << outneighbor << " consec[x.source] " << consec[x->source] << " consec[outneighbor] " << consec[outneighbor] << endl;

				//cout<<"C_i\(s,u) graph:";
				//cur_C_except_edge.printGraph();

				//these SCCs will be annotated with source and target nodes belonging only to the nodes of this graph -- maybe wrong
				/*for(auto &itz : consec)
				{
					cout<<"Map: "<<itz.first<<"->"<<itz.second<<endl;
				}*/

				cur_C_except_edge.SCC_map_iter(consec);

				/*for (int i = 0; i < cur_C_except_edge.sccs.size(); ++i)
				{
					cout << "scc has nodes: " << cur_C_except_edge.sccs[i]->nodes.size() << " and edges: " << cur_C_except_edge.sccs[i]->edges.size() << endl;
				}*/

				// ==================================================
				// making a CHILD NODE of the SUM parent node
				// ==================================================

				// NOW we made a possible step (one element of the father SUM xtnode)
				// this may have multiple SCCs, which go in a PROD xtnode, but if just one, then it can be a leaf

				if (cur_C_except_edge.sccs.size() == 1) { //////////////////////// ONLY 1 SCC, this is a LEAF node

					// scc * childSCC = (cur_C_except_edge.sccs[0]).clone();
					scc_vec * childSCC = cur_C_except_edge.sccs[0];

					if (childSCC->nodes.size() <= 1) { // trivial SCC, just add 1 to lb of parent (which is a SUM node)
						node_of_x->lb++;
						// TODO erase childSCC?
						delete childSCC;
					} else { // proper SCC
						childSCC->id = currSccID;
						xtnode * childNode = new xtnode(LEAF, currSccID, childSCC->lb, node_of_x);
						childSCC->xn = childNode;

						SCCS.push_back(childSCC); //INSERT IN VEC
						currSccID++;

						node_of_x->children.push_back(childNode); // insert in TREE
						node_of_x->lb += childSCC->lb;

					}
				}
				else if (cur_C_except_edge.sccs.size() > 1) { ///////////////////////////// MULTIPLE SCCs, PROD node in the "middle" with each SCC as child

					xtnode * middlePNode = new xtnode(PROD, 0, 1, node_of_x);
					node_of_x->children.push_back(middlePNode);
					if (DEBUG) cout << " --- FOR on " << cur_C_except_edge.sccs.size() << " SCCs\n";

					scc_vec * childSCC = NULL;

					for (int i = 0; i < cur_C_except_edge.sccs.size(); i++)
					{
						// childSCC = (cur_C_except_edge.sccs[i]).clone();
						childSCC = cur_C_except_edge.sccs[i];

						if (childSCC->nodes.size() <= 1) { // trivial SCC, do nothing since parent is a PROD node
							// TODO erase childSCC?
							delete childSCC;
							if (DEBUG) cout << "-trivial SCC\n";
						} else { // proper SCC
							childSCC->id = currSccID;
							xtnode * childNode = new xtnode(LEAF, currSccID, childSCC->lb, middlePNode);
							childSCC->xn = childNode;

							SCCS.push_back(childSCC);
							currSccID++;

							middlePNode->children.push_back(childNode); // insert in TREE under middlePNode
							middlePNode->lb *= childSCC->lb;

						}
					}

					node_of_x->lb += middlePNode->lb; //update LB of node_of_x
				}

				// TODO check if nodes generated were real or should be "CLOSED";

				// VVV TREE REMOVAL
				// 					xtnode * parent_of_x = node_of_x->parent;
				// if(parent_of_x == NULL){
				// 	cout << "[FATAL] chosen x to be removed is the root of XTREE (has NULL parent)";
				// 	return -1;
				// }
				// 	remove(parent_of_x->children.begin(),)

				// ==================================================

				//cout<<"I add: "<<outneighbor<<endl;
				seen_outneighbors.insert(outneighbor);   //add the out-neighbor to seen
				//cout<<"seen_outneighbors.size():"<<seen_outneighbors.size()<<endl;

			}// end of  if(it2==seen_outneighbors.end())
			//}  // end of if(g1.node_scc_assoc[*it]==g1.node_scc_assoc[x->source])
		} //end of for(vector<pair<int,int> >::const_iterator it=x->edges.begin(); it!=x->edges.end(); ++it)

		// ============= NEW STATS CODE ===============  checking whether iteration generated at least 2 nodes
		if(DEBUG) cout << "Iteration generated " << node_of_x->children.size() << "nodes\n";
		if( node_of_x->children.size() <= 1) unaryIterations++;
		// ============= END NEW STATS CODE ===============  
		
		// TODO check case where internal node remains without leaves
		// TODO also I have seen SUM nodes with lb=0

		if (DEBUG) cout << "Time to update TREE LB! \n";
		// ==================================================
		// UPDATE LB! node_of_x has the right LB, needs to be propagated up until root via PROD and SUM
		// ==================================================

		if (DEBUG) printXTREE(TREE);

		xtnode * tempNode = node_of_x;
		unsigned old_bound = old_x_lb;
		unsigned pb, old_pb;

		if (DEBUG) cout << "TYPES: SUM " << SUM << " PROD " << PROD << " LEAF " << LEAF << endl;
		if (DEBUG) cout << "node_of_x id: " << node_of_x->id << " type: " << node_of_x->type << endl;

		if (DEBUG) {
			cout << "=================! children SIZE: ";
			cout << node_of_x->children.size() << endl;
			cout << "children [id, type]: ";
			for (xtnode* ch : node_of_x->children) { // ch : xtnode*
				cout << "[" << ch->id << "," << ch->type << "] ";
			}
			cout << endl;
		}
		// int dbgcount = 0;

		while (tempNode->parent != NULL) { // update the bound of parent

			// dbgcount++;
			// if(dbgcount > 10) return 0;
			xtnode * parentNode = tempNode->parent;
			old_pb = parentNode->lb;

			if (DEBUG) cout << "id: " << tempNode->id << " type: " << tempNode->type << " | parent id: " << parentNode->id << " type: " << parentNode->type << endl;
			if (parentNode->type == PROD) {
				pb = (old_pb / old_bound) * tempNode->lb;
				parentNode->lb = pb;
				old_bound = old_pb;
				tempNode = parentNode;
			} else if (parentNode->type == SUM) {
				pb = (old_pb - old_bound) + tempNode->lb;
				parentNode->lb = pb;
				old_bound = old_pb;
				tempNode = parentNode;
			} else { // ERROR!!
				cout << "== [FATAL] INTERNAL NODE MARKED LEAF, aborting ==\n";
				cout << "ID " << parentNode->id << " TYPE " << parentNode->type << endl;
				return 1;
			}
		}
		if (DEBUG) cout << "=====================================\n";
		if (DEBUG) cout << "CURRENT BOUND : " << TREE->lb << endl;
		if (DEBUG) cout << "=====================================\n";

		/////////////TODO MUST delete x (old scc of expanded node) from memory
		delete x;


		if (DEBUG) {
			printXTREE(TREE);

		}

		bound_vec.push_back(TREE->lb);
	} // end of while(lb < z)

	//GREG
    double mean=0.0;
    auto p=getMeanVariance(bound_vec);
    cout<<"bound mean:"<<p.first<<" stdev: "<<sqrt(p.second)<<" min:"<<min_bound<<" max:"<<TREE->lb<<endl;

	if (while_count == 0) while_count++;

	uint64_t t1 = timeMs(); // end of algorithm


	bool answer = false;

	// ============= NEW STATS CODE ===============  updated prints
	cout << "ANSWER, LB, Z, ITERATIONS, LB/ITERATIONS, time (ms), UNARY ITERATIONS\n";

	if (MAP_is_empty)
	{	cout << "\n\nNO, " << TREE->lb << ", " << z << ", " << while_count << ", " << TREE->lb / while_count << ", " << (t1 - t0) << "," << unaryIterations << "\n";

	}
	else
	{	cout << "\n\nYES, " << TREE->lb << ", " << z << ", " << while_count << ", " << TREE->lb / while_count << ", " << (t1 - t0) << ", " << unaryIterations << "\n";
		// cout<<"\n\n YES! (current bound: "<< TREE->lb <<") \n\n";
		answer = true;
	}
	// ============= END NEW STATS CODE ===============  


	////////////////////// clear everything
	for (auto iii : SCCS) {
		delete iii;
	}

	if (SCCS.empty()) cout << "VEC IS EMPTY" << endl;

	delete_tree(TREE);

	if (answer == true)
		return true;
	else
		return false;

	/*	cout<<"Copy: \n";
		Graph g2(g1.getV(),g1.getSource(),g1.getTarget());
		g2.copy(g1);
		g2.printGraph();
	*/



	return false;

}

void euler_path(string &x, int d, string output, INT &line_cnt)
{
	INT n = x.length();

	vector<string> B_vec;
	unordered_map<string, int> B_map;


	INT id = 0;

	for (int m = 0; m <= n - d; m++)
	{
		string pref = x.substr(m, d - 1);
		string suff = x.substr(m + 1, d - 1);


		B_vec.push_back(pref);
		B_vec.push_back(suff);


		unordered_map<string, int>::const_iterator it = B_map.find(pref);
		if (it == B_map.end())
		{
			B_map.insert(std::pair<string, int>(pref, id));
			++id;
		}

		it = B_map.find(suff);
		if (it == B_map.end())
		{
			B_map.insert(std::pair<string, int>(suff, id));
			++id;
		}

	}



	vector< vector<int> > adj;

	adj.resize(id);

	unordered_map<int, string> reverse_B_map;

	for (vector<string>::const_iterator it = B_vec.begin(); it != B_vec.end(); ++it)
	{
		auto it2 = it;
		++it2;
		if (it2 != B_vec.end())
		{
			adj[B_map[*it]].push_back(B_map[*it2]);

			reverse_B_map.insert(std::pair<int, string>(B_map[*it], *it));
			reverse_B_map.insert(std::pair<int, string>(B_map[*it2], *it2));
			++it;
		}
	}
	cout << endl;


	printCircuit(adj, reverse_B_map, output, line_cnt);
}


void printCircuit(vector< vector<int> > adj, unordered_map<int, string> reverse_B_map, string output, INT &line_cnt)
{

	unordered_map<int, int> edge_count;

	for (int i = 0; i < adj.size(); i++)
	{

		edge_count[i] = adj[i].size();
	}

	if (!adj.size())
		return;
	stack<int> curr_path;

	vector<int> circuit;

	curr_path.push(0);
	int curr_v = 0;

	while (!curr_path.empty())
	{

		if (edge_count[curr_v])
		{
			curr_path.push(curr_v);
			int next_v = adj[curr_v].back();

			edge_count[curr_v]--;
			adj[curr_v].pop_back();
			curr_v = next_v;
		}


		else
		{
			circuit.push_back(curr_v);
			curr_v = curr_path.top();
			curr_path.pop();
		}
	}


	ofstream output_file;
	output_file.open(output, std::ios_base::app);

	output_file << ">" << line_cnt << endl;

	for (int i = circuit.size() - 1; i >= 0; i--)
	{
		if (i != circuit.size() - 1)
		{

			output_file << reverse_B_map[circuit[i]].substr(reverse_B_map[circuit[i]].length() - 1, 1);
		}
		else
		{
			output_file << reverse_B_map[circuit[i]];
		}

	}

	output_file << endl;
	output_file.close();
}


