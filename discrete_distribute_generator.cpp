#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
using namespace std;
//preprocessing:
//sort input
//build n(leaves) + n-1(internal) length array

/*algorithm:
1. firstly, take two elements in array, insert to nth element,
2. compare (nth element and 3th element) and (3th, 4th element),insert the smaller one into n+1th element
3. general case, compare (first two in frequency marker), (first one in frequency marker and first one in tree marker) and ( first two in tree marker)
4. when frequency marker == n, leaves are finished
5.

*/

/*
data structure
1. number
2. strings or positions of children
*/
typedef unsigned int indx;
struct node
{
	unsigned int number;
	union{
		// unsigned int children[2];//left and right
		struct{indx left;indx right;} children;
		char kmer[20];
	} u;
};
void make_node(node nodes[], indx* next_leaves, indx* next_internal, indx* last_internal){
	if(next_internal - last_internal <= 1){
		unsigned int comp1 = nodes[*last_internal].number + nodes[*next_leaves].number;
		unsigned int comp2 = nodes[*next_leaves].number + nodes[*next_leaves + 1].number;	
		if(comp1>comp2){
			//transfer the struct
			node temp;
			temp.number = comp2;
			temp.u.children.left = *next_leaves;
			temp.u.children.right = *next_leaves + 1;
			nodes[*next_internal] = temp;
			// nodes[*next_internal].number = comp2;
			// nodes[*next_internal].u.children[0] = *next_leaves;
			// nodes[*next_internal].u.children[1] = *next_leaves + 1;
			(*next_internal)++;
			*next_leaves = *next_leaves + 2;

		}	
		else{
			nodes[*next_internal].number = comp1;
			nodes[*next_internal].u.children.left = *next_leaves;
			nodes[*next_internal].u.children.right = *last_internal;
			(*next_internal)++;
			(*next_leaves)++;
			(*last_internal)++;		
		}
	}
	else{
		unsigned int comp1 = nodes[*next_leaves].number + nodes[*next_leaves + 1].number;
		unsigned int comp2 = nodes[*next_leaves].number + nodes[*last_internal].number;
		unsigned int comp3 = nodes[*last_internal].number + nodes[*last_internal + 1].number;	
		if(comp1 > comp2){
			if(comp2 > comp3){
					//take comp3
				nodes[*next_internal].number = comp3;
				nodes[*next_internal].u.children.left = *last_internal;
				nodes[*next_internal].u.children.right = *last_internal + 1;
				*last_internal = *last_internal + 2;
				(*next_internal)++;
			}
			else{
				//take comp2
				nodes[*next_internal].number = comp2;
				nodes[*next_internal].u.children.left = *next_leaves;
				nodes[*next_internal].u.children.right = *last_internal;
				(*last_internal)++;
				(*next_leaves)++;
				(*next_internal)++;
			}
		}
		else if(comp1 < comp3){
			//take comp1
			nodes[*last_internal].number = comp1;
			nodes[*last_internal].u.children.left = *next_leaves;
			nodes[*last_internal].u.children.right = *next_leaves + 1;
			*next_leaves = *next_leaves + 2;
			(*next_internal)++;
		}
		else{
			//take comp3
			nodes[*next_internal].number = comp3;
			nodes[*next_internal].u.children.left = *last_internal;
			nodes[*next_internal].u.children.right = *last_internal + 1;
			*last_internal = *last_internal + 2;
			(*next_internal)++;				
		}
	}
}
void make_node(node nodes[], indx* next_internal, indx* last_internal){
	nodes[*next_internal].number = nodes[*last_internal].number + nodes[*last_internal+1].number;
	nodes[*next_internal].u.children.left = *last_internal;
	nodes[*next_internal].u.children.right = *last_internal + 1;
	*last_internal = *last_internal + 2;
	(*next_internal)++;
}
void huffman_tree(node nodes[],int n){
	//n is the number of leavs node
	//first case
	indx next_leaves = 0;
	indx last_internal = n;
	indx next_internal = n;
	nodes[n].number = nodes[0].number + nodes[1].number;
	nodes[n].u.children.left = 0;
	nodes[n].u.children.right = 1;
	next_leaves = 2;
	next_internal = n + 1;
	indx* next_l = &next_leaves;
	indx* next_i = &next_internal;
	indx* last_i = &last_internal;
	//last and next
	//general case
	while((*next_l)<n-1){//when it is necessary to examine leaves part
		make_node(nodes, next_l, next_i, last_i);
	}
	while((*next_l) == n-1){

		if(*next_i - *last_i > 1){
			unsigned int comp1 = nodes[n-1].number + nodes[*last_i].number;
			unsigned int comp2 = nodes[*last_i].number + nodes[*last_i+1].number;
			if(comp1 > comp2){
				//take comp2
				nodes[*next_i].number = comp2;
				nodes[*next_i].u.children.left = last_internal;
				nodes[*next_i].u.children.right = last_internal + 1;
				*last_i = *last_i + 2;
				(*next_i)++;
			}
			else{
				//take comp1
				nodes[*next_i].number = comp1;
				nodes[*next_i].u.children.left = n-1;
				nodes[*next_i].u.children.right = *last_i;
				(*last_i)++;
				(*next_i)++;
				(*next_l)++;
			}
		}
		else{
			nodes[*next_i].number = nodes[*next_l].number + nodes[*last_i].number;
			nodes[*next_i].u.children.left = *last_i;
			nodes[*next_i].u.children.right = *next_l;
			(*next_l)++;			
			(*next_i)++;
			(*last_i)++;
		}

	}
	while((next_internal < (2*n - 1)) && (next_internal - last_internal > 1)){
		make_node(nodes, next_i,last_i);
	}
}
indx generator(unsigned int random, node nodes[], indx x, indx y, indx n){
	if(random<=nodes[x].number){
		if(x < n){
			return x;
		}
		else{
			return generator(random, nodes, nodes[x].u.children.left, nodes[x].u.children.right,n);
		}
	}
	else{
		if(y < n){
			return y;
		}
		else{
			return generator(random - nodes[x].number, nodes, nodes[y].u.children.left, nodes[y].u.children.right,n);
		}
	}
}
int main(){ // so far, main is used for testing only
	int number_of_lines = 0;
	string line;
	ifstream myfile("/Users/jianhualee/Documents/Project/bioproject/toy_data-counts_and frequency_counts/reads-sister_01-sorted_by_count.txt");
	if(myfile.is_open()){
		while(!myfile.eof()){
			getline(myfile,line);
			number_of_lines++;
		}
	}
	cout << "number of line is: " << number_of_lines << endl;

	indx test_number = 6;
	node* nodes = new node[2*test_number - 1];//1165956
	int i = 0;
	string line1;
	ifstream myfile1("/Users/jianhualee/Documents/Project/bioproject/test.txt");
	// ifstream myfile1("/Users/jianhualee/Documents/Project/bioproject/toy_data-counts_and frequency_counts/reads-sister_01-sorted_by_count.txt");
	if(myfile1.is_open()){
		while(i<test_number){
			getline(myfile1,line1);

			nodes[i].number =  stoi(line1.substr(21,3),0,0);
			for(int j=0;j<test_number;j++){
				nodes[i].u.kmer[j] = line1[j];
			}
			i++;
		}
	}
	try{
		for(int i=0;i<test_number-1;i++){
			if(nodes[i].number > nodes[i+1].number){
				throw invalid_argument("not sorted");
			}
		}//raise exception here
	}
	catch(const invalid_argument& e){
		cout << "not sorted" << endl;
		return 0;
	}
	huffman_tree(nodes,test_number);
	for(int i=0;i<(2*test_number-1);i++){
		cout << nodes[i].number << endl;
	}
	cout << "generate sample" << endl;
	indx result;
	unsigned int w = 7;
	cout << nodes[2*test_number-3].u.children.left << " " << nodes[2*test_number-3].u.children.right << endl;
	result = generator(w, nodes, 2*test_number-3, 2*test_number-4, test_number);
	cout << result << endl;
	return 0;
}

