#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"tripod.h"

int main(){
	struct bfs_struct b;
	struct rmq_struct r;
	struct tripod_decomposition_struct t;

	init(&b, &r, &t);

	decompose(&b, &r, &t);

	/*
	printf("vertex_tripod_assign = [ ");
	for (int k = 0; k < ((&b)->v); k++){
		printf("%d ", ((&t)->vertex_tripod_assign)[k]);
	}
	printf("]\n");
	*/
	//instead of printing this to the screen, we should print them to files with relevant names
	//we could output:
	//assignment of vertices to tripods (vertex_tripod_assign) -> 1 int per line
	//the relations between tripods (tripod_adjacency_list) -> 3 ints per line
	//bfs tree -> 1 int per line
	//tripods -> groups of three lines with many ints (vertices) per line

	for (int k = 0; k < ((&b)->v); k++){
		if (((&t)->vertex_tripod_assign)[k] == -1){
			printf("error. some vertices are not labeled as belonging to a tripod.\n");
			exit(0);
		}
	}

	if (!three_tree_test_pt2(&b, &t)) printf("three tree test unsuccessful\n");
	
	tripod_free(&b, &t);
	LCA_free(&r);
	BFS_free(&b);
}