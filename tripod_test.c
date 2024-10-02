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

	printf("vertex_tripod_assign = [ ");
	for (int k = 0; k < ((&b)->v); k++){
		printf("%d ", ((&t)->vertex_tripod_assign)[k]);
	}
	printf("]\n");

	BFS_free(&b);
	LCA_free(&r);
	tripod_free(&t);
}