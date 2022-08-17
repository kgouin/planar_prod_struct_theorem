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

	BFS_free(&b);
	LCA_free(&r);
}