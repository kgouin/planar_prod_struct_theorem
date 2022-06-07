#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"lca.h"

int main(){
	//complete binary tree creation (where the number of nodes equals 2^height-1)
	int h = 20;
	int m = ((1<<h)-1);

	clock_t start;

	printf("Creating %d-node adjacency list representation of binary tree...", m);
	fflush(stdout);
    start = clock();
	int** a = (int**)calloc(m,sizeof(int*));
	for (int i = 0; i < m; i++){
		a[i] = (int*)calloc(3,sizeof(int));
	}
	//fill in parents
	a[0][0] = -1;
	for (int i = 1; i < m; i++){
		a[i][0] = (int)ceil((double)i/2) - 1;
	}
	//fill in children
	for (int i = 0; i < m; i++){
		if (i < ((m-1)/2)){
			a[i][1] = i+i+1;
			a[i][2] = i+i+1+1;
		}
		else{
			a[i][1] = -1;
			a[i][2] = -1;
		}
	}
	printf("done (%.4fs)\n", ((double)clock()-start)/CLOCKS_PER_SEC);

	/*for (int i = 0; i < m; i++){
		for (int j = 0; j < 3; j++){
			if (a[i][j] == -1) printf("%d  ", a[i][j]);
			else printf("%d   ", a[i][j]);
		}
		printf("\n");
	}*/

	struct rmq_struct rs;
	printf("Creating RMQ structure for LCA...");
	fflush(stdout);
	LCA_init(&rs, a, ((2*(m))-1));
	printf("done (%.4fs)\n", ((double)clock()-start)/CLOCKS_PER_SEC);
	srand(time(NULL));

	printf("Performing %d queries...", m);
	fflush(stdout);
	int k = 0;
	while (k < m) {
		int i = (int)(((double)m/RAND_MAX) * rand());
		int j = (int)(((double)m/RAND_MAX) * rand());

		LCA_query(&rs, i, j);
		k++;
	}
	printf("done (%.4fs)\n", ((double)clock()-start)/CLOCKS_PER_SEC);

	//cleanup
	LCA_free(&rs);
	for (int i = 0; i < m; i++){
		free(a[i]);
	}
	free(a);
}























