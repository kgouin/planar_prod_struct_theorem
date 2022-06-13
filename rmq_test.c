#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"rmq.h"

int main(){
	int k = 10000000;
	int l = 0;
	struct rmq_struct s1;
	clock_t start;
	double elapsed;

	printf("Creating input %d-element input array...", k);
	fflush(stdout);
	start = clock();
	srand(time(NULL));
	s1.n = k;
	s1.d = (int*)malloc(k * sizeof(int));
	s1.d[0] = 0;
	for (int i = 1; i < k; i++){
		if ((rand() % 2) == 0) s1.d[i] = (s1.d[i-1]) + 1;
		else s1.d[i] = (s1.d[i-1]) - 1;
	}
	printf("done (%.4fs)\n", ((double)clock()-start)/CLOCKS_PER_SEC);

	printf("Creating RMQ structure...\n");
	fflush(stdout);
	RMQ_init(&s1);
	printf("done (%.4fs)\n", ((double)clock()-start)/CLOCKS_PER_SEC);

	printf("Block size = %d\n", s1.b);

	printf("Performing %d queries...\n", k);
	fflush(stdout);
	while (l < k) {
		int j = (int)(((double)k/RAND_MAX) * rand());
		int i = (int)(((double)(j)/RAND_MAX) * rand());
		
		/*if (RMQ_query(&s1, i, j) != RMQ_simple(&s1, i, j)){
			printf("ERROR\n");
			printf("RMQ_query = %d, RMQ_simple = %d\n", RMQ_query(&s1, i, j), RMQ_simple(&s1, i, j));
			printf("i = %d, j = %d\n", i, j);
			break;
		}*/

		RMQ_query(&s1, i, j);
		l++;
	}
	printf("done (%.4fs)\n", ((double)clock()-start)/CLOCKS_PER_SEC);
	RMQ_free(&s1);
}
