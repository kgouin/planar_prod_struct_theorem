#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"bfs.h"

int main(){
	struct bfs_struct s;
	BFS_init(&s);
	BFS(&s);
	printf("\n");
	printf("bt:\n");
	for (int i = 0; i < (s.v); i++){
		printf("%d ", s.bt[i]);
	}
	printf("\n\n");
	printf("ct:\n");
	for (int i = 0; i < (s.f); i++){
		for (int j = 0; j < 3; j++){
			if (s.ct[i][j] == -1) printf("%d ", s.ct[i][j]);
			else printf("%d  ", s.ct[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	printf("\n");
	printf("sim:\n");
	for (int i = 0; i < (s.f); i++){
		for (int j = 0; j < 3; j++){
			if (s.sim[i][j] == -1) printf("%d ", s.sim[i][j]);
			else printf("%d  ", s.sim[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	BFS_free(&s);
}
