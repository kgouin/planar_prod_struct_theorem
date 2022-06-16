#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"bfs.h"

int main(){
	struct bfs_struct s;
	BFS_init(&s);
	BFS(&s);
	for (int i = 0; i < (s.v); i++){
		printf("%d ", s.bfs[i]);
	}
	printf("\n");
	BFS_free(&s);
}
