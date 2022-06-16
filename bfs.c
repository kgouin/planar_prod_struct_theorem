#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"bfs.h"

void BFS_init(struct bfs_struct* s){
	FILE *fd;
	fd = fopen("adjacencies.txt", "r");
	if (!fd) exit(1); //could not open text file

	fscanf(fd, "%d", &(s->v));
	//printf("%d\n", (s->f));
	s->a = malloc((s->v)*sizeof(int*));
	s->n = malloc((s->v)*sizeof(int*));
	for (int i = 0; i < (s->v); i++){
		fscanf(fd, "%d", &s->n[i]);
		s->a[i] = malloc(s->n[i]*sizeof(int));
		if (i == 0){
			s->r = malloc(s->n[i]*sizeof(int));
			for (int j = 0; j < s->n[i]; j++){
				fscanf(fd, "%d", &s->a[i][j]);
				s->r[j] = s->a[i][j];
				//printf("%d ", s->r[j]);
				//printf("%d ", s->a[i][j]);
			}
			//printf("\n");
		}
		else {
			for (int j = 0; j < s->n[i]; j++){
				fscanf(fd, "%d", &s->a[i][j]);
				//printf("%d ", s->a[i][j]);
			}
			//printf("\n");
		}
	}

	fclose(fd);
}

int* BFS(struct bfs_struct* s){
	s->bfs = malloc((s->v)*sizeof(int));
	//for (int i = 0; i < s->v; i++){ //for testing purposes
	//	bfs_array[i] = -1;
	//}
	int* seen = calloc((s->v),sizeof(int));

	int c = 0;
	s->bfs[0] = 0;
	seen[0] = 1;
	for (int i = 0; i < s->v; i++){
		for (int j = 0; j < s->n[i]; j++){
			if (!seen[s->a[i][j]]){
				s->bfs[++c] = s->a[i][j];
				seen[s->a[i][j]] = 1;
			}
		}
	}

	free(seen);

	return s->bfs;
}

void BFS_free(struct bfs_struct* s){
	for (int i = 0; i < (s->v); i++){
		free(s->a[i]);
	}
	free(s->a);
	free(s->n);
	free(s->r);
	free(s->bfs);
}
