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
	int* q = malloc((s->v)*sizeof(int));
	int* seen = calloc((s->v),sizeof(int));

	s->bfs[0] = 0;
	q[0] = 0;
	seen[0] = 1;

	int start = 0;
	int end = 1;
	//we use start and end markers to implement a queue, or a circular array
	int temp;
	while (start < end){
		//'remove' elt at start of q (aka. increase start marker)
		temp = q[start];
		start++;
		//add all neighbours of 'removed' elt to q and increase end marker
		for (int k = 0; k < s->n[temp]; k++){
			if (!seen[s->a[temp][k]]){
				s->bfs[end] = s->a[temp][k];
				q[end] = s->a[temp][k];
				seen[s->a[temp][k]] = 1;
				end++;
			}
		}
	}

	free(q);
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
