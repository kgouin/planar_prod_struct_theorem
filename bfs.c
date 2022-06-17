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
			for (int j = 0; j < s->n[i]; j++){
				fscanf(fd, "%d", &s->a[i][j]);
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
	int* q = malloc((s->v)*sizeof(int));
	s->p = malloc((s->v)*sizeof(int));
	s->p[0] = -2;
	for (int i = 1; i < (s->v); i++){
		s->p[i] = -1;
	}

	q[0] = 0;

	int index;
	int start = 0;
	int end = 1;
	//we use start and end markers to implement a queue
	int temp;
	while (start < end){
		//'remove' elt at start of q (aka. increase start marker)
		temp = q[start];
		start++;
		//add all neighbours of 'removed' elt to q and increase end marker
		for (int k = 0; k < s->n[temp]; k++){
			index = s->a[temp][k];
			if (s->p[index] == -1){
				s->p[index] = temp;
				q[end] = index;
				end++;
			}
		}
	}

	free(q);

	return s->p;
}

void BFS_free(struct bfs_struct* s){
	for (int i = 0; i < (s->v); i++){
		free(s->a[i]);
	}
	free(s->a);
	free(s->n);
	free(s->p);
}
