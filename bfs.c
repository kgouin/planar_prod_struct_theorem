#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"bfs.h"

void BFS_init(struct bfs_struct* s){
	FILE *fd;
	fd = fopen("adjacencies.txt", "r");
	if (!fd) exit(1);

	fscanf(fd, "%d", &(s->v));
	s->al = malloc((s->v)*sizeof(int*));
	s->il = malloc((s->v)*sizeof(int*));
	s->n = malloc((s->v)*sizeof(int*));

	for (int i = 0; i < (s->v); i++){
		fscanf(fd, "%d", &s->n[i]);
		s->al[i] = malloc(s->n[i]*sizeof(int));
		s->il[i] = malloc(s->n[i]*sizeof(int));
		for (int j = 0; j < 2*(s->n[i]); j++){
			if (j < (s->n[i])) fscanf(fd, "%d", &s->al[i][j]);
			else fscanf(fd, "%d", &s->il[i][j-(s->n[i])]);
		}
	}

	fclose(fd);
}

int* BFS(struct bfs_struct* s){
	//initialization
	int* q = malloc((s->v)*sizeof(int));
	s->bt = malloc((s->v)*sizeof(int));
	s->ct = malloc(((2*(s->v))-4)*sizeof(int*));
	for (int k = 0; k < ((2*(s->v))-4); k++){
		s->ct[k] = malloc(3*sizeof(int));
		for (int j = 0; j < 3; j++){
			s->ct[k][j] = -1;
		}
	}
	for (int k = 0; k < s->v; k++){
		s->bt[k] = -1;
	}
	s->bt[0] = -2;

	//bfs spanning tree and cotree construction
	q[0] = 0;
	int start = 0;
	int end = 1;
	int temp;
	int x;
	int unique;
	int y;
	
	while (start < end){
		//'remove' elt at start of q (aka. increase start marker)
		temp = q[start];
		start++;
		//add all neighbours of 'removed' elt to q and increase end marker
		for (int k = 0; k < s->n[temp]; k++){
			if (s->bt[s->al[temp][k]] == -1){ //construct bfs tree
				s->bt[s->al[temp][k]] = temp;
				q[end] = s->al[temp][k];
				end++;
			}
			else if (s->bt[temp] != s->al[temp][k]){ //contruct cotree
				(k-1 < 0) ? (y = k-1+(s->n[k])) : (y = k-1);
				x = 0;
				unique = 1;
				while (s->ct[s->il[temp][k]][x] != -1 && x < 2){
					if (s->ct[s->il[temp][k]][x] == s->il[temp][y]) unique = 0;
					x++;
				}
				if (unique) s->ct[s->il[temp][k]][x] = s->il[temp][y];
				x = 0;
				unique = 1;
				while (s->ct[s->il[temp][y]][x] != -1 && x < 2){
					if (s->ct[s->il[temp][y]][x] == s->il[temp][k]) unique = 0;
					x++;
				}
				if (unique) s->ct[s->il[temp][y]][x] = s->il[temp][k];
			}
		}
	}

	//build cotree adjacency list with proper orientation and face 0 as the root
	printf("\n");
	printf("ct\n");
	for (int i = 0; i < ((2*(s->v))-4); i++){
		for (int j = 0; j < 3; j++){
			if (s->ct[i][j] == -1) printf("%d ", s->ct[i][j]);
			else printf("%d  ", s->ct[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	s->ct[0][1] = s->ct[0][0];
	s->ct[0][0] = -1;
	s->ct[0][2] = -1;
	for (int k = 0; k < ((2*(s->v))-4); k++){
		for (int m = 1; m < 3; m++){
			printf("k = %d, m = %d\n", k, m);
			if (s->ct[k][m] != -1 && s->ct[s->ct[k][m]][1] == k){
				s->ct[s->ct[k][m]][1] = s->ct[s->ct[k][m]][0];
				s->ct[s->ct[k][m]][0] = k;
			}
			else if (s->ct[k][m] != -1 && s->ct[s->ct[k][m]][2] == k){
				s->ct[s->ct[k][m]][2] = s->ct[s->ct[k][m]][0];
				s->ct[s->ct[k][m]][0] = k;
			}
		}
	}
	//then we have to decide if a specific child is a left child or a right child

	printf("\n");
	printf("ct'\n");
	for (int i = 0; i < ((2*(s->v))-4); i++){
		for (int j = 0; j < 3; j++){
			if (s->ct[i][j] == -1) printf("%d ", s->ct[i][j]);
			else printf("%d  ", s->ct[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	free(q);

	return s->bt;
}

void BFS_free(struct bfs_struct* s){
	free(s->n);
	for (int k = 0; k < (s->v); k++){
		free(s->al[k]);
		free(s->il[k]);
	}
	for (int k = 0; k < ((2*(s->v))-4); k++){
		free(s->ct[k]);
	}
	free(s->al);
	free(s->il);
	free(s->ct);
	free(s->bt);
}
