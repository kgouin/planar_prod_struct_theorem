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

	for (int k = 0; k < (s->v); k++){
		fscanf(fd, "%d", &s->n[k]);
		s->al[k] = malloc(s->n[k]*sizeof(int));
		s->il[k] = malloc(s->n[k]*sizeof(int));
		for (int j = 0; j < 2*(s->n[k]); j++){
			if (j < (s->n[k])) fscanf(fd, "%d", &s->al[k][j]);
			else fscanf(fd, "%d", &s->il[k][j-(s->n[k])]);
		}
	}

	fclose(fd);

	fd = fopen("triangles.txt", "r");
	if (!fd) exit(1);

	fscanf(fd, "%d", &(s->f));
	s->tri = malloc((s->f)*sizeof(int*));
	for (int k = 0; k < (s->f); k++){
		s->tri[k] = malloc(3*sizeof(int));
	}

	for (int k = 0; k < (s->f); k++){
		for (int j = 0; j < 3; j++){
			fscanf(fd, "%d", &s->tri[k][j]);
		}
	}

	fclose(fd);
	
	fd = fopen("simplices.txt", "r");
	if (!fd) exit(1);
	
	fscanf(fd, "%d", &(s->f));
	s->sim = malloc((s->f)*sizeof(int*));
	for (int k = 0; k < (s->f); k++){
		s->sim[k] = malloc(3*sizeof(int));
	}
	
	for (int k = 0; k < (s->f); k++){
		for (int j = 0; j < 3; j++){
			fscanf(fd, "%d", &s->sim[k][j]);
		}
	}

	fclose(fd);
}

int* BFS(struct bfs_struct* s){
	//initialization
	int* q = malloc((s->f)*sizeof(int));
	s->bt = malloc((s->v)*sizeof(int));
	s->pin = malloc((s->v)*sizeof(int));
	s->pin[0] = -2;
	s->ct = malloc((s->f)*sizeof(int*));
	for (int k = 0; k < (s->f); k++){
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
				s->bt[s->al[temp][k]] = temp; //here we set the parent
				q[end] = s->al[temp][k];
				end++;
				for (int j = 0; j < s->n[s->al[temp][k]]; j++){
					if (s->al[s->al[temp][k]][j] == temp) s->pin[s->al[temp][k]] = j; //here we set the parent index within al
				}
			}
			else if (s->bt[temp] != s->al[temp][k] && s->bt[k] != -2){ //contruct cotree
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
	s->ct[0][1] = s->ct[0][0];
	s->ct[0][0] = -1;
	s->ct[0][2] = -1;

	q[0] = 0;
	start = 0;
	end = 1;

	while (start < end){
		temp = q[start];
		start++;
		if (s->ct[temp][1] != -1){ //look at left child
			//make sure parent is correct
			if (s->ct[s->ct[temp][1]][1] == temp){
				s->ct[s->ct[temp][1]][1] = s->ct[s->ct[temp][1]][0];
				s->ct[s->ct[temp][1]][0] = temp;
			}
			else if (s->ct[s->ct[temp][1]][2] == temp){
				s->ct[s->ct[temp][1]][2] = s->ct[s->ct[temp][1]][0];
				s->ct[s->ct[temp][1]][0] = temp;
			}
			//else parent is in the correct position

			q[end] = s->ct[temp][1];
			end++;
		}
		if (s->ct[temp][2] != -1){ //look at right child
			//make sure parent is correct
			if (s->ct[s->ct[temp][2]][1] == temp){
				s->ct[s->ct[temp][2]][1] = s->ct[s->ct[temp][2]][0];
				s->ct[s->ct[temp][2]][0] = temp;
			}
			else if (s->ct[s->ct[temp][2]][2] == temp){
				s->ct[s->ct[temp][2]][2] = s->ct[s->ct[temp][2]][0];
				s->ct[s->ct[temp][2]][0] = temp;
			}
			//else parent is in the correct position

			q[end] = s->ct[temp][2];
			end++;
		}
	}

	//the position the children might not matter
	int lc; //will hold the left child
	for (int k = 1; k < (s->f); k++){
		if (s->ct[k][1] != -1){ //if a node has children
			for (int j = 0; j < 3; j++){
				if (s->ct[k][0] == s->tri[k][j]){ //searching for parent in tri
					(j < 2) ? (lc = s->tri[k][j+1]) : (lc = s->tri[k][0]);
					if (s->ct[k][1] != lc && s->ct[k][2] != lc) lc = -1; //special case when child not part of cotree
					break;
				}
			}
			if (s->ct[k][2] == lc){
				s->ct[k][2] = s->ct[k][1];
				s->ct[k][1] = lc;
			}
		}
	}

	free(q);

	return s->bt;
}

void BFS_free(struct bfs_struct* s){
	free(s->n);
	for (int k = 0; k < (s->v); k++){
		free(s->al[k]);
		free(s->il[k]);
	}
	for (int k = 0; k < (s->f); k++){
		free(s->ct[k]);
		free(s->tri[k]);
		free(s->sim[k]);
	}
	free(s->al);
	free(s->il);
	free(s->ct);
	free(s->tri);
	free(s->sim);
	free(s->bt);
	free(s->pin);
}
