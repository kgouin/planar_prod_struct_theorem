#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"tripod.h"

void tripod_init(struct bfs_struct* b, struct rmq_struct* r){
	//initialize bfs_struct & rmq_struct
	BFS_init(b);
	BFS(b);
	LCA_init(r, b->ct, (2*(b->f))-1);
	printf("initialization complete\n");

	int* acc = malloc((b->v)*sizeof(int)); //memory leak

	for (int k = 0; k < (b->v); k++){
		acc[k] = -1; //vertices labeled with -1 do not yet belong to a tripod
	}
	//label the three edges incident to face 0 as belonging to three exterior tripods
	acc[b->sim[0][0]] = (b->f);
	acc[b->sim[0][1]] = (b->f)+1;
	acc[b->sim[0][2]] = (b->f)+2;

	//start decomposition with the three triangles adjacent to outer face
	tripod(b, r, b->tri[0][0], b->tri[0][1], b->tri[0][2], acc);
	acc[5] = 2;
	tripod(b, r, 3, 8, 5, acc);
}

int* tripod(struct bfs_struct* b, struct rmq_struct* r, int f1, int f2, int f3, int* acc){ //recursive function

	int sp; //this will hold our sperner triangle
	int u;
	int y;
	int v;

	if (f1 == f2 && f1 == f3 && f2 == f3){
		return acc;
	}

	else if (f1 != f2 && f1 != f3 && f2 != f3){
		if (LCA_query(r, f1, f2) == LCA_query(r, f1, f3) && LCA_query(r, f1, f3) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f2);
		else {
			if (LCA_query(r, f1, f2) == LCA_query(r, f1, f3)) sp = LCA_query(r, f2, f3);
			else if (LCA_query(r, f1, f2) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f3);
			else if (LCA_query(r, f1, f3) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f2);
		}

		//store tripod (aka. for each vertex incident to sp, walk up bfs tree until you hit an edge belonging to another tripod)
		//b->sim[sp][0], b->sim[sp][1], b->sim[sp][2] are the vertices incident to our sperner triangle
		u = b->sim[sp][0];
		while (acc[u] == -1){
			acc[u] = sp;
			u = b->bt[b->sim[sp][0]];
		}
		u = b->sim[sp][1];
		while (acc[u] == -1){
			acc[u] = sp;
			u = b->bt[b->sim[sp][1]];
		}
		u = b->sim[sp][2];
		while (acc[u] == -1){
			acc[u] = sp;
			u = b->bt[b->sim[sp][2]];
		}

		//recurse
		
		//printing the appropriate faces of vertices which belong to sp by walking up the bfs tree
		printf("\n");

		//take the last vertex we see before we hit either the root or another tripod when walking up the bfs
		//using that vertex, find left and right faces
		v = b->sim[sp][0];
		while (b->bt[b->bt[v]] != -2 && (acc[b->bt[v]] == -1 || acc[b->bt[v]] == sp)){
			v = b->bt[v];
		}
		printf("sp vertex a = %d\n", v);
		printf("left face = %d\n", b->il[v][b->pin[v]]);
		(b->pin[v] == 0) ? (y = (b->n[v])-1) : (y = b->pin[v]-1);
		printf("right face = %d\n", b->il[v][y]);
		printf("opposite face = %d\n", b->tri[sp][0]);

		v = b->sim[sp][1];
		while (b->bt[b->bt[v]] != -2 && (acc[b->bt[v]] == -1 || acc[b->bt[v]] == sp)){
			v = b->bt[v];
		}
		printf("sp vertex b = %d\n", v);
		printf("left face = %d\n", b->il[v][b->pin[v]]);
		(b->pin[v] == 0) ? (y = (b->n[v])-1) : (y = b->pin[v]-1);
		printf("right face = %d\n", b->il[v][y]);
		printf("opposite face = %d\n", b->tri[sp][1]);

		v = b->sim[sp][2];
		while (b->bt[b->bt[v]] != -2 && (acc[b->bt[v]] == -1 || acc[b->bt[v]] == sp)){
			v = b->bt[v];
		}
		printf("sp vertex c = %d\n", v);
		printf("left face = %d\n", b->il[v][b->pin[v]]);
		(b->pin[v] == 0) ? (y = (b->n[v])-1) : (y = b->pin[v]-1);
		printf("right face = %d\n", b->il[v][y]);
		printf("opposite face = %d\n", b->tri[sp][2]);

		printf("\n");

		//we need to associate the faces encountered to the appropriate subproblems
	}

	else { //two of the three faces are the same, the other is different
		//do not perform lca queries. simply take one of the three input faces as your sperner triangle
		sp = f1;

		//store tripod (aka. for each vertex incident to sp, walk up bfs tree until you hit an edge belonging to another tripod)
		//b->sim[sp][0], b->sim[sp][1], b->sim[sp][2] are the vertices incident to our sperner triangle
		u = b->sim[sp][0];
		while (acc[u] == -1){
			acc[u] = sp;
			u = b->bt[b->sim[sp][0]];
		}
		u = b->sim[sp][1];
		while (acc[u] == -1){
			acc[u] = sp;
			u = b->bt[b->sim[sp][1]];
		}
		u = b->sim[sp][2];
		while (acc[u] == -1){
			acc[u] = sp;
			u = b->bt[b->sim[sp][2]];
		}

		//recurse
	}

	printf("%d is one of our sperner triangles\n", sp); //sperner triangles are correctly identified
	
	printf("acc = ");
	for (int k = 0; k < (b->v); k++){
		printf("%d ", acc[k]);
	}
	printf("\n");
	printf("\n--------------------------------\n");

	return 0;
}
