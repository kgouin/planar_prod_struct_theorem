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

	int* acc = malloc((((r->n)/2)+1)*sizeof(int)); //memory leak

	for (int k = 0; k < (((r->n)/2)+1); k++){
		acc[k] = -1; //vertices labeled with -1 do not yet belong to a tripod
	}
	//label the three edges incident to face 0 as belonging to three exterior tripods
	acc[b->sim[0][0]] = -2;
	acc[b->sim[0][1]] = -3;
	acc[b->sim[0][2]] = -4;

	//start decomposition with the three triangles adjacent to outer face
	tripod(b, r, b->tri[0][0], b->tri[0][1], b->tri[0][2], acc);
}

int* tripod(struct bfs_struct* b, struct rmq_struct* r, int f1, int f2, int f3, int* acc){ //recursive function

	int sp; //this will hold our sperner triangle
	int u; //used to properly store tripods

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
		u = acc[b->sim[sp][0]];
		while (u == -1){
			acc[b->sim[sp][0]] = sp;
			u = acc[b->bt[b->sim[sp][0]]];
		}
		u = acc[b->sim[sp][1]];
		while (u == -1){
			acc[b->sim[sp][1]] = sp;
			u = acc[b->bt[b->sim[sp][1]]];
		}
		u = acc[b->sim[sp][2]];
		while (u == -1){
			acc[b->sim[sp][2]] = sp;
			u = acc[b->bt[b->sim[sp][2]]];
		}

		//recurse (we'll need to handle the case where an input face is right on a boundary)
	}

	else { //two of the three faces are the same, the other is different
		//do not perform lca queries. simply take one of the three input faces as your sperner triangle
		sp = f1;

		//store tripod (aka. for each vertex incident to sp, walk up bfs tree until you hit an edge belonging to another tripod)
		//b->sim[sp][0], b->sim[sp][1], b->sim[sp][2] are the vertices incident to our sperner triangle
		u = acc[b->sim[sp][0]];
		while (u == -1){
			acc[b->sim[sp][0]] = sp;
			u = acc[b->bt[b->sim[sp][0]]];
		}
		u = acc[b->sim[sp][1]];
		while (u == -1){
			acc[b->sim[sp][1]] = sp;
			u = acc[b->bt[b->sim[sp][1]]];
		}
		u = acc[b->sim[sp][2]];
		while (u == -1){
			acc[b->sim[sp][2]] = sp;
			u = acc[b->bt[b->sim[sp][2]]];
		}

		//recurse (we'll need to handle the case where an input face is right on a boundary)
	}

	printf("%d is one of our sperner triangles\n", sp); //sperner triangle is correctly identified
	
	printf("acc = ");
	for (int k = 0; k < ((r->n)/2)+1; k++){
		printf("%d ", acc[k]);
	}
	printf("\n");

	return 0;
}







