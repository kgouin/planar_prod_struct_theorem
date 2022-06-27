#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"tripod.h"

void tripod_init(struct bfs_struct* b, struct rmq_struct* r){
	BFS_init(b);
	BFS(b);
	LCA_init(r, b->ct, b->f); //preprocess the cotree for lca queries
	printf("initialization complete\n");

	tripod(r, b->tri[0][0], b->tri[0][1], b->tri[0][2]); //start the decomposition with the three triangles adjacent to the outer face
	//repeat triangles if necessary
	//in tripod() we'll check whether some (if not all) of these triangles are the same
}

int tripod(struct rmq_struct* r, int t1, int t2, int t3){ //recursive function

	if (t1 == t2 && t1 == t3 && t2 == t3){ //the three triangles are all the same
		//base case
	}
	else if (t1 != t2 && t1 != t3 && t2 != t3){ //the three triangles are all different
		//perform three lca queries
		//if one of the three triangles is an ancestor of the other two but there is a lower common ancestor,
			//then the lower common ancestor is our sperner triangle
		//if one of the three triangles is the lowest common ancestor (aka. two or three of the queries gives us the same answer),
			//then that answer is our sperner triangle
		//label & store the tripod (aka. walk along the three legs and label those vertices as belonging to the tripod)
		//recurse
	}
	else { //two of the three triangles are the same, the other is different
		//do not perform lca queries
		//simply take one of the three input triangles as your sperner triangle
		//label & store the tripod (aka. walk along the three legs and label those vertices as belonging to the tripod)
		//recurse
	}
	
	//we'll need to handle the case where an input triangle is right on a boundary
	
	//check: does the sperner triangle lead back to the three input faces when you walk up along the three legs of the tripod?

	return 0;
}







