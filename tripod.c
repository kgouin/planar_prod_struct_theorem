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
	printf("bfs tree = ");
	for (int k = 0; k < b->v; k++){
		printf("%d ", b->bt[k]);
	}
	printf("\n");

	int* acc = malloc((b->v)*sizeof(int)); //memory leak

	for (int k = 0; k < (b->v); k++){
		acc[k] = -1; //vertices labeled with -1 do not yet belong to a tripod
	}
	//label the three edges incident to face 0 as belonging to three exterior tripods
	acc[b->sim[0][0]] = (b->f);
	acc[b->sim[0][1]] = (b->f)+1;
	acc[b->sim[0][2]] = (b->f)+2;

	//start decomposition with the three triangles adjacent to outer face
	tripod(b, r, b->tri[0][0], b->tri[0][2], b->tri[0][1], acc);
}

int* tripod(struct bfs_struct* b, struct rmq_struct* r, int f1, int f2, int f3, int* acc){ //recursive function
	//definitions
	int sp; //sperner triangle
	int u; //used to label vertices as belonging to tripods
	int v_a; //the last vertex before we hit another tripod
	int v_a_next; //the vertex at which we hit another tripod
	int v_a_l; //the face to the left of the edge (v_a, v_a_next) when going up the bfs tree
	int v_a_r; //the face to the right of the edge (v_a, v_a_next) when going up the bfs tree
	int v_a_op = -1; //the triangle with bichromatic edge on the cycle defining the subproblem
	int v_a_mirror; //the triangle adjacent to newly found sp
	int y_a;
	int v_b;
	int v_b_next;
	int v_b_l;
	int v_b_r;
	int v_b_op = -1;
	int v_b_mirror;
	int y_b;
	int v_c;
	int v_c_next;
	int v_c_l;
	int v_c_r;
	int v_c_op = -1;
	int v_c_mirror;
	int y_c;

	if (f1 == f2 && f1 == f3 && f2 == f3) return 0;

	else if (f1 != f2 && f1 != f3 && f2 != f3){
		//find sperner triangle
		if (LCA_query(r, f1, f2) == LCA_query(r, f1, f3) && LCA_query(r, f1, f3) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f2);
		else {
			if (LCA_query(r, f1, f2) == LCA_query(r, f1, f3)) sp = LCA_query(r, f2, f3);
			else if (LCA_query(r, f1, f2) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f3);
			else if (LCA_query(r, f1, f3) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f2);
		}

		//store tripod
		u = b->sim[sp][0];
		v_a = u;
		v_a_next = -1;
		while (acc[u] == -1){
			acc[u] = sp;
			v_a = u;
			u = b->bt[u];
			v_a_next = u;
		}
		v_a_l = b->il[v_a][b->pin[v_a]];
		(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
		v_a_r = b->il[v_a][y_a];

		u = b->sim[sp][1];
		v_b = u;
		v_b_next = -1;
		while (acc[u] == -1){
			acc[u] = sp;
			v_b = u;
			u = b->bt[u];
			v_b_next = u;
		}
		v_b_l = b->il[v_b][b->pin[v_b]];
		(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
		v_b_r = b->il[v_b][y_b];

		u = b->sim[sp][2];
		v_c = u;
		v_c_next = -1;
		while (acc[u] == -1){
			acc[u] = sp;
			v_c = u;
			u = b->bt[u];
			v_c_next = u;
		}
		v_c_l = b->il[v_c][b->pin[v_c]];
		(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
		v_c_r = b->il[v_c][y_c];

		printf("\n");
		printf("%d is one of our sperner triangles\n", sp); //sperner triangles are correctly identified
		printf("acc = ");
		for (int k = 0; k < (b->v); k++){
			printf("%d ", acc[k]);
		}
		printf("\n");

		//recurse

		v_a_op = f1;
		v_b_op = f2;
		v_c_op = f3;

		v_a_mirror = b->tri[sp][0];
		v_b_mirror = b->tri[sp][1];
		v_c_mirror = b->tri[sp][2];

		int partial_match;
		int double_match;

		if (v_a_next == -1 && v_b_next == -1){ //check to see where we hit the next tripod
			//acc[v_a] gives us the colour of the tripod that vertex a touches
			//acc[v_b] gives us the colour of the tripod that vertex b touches
			
			partial_match = 0;
			double_match = 0;
			for (int i = 0; i < 3; i++){ //look at all three vertices of our v_a_op triangle
				if ((acc[b->sim[v_a_op][i]] == acc[v_a] || acc[b->sim[v_a_op][i]] == acc[v_b]) && !partial_match){
					partial_match = 1;
				}
				else if ((acc[b->sim[v_a_op][i]] == acc[v_a] || acc[b->sim[v_a_op][i]] == acc[v_b]) && partial_match){
					double_match = 1;
				}
			}
			if (double_match){
				printf("NO SWITCH\n");
				//v_a_op is part of the subproblem between vertex v_a and vertex v_b
				//v_b_op is part of the subproblem between vertex v_b and vertex v_c
				//v_c_op is part of the subproblem between vertex v_c and vertex v_a

				//this is already what we have. no need to change things
			}
			else { //no double match. reset values and try with v_b_op
				partial_match = 0;
				double_match = 0;
				for (int i = 0; i < 3; i++){ //look at all three vertices of our v_b_op triangle
					if ((acc[b->sim[v_b_op][i]] == acc[v_a] || acc[b->sim[v_b_op][i]] == acc[v_b]) && !partial_match){
						partial_match = 1;
					}
					else if ((acc[b->sim[v_b_op][i]] == acc[v_a] || acc[b->sim[v_b_op][i]] == acc[v_b]) && partial_match){
						double_match = 1;
					}
				}
				if (double_match){
					printf("SWITCH 1\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_a;
					//int v_temp_next = v_a_next;
					//int v_temp_l = v_a_l;
					//int v_temp_r = v_a_r;
					int v_temp_op = v_a_op;
					//int v_temp_mirror = v_a_mirror;

					//v_a = v_b;
					//v_a_next = v_b_next;
					//v_a_l = v_b_l;
					//v_a_r = v_b_r;
					v_a_op = v_b_op;
					//v_a_mirror = v_b_mirror;

					//v_b = v_c;
					//v_b_next = v_c_next;
					//v_b_l = v_c_l;
					//v_b_r = v_c_r;
					v_b_op = v_c_op;
					//v_b_mirror = v_c_mirror;

					//v_c = v_temp;
					//v_c_next = v_temp_next;
					//v_c_l = v_temp_l;
					//v_c_r = v_temp_r;
					v_c_op = v_temp_op;
					//v_c_mirror = v_temp_mirror;
				}
				else { //no double match. v_c_op must match.
					printf("SWITCH 2\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_c;
					//int v_temp_next = v_c_next;
					//int v_temp_l = v_c_l;
					//int v_temp_r = v_c_r;
					int v_temp_op = v_c_op;
					//int v_temp_mirror = v_c_mirror;

					//v_c = v_b;
					//v_c_next = v_b_next;
					//v_c_l = v_b_l;
					//v_c_r = v_b_r;
					v_c_op = v_b_op;
					//v_c_mirror = v_b_mirror;

					//v_b = v_a;
					//v_b_next = v_a_next;
					//v_b_l = v_a_l;
					//v_b_r = v_a_r;
					v_b_op = v_a_op;
					//v_b_mirror = v_a_mirror;

					//v_a = v_temp;
					//v_a_next = v_temp_next;
					//v_a_l = v_temp_l;
					//v_a_r = v_temp_r;
					v_a_op = v_temp_op;
					//v_a_mirror = v_temp_mirror;
				}
			}
		}
		else if (v_a_next == -1 && v_b_next != -1){ //check to see where we hit the next tripod
			//acc[v_a] gives us the colour of the tripod that vertex a touches
			//acc[v_b_next] gives us the colour of the tripod that vertex b touches
			
			partial_match = 0;
			double_match = 0;
			for (int i = 0; i < 3; i++){ //look at all three vertices of our v_a_op triangle
				if ((acc[b->sim[v_a_op][i]] == acc[v_a] || acc[b->sim[v_a_op][i]] == acc[v_b_next]) && !partial_match){
					partial_match = 1;
				}
				else if ((acc[b->sim[v_a_op][i]] == acc[v_a] || acc[b->sim[v_a_op][i]] == acc[v_b_next]) && partial_match){
					double_match = 1;
				}
			}
			if (double_match){
				printf("NO SWITCH\n");
				//v_a_op is part of the subproblem between vertex v_a and vertex v_b
				//v_b_op is part of the subproblem between vertex v_b and vertex v_c
				//v_c_op is part of the subproblem between vertex v_c and vertex v_a

				//this is already what we have. no need to change things
			}
			else { //no double match. reset values and try with v_b_op
				partial_match = 0;
				double_match = 0;
				for (int i = 0; i < 3; i++){ //look at all three vertices of our v_b_op triangle
					if ((acc[b->sim[v_b_op][i]] == acc[v_a] || acc[b->sim[v_b_op][i]] == acc[v_b_next]) && !partial_match){
						partial_match = 1;
					}
					else if ((acc[b->sim[v_b_op][i]] == acc[v_a] || acc[b->sim[v_b_op][i]] == acc[v_b_next]) && partial_match){
						double_match = 1;
					}
				}
				if (double_match){
					printf("SWITCH 1\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_a;
					//int v_temp_next = v_a_next;
					//int v_temp_l = v_a_l;
					//int v_temp_r = v_a_r;
					int v_temp_op = v_a_op;
					//int v_temp_mirror = v_a_mirror;

					//v_a = v_b;
					//v_a_next = v_b_next;
					//v_a_l = v_b_l;
					//v_a_r = v_b_r;
					v_a_op = v_b_op;
					//v_a_mirror = v_b_mirror;

					//v_b = v_c;
					//v_b_next = v_c_next;
					//v_b_l = v_c_l;
					//v_b_r = v_c_r;
					v_b_op = v_c_op;
					//v_b_mirror = v_c_mirror;

					//v_c = v_temp;
					//v_c_next = v_temp_next;
					//v_c_l = v_temp_l;
					//v_c_r = v_temp_r;
					v_c_op = v_temp_op;
					//v_c_mirror = v_temp_mirror;
				}
				else { //no double match. v_c_op must match.
					printf("SWITCH 2\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_c;
					//int v_temp_next = v_c_next;
					//int v_temp_l = v_c_l;
					//int v_temp_r = v_c_r;
					int v_temp_op = v_c_op;
					//int v_temp_mirror = v_c_mirror;

					//v_c = v_b;
					//v_c_next = v_b_next;
					//v_c_l = v_b_l;
					//v_c_r = v_b_r;
					v_c_op = v_b_op;
					//v_c_mirror = v_b_mirror;

					//v_b = v_a;
					//v_b_next = v_a_next;
					//v_b_l = v_a_l;
					//v_b_r = v_a_r;
					v_b_op = v_a_op;
					//v_b_mirror = v_a_mirror;

					//v_a = v_temp;
					//v_a_next = v_temp_next;
					//v_a_l = v_temp_l;
					//v_a_r = v_temp_r;
					v_a_op = v_temp_op;
					//v_a_mirror = v_temp_mirror;
				}
			}
		}
		else if (v_a_next != -1 && v_b_next == -1){ //check to see where we hit the next tripod
			//acc[v_a_next] gives us the colour of the tripod that vertex a touches
			//acc[v_b] gives us the colour of the tripod that vertex b touches
			
			partial_match = 0;
			double_match = 0;
			for (int i = 0; i < 3; i++){ //look at all three vertices of our v_a_op triangle
				if ((acc[b->sim[v_a_op][i]] == acc[v_a_next] || acc[b->sim[v_a_op][i]] == acc[v_b]) && !partial_match){
					partial_match = 1;
				}
				else if ((acc[b->sim[v_a_op][i]] == acc[v_a_next] || acc[b->sim[v_a_op][i]] == acc[v_b]) && partial_match){
					double_match = 1;
				}
			}
			if (double_match){
				printf("NO SWITCH\n");
				//v_a_op is part of the subproblem between vertex v_a and vertex v_b
				//v_b_op is part of the subproblem between vertex v_b and vertex v_c
				//v_c_op is part of the subproblem between vertex v_c and vertex v_a

				//this is already what we have. no need to change things
			}
			else { //no double match. reset values and try with v_b_op
				partial_match = 0;
				double_match = 0;
				for (int i = 0; i < 3; i++){ //look at all three vertices of our v_b_op triangle
					if ((acc[b->sim[v_b_op][i]] == acc[v_a_next] || acc[b->sim[v_b_op][i]] == acc[v_b]) && !partial_match){
						partial_match = 1;
					}
					else if ((acc[b->sim[v_b_op][i]] == acc[v_a_next] || acc[b->sim[v_b_op][i]] == acc[v_b]) && partial_match){
						double_match = 1;
					}
				}
				if (double_match){
					printf("SWITCH 1\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_a;
					//int v_temp_next = v_a_next;
					//int v_temp_l = v_a_l;
					//int v_temp_r = v_a_r;
					int v_temp_op = v_a_op;
					//int v_temp_mirror = v_a_mirror;

					//v_a = v_b;
					//v_a_next = v_b_next;
					//v_a_l = v_b_l;
					//v_a_r = v_b_r;
					v_a_op = v_b_op;
					//v_a_mirror = v_b_mirror;

					//v_b = v_c;
					//v_b_next = v_c_next;
					//v_b_l = v_c_l;
					//v_b_r = v_c_r;
					v_b_op = v_c_op;
					//v_b_mirror = v_c_mirror;

					//v_c = v_temp;
					//v_c_next = v_temp_next;
					//v_c_l = v_temp_l;
					//v_c_r = v_temp_r;
					v_c_op = v_temp_op;
					//v_c_mirror = v_temp_mirror;
				}
				else { //no double match. v_c_op must match.
					printf("SWITCH 2\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_c;
					//int v_temp_next = v_c_next;
					//int v_temp_l = v_c_l;
					//int v_temp_r = v_c_r;
					int v_temp_op = v_c_op;
					//int v_temp_mirror = v_c_mirror;

					//v_c = v_b;
					//v_c_next = v_b_next;
					//v_c_l = v_b_l;
					//v_c_r = v_b_r;
					v_c_op = v_b_op;
					//v_c_mirror = v_b_mirror;

					//v_b = v_a;
					//v_b_next = v_a_next;
					//v_b_l = v_a_l;
					//v_b_r = v_a_r;
					v_b_op = v_a_op;
					//v_b_mirror = v_a_mirror;

					//v_a = v_temp;
					//v_a_next = v_temp_next;
					//v_a_l = v_temp_l;
					//v_a_r = v_temp_r;
					v_a_op = v_temp_op;
					//v_a_mirror = v_temp_mirror;
				}
			}
		}
		else if (v_a_next != -1 && v_b_next != -1){ //check to see where we hit the next tripod
			//acc[v_a_next] gives us the colour of the tripod that vertex a touches
			//acc[v_b_next] gives us the colour of the tripod that vertex b touches
			
			partial_match = 0;
			double_match = 0;
			for (int i = 0; i < 3; i++){ //look at all three vertices of our v_a_op triangle
				if ((acc[b->sim[v_a_op][i]] == acc[v_a_next] || acc[b->sim[v_a_op][i]] == acc[v_b_next]) && !partial_match){
					partial_match = 1;
				}
				else if ((acc[b->sim[v_a_op][i]] == acc[v_a_next] || acc[b->sim[v_a_op][i]] == acc[v_b_next]) && partial_match){
					double_match = 1;
				}
			}
			if (double_match){
				printf("NO SWITCH\n");
				//v_a_op is part of the subproblem between vertex v_a and vertex v_b
				//v_b_op is part of the subproblem between vertex v_b and vertex v_c
				//v_c_op is part of the subproblem between vertex v_c and vertex v_a

				//this is already what we have. no need to change things
			}
			else { //no double match. reset values and try with v_b_op
				partial_match = 0;
				double_match = 0;
				for (int i = 0; i < 3; i++){ //look at all three vertices of our v_b_op triangle
					if ((acc[b->sim[v_b_op][i]] == acc[v_a_next] || acc[b->sim[v_b_op][i]] == acc[v_b_next]) && !partial_match){
						partial_match = 1;
					}
					else if ((acc[b->sim[v_b_op][i]] == acc[v_a_next] || acc[b->sim[v_b_op][i]] == acc[v_b_next]) && partial_match){
						double_match = 1;
					}
				}
				if (double_match){
					printf("SWITCH 1\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_a;
					//int v_temp_next = v_a_next;
					//int v_temp_l = v_a_l;
					//int v_temp_r = v_a_r;
					int v_temp_op = v_a_op;
					//int v_temp_mirror = v_a_mirror;

					//v_a = v_b;
					//v_a_next = v_b_next;
					//v_a_l = v_b_l;
					//v_a_r = v_b_r;
					v_a_op = v_b_op;
					//v_a_mirror = v_b_mirror;

					//v_b = v_c;
					//v_b_next = v_c_next;
					//v_b_l = v_c_l;
					//v_b_r = v_c_r;
					v_b_op = v_c_op;
					//v_b_mirror = v_c_mirror;

					//v_c = v_temp;
					//v_c_next = v_temp_next;
					//v_c_l = v_temp_l;
					//v_c_r = v_temp_r;
					v_c_op = v_temp_op;
					//v_c_mirror = v_temp_mirror;
				}
				else { //no double match. v_c_op must match.
					printf("SWITCH 2\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_c;
					//int v_temp_next = v_c_next;
					//int v_temp_l = v_c_l;
					//int v_temp_r = v_c_r;
					int v_temp_op = v_c_op;
					//int v_temp_mirror = v_c_mirror;

					//v_c = v_b;
					//v_c_next = v_b_next;
					//v_c_l = v_b_l;
					//v_c_r = v_b_r;
					v_c_op = v_b_op;
					//v_c_mirror = v_b_mirror;

					//v_b = v_a;
					//v_b_next = v_a_next;
					//v_b_l = v_a_l;
					//v_b_r = v_a_r;
					v_b_op = v_a_op;
					//v_b_mirror = v_a_mirror;

					//v_a = v_temp;
					//v_a_next = v_temp_next;
					//v_a_l = v_temp_l;
					//v_a_r = v_temp_r;
					v_a_op = v_temp_op;
					//v_a_mirror = v_temp_mirror;
				}
			}
		}

		printf("v_a = %d\n", v_a);
		printf("v_b = %d\n", v_b);
		printf("v_c = %d\n", v_c);
		printf("v_a_mirror = %d\n", v_a_mirror);
		printf("v_b_mirror = %d\n", v_b_mirror);
		printf("v_c_mirror = %d\n", v_c_mirror);
		printf("v_a_op = %d\n", v_a_op);
		printf("v_b_op = %d\n", v_b_op);
		printf("v_c_op = %d\n", v_c_op);

		if (!(acc[v_a] != sp && acc[v_b] != sp)){ //if the subproblem exists
			if (acc[v_a] == sp && acc[v_b] == sp){ //leg a is non-empty && leg b is non-empty
				printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_b_r, v_a_l);
				tripod(b, r, v_a_op, v_b_r, v_a_l, acc);
			}
			else if (acc[v_a] == sp && acc[v_b] != sp){ //leg a is non-empty && leg b is empty
				printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_a_mirror, v_a_l);
				tripod(b, r, v_a_op, v_a_mirror, v_a_l, acc);
			}
			else if (acc[v_a] != sp && acc[v_b] == sp){ //leg a is empty && leg b is non-empty
				printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_b_r, v_a_mirror);
				tripod(b, r, v_a_op, v_b_r, v_a_mirror, acc);
			}
			/*if (acc[v_a] != sp && acc[v_b] != sp){ //leg a is empty && leg b is empty
				printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_a_mirror, v_a_mirror);
				tripod(b, r, v_a_op, v_a_mirror, v_a_mirror, acc);
			}*/
		}

		if (!(acc[v_b] != sp && acc[v_c] != sp)){ //if the subproblem exists
			if (acc[v_b] == sp && acc[v_c] == sp){ //leg b is non-empty && leg c is non-empty
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_c_r, v_b_l);
				tripod(b, r, v_b_op, v_c_r, v_b_l, acc);
			}
			else if (acc[v_b] == sp && acc[v_c] != sp){ //leg b is non-empty && leg c is empty
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_b_mirror, v_b_l);
				tripod(b, r, v_b_op, v_b_mirror, v_b_l, acc);
			}
			else if (acc[v_b] != sp && acc[v_c] == sp){ //leg b is empty && leg c is non-empty
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_c_r, v_b_mirror);
				tripod(b, r, v_b_op, v_c_r, v_b_mirror, acc);
			}
			/*if (acc[v_b] != sp && acc[v_c] != sp){ //leg b is empty && leg c is empty
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_b_mirror, v_b_mirror);
				tripod(b, r, v_b_op, v_b_mirror, v_b_mirror, acc);
			}*/
		}

		if (!(acc[v_c] != sp && acc[v_a] != sp)){ //if the subproblem exists
			if (acc[v_c] == sp && acc[v_a] == sp){ //leg c is non-empty && leg a is non-empty
				printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_a_r, v_c_l);
				tripod(b, r, v_c_op, v_a_r, v_c_l, acc);
			}
			else if (acc[v_c] == sp && acc[v_a] != sp){ //leg c is non-empty && leg a is empty
				printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_c_mirror, v_c_l);
				tripod(b, r, v_c_op, v_c_mirror, v_c_l, acc);
			}
			else if (acc[v_c] != sp && acc[v_a] == sp){ //leg c is empty && leg a is non-empty
				printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_a_r, v_c_mirror);
				tripod(b, r, v_c_op, v_a_r, v_c_mirror, acc);
			}
			/*if (acc[v_c] != sp && acc[v_a] != sp){ //leg c is empty && leg a is empty
				printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_c_mirror, v_c_mirror);
				tripod(b, r, v_c_op, v_c_mirror, v_c_mirror, acc);
			}*/
		}
	}

	else {
		//find sperner triangle
		sp = f1; //is this really what we want? no LCA query?

		//store tripod
		v_a_mirror = b->tri[sp][0];
		v_b_mirror = b->tri[sp][1];
		v_c_mirror = b->tri[sp][2];

		u = b->sim[sp][0];
		v_a = u;
		v_a_next = 0;
		while (acc[u] == -1){
			acc[u] = sp;
			v_a = u;
			u = b->bt[u];
			v_a_next = u;
		}
		v_a_l = b->il[v_a][b->pin[v_a]];
		(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
		v_a_r = b->il[v_a][y_a];

		u = b->sim[sp][1];
		v_b = u;
		v_b_next = 0;
		while (acc[u] == -1){
			acc[u] = sp;
			v_b = u;
			u = b->bt[u];
			v_b_next = u;
		}
		v_b_l = b->il[v_b][b->pin[v_b]];
		(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
		v_b_r = b->il[v_b][y_b];

		u = b->sim[sp][2];
		v_c = u;
		v_c_next = 0;
		while (acc[u] == -1){
			acc[u] = sp;
			v_c = u;
			u = b->bt[u];
			v_c_next = u;
		}
		v_c_l = b->il[v_c][b->pin[v_c]];
		(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
		v_c_r = b->il[v_c][y_c];

		printf("\n");
		printf("%d is one of our sperner triangles\n", sp); //sperner triangles are correctly identified
		printf("acc = ");
		for (int k = 0; k < (b->v); k++){
			printf("%d ", acc[k]);
		}
		printf("\n");

		//recurse

		v_a_op = f1;
		v_b_op = f2;
		v_c_op = f3;

		v_a_mirror = b->tri[sp][0];
		v_b_mirror = b->tri[sp][1];
		v_c_mirror = b->tri[sp][2];

		int partial_match;
		int double_match;

		if (v_a_next == -1 && v_b_next == -1){ //check to see where we hit the next tripod
			//acc[v_a] gives us the colour of the tripod that vertex a touches
			//acc[v_b] gives us the colour of the tripod that vertex b touches
			
			partial_match = 0;
			double_match = 0;
			for (int i = 0; i < 3; i++){ //look at all three vertices of our v_a_op triangle
				if ((acc[b->sim[v_a_op][i]] == acc[v_a] || acc[b->sim[v_a_op][i]] == acc[v_b]) && !partial_match){
					partial_match = 1;
				}
				else if ((acc[b->sim[v_a_op][i]] == acc[v_a] || acc[b->sim[v_a_op][i]] == acc[v_b]) && partial_match){
					double_match = 1;
				}
			}
			if (double_match){
				printf("NO SWITCH\n");
				//v_a_op is part of the subproblem between vertex v_a and vertex v_b
				//v_b_op is part of the subproblem between vertex v_b and vertex v_c
				//v_c_op is part of the subproblem between vertex v_c and vertex v_a

				//this is already what we have. no need to change things
			}
			else { //no double match. reset values and try with v_b_op
				partial_match = 0;
				double_match = 0;
				for (int i = 0; i < 3; i++){ //look at all three vertices of our v_b_op triangle
					if ((acc[b->sim[v_b_op][i]] == acc[v_a] || acc[b->sim[v_b_op][i]] == acc[v_b]) && !partial_match){
						partial_match = 1;
					}
					else if ((acc[b->sim[v_b_op][i]] == acc[v_a] || acc[b->sim[v_b_op][i]] == acc[v_b]) && partial_match){
						double_match = 1;
					}
				}
				if (double_match){
					printf("SWITCH 1\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_a;
					//int v_temp_next = v_a_next;
					//int v_temp_l = v_a_l;
					//int v_temp_r = v_a_r;
					int v_temp_op = v_a_op;
					//int v_temp_mirror = v_a_mirror;

					//v_a = v_b;
					//v_a_next = v_b_next;
					//v_a_l = v_b_l;
					//v_a_r = v_b_r;
					v_a_op = v_b_op;
					//v_a_mirror = v_b_mirror;

					//v_b = v_c;
					//v_b_next = v_c_next;
					//v_b_l = v_c_l;
					//v_b_r = v_c_r;
					v_b_op = v_c_op;
					//v_b_mirror = v_c_mirror;

					//v_c = v_temp;
					//v_c_next = v_temp_next;
					//v_c_l = v_temp_l;
					//v_c_r = v_temp_r;
					v_c_op = v_temp_op;
					//v_c_mirror = v_temp_mirror;
				}
				else { //no double match. v_c_op must match.
					printf("SWITCH 2\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_c;
					//int v_temp_next = v_c_next;
					//int v_temp_l = v_c_l;
					//int v_temp_r = v_c_r;
					int v_temp_op = v_c_op;
					//int v_temp_mirror = v_c_mirror;

					//v_c = v_b;
					//v_c_next = v_b_next;
					//v_c_l = v_b_l;
					//v_c_r = v_b_r;
					v_c_op = v_b_op;
					//v_c_mirror = v_b_mirror;

					//v_b = v_a;
					//v_b_next = v_a_next;
					//v_b_l = v_a_l;
					//v_b_r = v_a_r;
					v_b_op = v_a_op;
					//v_b_mirror = v_a_mirror;

					//v_a = v_temp;
					//v_a_next = v_temp_next;
					//v_a_l = v_temp_l;
					//v_a_r = v_temp_r;
					v_a_op = v_temp_op;
					//v_a_mirror = v_temp_mirror;
				}
			}
		}
		else if (v_a_next == -1 && v_b_next != -1){ //check to see where we hit the next tripod
			//acc[v_a] gives us the colour of the tripod that vertex a touches
			//acc[v_b_next] gives us the colour of the tripod that vertex b touches
			
			partial_match = 0;
			double_match = 0;
			for (int i = 0; i < 3; i++){ //look at all three vertices of our v_a_op triangle
				if ((acc[b->sim[v_a_op][i]] == acc[v_a] || acc[b->sim[v_a_op][i]] == acc[v_b_next]) && !partial_match){
					partial_match = 1;
				}
				else if ((acc[b->sim[v_a_op][i]] == acc[v_a] || acc[b->sim[v_a_op][i]] == acc[v_b_next]) && partial_match){
					double_match = 1;
				}
			}
			if (double_match){
				printf("NO SWITCH\n");
				//v_a_op is part of the subproblem between vertex v_a and vertex v_b
				//v_b_op is part of the subproblem between vertex v_b and vertex v_c
				//v_c_op is part of the subproblem between vertex v_c and vertex v_a

				//this is already what we have. no need to change things
			}
			else { //no double match. reset values and try with v_b_op
				partial_match = 0;
				double_match = 0;
				for (int i = 0; i < 3; i++){ //look at all three vertices of our v_b_op triangle
					if ((acc[b->sim[v_b_op][i]] == acc[v_a] || acc[b->sim[v_b_op][i]] == acc[v_b_next]) && !partial_match){
						partial_match = 1;
					}
					else if ((acc[b->sim[v_b_op][i]] == acc[v_a] || acc[b->sim[v_b_op][i]] == acc[v_b_next]) && partial_match){
						double_match = 1;
					}
				}
				if (double_match){
					printf("SWITCH 1\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_a;
					//int v_temp_next = v_a_next;
					//int v_temp_l = v_a_l;
					//int v_temp_r = v_a_r;
					int v_temp_op = v_a_op;
					//int v_temp_mirror = v_a_mirror;

					//v_a = v_b;
					//v_a_next = v_b_next;
					//v_a_l = v_b_l;
					//v_a_r = v_b_r;
					v_a_op = v_b_op;
					//v_a_mirror = v_b_mirror;

					//v_b = v_c;
					//v_b_next = v_c_next;
					//v_b_l = v_c_l;
					//v_b_r = v_c_r;
					v_b_op = v_c_op;
					//v_b_mirror = v_c_mirror;

					//v_c = v_temp;
					//v_c_next = v_temp_next;
					//v_c_l = v_temp_l;
					//v_c_r = v_temp_r;
					v_c_op = v_temp_op;
					//v_c_mirror = v_temp_mirror;
				}
				else { //no double match. v_c_op must match.
					printf("SWITCH 2\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_c;
					//int v_temp_next = v_c_next;
					//int v_temp_l = v_c_l;
					//int v_temp_r = v_c_r;
					int v_temp_op = v_c_op;
					//int v_temp_mirror = v_c_mirror;

					//v_c = v_b;
					//v_c_next = v_b_next;
					//v_c_l = v_b_l;
					//v_c_r = v_b_r;
					v_c_op = v_b_op;
					//v_c_mirror = v_b_mirror;

					//v_b = v_a;
					//v_b_next = v_a_next;
					//v_b_l = v_a_l;
					//v_b_r = v_a_r;
					v_b_op = v_a_op;
					//v_b_mirror = v_a_mirror;

					//v_a = v_temp;
					//v_a_next = v_temp_next;
					//v_a_l = v_temp_l;
					//v_a_r = v_temp_r;
					v_a_op = v_temp_op;
					//v_a_mirror = v_temp_mirror;
				}
			}
		}
		else if (v_a_next != -1 && v_b_next == -1){ //check to see where we hit the next tripod
			//acc[v_a_next] gives us the colour of the tripod that vertex a touches
			//acc[v_b] gives us the colour of the tripod that vertex b touches
			
			partial_match = 0;
			double_match = 0;
			for (int i = 0; i < 3; i++){ //look at all three vertices of our v_a_op triangle
				if ((acc[b->sim[v_a_op][i]] == acc[v_a_next] || acc[b->sim[v_a_op][i]] == acc[v_b]) && !partial_match){
					partial_match = 1;
				}
				else if ((acc[b->sim[v_a_op][i]] == acc[v_a_next] || acc[b->sim[v_a_op][i]] == acc[v_b]) && partial_match){
					double_match = 1;
				}
			}
			if (double_match){
				printf("NO SWITCH\n");
				//v_a_op is part of the subproblem between vertex v_a and vertex v_b
				//v_b_op is part of the subproblem between vertex v_b and vertex v_c
				//v_c_op is part of the subproblem between vertex v_c and vertex v_a

				//this is already what we have. no need to change things
			}
			else { //no double match. reset values and try with v_b_op
				partial_match = 0;
				double_match = 0;
				for (int i = 0; i < 3; i++){ //look at all three vertices of our v_b_op triangle
					if ((acc[b->sim[v_b_op][i]] == acc[v_a_next] || acc[b->sim[v_b_op][i]] == acc[v_b]) && !partial_match){
						partial_match = 1;
					}
					else if ((acc[b->sim[v_b_op][i]] == acc[v_a_next] || acc[b->sim[v_b_op][i]] == acc[v_b]) && partial_match){
						double_match = 1;
					}
				}
				if (double_match){
					printf("SWITCH 1\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_a;
					//int v_temp_next = v_a_next;
					//int v_temp_l = v_a_l;
					//int v_temp_r = v_a_r;
					int v_temp_op = v_a_op;
					//int v_temp_mirror = v_a_mirror;

					//v_a = v_b;
					//v_a_next = v_b_next;
					//v_a_l = v_b_l;
					//v_a_r = v_b_r;
					v_a_op = v_b_op;
					//v_a_mirror = v_b_mirror;

					//v_b = v_c;
					//v_b_next = v_c_next;
					//v_b_l = v_c_l;
					//v_b_r = v_c_r;
					v_b_op = v_c_op;
					//v_b_mirror = v_c_mirror;

					//v_c = v_temp;
					//v_c_next = v_temp_next;
					//v_c_l = v_temp_l;
					//v_c_r = v_temp_r;
					v_c_op = v_temp_op;
					//v_c_mirror = v_temp_mirror;
				}
				else { //no double match. v_c_op must match.
					printf("SWITCH 2\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_c;
					//int v_temp_next = v_c_next;
					//int v_temp_l = v_c_l;
					//int v_temp_r = v_c_r;
					int v_temp_op = v_c_op;
					//int v_temp_mirror = v_c_mirror;

					//v_c = v_b;
					//v_c_next = v_b_next;
					//v_c_l = v_b_l;
					//v_c_r = v_b_r;
					v_c_op = v_b_op;
					//v_c_mirror = v_b_mirror;

					//v_b = v_a;
					//v_b_next = v_a_next;
					//v_b_l = v_a_l;
					//v_b_r = v_a_r;
					v_b_op = v_a_op;
					//v_b_mirror = v_a_mirror;

					//v_a = v_temp;
					//v_a_next = v_temp_next;
					//v_a_l = v_temp_l;
					//v_a_r = v_temp_r;
					v_a_op = v_temp_op;
					//v_a_mirror = v_temp_mirror;
				}
			}
		}
		else if (v_a_next != -1 && v_b_next != -1){ //check to see where we hit the next tripod
			//acc[v_a_next] gives us the colour of the tripod that vertex a touches
			//acc[v_b_next] gives us the colour of the tripod that vertex b touches
			
			partial_match = 0;
			double_match = 0;
			for (int i = 0; i < 3; i++){ //look at all three vertices of our v_a_op triangle
				if ((acc[b->sim[v_a_op][i]] == acc[v_a_next] || acc[b->sim[v_a_op][i]] == acc[v_b_next]) && !partial_match){
					partial_match = 1;
				}
				else if ((acc[b->sim[v_a_op][i]] == acc[v_a_next] || acc[b->sim[v_a_op][i]] == acc[v_b_next]) && partial_match){
					double_match = 1;
				}
			}
			if (double_match){
				printf("NO SWITCH\n");
				//v_a_op is part of the subproblem between vertex v_a and vertex v_b
				//v_b_op is part of the subproblem between vertex v_b and vertex v_c
				//v_c_op is part of the subproblem between vertex v_c and vertex v_a

				//this is already what we have. no need to change things
			}
			else { //no double match. reset values and try with v_b_op
				partial_match = 0;
				double_match = 0;
				for (int i = 0; i < 3; i++){ //look at all three vertices of our v_b_op triangle
					if ((acc[b->sim[v_b_op][i]] == acc[v_a_next] || acc[b->sim[v_b_op][i]] == acc[v_b_next]) && !partial_match){
						partial_match = 1;
					}
					else if ((acc[b->sim[v_b_op][i]] == acc[v_a_next] || acc[b->sim[v_b_op][i]] == acc[v_b_next]) && partial_match){
						double_match = 1;
					}
				}
				if (double_match){
					printf("SWITCH 1\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_a;
					//int v_temp_next = v_a_next;
					//int v_temp_l = v_a_l;
					//int v_temp_r = v_a_r;
					int v_temp_op = v_a_op;
					//int v_temp_mirror = v_a_mirror;

					//v_a = v_b;
					//v_a_next = v_b_next;
					//v_a_l = v_b_l;
					//v_a_r = v_b_r;
					v_a_op = v_b_op;
					//v_a_mirror = v_b_mirror;

					//v_b = v_c;
					//v_b_next = v_c_next;
					//v_b_l = v_c_l;
					//v_b_r = v_c_r;
					v_b_op = v_c_op;
					//v_b_mirror = v_c_mirror;

					//v_c = v_temp;
					//v_c_next = v_temp_next;
					//v_c_l = v_temp_l;
					//v_c_r = v_temp_r;
					v_c_op = v_temp_op;
					//v_c_mirror = v_temp_mirror;
				}
				else { //no double match. v_c_op must match.
					printf("SWITCH 2\n");

					printf("v_a = %d\n", v_a);
					printf("v_b = %d\n", v_b);
					printf("v_c = %d\n", v_c);
					printf("v_a_mirror = %d\n", v_a_mirror);
					printf("v_b_mirror = %d\n", v_b_mirror);
					printf("v_c_mirror = %d\n", v_c_mirror);
					printf("v_a_op = %d\n", v_a_op);
					printf("v_b_op = %d\n", v_b_op);
					printf("v_c_op = %d\n", v_c_op);
					printf("----\n");

					//int v_temp = v_c;
					//int v_temp_next = v_c_next;
					//int v_temp_l = v_c_l;
					//int v_temp_r = v_c_r;
					int v_temp_op = v_c_op;
					//int v_temp_mirror = v_c_mirror;

					//v_c = v_b;
					//v_c_next = v_b_next;
					//v_c_l = v_b_l;
					//v_c_r = v_b_r;
					v_c_op = v_b_op;
					//v_c_mirror = v_b_mirror;

					//v_b = v_a;
					//v_b_next = v_a_next;
					//v_b_l = v_a_l;
					//v_b_r = v_a_r;
					v_b_op = v_a_op;
					//v_b_mirror = v_a_mirror;

					//v_a = v_temp;
					//v_a_next = v_temp_next;
					//v_a_l = v_temp_l;
					//v_a_r = v_temp_r;
					v_a_op = v_temp_op;
					//v_a_mirror = v_temp_mirror;
				}
			}
		}

		printf("v_a = %d\n", v_a);
		printf("v_b = %d\n", v_b);
		printf("v_c = %d\n", v_c);
		printf("v_a_mirror = %d\n", v_a_mirror);
		printf("v_b_mirror = %d\n", v_b_mirror);
		printf("v_c_mirror = %d\n", v_c_mirror);
		printf("v_a_op = %d\n", v_a_op);
		printf("v_b_op = %d\n", v_b_op);
		printf("v_c_op = %d\n", v_c_op);

		if (!(acc[v_a] != sp && acc[v_b] != sp)){ //if the subproblem exists //this might be too strong of an if statement
			if (acc[v_a] == sp && acc[v_b] == sp){ //leg a is non-empty && leg b is non-empty
				printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_b_r, v_a_l);
				tripod(b, r, v_a_op, v_b_r, v_a_l, acc);
			}
			else if (acc[v_a] == sp && acc[v_b] != sp){ //leg a is non-empty && leg b is empty
				printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_a_mirror, v_a_l);
				tripod(b, r, v_a_op, v_a_mirror, v_a_l, acc);
			}
			else if (acc[v_a] != sp && acc[v_b] == sp){ //leg a is empty && leg b is non-empty
				printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_b_r, v_a_mirror);
				tripod(b, r, v_a_op, v_b_r, v_a_mirror, acc);
			}
			/*if (acc[v_a] != sp && acc[v_b] != sp){ //leg a is empty && leg b is empty
				printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_a_mirror, v_a_mirror);
				tripod(b, r, v_a_op, v_a_mirror, v_a_mirror, acc);
			}*/
		}

		if (!(acc[v_b] != sp && acc[v_c] != sp)){ //if the subproblem exists
			if (acc[v_b] == sp && acc[v_c] == sp){ //leg b is non-empty && leg c is non-empty
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_c_r, v_b_l);
				tripod(b, r, v_b_op, v_c_r, v_b_l, acc);
			}
			else if (acc[v_b] == sp && acc[v_c] != sp){ //leg b is non-empty && leg c is empty
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_b_mirror, v_b_l);
				tripod(b, r, v_b_op, v_b_mirror, v_b_l, acc);
			}
			else if (acc[v_b] != sp && acc[v_c] == sp){ //leg b is empty && leg c is non-empty
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_c_r, v_b_mirror);
				tripod(b, r, v_b_op, v_c_r, v_b_mirror, acc);
			}
			/*if (acc[v_b] != sp && acc[v_c] != sp){ //leg b is empty && leg c is empty
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_b_mirror, v_b_mirror);
				tripod(b, r, v_b_op, v_b_mirror, v_b_mirror, acc);
			}*/
		}

		if (!(acc[v_c] != sp && acc[v_a] != sp)){ //if the subproblem exists
			if (acc[v_c] == sp && acc[v_a] == sp){ //leg c is non-empty && leg a is non-empty
				printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_a_r, v_c_l);
				tripod(b, r, v_c_op, v_a_r, v_c_l, acc);
			}
			else if (acc[v_c] == sp && acc[v_a] != sp){ //leg c is non-empty && leg a is empty
				printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_c_mirror, v_c_l);
				tripod(b, r, v_c_op, v_c_mirror, v_c_l, acc);
			}
			else if (acc[v_c] != sp && acc[v_a] == sp){ //leg c is empty && leg a is non-empty
				printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_a_r, v_c_mirror);
				tripod(b, r, v_c_op, v_a_r, v_c_mirror, acc);
			}
			/*if (acc[v_c] != sp && acc[v_a] != sp){ //leg c is empty && leg a is empty
				printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_c_mirror, v_c_mirror);
				tripod(b, r, v_c_op, v_c_mirror, v_c_mirror, acc);
			}*/
		}
	}

	printf("\n");
	return 0;
}
