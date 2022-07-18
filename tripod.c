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
	tripod(b, r, acc[b->sim[0][2]], b->tri[0][0], b->tri[0][1], b->tri[0][2], acc);
}

int* tripod(struct bfs_struct* b, struct rmq_struct* r, int old_sp, int f1, int f2, int f3, int* acc){ //recursive function
	//definitions
	int sp;
	int u;
	int v_a;
	int v_a_next;
	int v_a_l;
	int v_a_r;
	int v_a_op;
	int y_a;
	int v_b;
	int v_b_next;
	int v_b_l;
	int v_b_r;
	int v_b_op;
	int y_b;
	int v_c;
	int v_c_next;
	int v_c_l;
	int v_c_r;
	int v_c_op;
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
		v_a_op = b->tri[sp][0];
		v_b_op = b->tri[sp][1];
		v_c_op = b->tri[sp][2];

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
		if (acc[v_a] == sp || acc[v_b] == sp){
				if (acc[v_c] == old_sp && acc[v_c_next] == old_sp){
					printf("sub-problem 1 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_a, v_a_op, v_a_l, f1);
					tripod(b, r, sp, v_a_op, v_a_l, f1, acc);
				}
				else if (acc[v_a] == old_sp && acc[v_a_next] == old_sp){
					printf("sub-problem 1 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_a, v_a_op, v_a_l, f2);
					tripod(b, r, sp, v_a_op, v_a_l, f2, acc);
				}
				else if (acc[v_b] == old_sp && acc[v_b_next] == old_sp){
					printf("sub-problem 1 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_a, v_a_op, v_a_l, f3);
					tripod(b, r, sp, v_a_op, v_a_l, f3, acc);
				}
			else {
				printf("sub-problem 1 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_a, v_a_op, v_a_l, v_b_r);
				tripod(b, r, sp, v_a_op, v_a_l, v_b_r, acc);
			}
		}

		if (acc[v_b] == sp || acc[v_c] == sp){
				if (acc[v_b] == old_sp && acc[v_b_next] == old_sp){
					printf("sub-problem 2 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_b, v_b_op, v_b_l, f1);
					tripod(b, r, sp, v_b_op, v_b_l, f1, acc);
				}
				else if (acc[v_c] == old_sp && acc[v_c_next] == old_sp){
					printf("sub-problem 2 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_b, v_b_op, v_b_l, f2);
					tripod(b, r, sp, v_b_op, v_b_l, f2, acc);
				}
				else if (acc[v_a] == old_sp && acc[v_a_next] == old_sp){
					printf("sub-problem 2 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_b, v_b_op, v_b_l, f3);
					tripod(b, r, sp, v_b_op, v_b_l, f3, acc);
				}
			else {
				printf("sub-problem 2 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_b, v_b_op, v_b_l, v_c_r);
				tripod(b, r, sp, v_b_op, v_b_l, v_c_r, acc);
			}
		}

		if (acc[v_c] == sp || acc[v_a] == sp){ //check that the subproblem we're looking at est pas colle contre d'autres tripods
											   //aka. that there is indeed a subproblem here
				if (acc[v_a] == old_sp && acc[v_a_next] == old_sp){
					printf("sub-problem 3 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_c, v_c_op, v_c_l, f1);
					tripod(b, r, sp, v_c_op, v_c_l, f1, acc);
				}
				else if (acc[v_b] == old_sp && acc[v_b_next] == old_sp){
					printf("sub-problem 3 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_c, v_c_op, v_c_l, f2);
					tripod(b, r, sp, v_c_op, v_c_l, f2, acc);
				}
				else if (acc[v_c] == old_sp && acc[v_c_next] == old_sp){
					printf("sub-problem 3 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_c, v_c_op, v_c_l, f3);
					tripod(b, r, sp, v_c_op, v_c_l, f3, acc);
				}
			else {
				printf("sub-problem 3 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_c, v_c_op, v_c_l, v_a_r);
				tripod(b, r, sp, v_c_op, v_c_l, v_a_r, acc);
			}
		}
	}

	else {
		//find sperner triangle
		sp = f1;

		//store tripod
		v_a_op = b->tri[sp][0];
		v_b_op = b->tri[sp][1];
		v_c_op = b->tri[sp][2];

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
		if (acc[v_a] == sp || acc[v_b] == sp){
				if (acc[v_c] == old_sp && acc[v_c_next] == old_sp){
					printf("sub-problem 1 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_a, v_a_op, v_a_l, f1);
					tripod(b, r, sp, v_a_op, v_a_l, f1, acc);
				}
				else if (acc[v_a] == old_sp && acc[v_a_next] == old_sp){
					printf("sub-problem 1 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_a, v_a_op, v_a_l, f2);
					tripod(b, r, sp, v_a_op, v_a_l, f2, acc);
				}
				else if (acc[v_b] == old_sp && acc[v_b_next] == old_sp){
					printf("sub-problem 1 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_a, v_a_op, v_a_l, f3);
					tripod(b, r, sp, v_a_op, v_a_l, f3, acc);
				}
			else {
				printf("sub-problem 1 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_a, v_a_op, v_a_l, v_b_r);
				tripod(b, r, sp, v_a_op, v_a_l, v_b_r, acc);
			}
		}

		if (acc[v_b] == sp || acc[v_c] == sp){
				if (acc[v_b] == old_sp && acc[v_b_next] == old_sp){
					printf("sub-problem 2 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_b, v_b_op, v_b_l, f1);
					tripod(b, r, sp, v_b_op, v_b_l, f1, acc);
				}
				else if (acc[v_c] == old_sp && acc[v_c_next] == old_sp){
					printf("sub-problem 2 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_b, v_b_op, v_b_l, f2);
					tripod(b, r, sp, v_b_op, v_b_l, f2, acc);
				}
				else if (acc[v_a] == old_sp && acc[v_a_next] == old_sp){
					printf("sub-problem 2 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_b, v_b_op, v_b_l, f3);
					tripod(b, r, sp, v_b_op, v_b_l, f3, acc);
				}
			else {
				printf("sub-problem 2 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_b, v_b_op, v_b_l, v_c_r);
				tripod(b, r, sp, v_b_op, v_b_l, v_c_r, acc);
			}
		}

		if (acc[v_c] == sp || acc[v_a] == sp){ //check that the subproblem we're looking at est pas colle contre d'autres tripods
											   //aka. that there is indeed a subproblem here
				if (acc[v_a] == old_sp && acc[v_a_next] == old_sp){
					printf("sub-problem 3 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_c, v_c_op, v_c_l, f1);
					tripod(b, r, sp, v_c_op, v_c_l, f1, acc);
				}
				else if (acc[v_b] == old_sp && acc[v_b_next] == old_sp){
					printf("sub-problem 3 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_c, v_c_op, v_c_l, f2);
					tripod(b, r, sp, v_c_op, v_c_l, f2, acc);
				}
				else if (acc[v_c] == old_sp && acc[v_c_next] == old_sp){
					printf("sub-problem 3 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_c, v_c_op, v_c_l, f3);
					tripod(b, r, sp, v_c_op, v_c_l, f3, acc);
				}
			else {
				printf("sub-problem 3 for sp %d at vertex %d will be on faces %d, %d, %d\n", sp, v_c, v_c_op, v_c_l, v_a_r);
				tripod(b, r, sp, v_c_op, v_c_l, v_a_r, acc);
			}
		}
	}

	printf("\n");
	return 0;
}
