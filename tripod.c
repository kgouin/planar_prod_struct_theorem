#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"tripod.h"

void init(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t){
	//initialize bfs_struct & rmq_struct
	BFS_init(b);
	BFS(b);
	LCA_init(r, b->ct, (2*(b->f))-1);
	printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	printf("initialization complete\n");

	printf("bfs tree = "); //testing
	for (int k = 0; k < b->v; k++){
		printf("%d ", b->bt[k]);
	}
	printf("\n");

	int* acc = malloc((b->v)*sizeof(int)); //memory leak
	int* acc2 = malloc((b->f)*sizeof(int)); //memory leak
	for (int k = 0; k < (b->v); k++){
		acc[k] = -1; //vertices labeled with -1 do not yet belong to a tripod
	}
	for (int k = 0; k < (b->f); k++){
		acc2[k] = -1; //faces labeled with -1 have not been identified as sperner triangles
	}

	//label the three edges incident to face 0 as belonging to three exterior tripods
	acc[b->sim[0][0]] = (b->f);
	acc[b->sim[0][1]] = (b->f)+1;
	acc[b->sim[0][2]] = (b->f)+2;

	//start decomposition with the three triangles adjacent to outer face
	trichromatic_tripod(b, r, t, b->tri[0][0], b->tri[0][2], b->tri[0][1], acc, acc2);
}

/******************************************************** trichromatic ********************************************************/

int* trichromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int f2, int f3, int* acc, int* acc2){
	//sperner triangle identification
	int sp;
	if (f1 == f2 && f1 == f3) return 0;
	else if (f1 != f2 && f1 != f3 && f2 != f3){
		//find sperner triangle
		if (LCA_query(r, f1, f2) == LCA_query(r, f1, f3) && LCA_query(r, f1, f3) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f2);
		else {
			if (LCA_query(r, f1, f2) == LCA_query(r, f1, f3)) sp = LCA_query(r, f2, f3);
			else if (LCA_query(r, f1, f2) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f3);
			else if (LCA_query(r, f1, f3) == LCA_query(r, f2, f3)) sp = LCA_query(r, f1, f2);
		}
	}
	else sp = f1;

	trichromatic_store_tripod(b, t, sp, f1, f2, f3, acc);

	printf("\n%d is one of our sperner triangles\nacc = ", sp); //sperner triangles are correctly identified
	for (int k = 0; k < (b->v); k++){
		printf("%d ", acc[k]);
	}
	printf("\n");

	acc2[sp] = sp; //keep track of sperner triangles

	trichromatic_orient_subproblems(b, t, sp, f1, f2, f3, acc);

	tprint(t);

	trichromatic_decompose(b, r, t, sp, f1, f2, f3, acc, acc2);

	return 0;
}

void trichromatic_store_tripod(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3, int* acc){
	int u;
	u = b->sim[sp][0];
	t->v_a = u;
	t->v_a_next = -1;
	while (acc[u] == -1){
		acc[u] = sp;
		t->v_a = u;
		u = b->bt[u];
		t->v_a_next = u;
	}

	u = b->sim[sp][1];
	t->v_b = u;
	t->v_b_next = -1;
	while (acc[u] == -1){
		acc[u] = sp;
		t->v_b = u;
		u = b->bt[u];
		t->v_b_next = u;
	}

	u = b->sim[sp][2];
	t->v_c = u;
	t->v_c_next = -1;
	while (acc[u] == -1){
		acc[u] = sp;
		t->v_c = u;
		u = b->bt[u];
		t->v_c_next = u;
	}

	t->v_a_op = f1;
	t->v_b_op = f2;
	t->v_c_op = f3;

	t->v_a_mirror = b->tri[sp][0];
	t->v_b_mirror = b->tri[sp][1];
	t->v_c_mirror = b->tri[sp][2];
}

void trichromatic_orient_subproblems(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3, int* acc){
	int double_match;

	if (t->v_a_next == -1 && t->v_b_next == -1){ //check to see where we hit the next tripod
		//acc[v_a] gives us the colour of the tripod that vertex a touches
		//acc[v_b] gives us the colour of the tripod that vertex b touches

		double_match = 0;
		for (int i = 0; i < 3; i++){ //look at all three vertices of our v_a_op triangle
			if (acc[b->sim[t->v_a_op][i]] == acc[t->v_a]){
				for (int j = 0; j < 3; j++){
					if (acc[b->sim[t->v_a_op][j]] == acc[t->v_b]){
						double_match = 1;
					}
				}
			}
		}
		if (double_match){
			printf("NO SWITCH\n");
		}
		if (!double_match) { //no double match. reset values and try with v_b_op
			double_match = 0;
			for (int i = 0; i < 3; i++){ //look at all three vertices of our v_b_op triangle
				if (acc[b->sim[t->v_b_op][i]] == acc[t->v_a]){
					for (int j = 0; j < 3; j++){
						if (acc[b->sim[t->v_b_op][j]] == acc[t->v_b]){
							double_match = 1;
						}
					}
				}
			}
			if (double_match){
				printf("SWITCH 1\n");
				tprint(t);

				int v_temp_op = t->v_a_op;
				t->v_a_op = t->v_b_op;
				t->v_b_op = t->v_c_op;
				t->v_c_op = v_temp_op;
			}
			else { //no double match. v_c_op must match.
				printf("SWITCH 2\n");
				tprint(t);

				int v_temp_op = t->v_c_op;
				t->v_c_op = t->v_b_op;
				t->v_b_op = t->v_a_op;
				t->v_a_op = v_temp_op;
			}
		}
	}

	else if (t->v_a_next == -1 && t->v_b_next != -1){ //check to see where we hit the next tripod
		//acc[v_a] gives us the colour of the tripod that vertex a touches
		//acc[v_b_next] gives us the colour of the tripod that vertex b touches

		double_match = 0;
		for (int i = 0; i < 3; i++){ //look at all three vertices of our v_a_op triangle
			if (acc[b->sim[t->v_a_op][i]] == acc[t->v_a]){
				for (int j = 0; j < 3; j++){
					if (acc[b->sim[t->v_a_op][j]] == acc[t->v_b_next]){
						double_match = 1;
					}
				}
			}
		}
		if (double_match){
			printf("NO SWITCH\n");
		}
		if (!double_match) { //no double match. reset values and try with v_b_op
			double_match = 0;
			for (int i = 0; i < 3; i++){ //look at all three vertices of our v_b_op triangle
				if (acc[b->sim[t->v_b_op][i]] == acc[t->v_a]){
					for (int j = 0; j < 3; j++){
						if (acc[b->sim[t->v_b_op][j]] == acc[t->v_b_next]){
							double_match = 1;
						}
					}
				}
			}
			if (double_match){
				printf("SWITCH 1\n");
				tprint(t);

				int v_temp_op = t->v_a_op;
				t->v_a_op = t->v_b_op;
				t->v_b_op = t->v_c_op;
				t->v_c_op = v_temp_op;
			}
			else { //no double match. v_c_op must match.
				printf("SWITCH 2\n");
				tprint(t);

				int v_temp_op = t->v_c_op;
				t->v_c_op = t->v_b_op;
				t->v_b_op = t->v_a_op;
				t->v_a_op = v_temp_op;
			}
		}
	}

	else if (t->v_a_next != -1 && t->v_b_next == -1){ //check to see where we hit the next tripod
		//acc[v_a_next] gives us the colour of the tripod that vertex a touches
		//acc[v_b] gives us the colour of the tripod that vertex b touches

		double_match = 0;
		for (int i = 0; i < 3; i++){ //look at all three vertices of our v_a_op triangle
			if (acc[b->sim[t->v_a_op][i]] == acc[t->v_a_next]){
				for (int j = 0; j < 3; j++){
					if (acc[b->sim[t->v_a_op][j]] == acc[t->v_b]){
						double_match = 1;
					}
				}
			}
		}
		if (double_match){
			printf("NO SWITCH\n");
		}
		if (!double_match) { //no double match. reset values and try with v_b_op
			double_match = 0;
			for (int i = 0; i < 3; i++){ //look at all three vertices of our v_b_op triangle
				if (acc[b->sim[t->v_b_op][i]] == acc[t->v_a_next]){
					for (int j = 0; j < 3; j++){
						if (acc[b->sim[t->v_b_op][j]] == acc[t->v_b]){
							double_match = 1;
						}
					}
				}
			}
			if (double_match){
				printf("SWITCH 1\n");
				tprint(t);

				int v_temp_op = t->v_a_op;
				t->v_a_op = t->v_b_op;
				t->v_b_op = t->v_c_op;
				t->v_c_op = v_temp_op;
			}
			else { //no double match. v_c_op must match.
				printf("SWITCH 2\n");
				tprint(t);

				int v_temp_op = t->v_c_op;
				t->v_c_op = t->v_b_op;
				t->v_b_op = t->v_a_op;
				t->v_a_op = v_temp_op;
			}
		}
	}

	else if (t->v_a_next != -1 && t->v_b_next != -1){ //check to see where we hit the next tripod
		//acc[v_a_next] gives us the colour of the tripod that vertex a touches
		//acc[v_b_next] gives us the colour of the tripod that vertex b touches

		double_match = 0;
		for (int i = 0; i < 3; i++){ //look at all three vertices of our v_a_op triangle
			if (acc[b->sim[t->v_a_op][i]] == acc[t->v_a_next]){
				for (int j = 0; j < 3; j++){
					if (acc[b->sim[t->v_a_op][j]] == acc[t->v_b_next]){
						double_match = 1;
					}
				}
			}
		}
		if (double_match){
			printf("NO SWITCH\n");
		}
		if (!double_match) { //no double match. reset values and try with v_b_op
			double_match = 0;
			for (int i = 0; i < 3; i++){ //look at all three vertices of our v_b_op triangle
				if (acc[b->sim[t->v_b_op][i]] == acc[t->v_a_next]){
					for (int j = 0; j < 3; j++){
						if (acc[b->sim[t->v_b_op][j]] == acc[t->v_b_next]){
							double_match = 1;
						}
					}
				}
			}
			if (double_match){
				printf("SWITCH 1\n");
				tprint(t);

				int v_temp_op = t->v_a_op;
				t->v_a_op = t->v_b_op;
				t->v_b_op = t->v_c_op;
				t->v_c_op = v_temp_op;
			}
			else { //no double match. v_c_op must match.
				printf("SWITCH 2\n");
				tprint(t);

				int v_temp_op = t->v_c_op;
				t->v_c_op = t->v_b_op;
				t->v_b_op = t->v_a_op;
				t->v_a_op = v_temp_op;
			}
		}
	}
}

void trichromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3, int* acc, int* acc2){
	//definitions
	int v_a = t->v_a;
	int v_a_next = t->v_a_next;
	int v_a_l = t->v_a_l;
	int v_a_r = t-> v_a_r;
	int v_a_op = t->v_a_op;
	int v_a_mirror = t->v_a_mirror;
	int y_a = t->y_a;
	int v_b = t->v_b;
	int v_b_next = t->v_b_next;
	int v_b_l = t->v_b_l;
	int v_b_r = t->v_b_r;
	int v_b_op = t->v_b_op;
	int v_b_mirror = t->v_b_mirror;
	int y_b = t->y_b;
	int v_c = t->v_c;
	int v_c_next = t->v_c_next;
	int v_c_l = t->v_c_l;
	int v_c_r = t->v_c_r;
	int v_c_op = t->v_c_op;
	int v_c_mirror = t->v_c_mirror;
	int y_c = t->y_c;

	int a_flag = 0;
	int b_flag = 0;
	int a_np = 0;
	int b_np = 0;

	//subproblem a
	if (acc[v_a] == sp && acc[v_b] == sp){ //leg a is non-empty && leg b is non-empty
		v_a_l = b->il[v_a][b->pin[v_a]];
		(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
		v_b_r = b->il[v_b][y_b];
		printf("subproblem a for sp %d is trichromatic\n", sp);
		printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_b_r, v_a_l);
		trichromatic_tripod(b, r, t, v_a_op, v_b_r, v_a_l, acc, acc2);
		a_flag = 1;
	}
	else if (acc[v_a] == sp && acc[v_b] != sp){ //leg a is non-empty && leg b is empty
		v_a_l = b->il[v_a][b->pin[v_a]];
		printf("subproblem a for sp %d is trichromatic\n", sp);
		printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_a_mirror, v_a_l);
		trichromatic_tripod(b, r, t, v_a_op, v_a_mirror, v_a_l, acc, acc2);
		a_flag = 1;
	}
	else if (acc[v_a] != sp && acc[v_b] == sp){ //leg a is empty && leg b is non-empty
		(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
		v_b_r = b->il[v_b][y_b];
		printf("subproblem a for sp %d is trichromatic\n", sp);
		printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_b_r, v_a_mirror);
		trichromatic_tripod(b, r, t, v_a_op, v_b_r, v_a_mirror, acc, acc2);
		a_flag = 1;
	}
	else if (acc[v_a] != sp && acc[v_b] != sp){ //leg a is empty && leg b is empty
		if (sp == v_a_op){
			printf("no subproblem a for sp %d\n", sp);
			a_np = 1;
		}
		else {
			printf("subproblem a for sp %d is bichromatic\n", sp);
			printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_a_op, v_a_mirror);
			bichromatic_tripod( b, r, t, v_a_op, v_a_mirror, acc, acc2);
			a_flag = 1;
		}
	}

	//subproblem b
	if (acc[v_b] == sp && acc[v_c] == sp){ //leg b is non-empty && leg c is non-empty
		v_b_l = b->il[v_b][b->pin[v_b]];
		(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
		v_c_r = b->il[v_c][y_c];
		printf("subproblem b for sp %d is trichromatic\n", sp);
		printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_c_r, v_b_l);
		trichromatic_tripod(b, r, t, v_b_op, v_c_r, v_b_l, acc, acc2);
		b_flag = 1;
	}
	else if (acc[v_b] == sp && acc[v_c] != sp){ //leg b is non-empty && leg c is empty
		v_b_l = b->il[v_b][b->pin[v_b]];
		printf("subproblem b for sp %d is trichromatic\n", sp);
		printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_b_mirror, v_b_l);
		trichromatic_tripod(b, r, t, v_b_op, v_b_mirror, v_b_l, acc, acc2);
		b_flag = 1;
	}
	else if (acc[v_b] != sp && acc[v_c] == sp){ //leg b is empty && leg c is non-empty
		(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
		v_c_r = b->il[v_c][y_c];
		printf("subproblem b for sp %d is trichromatic\n", sp);
		printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_c_r, v_b_mirror);
		trichromatic_tripod(b, r, t, v_b_op, v_c_r, v_b_mirror, acc, acc2);
		b_flag = 1;
	}
	else if (acc[v_b] != sp && acc[v_c] != sp){ //leg b is empty && leg c is empty
		if (sp == v_b_op){
			printf("no subproblem b for sp %d\n", sp);
			b_np = 1;
		}
		else {
			int value;
			(a_flag) ? (value = v_b_mirror) : (value = v_a_mirror);
			printf("subproblem b for sp %d is bichromatic\n", sp);
			printf("subproblem b for sp %d will be on faces %d, %d\n", sp, v_b_op, value);
			bichromatic_tripod( b, r, t, v_b_op, value, acc, acc2);
			b_flag = 1;
		}
	}

	//subproblem c
	if (acc[v_c] == sp && acc[v_a] == sp){ //leg c is non-empty && leg a is non-empty
		v_c_l = b->il[v_c][b->pin[v_c]];
		(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
		v_a_r = b->il[v_a][y_a];
		printf("subproblem c for sp %d is trichromatic\n", sp);
		printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_a_r, v_c_l);
		trichromatic_tripod(b, r, t, v_c_op, v_a_r, v_c_l, acc, acc2);
	}
	else if (acc[v_c] == sp && acc[v_a] != sp){ //leg c is non-empty && leg a is empty
		v_c_l = b->il[v_c][b->pin[v_c]];
		printf("subproblem c for sp %d is trichromatic\n", sp);
		printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_c_mirror, v_c_l);
		trichromatic_tripod(b, r, t, v_c_op, v_c_mirror, v_c_l, acc, acc2);
	}
	else if (acc[v_c] != sp && acc[v_a] == sp){ //leg c is empty && leg a is non-empty
		(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
		v_a_r = b->il[v_a][y_a];
		printf("subproblem c for sp %d is trichromatic\n", sp);
		printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_a_r, v_c_mirror);
		trichromatic_tripod(b, r, t, v_c_op, v_a_r, v_c_mirror, acc, acc2);
	}
	else if (acc[v_c] != sp && acc[v_a] != sp){ //leg c is empty && leg a is empty
		if (sp == v_c_op){
			printf("no subproblem c for sp %d\n", sp);
		}
		else {
			int value;
			//if (a_np && b_np) value = v_a_mirror;
			//else if ((a_np && !b_np) || (!a_np && b_np)) value = v_b_mirror;
			//else if (!a_np && !b_np) value = v_c_mirror;
			//not sure if the lines below are the way to approach the incorrect value I was getting...
			if (a_np && b_np) (acc2[v_a_mirror] == -1) ? (value = v_a_mirror) : ((acc2[v_b_mirror] == -1) ? (value = v_b_mirror) : (value = v_c_mirror));
			else if ((a_np && !b_np) || (!a_np && b_np)) (acc2[v_b_mirror] == -1) ? (value = v_b_mirror) : ((acc2[v_c_mirror] == -1) ? (value = v_c_mirror) : (value = v_a_mirror));
			else if (!a_np && !b_np) (acc2[v_c_mirror] == -1) ? (value = v_c_mirror) : ((acc2[v_a_mirror] == -1) ? (value = v_a_mirror) : (value = v_b_mirror));
			//something off here...
			printf("OVER HERE\n");
			printf("v_a = %d\n", v_a);
			printf("v_b = %d\n", v_b);
			printf("v_c = %d\n", v_c);
			printf("v_a_next = %d\n", v_a_next);
			printf("v_b_next = %d\n", v_b_next);
			printf("v_c_next = %d\n", v_c_next);
			printf("v_a_mirror = %d\n", v_a_mirror);
			printf("v_b_mirror = %d\n", v_b_mirror);
			printf("v_c_mirror = %d\n", v_c_mirror);
			printf("v_a_op = %d\n", v_a_op);
			printf("v_b_op = %d\n", v_b_op);
			printf("v_c_op = %d\n", v_c_op);
			//printf("v_a_mirror = %d, v_b_mirror = %d, v_c_mirror = %d\n", v_a_mirror, v_b_mirror, v_c_mirror);
			printf("subproblem c for sp %d is bichromatic\n", sp);
			printf("subproblem c for sp %d will be on faces %d, %d\n", sp, v_c_op, value);
			bichromatic_tripod( b, r, t, v_c_op, value, acc, acc2);
		}
	}
}

/******************************************************** bichromatic ********************************************************/

int* bichromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int f2, int* acc, int* acc2){
	int sp;
	
	//if (acc2[f1] != -1 || acc2[f2] != -1) return 0; //this should not be needed
	if (f1 == f2) return 0;
	sp = f1;

	bichromatic_store_tripod(b, t, sp, f1, f2, acc);

	printf("\n%d is one of our sperner triangles\nacc = ", sp); //sperner triangles are correctly identified
	for (int k = 0; k < (b->v); k++){
		printf("%d ", acc[k]);
	}
	printf("\n");

	acc2[sp] = sp; //keep track of sperner triangles

	printf("v_a = %d\n", t->v_a);
	printf("v_b = %d\n", t->v_b);
	printf("v_c = %d\n", t->v_c);
	printf("v_a_next = %d\n", t->v_a_next);
	printf("v_b_next = %d\n", t->v_b_next);
	printf("v_c_next = %d\n", t->v_c_next);
	printf("v_a_mirror = %d\n", t->v_a_mirror);
	printf("v_b_mirror = %d\n", t->v_b_mirror);
	printf("v_c_mirror = %d\n", t->v_c_mirror);

	//no need to orient

	bichromatic_decompose(b, r, t, sp, f1, f2, acc, acc2);

	return 0;
}

void bichromatic_store_tripod(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int* acc){
	int u;
	u = b->sim[sp][0];
	t->v_a = u;
	t->v_a_next = -1;
	while (acc[u] == -1){
		acc[u] = sp;
		t->v_a = u;
		u = b->bt[u];
		t->v_a_next = u;
	}

	u = b->sim[sp][1];
	t->v_b = u;
	t->v_b_next = -1;
	while (acc[u] == -1){
		acc[u] = sp;
		t->v_b = u;
		u = b->bt[u];
		t->v_b_next = u;
	}

	u = b->sim[sp][2];
	t->v_c = u;
	t->v_c_next = -1;
	while (acc[u] == -1){
		acc[u] = sp;
		t->v_c = u;
		u = b->bt[u];
		t->v_c_next = u;
	}

	t->v_a_mirror = b->tri[sp][0];
	t->v_b_mirror = b->tri[sp][1];
	t->v_c_mirror = b->tri[sp][2];
}

void bichromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int* acc, int* acc2){
	if (f2 == sp) exit(0); //something wrong. abort

	//definitions
	int v_a = t->v_a;
	int v_a_next = t->v_a_next;
	int v_a_l = t->v_a_l;
	int v_a_r = t-> v_a_r;
	int v_a_op = t->v_a_op;
	int v_a_mirror = t->v_a_mirror;
	int y_a = t->y_a;
	int v_b = t->v_b;
	int v_b_next = t->v_b_next;
	int v_b_l = t->v_b_l;
	int v_b_r = t->v_b_r;
	int v_b_op = t->v_b_op;
	int v_b_mirror = t->v_b_mirror;
	int y_b = t->y_b;
	int v_c = t->v_c;
	int v_c_next = t->v_c_next;
	int v_c_l = t->v_c_l;
	int v_c_r = t->v_c_r;
	int v_c_op = t->v_c_op;
	int v_c_mirror = t->v_c_mirror;
	int y_c = t->y_c;

	//all subproblems
	if (acc[v_a] == sp && acc[v_b] != sp){ //leg a is non-empty && leg b is empty
		v_a_l = b->il[v_a][b->pin[v_a]];
		(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
		v_a_r = b->il[v_a][y_a];
		if (!(v_a_next == v_b || v_a_next == v_c)){ //if the path up the bfs tree from v_a does NOT lead to v_b or v_c
			if (acc[v_a_next] == acc[v_b]){
				if (v_a_l != f2 && v_a_mirror != f2){ //new check
					printf("subproblem a for sp %d is bichromatic\n", sp);
					printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_a_l, v_a_mirror);
					bichromatic_tripod( b, r, t, v_a_l, v_a_mirror, acc, acc2);
				}
				else return;

				printf("subproblem b for sp %d is trichromatic\n", sp);
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, f2, v_a_r, v_c_mirror);
				trichromatic_tripod( b, r, t, f2, v_a_r, v_c_mirror, acc, acc2);
			}
			else {
				if (v_c_mirror != f2 && v_a_r != f2){ //new check
					printf("subproblem a for sp %d is bichromatic\n", sp);
					printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_c_mirror, v_a_r);
					bichromatic_tripod( b, r, t, v_c_mirror, v_a_r, acc, acc2);
				}
				else return;

				printf("subproblem b for sp %d is trichromatic\n", sp);
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, f2, v_a_mirror, v_a_l);
				trichromatic_tripod( b, r, t, f2, v_a_mirror, v_a_l, acc, acc2);
			}
		}
		else { //the path up the bfs tree from v_a leads to v_b or v_c
			printf("subproblem for sp %d is trichromatic\n", sp);
			printf("subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_a_mirror, v_c_mirror);
			trichromatic_tripod( b, r, t, f2, v_a_mirror, v_c_mirror, acc, acc2);
		}
	}
	else if (acc[v_a] != sp && acc[v_b] == sp){ //leg a is empty && leg b is non-empty
		v_b_l = b->il[v_b][b->pin[v_b]];
		(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
		v_b_r = b->il[v_b][y_b];
		if (!(v_b_next == v_c || v_b_next == v_a)){ //if the path up the bfs tree from v_b does NOT lead to v_c or v_a
			if (acc[v_b_next] == acc[v_c]){
				if (v_b_l != f2 && v_b_mirror != f2){ //new check
					printf("subproblem a for sp %d is bichromatic\n", sp);
					printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_b_l, v_b_mirror);
					//we're not recalculating v_b_l
					bichromatic_tripod( b, r, t, v_b_l, v_b_mirror, acc, acc2);
				}
				else return;

				printf("subproblem b for sp %d is trichromatic\n", sp);
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, f2, v_b_r, v_a_mirror);
				trichromatic_tripod( b, r, t, f2, v_b_r, v_a_mirror, acc, acc2);
				
			}
			else {
				if (v_a_mirror != f2 && v_b_r != f2){ //new check
					printf("subproblem a for sp %d is bichromatic\n", sp);
					printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_a_mirror, v_b_r);
					bichromatic_tripod( b, r, t, v_a_mirror, v_b_r, acc, acc2);
				}
				else return;

				printf("subproblem b for sp %d is trichromatic\n", sp);
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, f2, v_b_mirror, v_b_l);
				trichromatic_tripod( b, r, t, f2, v_b_mirror, v_b_l, acc, acc2);
			}
		}
		else { //the path up the bfs tree from v_b leads to v_c or v_a
			printf("subproblem for sp %d is trichromatic\n", sp);
			printf("subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_b_mirror, v_a_mirror);
			trichromatic_tripod( b, r, t, f2, v_b_mirror, v_a_mirror, acc, acc2);
		}
	}
	else if (acc[v_a] != sp && acc[v_b] != sp){ //leg a is empty && leg b is empty
		if (acc[v_c] == acc[v_a]){ //if v_c is the same colour as v_a
			if (v_a_mirror != f2 && v_b_mirror != f2 && v_c_mirror != f2){ //if v_X_mirror == f2, then stop
				printf("here1\n");
				if (acc[v_c] == sp){ //leg c is not empty //maybe add something like this in more places
					printf("subproblem a for sp %d is monochromatic\n", sp);
					printf("subproblem a for sp %d will be on face %d\n", sp, v_c_mirror);
					//monochromatic_tripod( b, r, t, v_c_mirror, acc, acc2);
				}

				printf("subproblem b for sp %d is bichromatic\n", sp);
				printf("subproblem b for sp %d will be on faces %d, %d\n", sp, f2, v_a_mirror);
				bichromatic_tripod( b, r, t, f2, v_a_mirror, acc, acc2);
			}
		}
		else if (acc[v_c] == acc[v_b]){ //if v_c is the same colour as v_b
			printf("here2\n");
			if (v_a_mirror != f2 && v_b_mirror != f2 && v_c_mirror != f2){
				if (acc[v_c] == sp){ //leg c is not empty
					printf("subproblem a for sp %d is monochromatic\n", sp);
					printf("subproblem a for sp %d will be on face %d\n", sp, v_b_mirror);
					//monochromatic_tripod( b, r, t, v_b_mirror, acc, acc2);
				}

				printf("subproblem b for sp %d is bichromatic\n", sp);
				printf("subproblem b for sp %d will be on faces %d, %d\n", sp, f2, v_c_mirror);
				bichromatic_tripod( b, r, t, f2, v_c_mirror, acc, acc2);
			}
		}
		else { //v_c does not touch the cycle defining its subproblem(s)
			if (!(v_c_next == v_a || v_c_next == v_b)){ //if the path up the bfs tree from v_c does NOT lead to v_a or v_b
				if (acc[v_c] == sp){
					if (acc[v_c_next] == acc[v_a]){ //invalid read
						if (v_c_l != f2 && v_c_mirror != f2){ //new check
							printf("subproblem a for sp %d is bichromatic\n", sp);
							printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_c_l, v_c_mirror);
							bichromatic_tripod( b, r, t, v_c_l, v_c_mirror, acc, acc2);
						}
						else return;

						printf("subproblem b for sp %d is trichromatic\n", sp);
						printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, f2, v_c_r, v_b_mirror);
						trichromatic_tripod( b, r, t, f2, v_c_r, v_b_mirror, acc, acc2);
					}
				}
				else {
					if (v_b_mirror != f2 && v_c_r != f2){ //new check
						printf("subproblem a for sp %d is bichromatic\n", sp);
						printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_b_mirror, v_c_r);
						bichromatic_tripod( b, r, t, v_b_mirror, v_c_r, acc, acc2);
					}
					else return;

					printf("subproblem b for sp %d is trichromatic\n", sp);
					printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, f2, v_c_mirror, v_c_l);
					trichromatic_tripod( b, r, t, f2, v_c_mirror, v_c_l, acc, acc2);
				}
			}
			else { //the path up the bfs tree from v_c leads to v_a or v_b
				printf("subproblem for sp %d is trichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_c_mirror, v_b_mirror);
				trichromatic_tripod( b, r, t, f2, v_c_mirror, v_b_mirror, acc, acc2);
			}
		}
	}
}

/******************************************************** monochromatic ********************************************************/
//

void tprint(struct tripod_decomposition_struct* t){
	printf("\n");
	printf("v_a = %d\n", t->v_a);
	printf("v_b = %d\n", t->v_b);
	printf("v_c = %d\n", t->v_c);
	//printf("v_a_r = %d\n", t->v_a_r);
	//printf("v_b_r = %d\n", t->v_b_r);
	//printf("v_c_r = %d\n", t->v_c_r);
	//printf("v_a_l = %d\n", t->v_a_l);
	//printf("v_b_l = %d\n", t->v_b_l);
	//printf("v_c_l = %d\n", t->v_c_l);
	printf("v_a_next = %d\n", t->v_a_next);
	printf("v_b_next = %d\n", t->v_b_next);
	printf("v_c_next = %d\n", t->v_c_next);
	printf("v_a_mirror = %d\n", t->v_a_mirror);
	printf("v_b_mirror = %d\n", t->v_b_mirror);
	printf("v_c_mirror = %d\n", t->v_c_mirror);
	printf("v_a_op = %d\n", t->v_a_op);
	printf("v_b_op = %d\n", t->v_b_op);
	printf("v_c_op = %d\n", t->v_c_op);
	printf("\n");
}
