#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"tripod.h"

void tripod_init(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t){
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
	tripod(b, r, t, b->tri[0][0], b->tri[0][2], b->tri[0][1], acc, acc2);
}

int* tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int f2, int f3, int* acc, int* acc2){ //recursive function
	//definitions
	int sp; //sperner triangle

	if (f1 == f2 && f1 == f3 && f2 == f3) return 0;

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

	store_tripod(b, t, sp, f1, f2, f3, acc);

	printf("\n%d is one of our sperner triangles\nacc = ", sp); //sperner triangles are correctly identified
	for (int k = 0; k < (b->v); k++){
		printf("%d ", acc[k]);
	}
	printf("\n");

	acc2[sp] = sp; //keep track of sperner triangles

	orient_subproblems(b, t, sp, f1, f2, f3, acc);

	tprint(t);

	decompose(b, r, t, sp, f1, f2, f3, acc, acc2);

	return 0;
}

void store_tripod(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3, int* acc){
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

void orient_subproblems(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3, int* acc){
	int partial_match;
	int double_match;

	if (t->v_a_next == -1 && t->v_b_next == -1){ //check to see where we hit the next tripod
		//acc[v_a] gives us the colour of the tripod that vertex a touches
		//acc[v_b] gives us the colour of the tripod that vertex b touches

		partial_match = 0;
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
			printf("01\n");
			//v_a_op is part of the subproblem between vertex v_a and vertex v_b
			//v_b_op is part of the subproblem between vertex v_b and vertex v_c
			//v_c_op is part of the subproblem between vertex v_c and vertex v_a
			//this is already what we have. no need to change things
		}
		if (!double_match) { //no double match. reset values and try with v_b_op
			partial_match = 0;
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
				printf("01\n");
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

		partial_match = 0;
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
			printf("02\n");
			//v_a_op is part of the subproblem between vertex v_a and vertex v_b
			//v_b_op is part of the subproblem between vertex v_b and vertex v_c
			//v_c_op is part of the subproblem between vertex v_c and vertex v_a
			//this is already what we have. no need to change things
		}
		if (!double_match) { //no double match. reset values and try with v_b_op
			partial_match = 0;
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
				printf("02\n");
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

		partial_match = 0;
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
			printf("03\n");
			//v_a_op is part of the subproblem between vertex v_a and vertex v_b
			//v_b_op is part of the subproblem between vertex v_b and vertex v_c
			//v_c_op is part of the subproblem between vertex v_c and vertex v_a
			//this is already what we have. no need to change things
		}
		if (!double_match) { //no double match. reset values and try with v_b_op
			partial_match = 0;
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
				printf("03\n");
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

		partial_match = 0;
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
			printf("04\n");
			//v_a_op is part of the subproblem between vertex v_a and vertex v_b
			//v_b_op is part of the subproblem between vertex v_b and vertex v_c
			//v_c_op is part of the subproblem between vertex v_c and vertex v_a
			//this is already what we have. no need to change things
		}
		if (!double_match) { //no double match. reset values and try with v_b_op
			partial_match = 0;
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
				printf("04\n");
				tprint(t);

				int v_temp_op = t->v_c_op;
				t->v_c_op = t->v_b_op;
				t->v_b_op = t->v_a_op;
				t->v_a_op = v_temp_op;
			}
		}
	}
}

void decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3, int* acc, int* acc2){
	int a_incomplete = 1;
	int b_incomplete = 1;
	int c_incomplete = 1;
	int stop = 0;

	//if (f1 == 34 && f2 == 34 && f3 == 35) exit(0);

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

	if (acc[v_a] == sp && acc[v_b] == sp){ //leg a is non-empty && leg b is non-empty
		a_incomplete = 0;
		v_a_l = b->il[v_a][b->pin[v_a]];
		(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
		v_b_r = b->il[v_b][y_b];
		//no need to check for bichromatic cycle
		printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_b_r, v_a_l);
		tripod(b, r, t, v_a_op, v_b_r, v_a_l, acc, acc2);
	}
	else if (acc[v_a] == sp && acc[v_b] != sp){ //leg a is non-empty && leg b is empty
		a_incomplete = 0;
		//conditional below checks whether cycle defining subproblem is bichromatic. aka. if there is only one subproblem
		if (v_a_next == v_b || v_a_next == v_c) { //if the path up the bfs tree from v_a leads to v_b or v_c
			(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
			v_a_r = b->il[v_a][y_a];
			printf("the only subproblem for sp %d will be on faces %d, %d, %d\n", sp, v_c_mirror, v_a_r, v_a_mirror);
			tripod(b, r, t, v_c_mirror, v_a_r, v_a_mirror, acc, acc2);
			stop = 1;
		}
		else {
			v_a_l = b->il[v_a][b->pin[v_a]];
			printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_a_mirror, v_a_l);
			tripod(b, r, t, v_a_op,v_a_mirror, v_a_l, acc, acc2);
		}
	}
	else if (acc[v_a] != sp && acc[v_b] == sp){ //leg a is empty && leg b is non-empty
		a_incomplete = 0;
		//conditional below checks whether cycle defining subproblem is bichromatic. aka. if there is only one subproblem
		if (v_b_next == v_c || v_b_next == v_a){ //if the path up the bfs tree from v_b leads to v_c or v_a
			(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
			v_b_r = b->il[v_b][y_b];
			printf("the only subproblem for sp %d will be on faces %d, %d, %d\n", sp, v_a_mirror, v_b_r, v_b_mirror);
			tripod(b, r, t, v_a_mirror, v_b_r, v_b_mirror, acc, acc2);
			stop = 1;
		}
		else {
			(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
			v_b_r = b->il[v_b][y_b];
			printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_b_r, v_a_mirror);
			tripod(b, r, t, v_a_op, v_b_r, v_a_mirror, acc, acc2);
		}
	}
	else if (acc[v_a] != sp && acc[v_b] != sp){ //leg a is empty && leg b is empty
		printf("special subproblem a\n");
		//int next;
		//int mirr;
		//(v_a_op == sp) ? ((v_a_op == v_b_op) ? (next = v_c_op) : (next = v_b_op)) : (next = v_a_op);
		/*if (acc2[v_a_op] != -1){ //if v_a_op is a sperner triangle
			printf("v_a_op == some sp\n");
			if (v_a_op == v_b_op){
				printf("v_a_op == v_b_op\n");
				next = v_c_op;
				mirr = v_c_mirror;
			}
			else {
				printf("v_a_op != v_b_op\n");
				next = v_b_op;
				mirr = v_b_mirror;
			}
		}
		else{
			printf("v_a_op != some sp\n");
			next = v_a_op;
			mirr = v_a_mirror;
		}*/
		/*if (acc2[v_a_op] != -1 || acc2[v_a_mirror] != -1){ //don't use a
			if (acc2[v_b_op] != -1 || acc2[v_b_mirror] != -1){ //don't use b
				if (acc2[v_c_op] != -1 || acc2[v_c_mirror] != -1){ //don't use c
					next = -1;
					mirr = -1;
				}
				else { //use c
					next = v_c_op;
					mirr = v_c_mirror;
				}
			}
			else { //use b
				next = v_b_op;
				mirr = v_b_mirror;
			}
		}
		else { //use a
			next = v_a_op;
			mirr = v_a_mirror;
		}

		printf("next = %d, mirr = %d\n", next, mirr);
		printf("special subproblem a for sp %d will be on faces %d, %d, %d\n", sp, next, mirr, mirr);
		tripod(b, r, t, next, mirr, mirr, acc, acc2);*/
		/*if (acc2[v_a_op] == -1){ //if op is not a sperner triangle (aka if we're not on a border)
			if (acc2[v_a_mirror] == -1 && acc2[v_a_op] == -1){ //make sure we're not recursing on a previous sp
				printf("special subproblem 1 for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_a_mirror, v_a_mirror);
				tripod(b, r, t, v_a_op, v_a_mirror, v_a_mirror, acc, acc2);
			}
		}*/
	}

	if (!stop){ //if there are more than one subproblems, carry on
		if (acc[v_b] == sp && acc[v_c] == sp){ //leg b is non-empty && leg c is non-empty
			b_incomplete = 0;
			v_b_l = b->il[v_b][b->pin[v_b]];
			(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
			v_c_r = b->il[v_c][y_c];
			//no need to check for bichromatic cycle
			printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_c_r, v_b_l);
			tripod(b, r, t, v_b_op, v_c_r, v_b_l, acc, acc2);
		}
		else if (acc[v_b] == sp && acc[v_c] != sp){ //leg b is non-empty && leg c is empty
			b_incomplete = 0;
			//conditional below checks whether cycle defining subproblem is bichromatic. aka. if there is only one subproblem
			if (v_b_next == v_c || v_b_next == v_a){ //if the path up the bfs tree from v_b leads to v_c or v_a
				(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
				v_b_r = b->il[v_b][y_b];
				printf("the only subproblem for sp %d will be on faces %d, %d, %d\n", sp, v_a_mirror, v_b_r, v_b_mirror);
				tripod(b, r, t, v_a_mirror, v_b_r, v_b_mirror, acc, acc2);
			}
			else {
				v_b_l = b->il[v_b][b->pin[v_b]];
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_b_mirror, v_b_l);
				tripod(b, r, t, v_b_op, v_b_mirror, v_b_l, acc, acc2);
			}
		}
		else if (acc[v_b] != sp && acc[v_c] == sp){ //leg b is empty && leg c is non-empty
			b_incomplete = 0;
			//conditional below checks whether cycle defining subproblem is bichromatic. aka. if there is only one subproblem
			if (v_c_next == v_a || v_c_next == v_b){ //if the path up the bfs tree from v_c leads to v_a or v_b
				(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
				v_c_r = b->il[v_c][y_c];
				printf("the only subproblem for sp %d will be on faces %d, %d, %d\n", sp, v_b_mirror, v_c_r, v_c_mirror);
				tripod(b, r, t, v_b_mirror, v_c_r, v_c_mirror, acc, acc2);
			}
			else {
				(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
				v_c_r = b->il[v_c][y_c];
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_c_r, v_b_mirror);
				tripod(b, r, t, v_b_op, v_c_r, v_b_mirror, acc, acc2);
			}
		}
		else if (acc[v_b] != sp && acc[v_c] != sp){ //leg b is empty && leg c is empty
			//printf("special case b\n");
			//int next;
			//(v_b_op == sp) ? ((v_b_op == v_c_op) ? (next = v_a_op) : (next = v_c_op)) : (next = v_b_op);
			//printf("next = %d, v_b_mirror = %d\n", next, v_b_mirror);
			/*if (acc2[v_b_op] == -1){ //if op is not a sperner triangle (aka if we're not on a border)
				if (acc2[v_b_mirror] == -1 && acc2[v_b_op] == -1){ //make sure we're not recursing on a previous sp
					printf("special subproblem 2 for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_b_mirror, v_b_mirror);
					tripod(b, r, t, v_b_op, v_b_mirror, v_b_mirror, acc, acc2);
				}
			}*/
		}

		if (acc[v_c] == sp && acc[v_a] == sp){ //leg c is non-empty && leg a is non-empty
			c_incomplete = 0;
			v_c_l = b->il[v_c][b->pin[v_c]];
			(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
			v_a_r = b->il[v_a][y_a];
			//no need to check for bichromatic cycle
			printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_a_r, v_c_l);
			tripod(b, r, t, v_c_op, v_a_r, v_c_l, acc, acc2);
		}
		else if (acc[v_c] == sp && acc[v_a] != sp){ //leg c is non-empty && leg a is empty
			c_incomplete = 0;
			//conditional below checks whether cycle defining subproblem is bichromatic. aka. if there is only one subproblem
			if (v_c_next == v_a || v_c_next == v_b){ //if the path up the bfs tree from v_c leads to v_a or v_b
				(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
				v_c_r = b->il[v_c][y_c];
				printf("the only subproblem for sp %d will be on faces %d, %d, %d\n", sp, v_b_mirror, v_c_r, v_c_mirror);
				tripod(b, r, t, v_b_mirror, v_c_r, v_c_mirror, acc, acc2);
			}
			else {
				v_c_l = b->il[v_c][b->pin[v_c]];
				printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_c_mirror, v_c_l);
				tripod(b, r, t, v_c_op, v_c_mirror, v_c_l, acc, acc2);
			}
		}
		else if (acc[v_c] != sp && acc[v_a] == sp){ //leg c is empty && leg a is non-empty
			c_incomplete = 0;
			//conditional below checks whether cycle defining subproblem is bichromatic. aka. if there is only one subproblem
			if (v_a_next == v_b || v_a_next == v_c){ //if the path up the bfs tree from v_a leads to v_b or v_c
				(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
				v_a_r = b->il[v_a][y_a];
				printf("the only subproblem for sp %d will be on faces %d, %d, %d\n", sp, v_c_mirror, v_a_r, v_a_mirror);
				tripod(b, r, t, v_c_mirror, v_a_r, v_a_mirror, acc, acc2);
			}
			else {
				(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
				v_a_r = b->il[v_a][y_a];
				printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_a_r, v_c_mirror);
				tripod(b, r, t, v_c_op, v_a_r, v_c_mirror, acc, acc2);
			}
		}
		else if (acc[v_c] != sp && acc[v_a] != sp){ //leg c is empty && leg a is empty
			//printf("special case c\n");
			//int next;
			//(v_c_op == sp) ? ((v_c_op == v_a_op) ? (next = v_b_op) : (next = v_a_op)) : (next = v_c_op);
			//printf("next = %d, v_c_mirror = %d\n", next, v_c_mirror);
			/*if (acc2[v_c_mirror] == -1 && acc2[v_c_op] == -1){ //make sure we're not recursing on a previous sp
				if (acc2[v_c_op] == -1){ //if op is not a sperner triangle (aka if we're not on a border)
					printf("special subproblem 3 for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_c_mirror, v_c_mirror);
					tripod(b, r, t, v_c_op, v_c_mirror, v_c_mirror, acc, acc2);
				}
			}*/
		}

		/*if (a_incomplete && b_incomplete && c_incomplete){
			printf("in special if\n");
			int next;

			(v_a_op == v_b_op) ? (next = v_c_op) : (next = v_b_op);
			printf("acc2[v_a_mirror] = acc2[%d] = %d\n", v_a_mirror, acc2[v_a_mirror]);
			printf("acc2[v_a_op] = acc2[%d] = %d\n", v_a_op, acc2[v_a_op]);
			if (acc2[v_a_mirror] == -1 && acc2[next] == -1){ //make sure we're not recursing on a previous sp
				printf("v_a_op = %d, v_b_op = %d, v_c_op = %d\n", v_a_op, v_b_op, v_c_op);
				printf("special subproblem for sp %d will be on faces %d, %d, %d\n", sp, next, v_a_mirror, v_a_mirror);
				tripod(b, r, t, next, v_a_mirror, v_a_mirror, acc, acc2);
			}

			if (acc2[v_a_mirror] == -1 && acc2[v_a_op] == -1){ //make sure we're not recursing on a previous sp
				printf("special subproblem 1 for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_a_mirror, v_a_mirror);
				tripod(b, r, t, v_a_op, v_a_mirror, v_a_mirror, acc, acc2);
			}
			if (acc2[v_b_mirror] == -1 && acc2[v_b_op] == -1){ //make sure we're not recursing on a previous sp
				printf("special subproblem 2 for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_b_mirror, v_b_mirror);
				tripod(b, r, t, v_b_op, v_b_mirror, v_b_mirror, acc, acc2);
			}
			if (acc2[v_c_mirror] == -1 && acc2[v_c_op] == -1){ //make sure we're not recursing on a previous sp
				printf("special subproblem 3 for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_c_mirror, v_c_mirror);
				tripod(b, r, t, v_c_op, v_c_mirror, v_c_mirror, acc, acc2);
			}
		}*/
	}
}

void tprint(struct tripod_decomposition_struct* t){
	printf("\n");
	printf("v_a = %d\n", t->v_a);
	printf("v_b = %d\n", t->v_b);
	printf("v_c = %d\n", t->v_c);
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
