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
	//for each vertex v, acc[v] represents the colour of the tripod v belongs to, -1 if it does not yet belong to a tripod
	//for each face f, acc2[f] is -1 if it does not belong to a tripod, and acc2[f]=f otherwise
	for (int k = 0; k < (b->v); k++){
		acc[k] = -1; //vertices labeled with -1 do not yet belong to a tripod
	}
	acc2[0] = 0;
	for (int k = 1; k < (b->f); k++){
		acc2[k] = -1; //faces labeled with -1 have not been identified as sperner triangles
	}

	//label the three edges incident to face 0 as belonging to three exterior tripods
	acc[b->sim[0][0]] = (b->f);
	acc[b->sim[0][1]] = (b->f)+1;
	acc[b->sim[0][2]] = (b->f)+2;

	//start decomposition with the three triangles adjacent to outer face
	//trichromatic_tripod(b, r, t, b->tri[0][0], b->tri[0][2], b->tri[0][1], acc, acc2);
	trichromatic_tripod(b, r, t, b->tri[0][0], b->tri[0][1], b->tri[0][2], acc, acc2); //new from aug.15
}

/******************************************************** trichromatic ********************************************************/

int* trichromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int f2, int f3, int* acc, int* acc2){
	//sperner triangle identification
	int sp;
	if (f1 == f2 && f1 == f3){
		acc2[f1] = f1; //keep track of sperner triangles
		return 0;
	}
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

	//store tripod
	store_tripod(b, t, sp, acc);

	printf("\n%d is one of our sperner triangles\nacc = ", sp); //sperner triangles are correctly identified
	for (int k = 0; k < (b->v); k++){
		printf("%d ", acc[k]);
	}
	printf("\nacc2 = ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", acc2[i]);
	}
	printf("\n");

	acc2[sp] = sp; //keep track of sperner triangles

	if (f1 == sp && f2 == sp && f3 == sp) return 0; //no subproblems exist

	trichromatic_orient_subproblems(b, t, sp, f1, f2, f3, acc);

	tprint(t);

	trichromatic_decompose(b, r, t, sp, f1, f2, f3, acc, acc2);

	return 0;
}

void trichromatic_orient_subproblems(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3, int* acc){
	int leg_colour_a;
	int leg_colour_b;
	int leg_colour_c;
	
	(t->v_a_next == -1) ? (leg_colour_a = acc[t->v_a]) : (leg_colour_a = acc[t->v_a_next]);
	(t->v_b_next == -1) ? (leg_colour_b = acc[t->v_b]) : (leg_colour_b = acc[t->v_b_next]);
	(t->v_c_next == -1) ? (leg_colour_c = acc[t->v_c]) : (leg_colour_c = acc[t->v_c_next]);

	trichromatic_orient_subproblems_pt2(b, t, sp, f1, f2, f3, leg_colour_a, leg_colour_b, leg_colour_c, acc);
}

void trichromatic_orient_subproblems_pt2(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3, int leg_colour_a, int leg_colour_b, int leg_colour_c, int* acc){
	int f;
	int next_f;
	int next_next_f;

	if (f1 == b->tri[0][0] && f2 == b->tri[0][1] && f3 == b->tri[0][2]){
		//if we're using the faces adjacent to the outer face, then the orientation is different from the rest
		if (f1 != sp){
			f = f1;
			printf("f = f1 = %d\n", f);
			next_f = f3;
			next_next_f = f2;
		}
		else if (f2 != sp){
			f = f2;
			printf("f = f2 = %d\n", f);
			next_f = f1;
			next_next_f = f3;
		}
		else {
			f = f3;
			printf("f = f3 = %d\n", f);
			next_f = f2;
			next_next_f = f1;
		}
	}
	else {
		if (f1 != sp){
			f = f1;
			printf("f = f1 = %d\n", f);
			next_f = f2;
			next_next_f = f3;
		}
		else if (f2 != sp){
			f = f2;
			printf("f = f2 = %d\n", f);
			next_f = f3;
			next_next_f = f1;
		}
		else {
			f = f3;
			printf("f = f3 = %d\n", f);
			next_f = f1;
			next_next_f = f2;
		}
	}

	trichromatic_orient_subproblems_pt3(b, t, sp, f1, f2, f3, leg_colour_a, leg_colour_b, leg_colour_c, f, next_f, next_next_f, acc);
}

void trichromatic_orient_subproblems_pt3(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3, int leg_colour_a, int leg_colour_b, int leg_colour_c, int f, int next_f, int next_next_f, int* acc){
	int double_match = 0;
	for (int i = 0; i < 3; i++){ //look at all three vertices of our f triangle
		if (acc[b->sim[f][i]] == leg_colour_a){
			for (int j = 0; j < 3; j++){
				if (acc[b->sim[f][j]] == leg_colour_b){
					double_match = 1;
				}
			}
		}
	}
	if (double_match){
		printf("ORIENTATION 1\n");

		t->v_a_op = f;
		t->v_b_op = next_f;
		t->v_c_op = next_next_f;
	}
	else { //no double match. reset values and try with leg_colour_b and leg_colour_c
		double_match = 0;
		for (int i = 0; i < 3; i++){ //look at all three vertices of our f triangle
			if (acc[b->sim[f][i]] == leg_colour_b){
				for (int j = 0; j < 3; j++){
					if (acc[b->sim[f][j]] == leg_colour_c){
						double_match = 1;
					}
				}
			}
		}
		if (double_match){
			printf("ORIENTATION 2\n");

			t->v_a_op = next_next_f;
			t->v_b_op = f;
			t->v_c_op = next_f;
		}
		else { //no double match. f must match with leg_colour_c and leg_colour_a.
			printf("ORIENTATION 3\n");

			t->v_a_op = next_f;
			t->v_b_op = next_next_f;
			t->v_c_op = f;
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

	//subproblem a
	if (acc[v_a] == sp && acc[v_b] == sp){ //leg a is non-empty && leg b is non-empty
		v_a_l = b->il[v_a][b->pin[v_a]];
		(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
		v_b_r = b->il[v_b][y_b];
		printf("subproblem a for sp %d is trichromatic\n", sp);
		printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_b_r, v_a_l);
		trichromatic_tripod(b, r, t, v_a_op, v_b_r, v_a_l, acc, acc2);
	}
	else if (acc[v_a] == sp && acc[v_b] != sp){ //leg a is non-empty && leg b is empty
		v_a_l = b->il[v_a][b->pin[v_a]];
		printf("subproblem a for sp %d is trichromatic\n", sp);
		printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_a_mirror, v_a_l);
		trichromatic_tripod(b, r, t, v_a_op, v_a_mirror, v_a_l, acc, acc2);
	}
	else if (acc[v_a] != sp && acc[v_b] == sp){ //leg a is empty && leg b is non-empty
		(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
		v_b_r = b->il[v_b][y_b];
		printf("subproblem a for sp %d is trichromatic\n", sp);
		printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_b_r, v_a_mirror);
		trichromatic_tripod(b, r, t, v_a_op, v_b_r, v_a_mirror, acc, acc2);
	}
	else if (acc[v_a] != sp && acc[v_b] != sp){ //leg a is empty && leg b is empty
		if (sp == v_a_op){
			printf("no subproblem a for sp %d\n", sp);
			v_a_mirror = -1; //set v_a_mirror to negative value to later check whether it is in our subproblem boundary
		}
		else {
			printf("subproblem a for sp %d is bichromatic\n", sp);
			printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_a_op, v_a_mirror);
			bichromatic_tripod( b, r, t, v_a_op, v_a_mirror, acc, acc2);
			printf("---- weird subproblem a for which we are testing ----\n");
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
	}
	else if (acc[v_b] == sp && acc[v_c] != sp){ //leg b is non-empty && leg c is empty
		v_b_l = b->il[v_b][b->pin[v_b]];
		printf("subproblem b for sp %d is trichromatic\n", sp);
		printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_b_mirror, v_b_l);
		trichromatic_tripod(b, r, t, v_b_op, v_b_mirror, v_b_l, acc, acc2);
	}
	else if (acc[v_b] != sp && acc[v_c] == sp){ //leg b is empty && leg c is non-empty
		(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
		v_c_r = b->il[v_c][y_c];
		printf("subproblem b for sp %d is trichromatic\n", sp);
		printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_c_r, v_b_mirror);
		trichromatic_tripod(b, r, t, v_b_op, v_c_r, v_b_mirror, acc, acc2);
	}
	else if (acc[v_b] != sp && acc[v_c] != sp){ //leg b is empty && leg c is empty
		if (sp == v_b_op){
			printf("no subproblem b for sp %d\n", sp);
			v_b_mirror = -1; //set v_b_mirror to negative value to later check whether it is in our subproblem boundary
		}
		else {
			printf("subproblem b for sp %d is bichromatic\n", sp);
			printf("subproblem b for sp %d will be on faces %d, %d\n", sp, v_b_op, v_b_mirror);
			bichromatic_tripod( b, r, t, v_b_op, v_b_mirror, acc, acc2);
			printf("---- weird subproblem b for which we are testing ----\n");
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
			v_c_mirror = -1; //set v_c_mirror to negative value to later check whether it is in our subproblem boundary
		}
		else {
			printf("subproblem c for sp %d is bichromatic\n", sp);
			printf("subproblem c for sp %d will be on faces %d, %d\n", sp, v_c_op, v_c_mirror);
			bichromatic_tripod( b, r, t, v_c_op, v_c_mirror, acc, acc2);
			printf("---- weird subproblem c for which we are testing ----\n");
		}
	}
}

/******************************************************** bichromatic ********************************************************/

int* bichromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int f2, int* acc, int* acc2){
	int sp;

	sp = f1;

	store_tripod(b, t, sp, acc);

	printf("\n%d is one of our sperner triangles\nacc = ", sp); //sperner triangles are correctly identified
	for (int k = 0; k < (b->v); k++){
		printf("%d ", acc[k]);
	}
	printf("\nacc2 = ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", acc2[i]);
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
	printf("v_a_op = %d\n", t->v_a_op);
	printf("v_b_op = %d\n", t->v_b_op);
	printf("v_c_op = %d\n", t->v_c_op);

	bichromatic_decompose(b, r, t, sp, f1, f2, acc, acc2);

	return 0;
}

void bichromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int* acc, int* acc2){
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

	if (acc[v_a] == sp || acc[v_b] == sp || acc[v_c] == sp){ //if one of {leg a, leg b, leg c} is non-empty
		//here we either have case 5.1 or 5.2
		//there can be at most one non-empty leg, since more empty legs lead to trichromatic subproblems
		//the vertices of the crotch of sp are one of exactly three colours
		if (acc[v_a] == sp){ //v_a is non-empty
			if (v_a_next != v_b && v_a_next != v_c){ //if the path up the bfs tree from v_a does NOT lead to v_b or v_c
				//here we have case 5.1
				printf("we have two subproblems for sp %d: trichromatic + bichromatic\n", sp);
				if (acc[v_a_next] == acc[v_c]){ //if v_a_next is the same colour as v_c
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_a_r, v_c_mirror);
					bichromatic_tripod( b, r, t, v_a_r, v_c_mirror, acc, acc2);
					printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_a_l, v_a_mirror);
					trichromatic_tripod(b, r, t, f2, v_a_l, v_a_mirror, acc, acc2);
				}
				else { //if v_a_next is the same colour as v_b
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_a_l, v_a_mirror);
					bichromatic_tripod( b, r, t, v_a_l, v_a_mirror, acc, acc2);
					printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_a_r, v_c_mirror);
					trichromatic_tripod(b, r, t, f2, v_a_r, v_c_mirror, acc, acc2);
				}
			}
			else {
				//here we have case 5.2
				printf("we have one subproblem for sp %d: trichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_c_mirror, v_a_mirror);
				trichromatic_tripod(b, r, t, f2, v_c_mirror, v_a_mirror, acc, acc2);
			}
		}
		else if (acc[v_b] == sp){ //v_b is non-empty
			if (v_b_next != v_a && v_b_next != v_c){ //if the path up the bfs tree from v_b does NOT lead to v_a or v_c
				//here we have case 5.1
				printf("we have two subproblems for sp %d: trichromatic + bichromatic\n", sp);
				if (acc[v_b_next] == acc[v_a]){ //if v_b_next is the same colour as v_a
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_b_r, v_a_mirror);
					bichromatic_tripod( b, r, t, v_b_r, v_a_mirror, acc, acc2);
					printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_b_l, v_b_mirror);
					trichromatic_tripod(b, r, t, f2, v_b_l, v_b_mirror, acc, acc2);
				}
				else { //if v_b_next is the same colour as v_c
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_b_l, v_b_mirror);
					bichromatic_tripod( b, r, t, v_b_l, v_b_mirror, acc, acc2);
					printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_b_r, v_a_mirror);
					trichromatic_tripod(b, r, t, f2, v_b_r, v_a_mirror, acc, acc2);
				}
			}
			else {
				//here we have case 5.2
				printf("we have one subproblem for sp %d: trichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_a_mirror, v_b_mirror);
				trichromatic_tripod(b, r, t, f2, v_a_mirror, v_b_mirror, acc, acc2);
			}
		}
		else { //v_c is non-empty
			if (v_c_next != v_a && v_c_next != v_b){ //if the path up the bfs tree from v_c does NOT lead to v_a or v_b
				//here we have case 5.1
				printf("we have two subproblems for sp %d: trichromatic + bichromatic\n", sp);
				if (acc[v_c_next] == acc[v_b]){ //if v_c_next is the same colour as v_b
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_c_r, v_b_mirror);
					bichromatic_tripod( b, r, t, v_c_r, v_b_mirror, acc, acc2);
					printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_c_l, v_c_mirror);
					trichromatic_tripod(b, r, t, f2, v_c_l, v_c_mirror, acc, acc2);
				}
				else { //if v_c_next is the same colour as v_a
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_c_l, v_c_mirror);
					bichromatic_tripod( b, r, t, v_c_l, v_c_mirror, acc, acc2);
					printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_c_r, v_b_mirror);
					trichromatic_tripod(b, r, t, f2, v_c_r, v_b_mirror, acc, acc2);
				}
			}
			else {
				//here we have case 5.2
				printf("we have one subproblem for sp %d: trichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_b_mirror, v_c_mirror);
				trichromatic_tripod(b, r, t, f2, v_b_mirror, v_c_mirror, acc, acc2);
			}
		}
	}

	else { //all legs are empty
		//here we either have case 5.3, 5.4, 5.5, or 5.6 (the empty case)
		//the vertices of the crotch of sp are one of exactly two colours

		if (acc[v_a] == acc[v_b]){ //if v_a and v_b are the same colour, with v_c a different colour
			//edge of interest is {v_a, v_b}
			//mirror triangle of interest is v_a_mirror
			
			//v_a_mirror will have two of its vertices the same colour as v_a and v_b
			//v_a_mirror's third vertex is either the same colour as v_a and v_b, or is uncoloured
			
			if (b->bt[v_b] == v_a || b->bt[v_a] == v_b){
				//if v_a is a bfs parent of v_b or v_b is a bfs parent of v_a, then we are in case 5.4
				//DIFFERENT CASES HERE, DEPENDING ON WHETHER V_A IS PARENT OF V_B OR V_B IS PARENT OF V_A !!!!
				printf("we have one subproblem for sp %d: bichromatic\n", sp);
				if (b->bt[v_b] == v_a){ //if v_a is a bfs parent of v_b
					if (acc2[v_b_mirror] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_b_mirror);
						bichromatic_tripod( b, r, t, f2, v_b_mirror, acc, acc2);
					}
				}
				else { //v_b is a bfs parent of v_a
					if (acc2[v_c_mirror] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_c_mirror);
						bichromatic_tripod( b, r, t, f2, v_c_mirror, acc, acc2);
					}
				}
			}
			//otherwise, we are in case 5.3 or 5.5
			else if (f1 == f2){
				if (acc2[v_a_mirror] != acc[v_a]){ //make sure subproblem exists with (acc2[v_a_mirror] != acc[v_a]) check
					//if f1 == f2 then we are in case 5.5
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_a_mirror);
					monochromatic_tripod( b, r, t, v_a_mirror, acc, acc2);
				}
			}
			else {
				//we are in case 5.3
				//for any coloured edge, either it's a tree edge or there's a sp on the other side
				if (acc2[v_a_mirror] != acc[v_a]){ //make sure monochromatic subproblem exists with (acc2[v_a_mirror] != acc[v_a]) check
					printf("we have two subproblems for sp %d: monochromatic + bichromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_a_mirror);
					monochromatic_tripod( b, r, t, v_a_mirror, acc, acc2);
				}
				else {
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
				}
				//DIFFERENT CASES HERE, DEPENDING ON SPECIFIC TRIANGLE ORIENTATION !!!!
				if ((b->bt[v_b] == v_c || b->bt[v_c] == v_b) || acc2[v_b_mirror] > -1){
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_c_mirror);
					bichromatic_tripod( b, r, t, f2, v_c_mirror, acc, acc2);
				}
				else {
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_b_mirror);
					bichromatic_tripod( b, r, t, f2, v_b_mirror, acc, acc2);
				}
			}
		}
		else if (acc[v_b] == acc[v_c]){ //if v_b and v_c are the same colour, with v_a a different colour
			//edge of interest is {v_b, v_c}
			//mirror triangle of interest is v_b_mirror
			
			//v_b_mirror will have two of its vertices the same colour as v_b and v_c
			//v_b_mirror's third vertex is either the same colour as v_b and v_c, or is uncoloured
			
			if (b->bt[v_c] == v_b || b->bt[v_b] == v_c){
				//if v_b is a bfs parent of v_c or v_c is a bfs parent of v_b, then we are in case 5.4
				//DIFFERENT CASES HERE, DEPENDING ON WHETHER V_B IS PARENT OF V_C OR V_C IS PARENT OF V_B !!!!
				printf("we have one subproblem for sp %d: bichromatic\n", sp);
				if (b->bt[v_c] == v_b){ //v_b is a bfs parent of v_c
					if (acc2[v_c_mirror] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_c_mirror);
						bichromatic_tripod( b, r, t, f2, v_c_mirror, acc, acc2);
					}
				}
				else { //v_c is a bfs parent of v_b
					if (acc2[v_a_mirror] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_a_mirror);
						bichromatic_tripod( b, r, t, f2, v_a_mirror, acc, acc2);
					}
				}
			}
			//otherwise, we are in case 5.3 or 5.5
			else if (f1 == f2){
				if (acc2[v_b_mirror] != acc[v_b]){ //make sure subproblem exists with (acc2[v_b_mirror] != acc[v_b]) check
					//if f1 == f2 then we are in case 5.5
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_b_mirror);
					monochromatic_tripod( b, r, t, v_b_mirror, acc, acc2);
				}
			}
			else {
				//we are in case 5.3
				//for any coloured edge, either it's a tree edge or there's a sp on the other side
				if (acc2[v_b_mirror] != acc[v_b]){ //make sure monochromatic subproblem exists with (acc2[v_b_mirror] != acc[v_b]) check
					printf("we have two subproblems for sp %d: bichromatic + monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_b_mirror);
					monochromatic_tripod( b, r, t, v_b_mirror, acc, acc2);
				}
				else {
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
				}
				//DIFFERENT CASES HERE, DEPENDING ON SPECIFIC TRIANGLE ORIENTATION !!!!
				if ((b->bt[v_c] == v_a || b->bt[v_a] == v_c) || acc2[v_c_mirror] > -1){
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_a_mirror);
					bichromatic_tripod( b, r, t, f2, v_a_mirror, acc, acc2);
				}
				else {
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_c_mirror);
					bichromatic_tripod( b, r, t, f2, v_c_mirror, acc, acc2);
				}
			}
		}
		else if (acc[v_c] == acc[v_a]){ //if v_c and v_a are the same colour, with v_b a different colour
			//edge of interest is {v_c, v_a}
			//mirror triangle of interest is v_c_mirror
			
			//v_c_mirror will have two of its vertices the same colour as v_c and v_a
			//v_c_mirror's third vertex is either the same colour as v_c and v_a, or is uncoloured
			
			if (b->bt[v_a] == v_c || b->bt[v_c] == v_a){
				//if v_c is a bfs parent of v_a or v_a is a bfs parent of v_c, then we are in case 5.4
				//DIFFERENT CASES HERE, DEPENDING ON WHETHER V_C IS PARENT OF V_A OR V_A IS PARENT OF V_C !!!!
				printf("we have one subproblem for sp %d: bichromatic\n", sp);
				if (b->bt[v_a] == v_c){ //if v_c is a bfs parent of v_a
					if (acc2[v_a_mirror] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_a_mirror);
						bichromatic_tripod( b, r, t, f2, v_a_mirror, acc, acc2);
					}
				}
				else { //v_a is a bfs parent of v_c
					if (acc2[v_b_mirror] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_b_mirror);
						bichromatic_tripod( b, r, t, f2, v_b_mirror, acc, acc2);
					}
				}
			}
			//otherwise, we are in case 5.3 or 5.5
			else if (f1 == f2){
				if (acc2[v_c_mirror] != acc[v_c]){ //make sure subproblem exists with (acc2[v_c_mirror] != acc[v_c]) check
					//if f1 == f2 then we are in case 5.5
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_c_mirror);
					monochromatic_tripod( b, r, t, v_c_mirror, acc, acc2);
				}
			}
			else {
				//we are in case 5.3
				//for any coloured edge, either it's a tree edge or there's a sp on the other side
				if (acc2[v_c_mirror] != acc[v_c]){ //make sure monochromatic subproblem exists with (acc2[v_c_mirror] != acc[v_c]) check
					printf("we have two subproblems for sp %d: bichromatic + monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_c_mirror);
					monochromatic_tripod( b, r, t, v_c_mirror, acc, acc2);
				}
				else {
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
				}
				//DIFFERENT CASES HERE, DEPENDING ON SPECIFIC TRIANGLE ORIENTATION !!!!
				if ((b->bt[v_a] == v_b || b->bt[v_b] == v_a) || acc2[v_a_mirror] > -1){
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_b_mirror);
					bichromatic_tripod( b, r, t, f2, v_b_mirror, acc, acc2);
				}
				else {
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_a_mirror);
					bichromatic_tripod( b, r, t, f2, v_a_mirror, acc, acc2);
				}
			}
		}
		else {
			printf("special case involving the boundary of our planar graph\n"); //not necessarily the boundary...
			printf("all vertices are different colours\n");
			printf("are we done?\n"); //looks like it
		}
	}
}

/******************************************************** monochromatic ********************************************************/
int* monochromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int* acc, int* acc2){
	int sp = f1;

	store_tripod(b, t, sp, acc);

	printf("\n%d is one of our sperner triangles\nacc = ", sp); //sperner triangles are correctly identified
	for (int k = 0; k < (b->v); k++){
		printf("%d ", acc[k]);
	}
	printf("\nacc2 = ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", acc2[i]);
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

	monochromatic_decompose(b, r, t, sp, f1, acc, acc2); //passing sp and f1 is redundant

	return 0;
}

void monochromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1, int* acc, int* acc2){
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
	if (acc[v_a] == sp || acc[v_b] == sp || acc[v_c] == sp){ //if one of {leg a, leg b, leg c} is non-empty
		//here we either have case 6.1 or 6.2
		//there can be at most one non-empty leg
		//the vertices of the crotch of sp are one of exactly three colours
		if (acc[v_a] == sp){ //v_a is non-empty
			if (v_a_next != v_b && v_a_next != v_c){ //if the path up the bfs tree from v_a does NOT lead to v_b or v_c
				//here we have case 6.1
				printf("we have two subproblems for sp %d: bichromatic + bichromatic\n", sp);
				printf("first bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_a_r, v_c_mirror);
				bichromatic_tripod( b, r, t, v_a_r, v_c_mirror, acc, acc2);
				printf("second bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_a_l, v_a_mirror);
				bichromatic_tripod( b, r, t, v_a_l, v_a_mirror, acc, acc2);

			}
			else {
				//here we have case 6.2
				printf("we have one subproblem for sp %d: bichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d\n", sp, v_c_mirror, v_a_mirror);
				bichromatic_tripod(b, r, t, v_c_mirror, v_a_mirror, acc, acc2);
			}
		}
		else if (acc[v_b] == sp){ //v_b is non-empty
			if (v_b_next != v_a && v_b_next != v_c){ //if the path up the bfs tree from v_b does NOT lead to v_a or v_c
				//here we have case 6.1
				printf("we have two subproblems for sp %d: bichromatic + bichromatic\n", sp);
				printf("first bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_b_r, v_a_mirror);
				bichromatic_tripod( b, r, t, v_b_r, v_a_mirror, acc, acc2);
				printf("second bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_b_l, v_b_mirror);
				bichromatic_tripod( b, r, t, v_b_l, v_b_mirror, acc, acc2);
			}
			else {
				//here we have case 6.2
				printf("we have one subproblem for sp %d: bichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d\n", sp, v_a_mirror, v_b_mirror);
				bichromatic_tripod(b, r, t, v_a_mirror, v_b_mirror, acc, acc2);
			}
		}
		else { //v_c is non-empty
			if (v_c_next != v_a && v_c_next != v_b){ //if the path up the bfs tree from v_c does NOT lead to v_a or v_b
				//here we have case 6.1
				printf("we have two subproblems for sp %d: bichromatic + bichromatic\n", sp);
				printf("first bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_c_r, v_b_mirror);
				bichromatic_tripod( b, r, t, v_c_r, v_b_mirror, acc, acc2);
				printf("second bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_c_l, v_c_mirror);
				bichromatic_tripod( b, r, t, v_c_l, v_c_mirror, acc, acc2);
			}
			else {
				//here we have case 6.2
				printf("we have one subproblem for sp %d: bichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d\n", sp, v_b_mirror, v_c_mirror);
				bichromatic_tripod(b, r, t, v_b_mirror, v_c_mirror, acc, acc2);
			}
		}
	}
	else { //all legs are empty
		//here we either have case 6.3, 6.4, or 6.5 (the empty case)
		//the vertices of sp are one colour
		if (b->bt[v_b] == v_a || b->bt[v_a] == v_b){ //if v_a is a bfs parent of v_b or v_b is a bfs parent of v_a
			if (b->bt[v_a] == v_b){ //if v_b is a bfs parent of v_a
				//here, if the subproblem is empty, a given edge with purple endpoints is either an edge of the purple tripod, or it is adjacent to a sperner triangle
				if (acc2[v_c_mirror] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_c_mirror);
					monochromatic_tripod( b, r, t, v_c_mirror, acc, acc2);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
			else { //v_a is a bfs parent of v_b
				if (acc2[v_b_mirror] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_b_mirror);
					monochromatic_tripod( b, r, t, v_b_mirror, acc, acc2);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
		}
		else if (b->bt[v_c] == v_b || b->bt[v_b] == v_c){ //if v_b is a bfs parent of v_c or v_c is a bfs parent of v_b
			if (b->bt[v_b] == v_c){ //v_c is a bfs parent of v_b
				if (acc2[v_a_mirror] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_a_mirror);
					monochromatic_tripod( b, r, t, v_a_mirror, acc, acc2);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
			else { //v_b is a bfs parent of v_c
				if (acc2[v_c_mirror] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_c_mirror);
					monochromatic_tripod( b, r, t, v_c_mirror, acc, acc2);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
		}
		else if (b->bt[v_a] == v_c || b->bt[v_c] == v_a){ //if v_c is a bfs parent of v_a or v_a is a bfs parent of v_c
			if (b->bt[v_c] == v_a){ //if v_a is a bfs parent of v_c
				if (acc2[v_b_mirror] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_b_mirror);
					monochromatic_tripod( b, r, t, v_b_mirror, acc, acc2);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
			else { //v_c is a bfs parent of v_a
				if (acc2[v_a_mirror] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_a_mirror);
					monochromatic_tripod( b, r, t, v_a_mirror, acc, acc2);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
		}
		else { //if none of v_a, v_b, v_c are parents to one another
			//here we are in case 6.3
			printf("we have two subproblems for sp %d: monochromatic + monochromatic\n", sp);
			if (acc2[v_a_mirror] > -1){ //if v_a_mirror is a previous sp
				//recurse on v_b_mirror and v_c_mirror
				printf("first monochromatic subproblem for sp %d will be on face %d\n", sp, v_b_mirror);
				monochromatic_tripod( b, r, t, v_b_mirror, acc, acc2);
				printf("second monochromatic subproblem for sp %d will be on face %d\n", sp, v_c_mirror);
				monochromatic_tripod( b, r, t, v_c_mirror, acc, acc2);
			}
			else if (acc2[v_b_mirror] > -1){ //if v_b_mirror is a previous sp
				//recurse on v_c_mirror and v_a_mirror
				printf("first monochromatic subproblem for sp %d will be on face %d\n", sp, v_c_mirror);
				monochromatic_tripod( b, r, t, v_c_mirror, acc, acc2);
				printf("second monochromatic subproblem for sp %d will be on face %d\n", sp, v_a_mirror);
				monochromatic_tripod( b, r, t, v_a_mirror, acc, acc2);
			}
			else { //v_c_mirror is a previous sp
				//recurse on v_a_mirror and v_b_mirror
				printf("first monochromatic subproblem for sp %d will be on face %d\n", sp, v_a_mirror);
				monochromatic_tripod( b, r, t, v_a_mirror, acc, acc2);
				printf("second monochromatic subproblem for sp %d will be on face %d\n", sp, v_b_mirror);
				monochromatic_tripod( b, r, t, v_b_mirror, acc, acc2);
			}
		}
	}
}

void store_tripod(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int* acc){
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
