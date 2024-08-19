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
	for (int k = 0; k < (b->f); k++){
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

	//store tripod
	store_tripod(b, t, sp, acc);

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

	if (f1 == sp && f2 == sp && f3 == sp){
		printf("no subproblems to orient\n");
		t->v_a_op = sp; //setting these to real values so I don't have an error later on in the decomposition
		t->v_b_op = sp;
		t->v_c_op = sp;
		return;
	}
	else {
		if (f1 == b->tri[0][0] && f2 == b->tri[0][1] && f3 == b->tri[0][2]){
			//if we're using the faces adjacent to the outer face, then the orientation is different from the rest
			printf("special case!\n");
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

	int a_flag = 0;
	int b_flag = 0;
	int c_flag = 0;
	int a_np = 0;
	int b_np = 0;
	int c_np = 0;
	int a_near_np = 0;
	int b_near_np = 0;
	int c_near_np = 0;

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
			v_a_mirror = -1; //set v_a_mirror to negative value to later check whether it is in our subproblem boundary
		}
		else {
			printf("subproblem a for sp %d is bichromatic\n", sp);
			printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_a_op, v_a_mirror);
			bichromatic_tripod( b, r, t, v_a_op, v_a_mirror, acc, acc2);
			a_flag = 1;
			a_near_np = 1;
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
			v_b_mirror = -1; //set v_b_mirror to negative value to later check whether it is in our subproblem boundary
		}
		else {
			printf("subproblem b for sp %d is bichromatic\n", sp);
			printf("subproblem b for sp %d will be on faces %d, %d\n", sp, v_b_op, v_b_mirror);
			bichromatic_tripod( b, r, t, v_b_op, v_b_mirror, acc, acc2);
			b_flag = 1;
			b_near_np = 1;
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
		c_flag = 1;
	}
	else if (acc[v_c] == sp && acc[v_a] != sp){ //leg c is non-empty && leg a is empty
		v_c_l = b->il[v_c][b->pin[v_c]];
		printf("subproblem c for sp %d is trichromatic\n", sp);
		printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_c_mirror, v_c_l);
		trichromatic_tripod(b, r, t, v_c_op, v_c_mirror, v_c_l, acc, acc2);
		c_flag = 1;
	}
	else if (acc[v_c] != sp && acc[v_a] == sp){ //leg c is empty && leg a is non-empty
		(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
		v_a_r = b->il[v_a][y_a];
		printf("subproblem c for sp %d is trichromatic\n", sp);
		printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_a_r, v_c_mirror);
		trichromatic_tripod(b, r, t, v_c_op, v_a_r, v_c_mirror, acc, acc2);
		c_flag = 1;
	}
	else if (acc[v_c] != sp && acc[v_a] != sp){ //leg c is empty && leg a is empty
		if (sp == v_c_op){
			printf("no subproblem c for sp %d\n", sp);
			c_np = 1;
			v_c_mirror = -1; //set v_c_mirror to negative value to later check whether it is in our subproblem boundary
		}
		else {
			printf("subproblem c for sp %d is bichromatic\n", sp);
			printf("subproblem c for sp %d will be on faces %d, %d\n", sp, v_c_op, v_c_mirror);
			bichromatic_tripod( b, r, t, v_c_op, v_c_mirror, acc, acc2);
			c_flag = 1;
			printf("---- weird subproblem c for which we are testing ----\n");
		}
	}
}

/******************************************************** bichromatic ********************************************************/

int* bichromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int f2, int* acc, int* acc2){
	int sp;

	if (f1 == f2) return 0;
	sp = f1;

	store_tripod(b, t, sp, acc);

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
	printf("v_a_op = %d\n", t->v_a_op);
	printf("v_b_op = %d\n", t->v_b_op);
	printf("v_c_op = %d\n", t->v_c_op);

	//no need to orient <-- I wouldn't be so sure...

	//bichromatic_decompose(b, r, t, sp, f1, f2, acc, acc2);
	printf("skipping bichromatic_decompose\n");

	return 0;
}

void bichromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int* acc, int* acc2){
	if (f2 == sp) exit(0); //something wrong. abort

	if (sp == 35){
		printf("-------------------------\nsp == 35. temporary exit.\n-------------------------\n");
		exit(0);
	}

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
					printf("over here 1\n");
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
					printf("over here 2\n");
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
					printf("over here 3\n");
					printf("subproblem a for sp %d is bichromatic\n", sp);
					printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_b_l, v_b_mirror);
					bichromatic_tripod( b, r, t, v_b_l, v_b_mirror, acc, acc2);
				}
				else return;

				printf("subproblem b for sp %d is trichromatic\n", sp);
				printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, f2, v_b_r, v_a_mirror);
				trichromatic_tripod( b, r, t, f2, v_b_r, v_a_mirror, acc, acc2);
				
			}
			else {
				if (v_a_mirror != f2 && v_b_r != f2){ //new check
					printf("over here 4\n");
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
		v_c_l = b->il[v_c][b->pin[v_c]];
		(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
		v_c_r = b->il[v_c][y_c];
		if (acc[v_c] == acc[v_a]){ //if v_c is the same colour as v_a
			if (v_a_mirror != f2 && v_b_mirror != f2 && v_c_mirror != f2){ //if v_X_mirror == f2, then stop
				printf("here1\n");
				if (acc[v_c] == sp){ //leg c is not empty //maybe add something like this in more places
					printf("over here 5\n");
					printf("subproblem a for sp %d is monochromatic\n", sp);
					printf("subproblem a for sp %d will be on face %d\n", sp, v_c_mirror);
					//monochromatic_tripod( b, r, t, v_c_mirror, acc, acc2); //maybe specify vertices which make up the exterior
						//for example, v_c & v_a
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
					printf("over here 6\n");
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
				if (acc[v_c] == sp){ //do we need an else for this if?
					if (acc[v_c_next] == acc[v_a]){
						if (v_c_l != f2 && v_c_mirror != f2){ //new check
							printf("over here 7\n");
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
					if (v_c_next == -1){
						printf("over here 8\n");
						printf("subproblem a for sp %d is bichromatic\n", sp);
						printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_c_mirror, f2);
						bichromatic_tripod( b, r, t, v_c_mirror, f2, acc, acc2);
					}
					else {
						if (v_b_mirror != f2 && v_c_r != f2){ //new check
							printf("over here 9\n");
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
int* monochromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int* acc, int* acc2){
	int sp = f1;

	store_tripod(b, t, sp, acc);

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

	//no need to orient // I wouldn't be so sure...

	//monochromatic_decompose(b, r, t, sp, f1, acc, acc2); //passing sp and f1 is redundant
	printf("skipping monochromatic_decompose\n");

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
	if (acc[v_a] == sp && acc[v_b] != sp){ //leg a is non-empty && leg b is empty
		if (v_a_next == v_b || v_a_next == v_c){ //if the path up the bfs tree from v_a leads to v_b or v_c
			if (acc[v_a] != acc[v_b]){ //if v_a does not touch the cycle defining the subproblem
				//we have one bichromatic problem
				printf("unique subproblem for sp %d is bichromatic\n", sp);
				printf("unique subproblem for sp %d will be on faces %d, %d\n", sp, v_c_mirror, v_a_mirror);
				bichromatic_tripod( b, r, t, v_c_mirror, v_a_mirror, acc, acc2);
			}
			else { //if v_a touches the cycle defining the subproblem //not sure if this can actually happen here
				//we have two monochromatic subproblems
			}
		}
		else {
			//we have two bichromatic problems
			printf("subproblem a for sp %d is bichromatic\n", sp);
			printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_a_r, v_c_mirror);
			bichromatic_tripod( b, r, t, v_a_l, v_a_mirror, acc, acc2);

			printf("subproblem b for sp %d is bichromatic\n", sp);
			printf("subproblem b for sp %d will be on faces %d, %d\n", sp, v_a_l, v_a_mirror);
			bichromatic_tripod( b, r, t, v_a_l, v_a_mirror, acc, acc2);
		}
	}
	else if (acc[v_a] != sp && acc[v_b] == sp){ //leg a is empty && leg b is non-empty
		//we either have two bichromatic problems, or one bichromatic problem
	}
	else if (acc[v_a] != sp && acc[v_b] != sp){ //leg a is empty && leg b is empty
		//we either have one monochromatic problem, two monochromatic problems, or zero problems (our base case)
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

	t->v_a_mirror = b->tri[sp][0]; //check whether that simplex has the correct vertices
	t->v_b_mirror = b->tri[sp][1]; //are these clockwise or counter-clockwise? counter-clockwise.
	t->v_c_mirror = b->tri[sp][2]; //maybe switch v_b_mirror and v_c_mirror?
	//does it matter that these don't switch if switching the v_x_op? (i don't think so, but double check)
	//check whether these are okay when identifing sp
	//the spot where we check for empty subproblems by checking if v_x_op == sp,
	//if we find one, then maybe set v_x_mirror to -1 so that we can easily check whether it is in the subproblem boundary
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
