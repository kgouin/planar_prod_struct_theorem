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

	//initialize tripod_decomposition_struct
	t->vertex_tripod_assign = malloc((b->v)*sizeof(int));
	//label the three vertices incident to face 0 as belonging to three exterior tripods
	t->vertex_tripod_assign[b->sim[0][0]] = (b->f);
	t->vertex_tripod_assign[b->sim[0][1]] = (b->f)+1;
	t->vertex_tripod_assign[b->sim[0][2]] = (b->f)+2;
	for (int k = 3; k < (b->v); k++){
		t->vertex_tripod_assign[k] = -1; //label remaining vertices with -1
	}
	t->face_tripod_assign = malloc((b->f)*sizeof(int));
	t->face_tripod_assign[0] = 0; //this is a special case, since 0 is the outer face
	for (int k = 1; k < (b->f); k++){
		t->face_tripod_assign[k] = -1;
	}
	t->tripod_assign_order = malloc(((b->f)+2)*sizeof(int));
	//label the first three entries as the three exterior tripods that constitute our base case
	t->tripod_assign_order[0] = (b->f);
	t->tripod_assign_order[1] = (b->f)+1;
	t->tripod_assign_order[2] = (b->f)+2;
	for (int k = 3; k < ((b->f)+2); k++){
		t->tripod_assign_order[k] = -1;
	}
	t->tripod_assign_order_index = 3;

	printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	printf("initialization complete\n");

	printf("bfs tree = "); //testing
	for (int k = 0; k < b->v; k++){
		printf("%d ", b->bt[k]);
	}
	printf("\n");
}

int* decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t){

	//start decomposition with the three triangles adjacent to outer face
	trichromatic_tripod(b, r, t, b->tri[0][0], b->tri[0][2], b->tri[0][1]); //these are in counterclockwise order

	printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	printf("decomposition complete\n");

	//t->face_tripod_assign and t->tripod_assign_order exist for testing purposes
	printf("face_tripod_assign = [ ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", t->face_tripod_assign[i]);
	}
	printf("]\n");
	printf("tripod_assign_order = [ ");
	for (int i = 0; i < t->tripod_assign_order_index; i++){
		printf("%d ", t->tripod_assign_order[i]);
	}
	printf("]\n");

	return t->vertex_tripod_assign;
}

/******************************************************** trichromatic ********************************************************/

int* trichromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int f2, int f3){
	//sperner triangle identification
	int sp;
	if (f1 == f2 && f1 == f3){
		sp = f1;
		t->face_tripod_assign[sp] = sp; //keep track of sperner triangles
		//don't add tripod to t->tripod_assign_order
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
	store_tripod(b, t, sp);

	if (t->vertex_tripod_assign[t->v_a] == sp || t->vertex_tripod_assign[t->v_b] == sp || t->vertex_tripod_assign[t->v_c] == sp){ //if one of {leg a, leg b, leg c} is non-empty
		//add tripod to t->tripod_assign_order at index t->tripod_assign_order_index
		t->tripod_assign_order[t->tripod_assign_order_index] = sp;
		t->tripod_assign_order_index++;
	}

	//testing for 3-tree property
	if (!three_tree_test(b, t, sp)){
		printf("3-tree test failed\n");
		exit(0);
	}

	printf("\n%d is one of our sperner triangles\nvertex_tripod_assign = ", sp); //sperner triangles are correctly identified
	for (int k = 0; k < (b->v); k++){
		printf("%d ", t->vertex_tripod_assign[k]);
	}
	printf("\nface_tripod_assign = ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", t->face_tripod_assign[i]);
	}
	printf("\n");

	t->face_tripod_assign[sp] = sp; //keep track of sperner triangles

	if (f1 == sp && f2 == sp && f3 == sp) return 0; //no subproblems exist

	trichromatic_orient_subproblems(b, t, sp, f1, f2, f3);

	tprint(t);

	trichromatic_decompose(b, r, t, sp, f1, f2, f3);

	return 0;
}

void trichromatic_orient_subproblems(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3){
	int leg_colour_a;
	int leg_colour_b;
	int leg_colour_c;
	
	(t->v_a_next == -1) ? (leg_colour_a = t->vertex_tripod_assign[t->v_a]) : (leg_colour_a = t->vertex_tripod_assign[t->v_a_next]);
	(t->v_b_next == -1) ? (leg_colour_b = t->vertex_tripod_assign[t->v_b]) : (leg_colour_b = t->vertex_tripod_assign[t->v_b_next]);
	(t->v_c_next == -1) ? (leg_colour_c = t->vertex_tripod_assign[t->v_c]) : (leg_colour_c = t->vertex_tripod_assign[t->v_c_next]);

	trichromatic_orient_subproblems_pt2(b, t, sp, f1, f2, f3, leg_colour_a, leg_colour_b, leg_colour_c);
}

void trichromatic_orient_subproblems_pt2(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3, int leg_colour_a, int leg_colour_b, int leg_colour_c){
	int f;
	int next_f;
	int next_next_f;

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

	trichromatic_orient_subproblems_pt3(b, t, sp, f1, f2, f3, leg_colour_a, leg_colour_b, leg_colour_c, f, next_f, next_next_f);
}

void trichromatic_orient_subproblems_pt3(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3, int leg_colour_a, int leg_colour_b, int leg_colour_c, int f, int next_f, int next_next_f){
	int double_match = 0;
	for (int i = 0; i < 3; i++){ //look at all three vertices of our f triangle
		if (t->vertex_tripod_assign[b->sim[f][i]] == leg_colour_a){
			for (int j = 0; j < 3; j++){
				if (t->vertex_tripod_assign[b->sim[f][j]] == leg_colour_b){
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
			if (t->vertex_tripod_assign[b->sim[f][i]] == leg_colour_b){
				for (int j = 0; j < 3; j++){
					if (t->vertex_tripod_assign[b->sim[f][j]] == leg_colour_c){
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

void trichromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3){
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
	if (t->vertex_tripod_assign[v_a] == sp && t->vertex_tripod_assign[v_b] == sp){ //leg a is non-empty && leg b is non-empty
		v_a_l = b->il[v_a][b->pin[v_a]];
		(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
		v_b_r = b->il[v_b][y_b];
		printf("subproblem a for sp %d is trichromatic\n", sp);
		printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_b_r, v_a_l);
		trichromatic_tripod(b, r, t, v_a_op, v_b_r, v_a_l);
	}
	else if (t->vertex_tripod_assign[v_a] == sp && t->vertex_tripod_assign[v_b] != sp){ //leg a is non-empty && leg b is empty
		v_a_l = b->il[v_a][b->pin[v_a]];
		printf("subproblem a for sp %d is trichromatic\n", sp);
		printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_a_mirror, v_a_l);
		trichromatic_tripod(b, r, t, v_a_op, v_a_mirror, v_a_l);
	}
	else if (t->vertex_tripod_assign[v_a] != sp && t->vertex_tripod_assign[v_b] == sp){ //leg a is empty && leg b is non-empty
		(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
		v_b_r = b->il[v_b][y_b];
		printf("subproblem a for sp %d is trichromatic\n", sp);
		printf("subproblem a for sp %d will be on faces %d, %d, %d\n", sp, v_a_op, v_b_r, v_a_mirror);
		trichromatic_tripod(b, r, t, v_a_op, v_b_r, v_a_mirror);
	}
	else if (t->vertex_tripod_assign[v_a] != sp && t->vertex_tripod_assign[v_b] != sp){ //leg a is empty && leg b is empty
		if (t->face_tripod_assign[v_a_op] > -1){
			printf("no subproblem a for sp %d\n", sp);
			v_a_mirror = -1; //set v_a_mirror to negative value to later check whether it is in our subproblem boundary
		}
		else {
			printf("subproblem a for sp %d is bichromatic\n", sp);
			printf("subproblem a for sp %d will be on faces %d, %d\n", sp, v_a_op, v_a_mirror);
			bichromatic_tripod( b, r, t, v_a_op, v_a_mirror);
		}
	}

	//subproblem b
	if (t->vertex_tripod_assign[v_b] == sp && t->vertex_tripod_assign[v_c] == sp){ //leg b is non-empty && leg c is non-empty
		v_b_l = b->il[v_b][b->pin[v_b]];
		(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
		v_c_r = b->il[v_c][y_c];
		printf("subproblem b for sp %d is trichromatic\n", sp);
		printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_c_r, v_b_l);
		trichromatic_tripod(b, r, t, v_b_op, v_c_r, v_b_l);
	}
	else if (t->vertex_tripod_assign[v_b] == sp && t->vertex_tripod_assign[v_c] != sp){ //leg b is non-empty && leg c is empty
		v_b_l = b->il[v_b][b->pin[v_b]];
		printf("subproblem b for sp %d is trichromatic\n", sp);
		printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_b_mirror, v_b_l);
		trichromatic_tripod(b, r, t, v_b_op, v_b_mirror, v_b_l);
	}
	else if (t->vertex_tripod_assign[v_b] != sp && t->vertex_tripod_assign[v_c] == sp){ //leg b is empty && leg c is non-empty
		(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
		v_c_r = b->il[v_c][y_c];
		printf("subproblem b for sp %d is trichromatic\n", sp);
		printf("subproblem b for sp %d will be on faces %d, %d, %d\n", sp, v_b_op, v_c_r, v_b_mirror);
		trichromatic_tripod(b, r, t, v_b_op, v_c_r, v_b_mirror);
	}
	else if (t->vertex_tripod_assign[v_b] != sp && t->vertex_tripod_assign[v_c] != sp){ //leg b is empty && leg c is empty
		if (t->face_tripod_assign[v_b_op] > -1){
			printf("no subproblem b for sp %d\n", sp);
			v_b_mirror = -1; //set v_b_mirror to negative value to later check whether it is in our subproblem boundary
		}
		else {
			printf("subproblem b for sp %d is bichromatic\n", sp);
			printf("subproblem b for sp %d will be on faces %d, %d\n", sp, v_b_op, v_b_mirror);
			bichromatic_tripod( b, r, t, v_b_op, v_b_mirror);
		}
	}

	//subproblem c
	if (t->vertex_tripod_assign[v_c] == sp && t->vertex_tripod_assign[v_a] == sp){ //leg c is non-empty && leg a is non-empty
		v_c_l = b->il[v_c][b->pin[v_c]];
		(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
		v_a_r = b->il[v_a][y_a];
		printf("subproblem c for sp %d is trichromatic\n", sp);
		printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_a_r, v_c_l);
		trichromatic_tripod(b, r, t, v_c_op, v_a_r, v_c_l);
	}
	else if (t->vertex_tripod_assign[v_c] == sp && t->vertex_tripod_assign[v_a] != sp){ //leg c is non-empty && leg a is empty
		v_c_l = b->il[v_c][b->pin[v_c]];
		printf("subproblem c for sp %d is trichromatic\n", sp);
		printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_c_mirror, v_c_l);
		trichromatic_tripod(b, r, t, v_c_op, v_c_mirror, v_c_l);
	}
	else if (t->vertex_tripod_assign[v_c] != sp && t->vertex_tripod_assign[v_a] == sp){ //leg c is empty && leg a is non-empty
		(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
		v_a_r = b->il[v_a][y_a];
		printf("subproblem c for sp %d is trichromatic\n", sp);
		printf("subproblem c for sp %d will be on faces %d, %d, %d\n", sp, v_c_op, v_a_r, v_c_mirror);
		trichromatic_tripod(b, r, t, v_c_op, v_a_r, v_c_mirror);
	}
	else if (t->vertex_tripod_assign[v_c] != sp && t->vertex_tripod_assign[v_a] != sp){ //leg c is empty && leg a is empty
		if (t->face_tripod_assign[v_c_op] > -1){
			printf("no subproblem c for sp %d\n", sp);
			v_c_mirror = -1; //set v_c_mirror to negative value to later check whether it is in our subproblem boundary
		}
		else {
			printf("subproblem c for sp %d is bichromatic\n", sp);
			printf("subproblem c for sp %d will be on faces %d, %d\n", sp, v_c_op, v_c_mirror);
			bichromatic_tripod( b, r, t, v_c_op, v_c_mirror);
		}
	}
}

/******************************************************** bichromatic ********************************************************/

int* bichromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1, int f2){
	int sp;

	sp = f1;

	store_tripod(b, t, sp);

	if (t->vertex_tripod_assign[t->v_a] == sp || t->vertex_tripod_assign[t->v_b] == sp || t->vertex_tripod_assign[t->v_c] == sp){ //if one of {leg a, leg b, leg c} is non-empty
		//add tripod to t->tripod_assign_order at index t->tripod_assign_order_index
		t->tripod_assign_order[t->tripod_assign_order_index] = sp;
		t->tripod_assign_order_index++;
	}

	//testing for 3-tree property
	if (!three_tree_test(b, t, sp)){
		printf("3-tree test failed\n");
		exit(0);
	}

	printf("\n%d is one of our sperner triangles\nvertex_tripod_assign = ", sp); //sperner triangles are correctly identified
	for (int k = 0; k < (b->v); k++){
		printf("%d ", t->vertex_tripod_assign[k]);
	}
	printf("\nface_tripod_assign = ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", t->face_tripod_assign[i]);
	}
	printf("\n");

	for (int i = 0; i < (b->f); i++){
		if (t->face_tripod_assign[sp] == sp){
			printf("sperner triangle already found. exiting.\n");
			exit(0);
		}
	}

	t->face_tripod_assign[sp] = sp; //keep track of sperner triangles

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

	bichromatic_decompose(b, r, t, sp, f1, f2);

	return 0;
}

void bichromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1, int f2){
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

	if (t->vertex_tripod_assign[v_a] == sp || t->vertex_tripod_assign[v_b] == sp || t->vertex_tripod_assign[v_c] == sp){ //if one of {leg a, leg b, leg c} is non-empty
		//here we either have case 5.1 or 5.2
		//there can be at most one non-empty leg, since more empty legs lead to trichromatic subproblems
		//the vertices of the crotch of sp are one of exactly three colours
		if (t->vertex_tripod_assign[v_a] == sp){ //v_a is non-empty
			if (v_a_next != v_b && v_a_next != v_c){ //if the path up the bfs tree from v_a does NOT lead to v_b or v_c
				//here we have case 5.1
				(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
				v_a_r = b->il[v_a][y_a];
				v_a_l = b->il[v_a][b->pin[v_a]];
				printf("we have two subproblems for sp %d: trichromatic + bichromatic\n", sp);
				if (t->vertex_tripod_assign[v_a_next] == t->vertex_tripod_assign[v_c]){ //if v_a_next is the same colour as v_c
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_a_r, v_c_mirror);
					bichromatic_tripod( b, r, t, v_a_r, v_c_mirror);
					printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_a_mirror, v_a_l);
					trichromatic_tripod(b, r, t, f2, v_a_mirror, v_a_l);
				}
				else { //if v_a_next is the same colour as v_b
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_a_l, v_a_mirror);
					bichromatic_tripod( b, r, t, v_a_l, v_a_mirror);
					printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_a_r, v_c_mirror);
					trichromatic_tripod(b, r, t, f2, v_a_r, v_c_mirror);
				}
			}
			else {
				//here we have case 5.2
				printf("we have one subproblem for sp %d: trichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d, %d\n", sp, v_a_mirror, v_c_mirror, f2);
				trichromatic_tripod(b, r, t, v_a_mirror, v_c_mirror, f2);
			}
		}
		else if (t->vertex_tripod_assign[v_b] == sp){ //v_b is non-empty
			if (v_b_next != v_a && v_b_next != v_c){ //if the path up the bfs tree from v_b does NOT lead to v_a or v_c
				//here we have case 5.1
				(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
				v_b_r = b->il[v_b][y_b];
				v_b_l = b->il[v_b][b->pin[v_b]];
				printf("we have two subproblems for sp %d: trichromatic + bichromatic\n", sp);
				if (t->vertex_tripod_assign[v_b_next] == t->vertex_tripod_assign[v_a]){ //if v_b_next is the same colour as v_a
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_b_r, v_a_mirror);
					bichromatic_tripod( b, r, t, v_b_r, v_a_mirror);
					printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_b_mirror, v_b_l);
					trichromatic_tripod(b, r, t, f2, v_b_mirror, v_b_l);
				}
				else { //if v_b_next is the same colour as v_c
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_b_l, v_b_mirror);
					bichromatic_tripod( b, r, t, v_b_l, v_b_mirror);
					printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_b_r, v_a_mirror);
					trichromatic_tripod(b, r, t, f2, v_b_r, v_a_mirror);
				}
			}
			else {
				//here we have case 5.2
				printf("we have one subproblem for sp %d: trichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d, %d\n", sp, v_b_mirror, v_a_mirror, f2);
				trichromatic_tripod(b, r, t, v_b_mirror, v_a_mirror, f2);
			}
		}
		else { //v_c is non-empty
			if (v_c_next != v_a && v_c_next != v_b){ //if the path up the bfs tree from v_c does NOT lead to v_a or v_b
				//here we have case 5.1
				(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
				v_c_r = b->il[v_c][y_c];
				v_c_l = b->il[v_c][b->pin[v_c]];
				printf("we have two subproblems for sp %d: trichromatic + bichromatic\n", sp);
				if (t->vertex_tripod_assign[v_c_next] == t->vertex_tripod_assign[v_b]){ //if v_c_next is the same colour as v_b
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_c_r, v_b_mirror);
					bichromatic_tripod( b, r, t, v_c_r, v_b_mirror);
					printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_c_mirror, v_c_l);
					trichromatic_tripod(b, r, t, f2, v_c_mirror, v_c_l);
				}
				else { //if v_c_next is the same colour as v_a
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_c_l, v_c_mirror);
					bichromatic_tripod( b, r, t, v_c_l, v_c_mirror);
					printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_c_r, v_b_mirror);
					trichromatic_tripod(b, r, t, f2, v_c_r, v_b_mirror);
				}
			}
			else {
				//here we have case 5.2
				printf("we have one subproblem for sp %d: trichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d, %d\n", sp, v_c_mirror, v_b_mirror, f2);
				trichromatic_tripod(b, r, t, v_c_mirror, v_b_mirror, f2);
			}
		}
	}

	else { //all legs are empty
		//here we either have case 5.3, 5.4, 5.5, or 5.6 (the empty case)
		//the vertices of the crotch of sp are one of exactly two colours

		if (t->vertex_tripod_assign[v_a] == t->vertex_tripod_assign[v_b]){ //if v_a and v_b are the same colour, with v_c a different colour
			//edge of interest is {v_a, v_b}
			//mirror triangle of interest is v_a_mirror
			
			//v_a_mirror will have two of its vertices the same colour as v_a and v_b
			//v_a_mirror's third vertex is either the same colour as v_a and v_b, or is uncoloured
			
			if ((b->bt[v_b] == v_a || b->bt[v_a] == v_b) && (f1 != f2)){
				//if v_a is a bfs parent of v_b or v_b is a bfs parent of v_a, then we are in case 5.4
				if (t->face_tripod_assign[v_b_mirror] == -1){
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_b_mirror);
					bichromatic_tripod( b, r, t, f2, v_b_mirror);
				}
				else if (t->face_tripod_assign[v_c_mirror] == -1){
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_c_mirror);
					bichromatic_tripod( b, r, t, f2, v_c_mirror);
				}
			}
			//otherwise, we are in case 5.3 or 5.5
			else if (f1 == f2){
				if (t->face_tripod_assign[v_a_mirror] == -1){ //make sure subproblem exists
					//if f1 == f2 then we are in case 5.5
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_a_mirror);
					monochromatic_tripod( b, r, t, v_a_mirror);
				}
			}
			else {
				//we are in case 5.3
				//for any coloured edge, either it's a tree edge or there's a sp on the other side
				if (t->face_tripod_assign[v_a_mirror] == -1){ //make sure monochromatic subproblem exists with (t->face_tripod_assign[v_a_mirror] != t->vertex_tripod_assign[v_a]) check
					printf("we have two subproblems for sp %d: monochromatic + bichromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_a_mirror);
					monochromatic_tripod( b, r, t, v_a_mirror);
				}
				else {
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
				}
				//DIFFERENT CASES HERE, DEPENDING ON SPECIFIC TRIANGLE ORIENTATION !!!!
				if ((b->bt[v_b] == v_c || b->bt[v_c] == v_b) || t->face_tripod_assign[v_b_mirror] > -1){
					if (t->face_tripod_assign[v_c_mirror] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_c_mirror);
						bichromatic_tripod( b, r, t, f2, v_c_mirror);
					}
				}
				else {
					if (t->face_tripod_assign[v_b_mirror] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_b_mirror);
						bichromatic_tripod( b, r, t, f2, v_b_mirror);
					}
				}
			}
		}
		else if (t->vertex_tripod_assign[v_b] == t->vertex_tripod_assign[v_c]){ //if v_b and v_c are the same colour, with v_a a different colour
			//edge of interest is {v_b, v_c}
			//mirror triangle of interest is v_b_mirror
			
			//v_b_mirror will have two of its vertices the same colour as v_b and v_c
			//v_b_mirror's third vertex is either the same colour as v_b and v_c, or is uncoloured
			
			if ((b->bt[v_c] == v_b || b->bt[v_b] == v_c) && (f1 != f2)){
				//if v_b is a bfs parent of v_c or v_c is a bfs parent of v_b, then we are in case 5.4
				if (t->face_tripod_assign[v_c_mirror] == -1){
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_c_mirror);
					bichromatic_tripod( b, r, t, f2, v_c_mirror);
				}
				else if (t->face_tripod_assign[v_a_mirror] == -1){
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_a_mirror);
					bichromatic_tripod( b, r, t, f2, v_a_mirror);
				}
			}
			//otherwise, we are in case 5.3 or 5.5
			else if (f1 == f2){
				if (t->face_tripod_assign[v_b_mirror] == -1){ //make sure subproblem exists
					//if f1 == f2 then we are in case 5.5
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_b_mirror);
					monochromatic_tripod( b, r, t, v_b_mirror);
				}
			}
			else {
				//we are in case 5.3
				//for any coloured edge, either it's a tree edge or there's a sp on the other side
				if (t->face_tripod_assign[v_b_mirror] == -1){ //make sure monochromatic subproblem exists with (t->face_tripod_assign[v_b_mirror] != t->vertex_tripod_assign[v_b]) check
					printf("we have two subproblems for sp %d: bichromatic + monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_b_mirror);
					monochromatic_tripod( b, r, t, v_b_mirror);
				}
				else {
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
				}
				//DIFFERENT CASES HERE, DEPENDING ON SPECIFIC TRIANGLE ORIENTATION !!!!
				if ((b->bt[v_c] == v_a || b->bt[v_a] == v_c) || t->face_tripod_assign[v_c_mirror] > -1){
					if (t->face_tripod_assign[v_a_mirror] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_a_mirror);
						bichromatic_tripod( b, r, t, f2, v_a_mirror);
					}
				}
				else {
					if (t->face_tripod_assign[v_c_mirror] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_c_mirror);
						bichromatic_tripod( b, r, t, f2, v_c_mirror);
					}
				}
			}
		}
		else if (t->vertex_tripod_assign[v_c] == t->vertex_tripod_assign[v_a]){ //if v_c and v_a are the same colour, with v_b a different colour
			//edge of interest is {v_c, v_a}
			//mirror triangle of interest is v_c_mirror
			
			//v_c_mirror will have two of its vertices the same colour as v_c and v_a
			//v_c_mirror's third vertex is either the same colour as v_c and v_a, or is uncoloured
			
			if ((b->bt[v_a] == v_c || b->bt[v_c] == v_a) && (f1 != f2)){
				//if v_c is a bfs parent of v_a or v_a is a bfs parent of v_c, then we are in case 5.4
				if (t->face_tripod_assign[v_a_mirror] == -1){
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_a_mirror);
					bichromatic_tripod( b, r, t, f2, v_a_mirror);
				}
				else if (t->face_tripod_assign[v_b_mirror] == -1){
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_b_mirror);
					bichromatic_tripod( b, r, t, f2, v_b_mirror);
				}
			}
			//otherwise, we are in case 5.3 or 5.5
			else if (f1 == f2){
				if (t->face_tripod_assign[v_c_mirror] == -1){ //make sure subproblem exists
					//if f1 == f2 then we are in case 5.5
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_c_mirror);
					monochromatic_tripod( b, r, t, v_c_mirror);
				}
			}
			else {
				//we are in case 5.3
				//for any coloured edge, either it's a tree edge or there's a sp on the other side
				if (t->face_tripod_assign[v_c_mirror] == -1){ //make sure monochromatic subproblem exists with (t->face_tripod_assign[v_c_mirror] != t->vertex_tripod_assign[v_c]) check
					printf("we have two subproblems for sp %d: bichromatic + monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_c_mirror);
					monochromatic_tripod( b, r, t, v_c_mirror);
				}
				else {
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
				}
				//DIFFERENT CASES HERE, DEPENDING ON SPECIFIC TRIANGLE ORIENTATION !!!!
				if ((b->bt[v_a] == v_b || b->bt[v_b] == v_a) || t->face_tripod_assign[v_a_mirror] > -1){
					if (t->face_tripod_assign[v_b_mirror] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_b_mirror);
						bichromatic_tripod( b, r, t, f2, v_b_mirror);
					}
				}
				else {
					if (t->face_tripod_assign[v_a_mirror] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_a_mirror);
						bichromatic_tripod( b, r, t, f2, v_a_mirror);
					}
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
int* monochromatic_tripod(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int f1){
	int sp = f1;

	store_tripod(b, t, sp);

	if (t->vertex_tripod_assign[t->v_a] == sp || t->vertex_tripod_assign[t->v_b] == sp || t->vertex_tripod_assign[t->v_c] == sp){ //if one of {leg a, leg b, leg c} is non-empty
		//add tripod to t->tripod_assign_order at index t->tripod_assign_order_index
		t->tripod_assign_order[t->tripod_assign_order_index] = sp;
		t->tripod_assign_order_index++;
	}

	//testing for 3-tree property
	if (!three_tree_test(b, t, sp)){
		printf("3-tree test failed\n");
		exit(0);
	}

	printf("\n%d is one of our sperner triangles\nvertex_tripod_assign = ", sp); //sperner triangles are correctly identified
	for (int k = 0; k < (b->v); k++){
		printf("%d ", t->vertex_tripod_assign[k]);
	}
	printf("\nface_tripod_assign = ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", t->face_tripod_assign[i]);
	}
	printf("\n");

	for (int i = 0; i < (b->f); i++){
		if (t->face_tripod_assign[sp] == sp){
			printf("sperner triangle already found. exiting.\n");
			exit(0);
		}
	}

	t->face_tripod_assign[sp] = sp; //keep track of sperner triangles

	printf("v_a = %d\n", t->v_a);
	printf("v_b = %d\n", t->v_b);
	printf("v_c = %d\n", t->v_c);
	printf("v_a_next = %d\n", t->v_a_next);
	printf("v_b_next = %d\n", t->v_b_next);
	printf("v_c_next = %d\n", t->v_c_next);
	printf("v_a_mirror = %d\n", t->v_a_mirror);
	printf("v_b_mirror = %d\n", t->v_b_mirror);
	printf("v_c_mirror = %d\n", t->v_c_mirror);

	monochromatic_decompose(b, r, t, sp, f1); //passing sp and f1 is redundant

	return 0;
}

void monochromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1){
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
	if (t->vertex_tripod_assign[v_a] == sp || t->vertex_tripod_assign[v_b] == sp || t->vertex_tripod_assign[v_c] == sp){ //if one of {leg a, leg b, leg c} is non-empty
		//here we either have case 6.1 or 6.2
		//there can be at most one non-empty leg
		//the vertices of the crotch of sp are one of exactly three colours
		if (t->vertex_tripod_assign[v_a] == sp){ //v_a is non-empty
			if (v_a_next != v_b && v_a_next != v_c){ //if the path up the bfs tree from v_a does NOT lead to v_b or v_c
				//here we have case 6.1
				(b->pin[v_a] == 0) ? (y_a = (b->n[v_a])-1) : (y_a = b->pin[v_a]-1);
				v_a_r = b->il[v_a][y_a];
				v_a_l = b->il[v_a][b->pin[v_a]];
				printf("we have two subproblems for sp %d: bichromatic + bichromatic\n", sp);
				printf("first bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_a_r, v_c_mirror);
				bichromatic_tripod( b, r, t, v_a_r, v_c_mirror);
				printf("second bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_a_l, v_a_mirror);
				bichromatic_tripod( b, r, t, v_a_l, v_a_mirror);

			}
			else {
				//here we have case 6.2
				printf("we have one subproblem for sp %d: bichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d\n", sp, v_c_mirror, v_a_mirror);
				bichromatic_tripod(b, r, t, v_c_mirror, v_a_mirror);
			}
		}
		else if (t->vertex_tripod_assign[v_b] == sp){ //v_b is non-empty
			if (v_b_next != v_a && v_b_next != v_c){ //if the path up the bfs tree from v_b does NOT lead to v_a or v_c
				//here we have case 6.1
				(b->pin[v_b] == 0) ? (y_b = (b->n[v_b])-1) : (y_b = b->pin[v_b]-1);
				v_b_r = b->il[v_b][y_b];
				v_b_l = b->il[v_b][b->pin[v_b]];
				printf("we have two subproblems for sp %d: bichromatic + bichromatic\n", sp);
				printf("first bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_b_r, v_a_mirror);
				bichromatic_tripod( b, r, t, v_b_r, v_a_mirror);
				printf("second bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_b_l, v_b_mirror);
				bichromatic_tripod( b, r, t, v_b_l, v_b_mirror);
			}
			else {
				//here we have case 6.2
				printf("we have one subproblem for sp %d: bichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d\n", sp, v_a_mirror, v_b_mirror);
				bichromatic_tripod(b, r, t, v_a_mirror, v_b_mirror);
			}
		}
		else { //v_c is non-empty
			if (v_c_next != v_a && v_c_next != v_b){ //if the path up the bfs tree from v_c does NOT lead to v_a or v_b
				//here we have case 6.1
				(b->pin[v_c] == 0) ? (y_c = (b->n[v_c])-1) : (y_c = b->pin[v_c]-1);
				v_c_r = b->il[v_c][y_c];
				v_c_l = b->il[v_c][b->pin[v_c]];
				printf("we have two subproblems for sp %d: bichromatic + bichromatic\n", sp);
				printf("first bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_c_r, v_b_mirror);
				bichromatic_tripod( b, r, t, v_c_r, v_b_mirror);
				printf("second bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_c_l, v_c_mirror);
				bichromatic_tripod( b, r, t, v_c_l, v_c_mirror);
			}
			else {
				//here we have case 6.2
				printf("we have one subproblem for sp %d: bichromatic\n", sp);
				printf("subproblem for sp %d will be on faces %d, %d\n", sp, v_b_mirror, v_c_mirror);
				bichromatic_tripod(b, r, t, v_b_mirror, v_c_mirror);
			}
		}
	}
	else { //all legs are empty
		//here we either have case 6.3, 6.4, or 6.5 (the empty case)
		//the vertices of sp are one colour
		if (b->bt[v_b] == v_a || b->bt[v_a] == v_b){ //if v_a is a bfs parent of v_b or v_b is a bfs parent of v_a
			if (b->bt[v_a] == v_b){ //if v_b is a bfs parent of v_a
				//here, if the subproblem is empty, a given edge with purple endpoints is either an edge of the purple tripod, or it is adjacent to a sperner triangle
				if (t->face_tripod_assign[v_c_mirror] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_c_mirror);
					monochromatic_tripod( b, r, t, v_c_mirror);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
			else { //v_a is a bfs parent of v_b
				if (t->face_tripod_assign[v_b_mirror] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_b_mirror);
					monochromatic_tripod( b, r, t, v_b_mirror);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
		}
		else if (b->bt[v_c] == v_b || b->bt[v_b] == v_c){ //if v_b is a bfs parent of v_c or v_c is a bfs parent of v_b
			if (b->bt[v_b] == v_c){ //v_c is a bfs parent of v_b
				if (t->face_tripod_assign[v_a_mirror] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_a_mirror);
					monochromatic_tripod( b, r, t, v_a_mirror);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
			else { //v_b is a bfs parent of v_c
				if (t->face_tripod_assign[v_c_mirror] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_c_mirror);
					monochromatic_tripod( b, r, t, v_c_mirror);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
		}
		else if (b->bt[v_a] == v_c || b->bt[v_c] == v_a){ //if v_c is a bfs parent of v_a or v_a is a bfs parent of v_c
			if (b->bt[v_c] == v_a){ //if v_a is a bfs parent of v_c
				if (t->face_tripod_assign[v_b_mirror] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_b_mirror);
					monochromatic_tripod( b, r, t, v_b_mirror);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
			else { //v_c is a bfs parent of v_a
				if (t->face_tripod_assign[v_a_mirror] == -1){ //if the triangle on which we're about to recurse is NOT a sperner triangle
					//we are in case 6.4
					printf("we have one subproblem for sp %d: monochromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_a_mirror);
					monochromatic_tripod( b, r, t, v_a_mirror);
				}
				//otherwise we are in case 6.5 (the empty case)
			}
		}
		else { //if none of v_a, v_b, v_c are parents to one another
			//here we are in case 6.3
			printf("we have two subproblems for sp %d: monochromatic + monochromatic\n", sp); //problem here
			if (t->face_tripod_assign[v_a_mirror] > -1){ //if v_a_mirror is a previous sp
				//recurse on v_b_mirror and v_c_mirror
				if (t->face_tripod_assign[v_b_mirror] == -1){ //make sure v_b_mirror is not a sp
					printf("first monochromatic subproblem for sp %d will be on face %d\n", sp, v_b_mirror);
					monochromatic_tripod( b, r, t, v_b_mirror);
				}
				if (t->face_tripod_assign[v_c_mirror] == -1){ //make sure v_c_mirror is not a sp
					printf("second monochromatic subproblem for sp %d will be on face %d\n", sp, v_c_mirror);
					monochromatic_tripod( b, r, t, v_c_mirror);
				}
			}
			else if (t->face_tripod_assign[v_b_mirror] > -1){ //if v_b_mirror is a previous sp
				//recurse on v_c_mirror and v_a_mirror
				if (t->face_tripod_assign[v_c_mirror] == -1){ //make sure v_c_mirror is not a sp
					printf("first monochromatic subproblem for sp %d will be on face %d\n", sp, v_c_mirror);
					monochromatic_tripod( b, r, t, v_c_mirror);
				}
				if (t->face_tripod_assign[v_a_mirror] == -1){ //make sure v_a_mirror is not a sp
					printf("second monochromatic subproblem for sp %d will be on face %d\n", sp, v_a_mirror);
					monochromatic_tripod( b, r, t, v_a_mirror);
				}
			}
			else { //v_c_mirror is a previous sp //should this be an else if?
				//recurse on v_a_mirror and v_b_mirror
				if (t->face_tripod_assign[v_a_mirror] == -1){ //make sure v_a_mirror is not a sp
					printf("first monochromatic subproblem for sp %d will be on face %d\n", sp, v_a_mirror);
					monochromatic_tripod( b, r, t, v_a_mirror);
				}
				if (t->face_tripod_assign[v_b_mirror] == -1){ //make sure v_b_mirror is not a sp
					printf("second monochromatic subproblem for sp %d will be on face %d\n", sp, v_b_mirror);
					monochromatic_tripod( b, r, t, v_b_mirror);
				}
			}
		}
	}
}

void store_tripod(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp){
	int u;
	u = b->sim[sp][0];
	t->v_a = u;
	t->v_a_next = -1;
	while (t->vertex_tripod_assign[u] == -1){
		t->vertex_tripod_assign[u] = sp;
		t->v_a = u;
		u = b->bt[u];
		t->v_a_next = u;
	}

	u = b->sim[sp][1];
	t->v_b = u;
	t->v_b_next = -1;
	while (t->vertex_tripod_assign[u] == -1){
		t->vertex_tripod_assign[u] = sp;
		t->v_b = u;
		u = b->bt[u];
		t->v_b_next = u;
	}

	u = b->sim[sp][2];
	t->v_c = u;
	t->v_c_next = -1;
	while (t->vertex_tripod_assign[u] == -1){
		t->vertex_tripod_assign[u] = sp;
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

void tripod_free(struct tripod_decomposition_struct* t){
	free(t->vertex_tripod_assign);
	free(t->face_tripod_assign);
	free(t->tripod_assign_order);
}

int three_tree_test(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp){
	//sanity check for 3-tree property
	//for each vertex v of sp adjacent to another tripod
	//look at the tripod v belongs to
	//should be same tripod as v, or one of v's tripod parents (tripods on the boundary when we found the tripod)
	//all parents have smaller index
	//everyone has at most three parents

	int colour_a;
	int colour_b;
	int colour_c;
	(t->v_a_next == -1) ? (colour_a = t->vertex_tripod_assign[t->v_a]) : (colour_a = t->vertex_tripod_assign[t->v_a_next]);
	(t->v_b_next == -1) ? (colour_b = t->vertex_tripod_assign[t->v_b]) : (colour_b = t->vertex_tripod_assign[t->v_b_next]);
	(t->v_c_next == -1) ? (colour_c = t->vertex_tripod_assign[t->v_c]) : (colour_c = t->vertex_tripod_assign[t->v_c_next]);

	int vertex_index = -1;
	for (int m = 0; m < t->tripod_assign_order_index; m++){ //for each previously found tripod
		if (colour_a == t->tripod_assign_order[m]){
			vertex_index = m;
		}
	}
	if (vertex_index == -1){ //this means the parent tripod doesn't exist
		return 0;
	}
	else if (vertex_index > t->tripod_assign_order_index){ //this shouldn't be possible
		return 0;
	}

	vertex_index = -1;
	for (int m = 0; m < t->tripod_assign_order_index; m++){ //for each previously found tripod
		if (colour_b == t->tripod_assign_order[m]){
			vertex_index = m;
		}
	}
	if (vertex_index == -1){ //this means the parent tripod doesn't exist
		return 0;
	}
	else if (vertex_index > t->tripod_assign_order_index){ //this shouldn't be possible
		return 0;
	}

	vertex_index = -1;
	for (int m = 0; m < t->tripod_assign_order_index; m++){ //for each previously found tripod
		if (colour_c == t->tripod_assign_order[m]){
			vertex_index = m;
		}
	}
	if (vertex_index == -1){ //this means the parent tripod doesn't exist
		return 0;
	}
	else if (vertex_index > t->tripod_assign_order_index){ //this shouldn't be possible
		return 0;
	}

	return 1;
}
