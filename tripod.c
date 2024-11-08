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
	for (int i = 3; i < (b->v); i++){
		t->vertex_tripod_assign[i] = -1; //label remaining vertices with -1
	}
	t->face_tripod_assign = malloc((b->f)*sizeof(int));
	t->face_tripod_assign[0] = 0; //this is a special case, since 0 is the outer face
	for (int i = 1; i < (b->f); i++){
		t->face_tripod_assign[i] = -1;
	}
	t->tripod_assign_order = malloc(((b->f)+2)*sizeof(int));
	//label the first three entries as the three exterior tripods that constitute our base case
	t->tripod_assign_order[0] = (b->f);
	t->tripod_assign_order[1] = (b->f)+1;
	t->tripod_assign_order[2] = (b->f)+2;
	for (int i = 3; i < ((b->f)+2); i++){
		t->tripod_assign_order[i] = -1;
	}
	t->tripod_assign_order_index = 3;

	t->v_x = malloc(3*sizeof(int));
	t->v_x_next = malloc(3*sizeof(int));
	t->v_x_op = malloc(3*sizeof(int));
	t->v_x_mirror = malloc(3*sizeof(int));

	printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	printf("--------------------------------------------------------------------------------------------------------------------------------\n");
	printf("initialization complete\n");

	printf("bfs tree = "); //testing
	for (int i = 0; i < b->v; i++){
		printf("%d ", b->bt[i]);
	}
	printf("\n");

	//write bfs to file for later visualization
	FILE *fd;
	fd = fopen("bfs.txt", "w");
	if (!fd) exit(1);
	fprintf(fd, "%d\n", b->v);
	for (int i = 0; i < b->v; i++){
		fprintf(fd, "%d\n", b->bt[i]);
	}
	fclose(fd);
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
	else { //if two of the three input triangles are the same
		//then they are (it is) the sperner triangle
		//performing LCA on f1, f2, and f3 may lead to an incorrect selection of sp
		if (f1 == f2) sp = f1;
		else sp = f3;
	}

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
	for (int i = 0; i < (b->v); i++){
		printf("%d ", t->vertex_tripod_assign[i]);
	}
	printf("\nface_tripod_assign = ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", t->face_tripod_assign[i]);
	}
	printf("\n");

	for (int i = 0; i < (b->f); i++){
		if (t->face_tripod_assign[f1] == f1 || t->face_tripod_assign[f2] == f2 || t->face_tripod_assign[f3] == f3){
			printf("one of our input triangles is a sperner triangle. exiting.\n");
			printf("tripod_assign_order = [ ");
			for (int i = 0; i < t->tripod_assign_order_index; i++){
				printf("%d ", (t->tripod_assign_order)[i]);
			}
			printf("]\n");
			exit(0);
		}
	}

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

	int* f = malloc(3*sizeof(int)); //memory leak
	f[0] = f1;
	f[1] = f2;
	f[2] = f3;

	int k = 0;
	while (f[k] == sp) k++;
	//f[k], f[(k+1)%3], f[(k+2)%3]
	//rotate f

	printf("f = f%d = %d\n", k+1, f[k]); //temporary, so that diff will show no difference between the original code and the new modular arithmetic code

	int double_match = 0;
	for (int i = 0; i < 3; i++){ //look at all three vertices of our f triangle
		if (t->vertex_tripod_assign[b->sim[f[k]][i]] == leg_colour_a){
			for (int j = 0; j < 3; j++){
				if (t->vertex_tripod_assign[b->sim[f[k]][j]] == leg_colour_b){
					double_match = 1;
				}
			}
		}
	}
	if (double_match){
		printf("ORIENTATION 1\n");

		t->v_a_op = f[k];
		t->v_b_op = f[(k+1)%3];
		t->v_c_op = f[(k+2)%3];

		t->v_x_op[0] = f[k];
		t->v_x_op[1] = f[(k+1)%3];
		t->v_x_op[2] = f[(k+2)%3];
	}
	else { //no double match. reset values and try with leg_colour_b and leg_colour_c
		double_match = 0;
		for (int i = 0; i < 3; i++){ //look at all three vertices of our f triangle
			if (t->vertex_tripod_assign[b->sim[f[k]][i]] == leg_colour_b){
				for (int j = 0; j < 3; j++){
					if (t->vertex_tripod_assign[b->sim[f[k]][j]] == leg_colour_c){
						double_match = 1;
					}
				}
			}
		}
		if (double_match){
			printf("ORIENTATION 2\n");

			t->v_a_op = f[(k+2)%3];
			t->v_b_op = f[k];
			t->v_c_op = f[(k+1)%3];

			t->v_x_op[0] = f[(k+2)%3];
			t->v_x_op[1] = f[k];
			t->v_x_op[2] = f[(k+1)%3];
		}
		else { //no double match. f must match with leg_colour_c and leg_colour_a.
			printf("ORIENTATION 3\n");

			t->v_a_op = f[(k+1)%3];
			t->v_b_op = f[(k+2)%3];
			t->v_c_op = f[k];

			t->v_x_op[0] = f[(k+1)%3];
			t->v_x_op[1] = f[(k+2)%3];
			t->v_x_op[2] = f[k];
		}
	}
}

void trichromatic_decompose(struct bfs_struct* b, struct rmq_struct* r, struct tripod_decomposition_struct* t, int sp, int f1, int f2, int f3){
	//definitions
	int* v_x = malloc(3*sizeof(int)); //mem leak
	int* v_x_next = malloc(3*sizeof(int)); //mem leak
	int* v_x_op = malloc(3*sizeof(int)); //mem leak
	int* v_x_mirror = malloc(3*sizeof(int)); //mem leak
	int* v_x_l = malloc(3*sizeof(int)); //mem leak
	int* v_x_r = malloc(3*sizeof(int)); //mem leak
	for (int i = 0; i < 3; i++){
		v_x[i] = t->v_x[i];
		v_x_next[i] = t->v_x_next[i];
		v_x_op[i] = t->v_x_op[i];
		v_x_mirror[i] = t->v_x_mirror[i];
		v_x_l[i] = b->il[v_x[i]][b->pin[v_x[i]]];
		(b->pin[v_x[i]] == 0) ? (v_x_r[i] = b->il[v_x[i]][(b->n[v_x[i]])-1]) : (v_x_r[i] = b->il[v_x[i]][b->pin[v_x[i]]-1]);
	}

	for (int i = 0; i < 3; i++){
		if (t->vertex_tripod_assign[v_x[i]] == sp && t->vertex_tripod_assign[v_x[(i+1)%3]] == sp){ //leg i is non-empty && leg (i+1)%3 is non-empty
			printf("subproblem x for sp %d is trichromatic\n", sp);
			printf("subproblem x for sp %d will be on faces %d, %d, %d\n", sp, v_x_op[i], v_x_r[(i+1)%3], v_x_l[i]);
			trichromatic_tripod(b, r, t, v_x_op[i], v_x_r[(i+1)%3], v_x_l[i]);
		}
		else if (t->vertex_tripod_assign[v_x[i]] == sp && t->vertex_tripod_assign[v_x[(i+1)%3]] != sp){ //leg i is non-empty && leg (i+1)%3 is empty
			printf("subproblem x for sp %d is trichromatic\n", sp);
			printf("subproblem x for sp %d will be on faces %d, %d, %d\n", sp, v_x_op[i], v_x_mirror[i], v_x_l[i]);
			trichromatic_tripod(b, r, t, v_x_op[i], v_x_mirror[i], v_x_l[i]);
		}
		else if (t->vertex_tripod_assign[v_x[i]] != sp && t->vertex_tripod_assign[v_x[(i+1)%3]] == sp){ //leg i is empty && leg (i+1)%3 is non-empty
			printf("subproblem x for sp %d is trichromatic\n", sp);
			printf("subproblem x for sp %d will be on faces %d, %d, %d\n", sp, v_x_op[i], v_x_r[(i+1)%3], v_x_mirror[i]);
			trichromatic_tripod(b, r, t, v_x_op[i], v_x_r[(i+1)%3], v_x_mirror[i]);
		}
		else if (t->vertex_tripod_assign[v_x[i]] != sp && t->vertex_tripod_assign[v_x[(i+1)%3]] != sp){ //leg i is empty && leg (i+1)%3 is empty
			if (t->face_tripod_assign[v_x_op[i]] > -1){
				printf("no subproblem x for sp %d\n", sp);
				v_x_mirror[i] = -1; //set v_x_mirror[i] to negative value to later check whether it is in our subproblem boundary
			}
			else {
				printf("subproblem x for sp %d is bichromatic\n", sp);
				printf("subproblem x for sp %d will be on faces %d, %d\n", sp, v_x_op[i], v_x_mirror[i]);
				bichromatic_tripod( b, r, t, v_x_op[i], v_x_mirror[i]);
			}
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
	for (int i = 0; i < (b->v); i++){
		printf("%d ", t->vertex_tripod_assign[i]);
	}
	printf("\nface_tripod_assign = ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", t->face_tripod_assign[i]);
	}
	printf("\n");

	for (int i = 0; i < (b->f); i++){
		if (t->face_tripod_assign[f1] == f1 || t->face_tripod_assign[f2] == f2){
			printf("one of our input triangles is a sperner triangle. exiting.\n");
			printf("tripod_assign_order = [ ");
			for (int i = 0; i < t->tripod_assign_order_index; i++){
				printf("%d ", (t->tripod_assign_order)[i]);
			}
			printf("]\n");
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
	int v_a_mirror = t->v_a_mirror;
	int v_b = t->v_b;
	int v_b_next = t->v_b_next;
	int v_b_mirror = t->v_b_mirror;
	int v_c = t->v_c;
	int v_c_next = t->v_c_next;
	int v_c_mirror = t->v_c_mirror;

	int* v_x = malloc(3*sizeof(int)); //mem leak
	int* v_x_next = malloc(3*sizeof(int)); //mem leak
	int* v_x_op = malloc(3*sizeof(int)); //mem leak
	int* v_x_mirror = malloc(3*sizeof(int)); //mem leak
	int* v_x_l = malloc(3*sizeof(int)); //mem leak
	int* v_x_r = malloc(3*sizeof(int)); //mem leak
	for (int i = 0; i < 3; i++){
		v_x[i] = t->v_x[i];
		v_x_next[i] = t->v_x_next[i];
		v_x_mirror[i] = t->v_x_mirror[i];
		v_x_l[i] = b->il[v_x[i]][b->pin[v_x[i]]];
		(b->pin[v_x[i]] == 0) ? (v_x_r[i] = b->il[v_x[i]][(b->n[v_x[i]])-1]) : (v_x_r[i] = b->il[v_x[i]][b->pin[v_x[i]]-1]);
	}

	int v_a_l = b->il[v_a][b->pin[v_a]];
	int v_b_l = b->il[v_b][b->pin[v_b]];
	int v_c_l = b->il[v_c][b->pin[v_c]];
	int v_a_r; int v_b_r; int v_c_r;
	(b->pin[v_a] == 0) ? (v_a_r = b->il[v_a][(b->n[v_a])-1]) : (v_a_r = b->il[v_a][b->pin[v_a]-1]);
	(b->pin[v_b] == 0) ? (v_b_r = b->il[v_b][(b->n[v_b])-1]) : (v_b_r = b->il[v_b][b->pin[v_b]-1]);
	(b->pin[v_c] == 0) ? (v_c_r = b->il[v_c][(b->n[v_c])-1]) : (v_c_r = b->il[v_c][b->pin[v_c]-1]);

	int k;
	for (k = 0; k < 3 && t->vertex_tripod_assign[v_x[k]] != sp; k++){
		//this loop checks who is empty, k == 3 if none
		//if k \in {0, 1, 2} then we're in the first if statement
	}
	//use: k, (k+1)%3, (k+2)%3

	if (k < 3){ //if one of our three legs is non-empty (there can be at most one non-empty leg)
		//here we either have case 5.1 or 5.2
		//there can be at most one non-empty leg, since more empty legs lead to trichromatic subproblems
		//the vertices of the crotch of sp are one of exactly three colours
		if (v_x_next[k] != v_x[(k+1)%3] && v_x_next[k] != v_x[(k+2)%3]){ //if the path up the bfs tree from v_x[k] does NOT lead to v_x[(k+1)%3 or v_x[(k+2)%3]
			//here we have case 5.1
			printf("we have two subproblems for sp %d: trichromatic + bichromatic\n", sp);
			if (t->vertex_tripod_assign[v_x_next[k]] == t->vertex_tripod_assign[v_x[(k+1)%3]]){ //if v_x_next[k] is the same colour as v_x[(k+1)%3]
				printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_x_l[k], v_x_mirror[k]);
				bichromatic_tripod( b, r, t, v_x_l[k], v_x_mirror[k]);
				printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_x_r[k], v_x_mirror[(k+2)%3]);
				trichromatic_tripod(b, r, t, f2, v_x_r[k], v_x_mirror[(k+2)%3]);
			}
			else { //if v_x_next[k] is the same colour as v_x[(k+2)%3]
				printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, v_x_r[k], v_x_mirror[(k+2)%3]);
				bichromatic_tripod( b, r, t, v_x_r[k], v_x_mirror[(k+2)%3]);
				printf("trichromatic subproblem for sp %d will be on faces %d, %d, %d\n", sp, f2, v_x_mirror[k], v_x_l[k]);
				trichromatic_tripod(b, r, t, f2, v_x_mirror[k], v_x_l[k]);
			}
		}
		else {
			//here we have case 5.2
			printf("we have one subproblem for sp %d: trichromatic\n", sp);
			printf("subproblem for sp %d will be on faces %d, %d, %d\n", sp, v_x_mirror[k], v_x_mirror[(k+2)%3], f2);
			trichromatic_tripod(b, r, t, v_x_mirror[k], v_x_mirror[(k+2)%3], f2);
		}
	}
	else { //all legs are empty
		//here we either have case 5.3, 5.4, 5.5, or 5.6 (the empty case)
		//the vertices of the crotch of sp are one of exactly two colours

		int m;
		for (m = 0; m < 3 && t->vertex_tripod_assign[v_x[m]] != t->vertex_tripod_assign[v_x[(k+1)%3]]; m++){
			//this loop checks who is of the same colour
		}
		//use: m, (m+1)%3, (m+2)%3, with m and (m+1)%3 of the same colour

		if (m < 3){ //if m and (m+1)%3 are the same colour, with (m+2)%3 a different colour
			//edge of interest is {v_x[m], v_x[(m+1)%3]}
			//mirror triangle of interest is v_x_mirror[m]
			
			//v_x_mirror[m] will have two of its vertices the same colour as v_x[m] and v_x[(m+1)%3]
			//v_x_mirror[m]'s third vertex is either the same colour as v_x[m] and v_x[(m+1)%3], or is uncoloured
			
			if ((b->bt[v_x[(m+1)%3]] == v_x[m] || b->bt[v_x[m]] == v_x[(m+1)%3]) && (f1 != f2)){
				//if v_x[m] is a bfs parent of v_x[(m+1)%3] or v_x[(m+1)%3] is a bfs parent of v_x[m], then we are in case 5.4
				//problem here. we're not always picking the right mirror triangle
				if (t->face_tripod_assign[v_x_mirror[(m+1)%3]] == -1){
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_x_mirror[(m+1)%3]);
					bichromatic_tripod( b, r, t, f2, v_x_mirror[(m+1)%3]);
				}
				else if (t->face_tripod_assign[v_x_mirror[(m+2)%3]] == -1){
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
					printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_x_mirror[(m+2)%3]);
					bichromatic_tripod( b, r, t, f2, v_x_mirror[(m+2)%3]);
				}
			}
			//otherwise, we are in case 5.3 or 5.5
			else if (f1 == f2){
				if (b->bt[v_x[(m+1)%3]] != v_x[m] && b->bt[v_x[m]] != v_x[(m+1)%3]){ //are we crossing a coloured edge? if not:
					if (t->face_tripod_assign[v_x_mirror[m]] == -1){ //are we about to select a sp on which to recurse? if not:
						//if f1 == f2 then we are in case 5.5
						printf("we have one subproblem for sp %d: monochromatic\n", sp);
						printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_x_mirror[m]);
						monochromatic_tripod( b, r, t, v_x_mirror[m]);
					}
				}
			}
			else {
				//we are in case 5.3
				//for any coloured edge, either it's a tree edge or there's a sp on the other side
				if (t->face_tripod_assign[v_x_mirror[m]] == -1){ //make sure monochromatic subproblem exists
					printf("we have two subproblems for sp %d: monochromatic + bichromatic\n", sp);
					printf("monochromatic subproblem for sp %d will be on face %d\n", sp, v_x_mirror[m]);
					monochromatic_tripod( b, r, t, v_x_mirror[m]);
				}
				else {
					printf("we have one subproblem for sp %d: bichromatic\n", sp);
				}
				//different cases here, depending on specific triangle orientation
				if ((b->bt[v_x[(m+1)%3]] == v_x[(m+2)%3] || b->bt[v_x[(m+2)%3]] == v_x[(m+1)%3]) || t->face_tripod_assign[v_x_mirror[(m+1)%3]] > -1){
					if (t->face_tripod_assign[v_x_mirror[(m+2)%3]] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_x_mirror[(m+2)%3]);
						bichromatic_tripod( b, r, t, f2, v_x_mirror[(m+2)%3]);
					}
				}
				else {
					if (t->face_tripod_assign[v_x_mirror[(m+1)%3]] == -1){
						printf("bichromatic subproblem for sp %d will be on faces %d, %d\n", sp, f2, v_x_mirror[(m+1)%3]);
						bichromatic_tripod( b, r, t, f2, v_x_mirror[(m+1)%3]);
					}
				}
			}
		}
		else { // if m == 3
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
	for (int i = 0; i < (b->v); i++){
		printf("%d ", t->vertex_tripod_assign[i]);
	}
	printf("\nface_tripod_assign = ");
	for (int i = 0; i < b->f; i++){
		printf("%d ", t->face_tripod_assign[i]);
	}
	printf("\n");

	for (int i = 0; i < (b->f); i++){
		if (t->face_tripod_assign[sp] == sp){
			printf("sperner triangle already found. exiting.\n");
			printf("tripod_assign_order = [ ");
			for (int i = 0; i < t->tripod_assign_order_index; i++){
				printf("%d ", (t->tripod_assign_order)[i]);
			}
			printf("]\n");
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
	int v_a_mirror = t->v_a_mirror;
	int v_b = t->v_b;
	int v_b_next = t->v_b_next;
	int v_b_mirror = t->v_b_mirror;
	int v_c = t->v_c;
	int v_c_next = t->v_c_next;
	int v_c_mirror = t->v_c_mirror;

	int* v_x = malloc(3*sizeof(int)); //mem leak
	int* v_x_next = malloc(3*sizeof(int)); //mem leak
	int* v_x_op = malloc(3*sizeof(int)); //mem leak
	int* v_x_mirror = malloc(3*sizeof(int)); //mem leak
	int* v_x_l = malloc(3*sizeof(int)); //mem leak
	int* v_x_r = malloc(3*sizeof(int)); //mem leak
	for (int i = 0; i < 3; i++){
		v_x[i] = t->v_x[i];
		v_x_next[i] = t->v_x_next[i];
		v_x_mirror[i] = t->v_x_mirror[i];
		v_x_l[i] = b->il[v_x[i]][b->pin[v_x[i]]];
		(b->pin[v_x[i]] == 0) ? (v_x_r[i] = b->il[v_x[i]][(b->n[v_x[i]])-1]) : (v_x_r[i] = b->il[v_x[i]][b->pin[v_x[i]]-1]);
	}

	int v_a_l = b->il[v_a][b->pin[v_a]];
	int v_b_l = b->il[v_b][b->pin[v_b]];
	int v_c_l = b->il[v_c][b->pin[v_c]];
	int v_a_r; int v_b_r; int v_c_r;
	(b->pin[v_a] == 0) ? (v_a_r = b->il[v_a][(b->n[v_a])-1]) : (v_a_r = b->il[v_a][b->pin[v_a]-1]);
	(b->pin[v_b] == 0) ? (v_b_r = b->il[v_b][(b->n[v_b])-1]) : (v_b_r = b->il[v_b][b->pin[v_b]-1]);
	(b->pin[v_c] == 0) ? (v_c_r = b->il[v_c][(b->n[v_c])-1]) : (v_c_r = b->il[v_c][b->pin[v_c]-1]);

	//all subproblems
	if (t->vertex_tripod_assign[v_a] == sp || t->vertex_tripod_assign[v_b] == sp || t->vertex_tripod_assign[v_c] == sp){ //if one of {leg a, leg b, leg c} is non-empty
		//here we either have case 6.1 or 6.2
		//there can be at most one non-empty leg
		//the vertices of the crotch of sp are one of exactly three colours
		if (t->vertex_tripod_assign[v_a] == sp){ //v_a is non-empty
			if (v_a_next != v_b && v_a_next != v_c){ //if the path up the bfs tree from v_a does NOT lead to v_b or v_c
				//here we have case 6.1
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
			else if (t->face_tripod_assign[v_c_mirror] > -1){ //v_c_mirror is a previous sp
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
			else {
				printf("else statement which we should never reach!!!!\n"); //assert k is not 3
			}
		}
	}
}

void store_tripod(struct bfs_struct* b, struct tripod_decomposition_struct* t, int sp){
	int u;

	/*
	for (int i = 0; i < 3; i++){
		u = b->sim[sp][i];
		t->v_x[i] = u;
		t->v_x_next[i] = -1;
		while (t->vertex_tripod_assign[u] == -1){
			t->vertex_tripod_assign[u] = sp;
			t->v_x[i] = u;
			u = b->bt[u];
			t->v_x_next[i] = u;
		}
		t->v_x_mirror[i] = b->tri[sp][i];
	}
	*///while there exists a portion of code using non modular arithmetic, this commented block has to remain commented
	//and the semi-duplicated code below has to remain as is

	u = b->sim[sp][0];
	t->v_a = u;
	t->v_x[0] = u;
	t->v_a_next = -1;
	t->v_x_next[0] = -1;
	while (t->vertex_tripod_assign[u] == -1){
		t->vertex_tripod_assign[u] = sp;
		t->v_a = u;
		t->v_x[0] = u;
		u = b->bt[u];
		t->v_a_next = u;
		t->v_x_next[0] = u;
	}

	u = b->sim[sp][1];
	t->v_b = u;
	t->v_x[1] = u;
	t->v_b_next = -1;
	t->v_x_next[1] = -1;
	while (t->vertex_tripod_assign[u] == -1){
		t->vertex_tripod_assign[u] = sp;
		t->v_b = u;
		t->v_x[1] = u;
		u = b->bt[u];
		t->v_b_next = u;
		t->v_x_next[1] = u;
	}

	u = b->sim[sp][2];
	t->v_c = u;
	t->v_x[2] = u;
	t->v_c_next = -1;
	t->v_x_next[2] = -1;
	while (t->vertex_tripod_assign[u] == -1){
		t->vertex_tripod_assign[u] = sp;
		t->v_c = u;
		t->v_x[2] = u;
		u = b->bt[u];
		t->v_c_next = u;
		t->v_x_next[2] = u;
	}

	t->v_a_mirror = b->tri[sp][0];
	t->v_b_mirror = b->tri[sp][1];
	t->v_c_mirror = b->tri[sp][2];
	t->v_x_mirror[0] = b->tri[sp][0];
	t->v_x_mirror[1] = b->tri[sp][1];
	t->v_x_mirror[2] = b->tri[sp][2];
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
	free(t->v_x);
	free(t->v_x_next);
	free(t->v_x_op);
	free(t->v_x_mirror);
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
	for (int i = 0; i < t->tripod_assign_order_index; i++){ //for each previously found tripod
		if (colour_a == t->tripod_assign_order[i]){
			vertex_index = i;
		}
	}
	if (vertex_index == -1){ //this means the parent tripod doesn't exist
		return 0;
	}
	else if (vertex_index > t->tripod_assign_order_index){ //this shouldn't be possible
		return 0;
	}

	vertex_index = -1;
	for (int i = 0; i < t->tripod_assign_order_index; i++){ //for each previously found tripod
		if (colour_b == t->tripod_assign_order[i]){
			vertex_index = i;
		}
	}
	if (vertex_index == -1){ //this means the parent tripod doesn't exist
		return 0;
	}
	else if (vertex_index > t->tripod_assign_order_index){ //this shouldn't be possible
		return 0;
	}

	vertex_index = -1;
	for (int i = 0; i < t->tripod_assign_order_index; i++){ //for each previously found tripod
		if (colour_c == t->tripod_assign_order[i]){
			vertex_index = i;
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

int three_tree_test_pt2(struct bfs_struct* b, struct tripod_decomposition_struct* t){
	return 1;
}

//something else to check (with a new function)
//for each vertex make sure the parents are okay, and that a "parent" is not in fact a "grand-parent"
//use the three_tree_test function to keep a list of vertices and their parents (like an ajdacency list but for ancestry)
//then, at the end of the program, when the recursion terminates, call the new function
//this new function is similar to three_tree_test, but it should catch an extra case where v has a parent beyond the boundary of its three defining tripods
//so, for each vertex v, look at the parent of tripod v belongs to
//this should be same tripod as v, or one of v's tripod parents (tripods on the boundary when we found the tripod)
//all parents have smaller index
//everyone has at most three parents